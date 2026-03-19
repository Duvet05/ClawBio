#!/usr/bin/env python3
"""clinical-trial-finder: Find clinical trials via ClinicalTrials.gov API v2.

Optionally enriches gene symbols with disease associations from OpenTargets
before querying, giving gene → disease → trial chains grounded in evidence.
"""

import argparse
import hashlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
DEMO_DATA = SKILL_DIR / "demo_input.txt"

CT_API = "https://clinicaltrials.gov/api/v2/studies"
CT_FIELDS = ",".join([
    "NCTId", "BriefTitle", "OverallStatus", "Phase",
    "StartDate", "CompletionDate", "StudyType", "BriefSummary",
    "Condition", "InterventionName", "ConditionMeshId", "ConditionMeshTerm",
])
DEFAULT_PAGE_SIZE = 20

DISCLAIMER = (
    "*ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions.*"
)

# FHIR R4 value set mappings — source: hl7.org/fhir/R4/valueset-research-study-status.html
FHIR_STATUS: dict[str, str] = {
    "RECRUITING":             "active",
    "ACTIVE_NOT_RECRUITING":  "active",
    "NOT_YET_RECRUITING":     "approved",
    "COMPLETED":              "completed",
    "TERMINATED":             "administratively-completed",
    "WITHDRAWN":              "withdrawn",
    "SUSPENDED":              "temporarily-closed-to-accrual",
    "UNKNOWN":                "unknown",
}

FHIR_PHASE: dict[str, str] = {
    "PHASE1":          "phase-1",
    "PHASE2":          "phase-2",
    "PHASE3":          "phase-3",
    "PHASE4":          "phase-4",
    "PHASE1 / PHASE2": "phase-1-phase-2",
    "PHASE2 / PHASE3": "phase-2-phase-3",
    "NA":              "n-a",
    "N/A":             "n-a",   # CT.gov API v2 returns "N/A" (with slash)
}

FHIR_PHASE_DISPLAY: dict[str, str] = {
    "PHASE1":          "Phase 1",
    "PHASE2":          "Phase 2",
    "PHASE3":          "Phase 3",
    "PHASE4":          "Phase 4",
    "PHASE1 / PHASE2": "Phase 1/Phase 2",
    "PHASE2 / PHASE3": "Phase 2/Phase 3",
    "NA":              "N/A",
    "N/A":             "N/A",
}

_MONTH_ABBR: dict[str, str] = {
    "January": "01", "February": "02", "March": "03", "April": "04",
    "May": "05", "June": "06", "July": "07", "August": "08",
    "September": "09", "October": "10", "November": "11", "December": "12",
}

ALL_STATUSES = [
    "RECRUITING", "ACTIVE_NOT_RECRUITING", "NOT_YET_RECRUITING",
    "COMPLETED", "TERMINATED", "WITHDRAWN", "SUSPENDED", "UNKNOWN",
]

STATUS_EMOJI: dict[str, str] = {
    "RECRUITING":            "🟢",
    "ACTIVE_NOT_RECRUITING": "🟡",
    "COMPLETED":             "✅",
    "NOT_YET_RECRUITING":    "⏳",
    "TERMINATED":            "🔴",
    "WITHDRAWN":             "⚫",
    "SUSPENDED":             "🟠",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_RE_DATE_ISO   = re.compile(r"^\d{4}-\d{2}(-\d{2})?$")
_RE_DATE_MONTH = re.compile(r"^(\w+)\s+(\d{4})$")


def _to_fhir_date(date_str: str) -> str:
    """Normalize a ClinicalTrials.gov date to FHIR-valid format.

    Handles: "2024-01-15", "2024-01", "January 2024" → YYYY-MM[-DD].
    Returns the original string unchanged if it cannot be parsed.
    """
    if not date_str:
        return ""
    if _RE_DATE_ISO.match(date_str):
        return date_str
    m = _RE_DATE_MONTH.match(date_str)
    if m:
        month_num = _MONTH_ABBR.get(m.group(1))
        if month_num:
            return f"{m.group(2)}-{month_num}"
    return date_str


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def parse_input(input_path: Path) -> dict:
    """Read a query file — one search term per line, # lines are comments."""
    lines = input_path.read_text().splitlines()
    terms = [line.strip() for line in lines if line.strip() and not line.startswith("#")]
    if not terms:
        raise ValueError(f"No search terms in {input_path}. Add at least one non-comment line.")
    return {"query": " ".join(terms), "terms": terms}


# ---------------------------------------------------------------------------
# ClinicalTrials.gov
# ---------------------------------------------------------------------------

def fetch_trials(query: str, max_results: int = DEFAULT_PAGE_SIZE) -> list[dict]:
    """Query ClinicalTrials.gov API v2. Returns normalised trial records."""
    params = urllib.parse.urlencode({
        "query.cond": query,
        "pageSize":   max_results,
        "fields":     CT_FIELDS,
        "format":     "json",
    })
    url = f"{CT_API}?{params}"

    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read().decode())
    except urllib.error.URLError as exc:
        raise RuntimeError(f"ClinicalTrials.gov unreachable: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"ClinicalTrials.gov returned malformed JSON: {exc}") from exc

    return [_normalise_trial(s) for s in data.get("studies", [])]


def _normalise_trial(study: dict) -> dict:
    """Extract flat record from a ClinicalTrials.gov study object."""
    proto      = study.get("protocolSection", {})
    id_mod     = proto.get("identificationModule", {})
    status_mod = proto.get("statusModule", {})
    design_mod = proto.get("designModule", {})
    desc_mod   = proto.get("descriptionModule", {})
    cond_mod   = proto.get("conditionsModule", {})
    interv_mod = proto.get("armsInterventionsModule", {})
    derived    = study.get("derivedSection", {})

    summary = desc_mod.get("briefSummary", "")
    interventions = [
        arm["name"] for arm in interv_mod.get("interventions", []) if arm.get("name")
    ]
    # MeSH-coded conditions from CT.gov's own derivedSection — no extra API call
    condition_meshes = derived.get("conditionBrowseModule", {}).get("meshes", [])

    return {
        "nct_id":            id_mod.get("nctId", ""),
        "title":             id_mod.get("briefTitle", ""),
        "status":            status_mod.get("overallStatus", "UNKNOWN"),
        "phase":             " / ".join(design_mod.get("phases", [])),
        "study_type":        design_mod.get("studyType", ""),
        "start_date":        _to_fhir_date(status_mod.get("startDateStruct", {}).get("date", "")),
        "completion_date":   _to_fhir_date(status_mod.get("completionDateStruct", {}).get("date", "")),
        "conditions":        cond_mod.get("conditions", []),
        "condition_meshes":  condition_meshes,
        "interventions":     interventions,
        "summary":           summary[:300] + "…" if len(summary) > 300 else summary,
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _count_recruiting(trials: list[dict]) -> int:
    return sum(1 for t in trials if t["status"] == "RECRUITING")


def write_report(
    query_info: dict,
    trials: list[dict],
    output_dir: Path,
    gene_context: dict | None = None,
) -> Path:
    """Write markdown report. gene_context is set when --gene was used."""
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "report.md"

    recruiting  = _count_recruiting(trials)
    phases      = sorted({t["phase"] for t in trials if t["phase"]})
    study_types = sorted({t["study_type"] for t in trials if t["study_type"]})
    timestamp   = datetime.now().isoformat(timespec="seconds")

    lines = ["# Clinical Trial Finder Report", ""]

    if gene_context:
        lines += [
            f"**Gene**: `{gene_context['symbol']}` — {gene_context['name']}  ",
            f"**Via**: OpenTargets Platform (association score ≥ {gene_context['min_score']})  ",
            f"**Associated diseases queried**: {', '.join(gene_context['diseases'])}  ",
        ]
    else:
        lines += [f"**Query**: `{query_info['query']}`  "]

    lines += [
        f"**Generated**: {timestamp}  ",
        f"**Source**: ClinicalTrials.gov API v2  ",
        f"**Trials found**: {len(trials)} ({recruiting} recruiting)",
        "",
        "---",
        "",
        "## Summary",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| Total trials | {len(trials)} |",
        f"| Recruiting now | {recruiting} |",
        f"| Phases | {', '.join(phases) or 'N/A'} |",
        f"| Study types | {', '.join(study_types) or 'N/A'} |",
        "",
        "## Trials",
        "",
    ]

    for t in trials:
        emoji   = STATUS_EMOJI.get(t["status"], "⬜")
        nct_url = f"https://clinicaltrials.gov/study/{t['nct_id']}"
        lines += [
            f"### {emoji} {t['title']}",
            "",
            f"- **NCT ID**: [{t['nct_id']}]({nct_url})",
            f"- **Status**: {t['status']}",
            f"- **Phase**: {t['phase'] or 'N/A'}",
            f"- **Type**: {t['study_type'] or 'N/A'}",
            f"- **Start**: {t['start_date'] or 'N/A'} | **Est. completion**: {t['completion_date'] or 'N/A'}",
        ]
        if t["conditions"]:
            lines.append(f"- **Conditions**: {', '.join(t['conditions'][:3])}")
        if t["interventions"]:
            lines.append(f"- **Interventions**: {', '.join(t['interventions'][:3])}")
        if t["summary"]:
            lines += ["", f"> {t['summary']}"]
        lines.append("")

    lines += [
        "---", "",
        "## Reproducibility", "",
        "- `commands.sh` -- exact command to reproduce",
        "- `checksums.sha256` -- SHA-256 of all outputs",
        "",
        "---", "",
        DISCLAIMER,
    ]
    report_path.write_text("\n".join(lines))
    return report_path


def write_summary(query_info: dict, trials: list[dict], output_dir: Path) -> Path:
    """Write machine-readable JSON summary."""
    payload = {
        "query":      query_info["query"],
        "timestamp":  datetime.now().isoformat(),
        "source":     "clinicaltrials.gov/api/v2",
        "total":      len(trials),
        "recruiting": _count_recruiting(trials),
        "trials":     trials,
    }
    path = output_dir / "summary.json"
    path.write_text(json.dumps(payload, indent=2))
    return path


def write_fhir_bundle(trials: list[dict], output_dir: Path) -> Path:
    """Write a FHIR R4 Bundle of ResearchStudy resources."""
    entries = [_trial_to_fhir(t) for t in trials]
    bundle = {
        "resourceType": "Bundle",
        "type":         "searchset",
        "timestamp":    datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00"),
        "total":        len(entries),
        "entry":        entries,
        "meta": {"tag": [{"system": "https://clawbio.ai", "code": "clinical-trial-finder"}]},
    }
    path = output_dir / "fhir_bundle.json"
    path.write_text(json.dumps(bundle, indent=2))
    return path


def _trial_to_fhir(t: dict) -> dict:
    """Map a normalised trial record to a FHIR R4 ResearchStudy Bundle entry."""
    resource: dict = {
        "resourceType": "ResearchStudy",
        "id":           t["nct_id"],
        "meta":         {"profile": ["http://hl7.org/fhir/StructureDefinition/ResearchStudy"]},
        "identifier":   [{"use": "official", "system": "https://clinicaltrials.gov", "value": t["nct_id"]}],
        "title":        t["title"],
        "status":       FHIR_STATUS.get(t["status"], "unknown"),
        "phase": {
            "coding": [{
                "system":  "http://terminology.hl7.org/CodeSystem/research-study-phase",
                "code":    FHIR_PHASE.get(t["phase"], "n-a"),
                "display": FHIR_PHASE_DISPLAY.get(t["phase"], ""),
            }]
        },
    }
    if t["study_type"]:
        resource["category"] = [{
            "coding": [{
                "system":  "http://clinicaltrials.gov/study-type",
                "code":    t["study_type"],
                "display": t["study_type"].replace("_", " ").title(),
            }]
        }]
    if t.get("condition_meshes"):
        # Prefer MeSH-coded conditions (authoritative, from CT.gov's own derivedSection)
        resource["condition"] = [
            {
                "coding": [{
                    "system":  "http://id.nlm.nih.gov/mesh/",
                    "code":    m["id"],
                    "display": m.get("term", ""),
                }],
                "text": m.get("term", m["id"]),
            }
            for m in t["condition_meshes"]
            if m.get("id")
        ]
    elif t["conditions"]:
        resource["condition"] = [{"text": c} for c in t["conditions"]]
    if t["interventions"]:
        resource["focus"] = [{"text": i} for i in t["interventions"]]
    if t["summary"]:
        resource["description"] = t["summary"]
    if t["start_date"]:
        resource["period"] = {"start": t["start_date"]}
        if t["completion_date"]:
            resource["period"]["end"] = t["completion_date"]

    return {
        "fullUrl":  f"https://clinicaltrials.gov/study/{t['nct_id']}",
        "resource": resource,
    }


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

_PHASE_DISPLAY: dict[str, str] = {
    "PHASE1": "Phase 1", "PHASE2": "Phase 2",
    "PHASE3": "Phase 3", "PHASE4": "Phase 4",
    "PHASE1 / PHASE2": "Phase 1/2", "PHASE2 / PHASE3": "Phase 2/3",
    "NA": "N/A",
    # Note: "N/A" (with slash) omitted — fallback `t["phase"] or "Unknown"` handles it
}

_STATUS_COLOR: dict[str, str] = {
    "RECRUITING":            "#2ecc71",
    "NOT_YET_RECRUITING":    "#3498db",
    "ACTIVE_NOT_RECRUITING": "#f1c40f",
    "SUSPENDED":             "#e67e22",
    "TERMINATED":            "#e74c3c",
    "WITHDRAWN":             "#7f8c8d",
    "COMPLETED":             "#95a5a6",
    "UNKNOWN":               "#bdc3c7",
}

_PHASE_ORDER = [
    "Phase 1", "Phase 1/2", "Phase 2", "Phase 2/3",
    "Phase 3", "Phase 4", "N/A", "Unknown",
]


def write_phase_chart(trials: list[dict], output_dir: Path, title: str = "") -> Path | None:
    """Write a stacked bar chart of trial counts by phase, coloured by status.

    Returns the Path to the PNG, or None if matplotlib is not installed.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    # Aggregate: phase → {status: count}
    phase_status: dict[str, Counter] = defaultdict(Counter)
    for t in trials:
        phase = _PHASE_DISPLAY.get(t["phase"], t["phase"] or "Unknown")
        phase_status[phase][t["status"]] += 1

    # Stable phase order
    phases = [p for p in _PHASE_ORDER if p in phase_status]
    phases += [p for p in phase_status if p not in phases]
    if not phases:
        return None

    # O(P) pass — collect which statuses actually appear, then filter in order
    present_statuses = {s for p in phases for s in phase_status[p]}
    all_statuses = [s for s in _STATUS_COLOR if s in present_statuses]

    fig, ax = plt.subplots(figsize=(max(7, len(phases) * 1.4), 5))
    bottoms = [0] * len(phases)

    for status in all_statuses:
        counts = [phase_status[p].get(status, 0) for p in phases]
        ax.bar(
            phases, counts, bottom=bottoms,
            color=_STATUS_COLOR[status],
            label=status.replace("_", " ").title(),
            edgecolor="white", linewidth=0.5,
        )
        bottoms = [b + c for b, c in zip(bottoms, counts)]

    # Total labels on top — proportional offset so they're visible at any scale
    for i, total in enumerate(bottoms):
        if total:
            ax.text(i, total * 1.02, str(total), ha="center", va="bottom",
                    fontsize=9, fontweight="bold")

    ax.set_xlabel("Phase", fontsize=11)
    ax.set_ylabel("Number of Trials", fontsize=11)
    ax.set_title(title or "Clinical Trial Phase Distribution", fontsize=13,
                 fontweight="bold", pad=14)
    ax.set_ylim(0, max(bottoms, default=1) * 1.25 + 1)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.85)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    path = fig_dir / "phase_distribution.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    """Return hex SHA-256 digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def write_commands(args: argparse.Namespace, output_dir: Path) -> Path:
    """Write commands.sh with the exact CLI invocation to reproduce the run."""
    parts = ["python skills/clinical-trial-finder/clinical_trial_finder.py"]
    if args.demo:
        parts.append("--demo")
    elif args.input:
        parts.append(f"--input {args.input}")
    elif args.query:
        parts.append(f"--query \"{args.query}\"")
    elif args.gene:
        parts.append(f"--gene {args.gene}")
    if args.status:
        parts.append(f"--status {args.status}")
    if args.max_results != DEFAULT_PAGE_SIZE:
        parts.append(f"--max-results {args.max_results}")
    if args.fhir:
        parts.append("--fhir")
    if args.ot_min_score != 0.6:
        parts.append(f"--ot-min-score {args.ot_min_score}")
    if args.ot_max_diseases != 5:
        parts.append(f"--ot-max-diseases {args.ot_max_diseases}")
    parts.append(f"--output {output_dir}")
    path = output_dir / "commands.sh"
    path.write_text(" \\\n  ".join(parts) + "\n")
    return path


def write_checksums(output_dir: Path) -> Path:
    """Write checksums.sha256 with SHA-256 digests of all generated outputs."""
    targets = [
        "report.md", "summary.json", "fhir_bundle.json",
        "figures/phase_distribution.png",
    ]
    lines = []
    for rel in targets:
        p = output_dir / rel
        if p.exists():
            lines.append(f"{_sha256(p)}  {rel}")
    path = output_dir / "checksums.sha256"
    path.write_text("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Clinical Trial Finder — ClinicalTrials.gov API v2 + OpenTargets"
    )
    src = p.add_mutually_exclusive_group()
    src.add_argument("--input", type=Path, help="Query file (one search term per line)")
    src.add_argument("--query", type=str,  help="Direct search query string")
    src.add_argument("--gene",  type=str,  help="Gene symbol (e.g. BRCA1) — enriched via OpenTargets")
    src.add_argument("--demo",  action="store_true", help="Run with built-in demo data (BRCA1)")

    p.add_argument("--output",      type=Path, default=Path("/tmp/clinical_trial_finder_output"))
    p.add_argument("--max-results", type=int,  default=DEFAULT_PAGE_SIZE)
    p.add_argument("--status",      type=str,  default=None, choices=ALL_STATUSES,
                   help="Filter trials by recruitment status (default: show all)")
    p.add_argument("--fhir",        action="store_true", help="Also write fhir_bundle.json")
    p.add_argument("--ot-min-score", type=float, default=0.6,
                   help="OpenTargets association score threshold (default: 0.6)")
    p.add_argument("--ot-max-diseases", type=int, default=5,
                   help="Max diseases to query per gene via OpenTargets (default: 5)")
    return p


def main() -> None:
    args = _build_parser().parse_args()

    gene_context: dict | None = None

    if args.demo:
        query_info = parse_input(DEMO_DATA)

    elif args.input:
        query_info = parse_input(args.input)

    elif args.query:
        query_info = {"query": args.query, "terms": [args.query]}

    elif args.gene:
        import opentargets  # local module — only imported when --gene is used
        symbol = args.gene.upper()
        print(f"Resolving {symbol!r} via OpenTargets…")
        ensembl_id, gene_name = opentargets.resolve_gene(symbol)
        diseases = opentargets.get_diseases(
            ensembl_id,
            min_score=args.ot_min_score,
            max_results=args.ot_max_diseases,
        )
        if not diseases:
            raise SystemExit(
                f"No disease associations found for {symbol} "
                f"with score ≥ {args.ot_min_score}. Try --ot-min-score 0.3"
            )
        disease_names = [d.name for d in diseases]
        print(f"✓ {symbol} → {len(diseases)} diseases: {', '.join(disease_names)}")

        # Query each disease individually, deduplicate by NCT ID
        per_disease = max(1, args.max_results // len(diseases))
        seen: set[str] = set()
        trials = []
        for disease in disease_names:
            print(f"  Querying trials for: {disease!r}")
            for t in fetch_trials(disease, max_results=per_disease):
                if t["nct_id"] not in seen:
                    seen.add(t["nct_id"])
                    trials.append(t)

        query_info = {"query": symbol, "terms": disease_names}
        gene_context = {
            "symbol":    symbol,
            "name":      gene_name,
            "ensembl":   ensembl_id,
            "diseases":  disease_names,
            "min_score": args.ot_min_score,
        }
        recruiting = _count_recruiting(trials)
        print(f"✓ Found {len(trials)} unique trials ({recruiting} recruiting)")

    else:
        _build_parser().error("Provide --input, --query, --gene, or --demo")

    if not args.gene:
        print(f"Querying ClinicalTrials.gov: {query_info['query']!r}")
        trials = fetch_trials(query_info["query"], max_results=args.max_results)
        recruiting = _count_recruiting(trials)
    print(f"✓ Found {len(trials)} trials ({recruiting} recruiting)")

    if args.status:
        trials = [t for t in trials if t["status"] == args.status]
        print(f"✓ Filtered by {args.status}: {len(trials)} trials remaining")

    args.output.mkdir(parents=True, exist_ok=True)

    report = write_report(query_info, trials, args.output, gene_context)
    print(f"✓ Report  → {report}")

    summary = write_summary(query_info, trials, args.output)
    print(f"✓ Summary → {summary}")

    if args.fhir:
        bundle = write_fhir_bundle(trials, args.output)
        print(f"✓ FHIR R4 → {bundle}")

    chart_title = (
        f"Phase Distribution — {gene_context['symbol']}" if gene_context
        else f"Phase Distribution — {query_info['query'][:40]}"
    )
    chart = write_phase_chart(trials, args.output, title=chart_title)
    if chart:
        print(f"✓ Chart   → {chart}")

    cmds = write_commands(args, args.output)
    print(f"✓ Repro   → {cmds}")

    checksums = write_checksums(args.output)
    print(f"✓ SHA-256 → {checksums}")


if __name__ == "__main__":
    main()
