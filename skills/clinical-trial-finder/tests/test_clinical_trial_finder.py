"""Tests for clinical-trial-finder — no network calls."""

import argparse
import importlib.util
import json
import sys
import urllib.error
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Load modules from skill directory without installing them as packages
SKILL_DIR = Path(__file__).resolve().parent.parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SKILL_DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ctf = _load("clinical_trial_finder")
ot  = _load("opentargets")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MOCK_TRIALS = [
    {
        "nct_id": "NCT00000001",
        "title": "Synthetic Trial A",
        "status": "RECRUITING",
        "phase": "PHASE2",
        "study_type": "INTERVENTIONAL",
        "start_date": "2024-01",
        "completion_date": "2026-12",
        "conditions": ["Breast Cancer"],
        "condition_meshes": [{"id": "D001943", "term": "Breast Neoplasms"}],
        "interventions": ["Drug X"],
        "summary": "Synthetic trial for testing.",
    },
    {
        "nct_id": "NCT00000002",
        "title": "Synthetic Trial B",
        "status": "COMPLETED",
        "phase": "PHASE3",
        "study_type": "INTERVENTIONAL",
        "start_date": "2020-01",
        "completion_date": "2023-12",
        "conditions": ["Breast Cancer", "BRCA1 Mutation"],
        "condition_meshes": [],   # no MeSH — falls back to text
        "interventions": [],
        "summary": "",
    },
    {
        "nct_id": "NCT00000003",
        "title": "Synthetic Trial C",
        "status": "TERMINATED",
        "phase": "",
        "study_type": "OBSERVATIONAL",
        "start_date": "",
        "completion_date": "",
        "conditions": [],
        "condition_meshes": [],
        "interventions": [],
        "summary": "",
    },
]

MOCK_QUERY = {"query": "BRCA1 breast cancer", "terms": ["BRCA1 breast cancer"]}


# ---------------------------------------------------------------------------
# parse_input
# ---------------------------------------------------------------------------

def test_parse_input_basic(tmp_path):
    f = tmp_path / "q.txt"
    f.write_text("# comment\nBRCA1\nbreast cancer\n")
    result = ctf.parse_input(f)
    assert result["query"] == "BRCA1 breast cancer"
    assert result["terms"] == ["BRCA1", "breast cancer"]


def test_parse_input_empty_raises(tmp_path):
    f = tmp_path / "empty.txt"
    f.write_text("# only comments\n\n")
    with pytest.raises(ValueError, match="No search terms"):
        ctf.parse_input(f)


def test_parse_input_single_term(tmp_path):
    f = tmp_path / "q.txt"
    f.write_text("EGFR\n")
    result = ctf.parse_input(f)
    assert result["query"] == "EGFR"


# ---------------------------------------------------------------------------
# _normalise_trial
# ---------------------------------------------------------------------------

def test_normalise_trial_truncates_summary():
    long_summary = "x" * 400
    raw = {
        "protocolSection": {
            "identificationModule": {"nctId": "NCT1", "briefTitle": "T"},
            "statusModule": {"overallStatus": "RECRUITING"},
            "designModule": {"phases": ["PHASE1"], "studyType": "INTERVENTIONAL"},
            "descriptionModule": {"briefSummary": long_summary},
            "conditionsModule": {},
            "armsInterventionsModule": {},
        }
    }
    trial = ctf._normalise_trial(raw)
    assert len(trial["summary"]) <= 304  # 300 chars + ellipsis
    assert trial["summary"].endswith("…")


def test_normalise_trial_missing_modules():
    """Graceful handling of empty study object."""
    trial = ctf._normalise_trial({"protocolSection": {}})
    assert trial["nct_id"] == ""
    assert trial["status"] == "UNKNOWN"
    assert trial["conditions"] == []


# ---------------------------------------------------------------------------
# fetch_trials
# ---------------------------------------------------------------------------

def _make_ct_response(trials_raw: list[dict]) -> bytes:
    return json.dumps({"studies": trials_raw}).encode()


def test_fetch_trials_network_error():
    with patch("urllib.request.urlopen", side_effect=urllib.error.URLError("timeout")):
        with pytest.raises(RuntimeError, match="ClinicalTrials.gov unreachable"):
            ctf.fetch_trials("BRCA1")


def test_fetch_trials_malformed_json():
    mock_resp = MagicMock()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    mock_resp.read.return_value = b"not-json{"
    with patch("urllib.request.urlopen", return_value=mock_resp):
        with pytest.raises(RuntimeError, match="malformed JSON"):
            ctf.fetch_trials("BRCA1")


def test_fetch_trials_empty_results():
    mock_resp = MagicMock()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    mock_resp.read.return_value = json.dumps({"studies": []}).encode()
    with patch("urllib.request.urlopen", return_value=mock_resp):
        result = ctf.fetch_trials("nothing-matches-this-query-xyz")
    assert result == []


# ---------------------------------------------------------------------------
# _count_recruiting
# ---------------------------------------------------------------------------

def test_count_recruiting():
    assert ctf._count_recruiting(MOCK_TRIALS) == 1


def test_count_recruiting_none():
    assert ctf._count_recruiting([MOCK_TRIALS[1], MOCK_TRIALS[2]]) == 0


# ---------------------------------------------------------------------------
# write_report
# ---------------------------------------------------------------------------

def test_write_report_creates_file(tmp_path):
    path = ctf.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    assert path.exists()


def test_write_report_contains_nct_ids(tmp_path):
    ctf.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.md").read_text()
    for t in MOCK_TRIALS:
        assert t["nct_id"] in content


def test_write_report_contains_disclaimer(tmp_path):
    ctf.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "research and educational tool" in content


def test_write_report_gene_context(tmp_path):
    gene_context = {
        "symbol": "BRCA1", "name": "BRCA1 DNA repair associated",
        "ensembl": "ENSG00000012048",
        "diseases": ["breast cancer", "ovarian cancer"],
        "min_score": 0.6,
    }
    ctf.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path, gene_context=gene_context)
    content = (tmp_path / "report.md").read_text()
    assert "BRCA1" in content
    assert "OpenTargets" in content


def test_write_report_terminated_trial_shown(tmp_path):
    """TERMINATED trials must not be hidden."""
    ctf.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "NCT00000003" in content


# ---------------------------------------------------------------------------
# write_summary
# ---------------------------------------------------------------------------

def test_write_summary(tmp_path):
    path = ctf.write_summary(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    data = json.loads(path.read_text())
    assert data["total"] == 3
    assert data["recruiting"] == 1
    assert len(data["trials"]) == 3


# ---------------------------------------------------------------------------
# write_fhir_bundle
# ---------------------------------------------------------------------------

def test_write_fhir_bundle_structure(tmp_path):
    path = ctf.write_fhir_bundle(MOCK_TRIALS, tmp_path)
    bundle = json.loads(path.read_text())
    assert bundle["resourceType"] == "Bundle"
    assert bundle["type"] == "searchset"
    assert bundle["total"] == 3


def test_fhir_status_mapping(tmp_path):
    ctf.write_fhir_bundle(MOCK_TRIALS, tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    by_id = {e["resource"]["id"]: e["resource"] for e in bundle["entry"]}
    assert by_id["NCT00000001"]["status"] == "active"
    assert by_id["NCT00000002"]["status"] == "completed"
    assert by_id["NCT00000003"]["status"] == "administratively-completed"


def test_fhir_phase_mapping(tmp_path):
    ctf.write_fhir_bundle(MOCK_TRIALS, tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    by_id = {e["resource"]["id"]: e["resource"] for e in bundle["entry"]}
    assert by_id["NCT00000001"]["phase"]["coding"][0]["code"] == "phase-2"


def test_fhir_all_statuses_covered():
    known = ["RECRUITING", "ACTIVE_NOT_RECRUITING", "NOT_YET_RECRUITING",
             "COMPLETED", "TERMINATED", "WITHDRAWN", "SUSPENDED", "UNKNOWN"]
    for s in known:
        assert s in ctf.FHIR_STATUS, f"Missing FHIR mapping for: {s}"


def test_fhir_trial_without_conditions(tmp_path):
    """Trial with no conditions should not have 'condition' key."""
    ctf.write_fhir_bundle([MOCK_TRIALS[2]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "condition" not in resource


def test_fhir_identifier_use_official(tmp_path):
    ctf.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    identifier = bundle["entry"][0]["resource"]["identifier"][0]
    assert identifier["use"] == "official"
    assert identifier["value"] == "NCT00000001"


def test_fhir_phase_display(tmp_path):
    ctf.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    coding = bundle["entry"][0]["resource"]["phase"]["coding"][0]
    assert coding["code"] == "phase-2"
    assert coding["display"] == "Phase 2"


def test_fhir_category_from_study_type(tmp_path):
    ctf.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "category" in resource
    assert resource["category"][0]["coding"][0]["code"] == "INTERVENTIONAL"


def test_fhir_focus_from_interventions(tmp_path):
    ctf.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "focus" in resource
    assert resource["focus"][0]["text"] == "Drug X"


def test_fhir_no_focus_when_no_interventions(tmp_path):
    ctf.write_fhir_bundle([MOCK_TRIALS[1]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "focus" not in resource


def test_fhir_condition_uses_mesh_when_available(tmp_path):
    """MeSH-coded conditions preferred over free-text when derivedSection provides them."""
    ctf.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    cond = bundle["entry"][0]["resource"]["condition"][0]
    assert "coding" in cond
    assert cond["coding"][0]["system"] == "http://id.nlm.nih.gov/mesh/"
    assert cond["coding"][0]["code"] == "D001943"
    assert cond["text"] == "Breast Neoplasms"


def test_fhir_condition_falls_back_to_text_without_mesh(tmp_path):
    """Trial without MeSH data falls back to text-only CodeableConcept."""
    ctf.write_fhir_bundle([MOCK_TRIALS[1]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    cond = bundle["entry"][0]["resource"]["condition"][0]
    assert "coding" not in cond
    assert cond["text"] == "Breast Cancer"


def test_normalise_trial_extracts_condition_meshes():
    raw = {
        "protocolSection": {
            "identificationModule": {"nctId": "NCT1", "briefTitle": "T"},
            "statusModule": {"overallStatus": "RECRUITING"},
            "designModule": {},
            "descriptionModule": {},
            "conditionsModule": {"conditions": ["Breast Cancer"]},
            "armsInterventionsModule": {},
        },
        "derivedSection": {
            "conditionBrowseModule": {
                "meshes": [{"id": "D001943", "term": "Breast Neoplasms"}]
            }
        },
    }
    trial = ctf._normalise_trial(raw)
    assert trial["condition_meshes"] == [{"id": "D001943", "term": "Breast Neoplasms"}]


def test_normalise_trial_empty_meshes_when_no_derived():
    raw = {
        "protocolSection": {
            "identificationModule": {},
            "statusModule": {},
            "designModule": {},
            "descriptionModule": {},
            "conditionsModule": {},
            "armsInterventionsModule": {},
        }
    }
    trial = ctf._normalise_trial(raw)
    assert trial["condition_meshes"] == []


# ---------------------------------------------------------------------------
# _to_fhir_date
# ---------------------------------------------------------------------------

def test_to_fhir_date_passthrough_iso():
    assert ctf._to_fhir_date("2024-01") == "2024-01"
    assert ctf._to_fhir_date("2024-01-15") == "2024-01-15"
    assert ctf._to_fhir_date("") == ""


def test_to_fhir_date_normalizes_month_year():
    assert ctf._to_fhir_date("January 2024") == "2024-01"
    assert ctf._to_fhir_date("December 2026") == "2026-12"
    assert ctf._to_fhir_date("March 2020") == "2020-03"


def test_to_fhir_date_unknown_format_passthrough():
    assert ctf._to_fhir_date("Q1 2024") == "Q1 2024"


# ---------------------------------------------------------------------------
# opentargets module
# ---------------------------------------------------------------------------

def _make_ot_response(data: dict) -> MagicMock:
    mock_resp = MagicMock()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    mock_resp.read.return_value = json.dumps({"data": data}).encode()
    return mock_resp


def test_ot_resolve_gene():
    resp = _make_ot_response({
        "search": {"hits": [{"id": "ENSG00000012048",
                             "object": {"approvedSymbol": "BRCA1",
                                        "approvedName": "BRCA1 DNA repair associated"}}]}
    })
    with patch("urllib.request.urlopen", return_value=resp):
        ensembl_id, name = ot.resolve_gene("BRCA1")
    assert ensembl_id == "ENSG00000012048"
    assert "BRCA1" in name


def test_ot_resolve_gene_not_found():
    resp = _make_ot_response({"search": {"hits": []}})
    with patch("urllib.request.urlopen", return_value=resp):
        with pytest.raises(ValueError, match="not found"):
            ot.resolve_gene("FAKEGENE999")


def test_ot_get_diseases_filters_by_score():
    resp = _make_ot_response({
        "target": {"associatedDiseases": {"rows": [
            {"disease": {"id": "MONDO_1", "name": "Breast Cancer"}, "score": 0.85},
            {"disease": {"id": "MONDO_2", "name": "Low Score Disease"}, "score": 0.3},
            {"disease": {"id": "MONDO_3", "name": "Ovarian Cancer"}, "score": 0.80},
        ]}}
    })
    with patch("urllib.request.urlopen", return_value=resp):
        diseases = ot.get_diseases("ENSG00000012048", min_score=0.6)
    names = [d.name for d in diseases]
    assert "Breast Cancer" in names
    assert "Ovarian Cancer" in names
    assert "Low Score Disease" not in names


def test_ot_network_error():
    with patch("urllib.request.urlopen", side_effect=urllib.error.URLError("refused")):
        with pytest.raises(RuntimeError, match="unreachable"):
            ot.resolve_gene("BRCA1")


def test_ot_graphql_error():
    mock_resp = MagicMock()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    mock_resp.read.return_value = json.dumps(
        {"errors": [{"message": "Field not found"}]}
    ).encode()
    with patch("urllib.request.urlopen", return_value=mock_resp):
        with pytest.raises(RuntimeError, match="GraphQL error"):
            ot.resolve_gene("BRCA1")


# ---------------------------------------------------------------------------
# --status filter
# ---------------------------------------------------------------------------

def test_status_filter_recruiting():
    filtered = [t for t in MOCK_TRIALS if t["status"] == "RECRUITING"]
    assert len(filtered) == 1
    assert filtered[0]["nct_id"] == "NCT00000001"


def test_status_filter_none_returns_all():
    filtered = [t for t in MOCK_TRIALS if True]  # no filter
    assert len(filtered) == len(MOCK_TRIALS)


# ---------------------------------------------------------------------------
# write_checksums
# ---------------------------------------------------------------------------

def test_write_checksums(tmp_path):
    (tmp_path / "report.md").write_text("# Test")
    (tmp_path / "summary.json").write_text("{}")
    path = ctf.write_checksums(tmp_path)
    assert path.exists()
    content = path.read_text()
    lines = [l for l in content.strip().split("\n") if l]
    assert len(lines) == 2
    for line in lines:
        digest, name = line.split("  ", 1)
        assert len(digest) == 64
        assert name in ("report.md", "summary.json")


def test_write_checksums_skips_missing(tmp_path):
    (tmp_path / "report.md").write_text("# Test")
    path = ctf.write_checksums(tmp_path)
    content = path.read_text().strip()
    assert "summary.json" not in content
    assert "report.md" in content


def test_write_checksums_includes_figures(tmp_path):
    (tmp_path / "report.md").write_text("x")
    fig_dir = tmp_path / "figures"
    fig_dir.mkdir()
    (fig_dir / "phase_distribution.png").write_bytes(b"\x89PNG")
    path = ctf.write_checksums(tmp_path)
    assert "figures/phase_distribution.png" in path.read_text()


# ---------------------------------------------------------------------------
# write_commands
# ---------------------------------------------------------------------------

def test_write_commands_demo(tmp_path):
    args = argparse.Namespace(
        demo=True, input=None, query=None, gene=None,
        status=None, max_results=20, fhir=False,
        ot_min_score=0.6, ot_max_diseases=5,
    )
    path = ctf.write_commands(args, tmp_path)
    assert path.exists()
    content = path.read_text()
    assert "--demo" in content
    assert f"--output {tmp_path}" in content


def test_write_commands_query_with_status(tmp_path):
    args = argparse.Namespace(
        demo=False, input=None, query="lung cancer", gene=None,
        status="RECRUITING", max_results=20, fhir=True,
        ot_min_score=0.6, ot_max_diseases=5,
    )
    path = ctf.write_commands(args, tmp_path)
    content = path.read_text()
    assert '--query "lung cancer"' in content
    assert "--status RECRUITING" in content
    assert "--fhir" in content


def test_write_commands_non_default_ot_params(tmp_path):
    args = argparse.Namespace(
        demo=False, input=None, query=None, gene="BRCA1",
        status=None, max_results=50, fhir=False,
        ot_min_score=0.3, ot_max_diseases=10,
    )
    path = ctf.write_commands(args, tmp_path)
    content = path.read_text()
    assert "--gene BRCA1" in content
    assert "--max-results 50" in content
    assert "--ot-min-score 0.3" in content
    assert "--ot-max-diseases 10" in content


# ---------------------------------------------------------------------------
# write_phase_chart — figures/ subdir
# ---------------------------------------------------------------------------

def test_chart_in_figures_subdir(tmp_path):
    chart = ctf.write_phase_chart(MOCK_TRIALS, tmp_path)
    if chart is None:
        pytest.skip("matplotlib not installed")
    assert chart.parent.name == "figures"
    assert chart.name == "phase_distribution.png"
    assert chart.exists()


# ---------------------------------------------------------------------------
# Integration sanity
# ---------------------------------------------------------------------------

def test_demo_data_exists_and_parseable():
    demo = SKILL_DIR / "demo_input.txt"
    assert demo.exists()
    result = ctf.parse_input(demo)
    assert result["query"]
