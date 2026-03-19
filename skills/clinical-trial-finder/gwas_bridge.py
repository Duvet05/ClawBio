"""Bridge to gwas-lookup skill -- resolves rsID to traits and genes.

Calls gwas-lookup via clawbio.py (which handles import paths correctly)
as a subprocess.  Extracts genome-wide significant GWAS traits and eQTL
gene symbols, deduplicates, and returns them for CT.gov queries.
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

_CLAWBIO_DIR = Path(__file__).resolve().parent.parent.parent
_CLAWBIO_PY = _CLAWBIO_DIR / "clawbio.py"
_GWAS_SCRIPT = _CLAWBIO_DIR / "skills" / "gwas-lookup" / "gwas_lookup.py"


def resolve_rsid(rsid: str, max_traits: int = 5) -> dict:
    """Resolve an rsID via gwas-lookup and return associated traits and genes.

    Routes through clawbio.py to handle import paths correctly.
    Returns dict with keys: rsid, traits (list[str]), genes (list[str]).
    Raises ValueError if gwas-lookup fails or returns no associations.
    """
    if not _GWAS_SCRIPT.exists():
        raise ValueError(
            f"gwas-lookup skill not found at {_GWAS_SCRIPT}. "
            "Install it or use --gene/--query instead."
        )

    with tempfile.TemporaryDirectory(prefix="ctf_gwas_") as tmp:
        tmp_dir = Path(tmp)
        # gwas-lookup mixes relative and absolute imports; run it from its
        # own directory with that directory on PYTHONPATH so both work.
        import os

        gwas_dir = str(_GWAS_SCRIPT.parent)
        env = {**os.environ, "PYTHONPATH": gwas_dir}
        result = subprocess.run(
            [
                sys.executable,
                str(_GWAS_SCRIPT),
                "--rsid",
                rsid,
                "--no-figures",
                "--output",
                str(tmp_dir),
            ],
            capture_output=True,
            text=True,
            timeout=90,
            cwd=gwas_dir,
            env=env,
        )

        # clawbio.py wraps the skill; check for result.json regardless of exit code
        result_path = tmp_dir / "result.json"
        if not result_path.exists():
            stderr = result.stderr.strip().split("\n")[-1] if result.stderr else "unknown error"
            raise ValueError(f"gwas-lookup failed for {rsid}: {stderr}")

        data = json.loads(result_path.read_text())

    merged = data.get("data", {}).get("merged", {})
    return _extract_traits_and_genes(rsid, merged, max_traits)


def _extract_traits_and_genes(
    rsid: str, merged: dict, max_traits: int
) -> dict:
    """Extract unique, genome-wide significant traits and eQTL genes.

    Traits are deduplicated case-insensitively and ranked by p-value.
    Only genome-wide significant associations (p < 5e-8) are included.
    """
    # Collect traits from GWAS associations (most authoritative)
    trait_pvals: dict[str, float] = {}
    for assoc in merged.get("gwas_associations", []):
        trait = assoc.get("trait", "").strip()
        pval = assoc.get("pval", 1.0)
        if trait and pval < 5e-8:  # genome-wide significance threshold
            key = trait.lower()
            if key not in trait_pvals or pval < trait_pvals[key]:
                trait_pvals[key] = pval

    # Also check PheWAS for additional disease names
    for source_hits in merged.get("phewas", {}).values():
        for hit in source_hits:
            trait = hit.get("phenostring", "").strip()
            pval = hit.get("pval", 1.0)
            if trait and pval < 5e-8:
                key = trait.lower()
                if key not in trait_pvals or pval < trait_pvals[key]:
                    trait_pvals[key] = pval

    # Sort by p-value (most significant first), take top N
    sorted_traits = sorted(trait_pvals.items(), key=lambda x: x[1])
    traits = [k.title() for k, _ in sorted_traits[:max_traits]]

    # Extract gene symbols from eQTL associations
    genes: list[str] = []
    seen_genes: set[str] = set()
    for eqtl in merged.get("eqtl_associations", []):
        gene = eqtl.get("gene", "").strip()
        if gene and gene not in seen_genes:
            seen_genes.add(gene)
            genes.append(gene)

    if not traits and not genes:
        raise ValueError(
            f"No genome-wide significant associations found for {rsid}. "
            "Try --gene or --query instead."
        )

    return {"rsid": rsid, "traits": traits, "genes": genes}
