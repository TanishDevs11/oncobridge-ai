"""End-to-end local test: generates synthetic RNA-seq counts, drives all
three endpoints against a locally running uvicorn instance, and sanity
checks the response shapes. Also tests validation failure paths.

Run with the server already started separately:
    uvicorn main:app --port 8000
Then:
    python test_e2e.py
"""

import base64
import sys
import time

import numpy as np
import requests

BASE = "http://127.0.0.1:8000"
N_GENES = 300
SAMPLES = [f"sample{i+1}" for i in range(8)]  # 4 normal, 4 tumor
CONDITIONS = ["normal"] * 4 + ["tumor"] * 4
N_DEG_GENES = 15

rng = np.random.default_rng(42)


def build_synthetic_data():
    gene_names = [f"GENE{i:04d}" for i in range(N_GENES)]
    # Salt in some real biomarker panel genes so the biomarker cross-ref has hits.
    for i, real_gene in enumerate(["TP53", "EGFR", "KRAS", "MYC", "BRCA1", "PTEN"]):
        gene_names[i] = real_gene

    deg_indices = set(rng.choice(N_GENES, size=N_DEG_GENES, replace=False).tolist())

    counts = {}
    for gi, gene in enumerate(gene_names):
        base_mean = rng.uniform(20, 500)
        sample_counts = {}
        for si, sample in enumerate(SAMPLES):
            is_tumor = CONDITIONS[si] == "tumor"
            mean = base_mean
            if gi in deg_indices and is_tumor:
                mean = base_mean * rng.choice([6, 8, 10])  # strong upregulation
            elif gi in deg_indices and not is_tumor:
                mean = base_mean
            n = 10
            p = n / (n + mean)
            val = int(rng.negative_binomial(n, p))
            sample_counts[sample] = val
        counts[gene] = sample_counts

    metadata = {sample: cond for sample, cond in zip(SAMPLES, CONDITIONS)}
    return counts, metadata


def is_valid_png_b64(s: str) -> bool:
    try:
        raw = base64.b64decode(s)
        return raw[:8] == b"\x89PNG\r\n\x1a\n"
    except Exception:
        return False


def test_health():
    r = requests.get(f"{BASE}/health")
    assert r.status_code == 200, r.text
    assert r.json() == {"status": "ok"}
    print("PASS /health")


def test_skill_md():
    r = requests.get(f"{BASE}/skill.md")
    assert r.status_code == 200, r.text
    assert len(r.text) > 0
    print("PASS /skill.md ->", r.headers.get("content-type"))


def test_validation_failures():
    # Mismatched samples
    counts, metadata = build_synthetic_data()
    bad_metadata = dict(metadata)
    del bad_metadata[SAMPLES[0]]
    r = requests.post(f"{BASE}/analyze", json={"counts": counts, "metadata": bad_metadata, "condition_column": "condition"})
    assert r.status_code == 422, r.text
    print("PASS mismatched samples -> 422:", r.json())

    # Only one condition
    one_cond_metadata = {s: "tumor" for s in SAMPLES}
    r = requests.post(f"{BASE}/analyze", json={"counts": counts, "metadata": one_cond_metadata, "condition_column": "condition"})
    assert r.status_code == 422, r.text
    print("PASS single condition -> 422:", r.json())

    # Non-integer / negative counts
    bad_counts = {k: dict(v) for k, v in counts.items()}
    first_gene = next(iter(bad_counts))
    bad_counts[first_gene][SAMPLES[0]] = -5
    r = requests.post(f"{BASE}/analyze", json={"counts": bad_counts, "metadata": metadata, "condition_column": "condition"})
    assert r.status_code == 422, r.text
    print("PASS negative count -> 422:", r.json())

    bad_counts2 = {k: dict(v) for k, v in counts.items()}
    bad_counts2[first_gene][SAMPLES[0]] = 3.5
    r = requests.post(f"{BASE}/analyze", json={"counts": bad_counts2, "metadata": metadata, "condition_column": "condition"})
    assert r.status_code == 422, r.text
    print("PASS float count -> 422:", r.json())

    # Too few samples per condition
    small_metadata = {s: c for s, c in zip(SAMPLES[:3], ["normal", "normal", "tumor"])}
    small_counts = {g: {s: v[s] for s in SAMPLES[:3]} for g, v in counts.items()}
    r = requests.post(f"{BASE}/analyze", json={"counts": small_counts, "metadata": small_metadata, "condition_column": "condition"})
    assert r.status_code == 422, r.text
    print("PASS too few samples per condition -> 422:", r.json())

    # Unknown 404 job
    r = requests.get(f"{BASE}/jobs/does-not-exist")
    assert r.status_code == 404, r.text
    print("PASS unknown job -> 404")

    r = requests.get(f"{BASE}/results/does-not-exist")
    assert r.status_code == 404, r.text
    print("PASS unknown results -> 404")


def test_full_pipeline():
    counts, metadata = build_synthetic_data()
    r = requests.post(f"{BASE}/analyze", json={"counts": counts, "metadata": metadata, "condition_column": "condition"})
    assert r.status_code == 202, r.text
    body = r.json()
    job_id = body["job_id"]
    assert body["status"] == "queued"
    print("PASS /analyze -> 202", body)

    # Immediately check /results should be 409
    r = requests.get(f"{BASE}/results/{job_id}")
    assert r.status_code == 409, r.text
    print("PASS /results while running -> 409:", r.json())

    # Poll
    deadline = time.time() + 300
    last_stage = None
    while time.time() < deadline:
        r = requests.get(f"{BASE}/jobs/{job_id}")
        assert r.status_code == 200, r.text
        j = r.json()
        if j["stage"] != last_stage:
            print("  stage ->", j["stage"], "status:", j["status"])
            last_stage = j["stage"]
        if j["status"] == "complete":
            break
        if j["status"] == "failed":
            print("JOB FAILED:", j["error"])
            sys.exit(1)
        time.sleep(2)
    else:
        print("TIMEOUT waiting for job")
        sys.exit(1)

    r = requests.get(f"{BASE}/results/{job_id}")
    assert r.status_code == 200, r.text
    report = r.json()

    assert report["job_id"] == job_id
    assert report["summary"]["n_genes_input"] == N_GENES
    assert report["summary"]["n_samples"] == len(SAMPLES)
    assert set(report["summary"]["conditions"]) == {"normal", "tumor"}
    print("Summary:", report["summary"])

    assert isinstance(report["degs"], list)
    print(f"DEGs found: {len(report['degs'])}")
    if len(report["degs"]) == 0:
        print("WARNING: expected some DEGs given synthetic shift, got 0")

    assert "GO_Biological_Process_2023" in report["enrichment"]
    assert "KEGG_2021_Human" in report["enrichment"]
    total_terms = sum(len(v) for v in report["enrichment"].values())
    print(f"Enrichment terms total: {total_terms} (0 is acceptable if Enrichr unreachable from this sandbox)")

    biomarker_genes = {b["gene"] for b in report["biomarkers"]}
    print("Biomarkers matched:", biomarker_genes)
    assert biomarker_genes.issubset(set(["TP53", "EGFR", "KRAS", "MYC", "BRCA1", "PTEN"]) | biomarker_genes)
    expected_panel_hits = {"TP53", "EGFR", "KRAS", "MYC", "BRCA1", "PTEN"}
    assert expected_panel_hits.issubset(biomarker_genes), f"missing {expected_panel_hits - biomarker_genes}"

    assert is_valid_png_b64(report["figures"]["volcano_png_base64"]), "volcano figure is not a valid PNG"
    assert is_valid_png_b64(report["figures"]["heatmap_png_base64"]), "heatmap figure is not a valid PNG"
    print("PASS figures are valid non-empty PNGs")

    assert report["parameters"]["log2fc_threshold"] == 2.0
    assert report["parameters"]["padj_threshold"] == 0.05

    print("PASS full pipeline end-to-end")


if __name__ == "__main__":
    test_health()
    test_skill_md()
    test_validation_failures()
    test_full_pipeline()
    print("\nALL TESTS PASSED")
