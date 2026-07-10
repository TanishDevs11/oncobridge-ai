# OncoBridge AI

OncoBridge AI runs a complete cancer transcriptomics workflow — differential expression, functional enrichment, cancer biomarker annotation, and visualization — on raw RNA-seq gene counts and returns a single structured JSON report.

Base URL: https://oncobridge-ai-production.up.railway.app

## Important constraint

This service only supports **two-group comparisons**. The `metadata` you send must contain **exactly two distinct condition labels** (e.g. `"normal"` vs `"tumor"`), with at least 2 samples per condition. Requests with one condition, three or more conditions, or fewer than 2 samples in either group are rejected with `422` before any job is created.

## Endpoints

### GET /health

Liveness check. Use this once as a warm-up call before anything else — free-tier hosting sleeps when idle and the first request after a quiet period can take 30-60 seconds to respond.

```
curl -s https://oncobridge-ai-production.up.railway.app/health
```

Real response:
```json
{"status":"ok"}
```

### GET /skill.md

Returns this file as plain text (`Content-Type: text/markdown`). Useful for an agent to re-fetch the latest version of these instructions at runtime.

```
curl -s https://oncobridge-ai-production.up.railway.app/skill.md
```

### POST /analyze

Submits a gene counts matrix and sample metadata for analysis. Returns immediately with a job ID; the actual pipeline (DESeq2 differential expression, Enrichr functional enrichment, biomarker annotation, figure generation) runs in the background. On a small panel (tens of genes, as in the example below) it typically finishes within about 5-20 seconds; larger gene panels take longer because DESeq2's model fit and the live Enrichr API call both scale with input size.

Request body fields:
- `counts` (object, required): `gene_symbol -> {sample_id: non-negative integer read count}`. Every gene must report counts for the exact same set of samples.
- `metadata` (object, required): `sample_id -> condition label`. Exactly 2 distinct labels required, at least 2 samples per label.
- `condition_column` (string, optional): informational label for the comparison column; the service always compares the two condition labels found in `metadata`.

Real example call (this exact payload was run against the live service):

```
curl -s -X POST https://oncobridge-ai-production.up.railway.app/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "counts": {
      "TP53":   {"sample1": 271, "sample2": 200, "sample3": 202, "sample4": 3039, "sample5": 1672, "sample6": 3189},
      "EGFR":   {"sample1": 83,  "sample2": 140, "sample3": 120, "sample4": 659,  "sample5": 1709, "sample6": 1055},
      "KRAS":   {"sample1": 205, "sample2": 161, "sample3": 360, "sample4": 2200, "sample5": 2838, "sample6": 1916},
      "MYC":    {"sample1": 289, "sample2": 178, "sample3": 125, "sample4": 1561, "sample5": 1784, "sample6": 2697},
      "BRCA1":  {"sample1": 262, "sample2": 171, "sample3": 83,  "sample4": 3026, "sample5": 2928, "sample6": 3261},
      "ACTB":   {"sample1": 249, "sample2": 324, "sample3": 176, "sample4": 295,  "sample5": 254,  "sample6": 312},
      "GAPDH":  {"sample1": 371, "sample2": 379, "sample3": 395, "sample4": 347,  "sample5": 307,  "sample6": 515}
    },
    "metadata": {
      "sample1": "normal", "sample2": "normal", "sample3": "normal",
      "sample4": "tumor",  "sample5": "tumor",  "sample6": "tumor"
    },
    "condition_column": "condition"
  }'
```

Real response (`202 Accepted`):
```json
{"job_id":"7a350bfa-bb86-4217-9449-e7342d43d557","status":"queued"}
```

Real validation-failure example — metadata with only one condition label:
```
curl -s -X POST https://oncobridge-ai-production.up.railway.app/analyze \
  -H "Content-Type: application/json" \
  -d '{"counts": {"TP53": {"s1": 100}}, "metadata": {"s1": "tumor"}, "condition_column": "condition"}'
```
Real response (`422 Unprocessable Entity`):
```json
{"detail":"Exactly 2 distinct condition labels are required (two-group comparison only), found 1: ['tumor']."}
```

### GET /jobs/{job_id}

Reports the current status and pipeline stage of a submitted job. Poll this every 5-10 seconds after submitting.

Stage progresses in this fixed order: `validating` -> `differential_expression` -> `filtering_degs` -> `enrichment` -> `biomarker_annotation` -> `visualization` -> `report_assembly`.

```
curl -s https://oncobridge-ai-production.up.railway.app/jobs/7a350bfa-bb86-4217-9449-e7342d43d557
```

Real response once finished:
```json
{"job_id":"7a350bfa-bb86-4217-9449-e7342d43d557","status":"complete","stage":"report_assembly","error":null}
```

While the job is still processing, `status` is `"running"` and `stage` names the current pipeline step (e.g. `"differential_expression"` or `"enrichment"`) instead of `"report_assembly"`. On small gene panels like the example above the whole pipeline can finish in just a few seconds, so you may only ever observe `queued` then `complete`.

`status` is one of `queued`, `running`, `complete`, `failed`. Unknown job IDs return `404`.

### GET /results/{job_id}

Returns the full analysis report. Only returns `200` once `status` is `complete`. Returns `409` with `{"status": "<current status>"}` while still queued/running, `404` for an unknown job ID, and `500` with the stored error message if the job failed.

```
curl -s https://oncobridge-ai-production.up.railway.app/results/7a350bfa-bb86-4217-9449-e7342d43d557
```

Real response (`200 OK`, trimmed here for length — `degs` had 10 entries, `enrichment` had 552 GO terms and 136 KEGG terms, `biomarkers` had 20 entries, and both figure fields are full base64 PNG strings tens of thousands of characters long):
```json
{
  "job_id": "7a350bfa-bb86-4217-9449-e7342d43d557",
  "summary": {"n_genes_input": 30, "n_samples": 6, "n_degs": 10, "conditions": ["normal", "tumor"]},
  "degs": [
    {"gene": "BRCA1", "log2FoldChange": 3.52375145570805, "pvalue": 1.8633719805316935e-18, "padj": 5.59011594159508e-17},
    {"gene": "ALK", "log2FoldChange": 2.9939917843151163, "pvalue": 3.412882384383047e-17, "padj": 5.11932357657457e-16},
    {"gene": "TP53", "log2FoldChange": 2.861365826951824, "pvalue": 1.4854683819820658e-15, "padj": 1.4854683819820658e-14}
  ],
  "enrichment": {
    "GO_Biological_Process_2023": [
      {"term": "Regulation Of Cell Population Proliferation (GO:0042127)", "p_value": 1.2765208625023462e-08, "adjusted_p_value": 7.04639516101295e-06, "genes": ["ALK", "MYC", "ERBB2", "PTEN", "KRAS", "TP53", "EGFR"]}
    ],
    "KEGG_2021_Human": [
      {"term": "Breast cancer", "p_value": 4.850917936031202e-19, "adjusted_p_value": 6.597248393002433e-17, "genes": ["RB1", "PIK3CA", "MYC", "ERBB2", "PTEN", "KRAS", "BRCA1", "TP53", "EGFR"]}
    ]
  },
  "biomarkers": [
    {"gene": "TP53", "role": "tumor_suppressor", "relevance": "Most frequently mutated gene in human cancer; guards genome integrity via cell-cycle arrest and apoptosis.", "in_degs": true},
    {"gene": "EGFR", "role": "oncogene", "relevance": "Receptor tyrosine kinase driving proliferation; activating mutations/amplification common in lung cancer and glioblastoma.", "in_degs": true},
    {"gene": "MET", "role": "proto-oncogene", "relevance": "Hepatocyte growth factor receptor; amplification/mutation drives proliferation and invasion, targetable in NSCLC.", "in_degs": false}
  ],
  "figures": {
    "volcano_png_base64": "iVBORw0KGgoAAAANSUhEUgAAAogAAAIYCAYAAADuLx35AAAAOnRFWHRTb2Z0...(36372 chars total, valid PNG)",
    "heatmap_png_base64": "iVBORw0KGgoAAAANSUhEUgAAArwAAAKhCAYAAACy8n0GAAAAOnRFWHRTb2Z0...(30604 chars total, valid PNG)"
  },
  "parameters": {"log2fc_threshold": 2.0, "padj_threshold": 0.05}
}
```

`degs` is filtered to `|log2FoldChange| >= 2.0` and `padj < 0.05` (see `parameters` in the response for the exact thresholds used). `biomarkers` lists every gene from a curated ~50-gene cancer panel that appears in your input `counts`, each flagged `in_degs: true/false` depending on whether it passed the DEG filter. `figures` contains full base64-encoded PNG images (volcano plot and a clustered heatmap of the top DEGs) — decode and save as `.png` to view them.

## How an agent should use this service

1. Build a `counts` object covering every gene you want analyzed, with a read count for every sample, and a `metadata` object giving each sample's condition. Confirm exactly 2 distinct condition labels and at least 2 samples per label.
2. `POST /analyze` with that JSON body. On success you get `202` and a `job_id`. On `422`, read `detail` — it names the exact problem (mismatched samples, wrong count type, wrong number of conditions, etc.) and the request should be corrected and resent, not retried as-is.
3. Poll `GET /jobs/{job_id}` every 5-10 seconds. Watch `status` move from `queued` to `running` (with `stage` advancing through the pipeline) to `complete` or `failed`. Real runs on a small gene panel (dozens of genes) complete in about 5-20 seconds; larger panels take longer because of the DESeq2 fit and the live Enrichr API call.
4. Once `status` is `complete`, `GET /results/{job_id}` for the full report. If `status` is `failed`, this returns `500` with the stored error message instead — do not keep polling a failed job.
5. Read `summary.n_degs` first to know whether any genes passed the significance filter before parsing the full `degs`/`enrichment`/`biomarkers` lists.
