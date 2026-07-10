# OncoBridge AI

A single FastAPI service that runs raw RNA-seq gene counts through a complete cancer transcriptomics workflow — differential expression, functional enrichment, cancer biomarker annotation, and visualization — and returns one structured JSON report. Built for NANDAHack 2026, Phase 2.

**Live service:** https://oncobridge-ai-production.up.railway.app
**Agent-facing docs:** [`SKILL.md`](./SKILL.md), also served live at [`/skill.md`](https://oncobridge-ai-production.up.railway.app/skill.md)

## What it does

Submit a gene x sample counts matrix and a condition label per sample. The service:

1. Validates the request synchronously (bad input fails fast with a specific `422`, not a background job that dies later).
2. Runs differential expression with [`pydeseq2`](https://pydeseq2.readthedocs.io/) (a pure-Python DESeq2 reimplementation — no R, no `rpy2`).
3. Filters to significant DEGs (`|log2FoldChange| >= 2.0` and `padj < 0.05` by default).
4. Runs functional enrichment against `GO_Biological_Process_2023` and `KEGG_2021_Human` via the live [Enrichr](https://maayanlab.cloud/Enrichr/) API using `gseapy`.
5. Cross-references every input gene against a curated ~50-gene cancer biomarker panel (oncogenes, tumor suppressors, proto-oncogenes) and flags which ones are DEGs.
6. Renders a volcano plot and a clustered heatmap of the top DEGs as base64-encoded PNGs.
7. Assembles everything into one JSON report.

**Constraint:** two-group comparisons only. `metadata` must contain exactly two distinct condition labels, with at least 2 samples per condition.

## API

| Method & path | Description |
|---|---|
| `POST /analyze` | Submit `counts` + `metadata`. Returns `202` + `job_id`, or `422` with a specific validation error. |
| `GET /jobs/{job_id}` | Poll job status and current pipeline stage. `404` if unknown. |
| `GET /results/{job_id}` | Full report once complete. `409` while running, `500` on failure, `404` if unknown. |
| `GET /health` | Liveness check. |
| `GET /skill.md` | This service's machine-readable API description, as `text/markdown`. |

Full request/response shapes, real `curl` examples, and polling guidance for an agent are in [`SKILL.md`](./SKILL.md) — every example there was captured from the live deployment, not hand-written.

## Pipeline stages

Tracked per job via `GET /jobs/{job_id}`, in order:

```
validating -> differential_expression -> filtering_degs -> enrichment -> biomarker_annotation -> visualization -> report_assembly
```

## Project layout

```
main.py         FastAPI app: routes, validation, and the background pipeline
biomarkers.py   Curated cancer gene panel (role + one-line relevance per gene)
test_e2e.py     Synthetic-data end-to-end test: drives all endpoints against a local server
requirements.txt   Pinned to the exact versions this was built and deployed against
Procfile        Railway/Heroku-style start command
.python-version Pins the runtime to Python 3.10
SKILL.md        Agent-facing API description, served live at /skill.md
```

## Running locally

Requires Python 3.10 (pydeseq2's dependency wheels are most reliable there).

```bash
python -m venv .venv
.venv/Scripts/activate   # or `source .venv/bin/activate` on macOS/Linux
pip install -r requirements.txt
uvicorn main:app --port 8000
```

Then, in another terminal:

```bash
python test_e2e.py
```

This generates a synthetic negative-binomial counts matrix, drives `/analyze` -> `/jobs/{id}` -> `/results/{id}` end to end, and exercises every `422` validation path.

Enrichment calls the live Enrichr API — it needs real internet access. In a network-restricted environment, `enrichment` degrades gracefully to empty lists instead of failing the job.

## Deployment notes

Deployed on Railway. Two things worth knowing if you redeploy this elsewhere:

- **Memory:** `pydeseq2`'s default `n_cpus` resolves to `os.cpu_count()`, and BLAS thread pools size themselves the same way — off the host's visible core count, not the container's memory limit. On a memory-constrained container with a host reporting many cores, this fans out enough worker processes to get OOM-killed even on a tiny gene panel. `main.py` pins `n_cpus=1`, `low_memory=True`, and single-threads BLAS via environment variables set before `numpy` is ever imported.
- **Cold starts:** free-tier hosting sleeps when idle. Send a `GET /health` warm-up request before testing other endpoints; the first request after a quiet period can take 30-60 seconds.

## License

Built for NANDAHack 2026, Phase 2.
