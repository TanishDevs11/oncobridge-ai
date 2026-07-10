"""OncoBridge AI: end-to-end cancer transcriptomics pipeline as a web service.

Pipeline: raw RNA-seq counts -> differential expression (pydeseq2) ->
DEG filtering -> functional enrichment (gseapy/Enrichr) -> biomarker
annotation -> visualization -> assembled JSON report.

Two-group comparisons only (exactly two condition labels in metadata).
"""

import os

# Must be set before numpy/scipy are imported anywhere in the process.
# Without this, BLAS libraries size their thread pools off the host's
# detected CPU count (not the container's cgroup memory limit), and
# pydeseq2's own multiprocessing defaults to one worker per core too --
# on a host reporting many cores, that blows well past a 1GB container
# limit even for a tiny gene panel. Force everything single-threaded.
for _env_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_env_var, "1")

import base64
import gc
import io
import json
import logging
import uuid
from typing import Dict, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from fastapi import BackgroundTasks, FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse, PlainTextResponse

from biomarkers import BIOMARKER_PANEL

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("oncobridge")

app = FastAPI(title="OncoBridge AI")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
LOG2FC_THRESHOLD = 2.0
PADJ_THRESHOLD = 0.05
MIN_SAMPLES_PER_CONDITION = 2
ENRICHR_GENE_SETS = ["GO_Biological_Process_2023", "KEGG_2021_Human"]
HEATMAP_MAX_GENES = 50
HEATMAP_MIN_GENES = 20

STAGES = [
    "validating",
    "differential_expression",
    "filtering_degs",
    "enrichment",
    "biomarker_annotation",
    "visualization",
    "report_assembly",
]

# ---------------------------------------------------------------------------
# In-memory job store: job_id -> {status, stage, error, result}
# ---------------------------------------------------------------------------
JOBS: Dict[str, dict] = {}


# ---------------------------------------------------------------------------
# Strict JSON parsing with duplicate-key detection
# ---------------------------------------------------------------------------
class DuplicateKeyError(ValueError):
    pass


def _dedupe_check_hook(pairs):
    seen = {}
    for k, v in pairs:
        if k in seen:
            raise DuplicateKeyError(f"Duplicate key '{k}' found in request JSON body.")
        seen[k] = v
    return seen


async def read_json_strict(request: Request) -> dict:
    raw = await request.body()
    try:
        return json.loads(raw, object_pairs_hook=_dedupe_check_hook)
    except DuplicateKeyError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except json.JSONDecodeError as e:
        raise HTTPException(status_code=422, detail=f"Request body is not valid JSON: {e}")


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def validate_analyze_payload(payload: dict):
    """Validate the /analyze request body. Raises HTTPException(422, ...) with
    a specific, actionable message on any failure. Returns (counts, metadata)
    cleaned and ready for the pipeline on success.
    """
    if not isinstance(payload, dict):
        raise HTTPException(422, detail="Request body must be a JSON object.")

    counts = payload.get("counts")
    metadata = payload.get("metadata")

    if counts is None or not isinstance(counts, dict):
        raise HTTPException(422, detail="Field 'counts' is required and must be an object mapping gene_id -> {sample_id: count}.")
    if metadata is None or not isinstance(metadata, dict):
        raise HTTPException(422, detail="Field 'metadata' is required and must be an object mapping sample_id -> condition label.")
    if len(counts) == 0:
        raise HTTPException(422, detail="Field 'counts' must contain at least one gene.")
    if len(metadata) == 0:
        raise HTTPException(422, detail="Field 'metadata' must contain at least one sample.")

    metadata_samples = set(metadata.keys())
    for s, label in metadata.items():
        if not isinstance(label, str) or not label.strip():
            raise HTTPException(422, detail=f"Metadata condition label for sample '{s}' must be a non-empty string.")

    # Validate each gene's sample set is non-negative integer counts and
    # matches metadata sample set exactly.
    expected_samples = None
    for gene, sample_counts in counts.items():
        if not isinstance(sample_counts, dict):
            raise HTTPException(422, detail=f"Gene '{gene}' must map to an object of {{sample_id: count}}.")
        gene_samples = set(sample_counts.keys())
        if expected_samples is None:
            expected_samples = gene_samples
        elif gene_samples != expected_samples:
            raise HTTPException(
                422,
                detail=(
                    f"Gene '{gene}' has sample set {sorted(gene_samples)}, which does not match the "
                    f"sample set of other genes {sorted(expected_samples)}. All genes must report counts "
                    f"for the exact same set of samples."
                ),
            )
        for sample_id, value in sample_counts.items():
            if isinstance(value, bool) or not isinstance(value, int):
                raise HTTPException(
                    422,
                    detail=f"Count for gene '{gene}', sample '{sample_id}' must be a non-negative integer, got {value!r}.",
                )
            if value < 0:
                raise HTTPException(
                    422,
                    detail=f"Count for gene '{gene}', sample '{sample_id}' must be non-negative, got {value}.",
                )

    if expected_samples != metadata_samples:
        missing_in_metadata = expected_samples - metadata_samples
        missing_in_counts = metadata_samples - expected_samples
        parts = []
        if missing_in_metadata:
            parts.append(f"samples in 'counts' with no entry in 'metadata': {sorted(missing_in_metadata)}")
        if missing_in_counts:
            parts.append(f"samples in 'metadata' with no entry in 'counts': {sorted(missing_in_counts)}")
        raise HTTPException(422, detail="Sample sets in 'counts' and 'metadata' must match exactly. " + "; ".join(parts) + ".")

    conditions = sorted(set(metadata.values()))
    if len(conditions) != 2:
        raise HTTPException(
            422,
            detail=(
                f"Exactly 2 distinct condition labels are required (two-group comparison only), "
                f"found {len(conditions)}: {conditions}."
            ),
        )

    for cond in conditions:
        n = sum(1 for v in metadata.values() if v == cond)
        if n < MIN_SAMPLES_PER_CONDITION:
            raise HTTPException(
                422,
                detail=f"Condition '{cond}' has only {n} sample(s); at least {MIN_SAMPLES_PER_CONDITION} are required per condition.",
            )

    return counts, metadata


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
def _fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=110)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def run_pipeline(job_id: str, counts: dict, metadata: dict):
    job = JOBS[job_id]
    try:
        job["status"] = "running"

        # ---- validating (structural setup) ----
        job["stage"] = "validating"
        genes = list(counts.keys())
        samples = sorted(metadata.keys())
        conditions = sorted(set(metadata.values()))
        cond_a, cond_b = conditions[0], conditions[1]

        counts_df = pd.DataFrame(
            {gene: [counts[gene][s] for s in samples] for gene in genes},
            index=samples,
        ).astype(int)
        metadata_df = pd.DataFrame({"condition": [metadata[s] for s in samples]}, index=samples)
        metadata_df["condition"] = metadata_df["condition"].astype("category")

        # ---- differential_expression ----
        job["stage"] = "differential_expression"
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata_df,
            design_factors="condition",
            refit_cooks=True,
            quiet=True,
            n_cpus=1,
            low_memory=True,
        )
        dds.deseq2()
        stat_res = DeseqStats(dds, contrast=["condition", cond_b, cond_a], quiet=True, n_cpus=1)
        stat_res.summary()
        res_df = stat_res.results_df.copy()
        res_df["gene"] = res_df.index
        res_df["padj"] = res_df["padj"].fillna(1.0)
        res_df["pvalue"] = res_df["pvalue"].fillna(1.0)
        del dds, stat_res
        gc.collect()

        # ---- filtering_degs ----
        job["stage"] = "filtering_degs"
        deg_mask = (res_df["padj"] < PADJ_THRESHOLD) & (res_df["log2FoldChange"].abs() >= LOG2FC_THRESHOLD)
        degs_df = res_df[deg_mask].sort_values("padj")
        deg_genes = set(degs_df["gene"])
        degs_list = [
            {
                "gene": row["gene"],
                "log2FoldChange": float(row["log2FoldChange"]),
                "pvalue": float(row["pvalue"]),
                "padj": float(row["padj"]),
            }
            for _, row in degs_df.iterrows()
        ]

        # ---- enrichment ----
        job["stage"] = "enrichment"
        enrichment = {gs: [] for gs in ENRICHR_GENE_SETS}
        if len(deg_genes) > 0:
            try:
                import gseapy as gp

                enr = gp.enrichr(
                    gene_list=list(deg_genes),
                    gene_sets=ENRICHR_GENE_SETS,
                    organism="human",
                    outdir=None,
                )
                results = enr.results
                for gs in ENRICHR_GENE_SETS:
                    subset = results[results["Gene_set"] == gs].sort_values("Adjusted P-value")
                    enrichment[gs] = [
                        {
                            "term": row["Term"],
                            "p_value": float(row["P-value"]),
                            "adjusted_p_value": float(row["Adjusted P-value"]),
                            "genes": row["Genes"].split(";") if row["Genes"] else [],
                        }
                        for _, row in subset.iterrows()
                    ]
            except Exception as e:
                logger.warning("Enrichr call failed for job %s, degrading gracefully: %s", job_id, e)
                enrichment = {gs: [] for gs in ENRICHR_GENE_SETS}
            gc.collect()

        # ---- biomarker_annotation ----
        job["stage"] = "biomarker_annotation"
        input_gene_set = set(genes)
        biomarkers = []
        for gene, info in BIOMARKER_PANEL.items():
            if gene in input_gene_set:
                biomarkers.append(
                    {
                        "gene": gene,
                        "role": info["role"],
                        "relevance": info["relevance"],
                        "in_degs": gene in deg_genes,
                    }
                )

        # ---- visualization ----
        job["stage"] = "visualization"
        volcano_b64 = _make_volcano(res_df, deg_mask)
        heatmap_b64 = _make_heatmap(counts_df, degs_df, metadata_df)

        # ---- report_assembly ----
        job["stage"] = "report_assembly"
        report = {
            "job_id": job_id,
            "summary": {
                "n_genes_input": len(genes),
                "n_samples": len(samples),
                "n_degs": len(degs_list),
                "conditions": [cond_a, cond_b],
            },
            "degs": degs_list,
            "enrichment": enrichment,
            "biomarkers": biomarkers,
            "figures": {
                "volcano_png_base64": volcano_b64,
                "heatmap_png_base64": heatmap_b64,
            },
            "parameters": {
                "log2fc_threshold": LOG2FC_THRESHOLD,
                "padj_threshold": PADJ_THRESHOLD,
            },
        }
        job["result"] = report
        job["status"] = "complete"

    except Exception as e:
        logger.exception("Job %s failed", job_id)
        job["status"] = "failed"
        job["error"] = str(e)


def _make_volcano(res_df: pd.DataFrame, deg_mask: pd.Series) -> str:
    x = res_df["log2FoldChange"].astype(float)
    y = -np.log10(res_df["padj"].astype(float).clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(x[~deg_mask], y[~deg_mask], s=8, c="#9aa0a6", alpha=0.6, label="Not significant")
    ax.scatter(x[deg_mask], y[deg_mask], s=10, c="#d62728", alpha=0.8, label="DEG")
    ax.axvline(LOG2FC_THRESHOLD, color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(-LOG2FC_THRESHOLD, color="grey", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(PADJ_THRESHOLD), color="grey", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(padj)")
    ax.set_title("Volcano plot")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    return _fig_to_base64(fig)


def _make_heatmap(counts_df: pd.DataFrame, degs_df: pd.DataFrame, metadata_df: pd.DataFrame) -> str:
    log_counts = np.log2(counts_df + 1)

    if len(degs_df) >= 1:
        n = min(HEATMAP_MAX_GENES, max(HEATMAP_MIN_GENES, len(degs_df)))
        n = min(n, len(degs_df))
        top_genes = degs_df.sort_values("padj").head(n)["gene"].tolist()
    else:
        variances = log_counts.var(axis=0).sort_values(ascending=False)
        top_genes = variances.head(min(HEATMAP_MIN_GENES, log_counts.shape[1])).index.tolist()

    mat = log_counts[top_genes].T  # genes x samples
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, 1), axis=0)

    categories = list(metadata_df["condition"].cat.categories)
    palette = dict(zip(categories, sns.color_palette("Set2", len(categories))))
    col_colors = pd.Series(
        [palette[c] for c in metadata_df["condition"].astype(str)],
        index=metadata_df.index,
        name="condition",
    )

    g = sns.clustermap(
        z,
        cmap="vlag",
        center=0,
        col_colors=col_colors,
        figsize=(max(6, 0.4 * z.shape[1] + 4), max(6, 0.25 * z.shape[0] + 2)),
        yticklabels=True,
        xticklabels=True,
    )
    g.ax_heatmap.set_xlabel("Sample")
    g.ax_heatmap.set_ylabel("Gene")
    buf = io.BytesIO()
    g.savefig(buf, format="png", bbox_inches="tight", dpi=110)
    plt.close(g.fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@app.post("/analyze", status_code=202)
async def analyze(request: Request, background_tasks: BackgroundTasks):
    payload = await read_json_strict(request)
    counts, metadata = validate_analyze_payload(payload)

    job_id = str(uuid.uuid4())
    JOBS[job_id] = {"job_id": job_id, "status": "queued", "stage": None, "error": None, "result": None}
    background_tasks.add_task(run_pipeline, job_id, counts, metadata)
    return JSONResponse(status_code=202, content={"job_id": job_id, "status": "queued"})


@app.get("/jobs/{job_id}")
async def get_job(job_id: str):
    job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"No job found with id '{job_id}'.")
    return {
        "job_id": job["job_id"],
        "status": job["status"],
        "stage": job["stage"],
        "error": job["error"],
    }


@app.get("/results/{job_id}")
async def get_results(job_id: str):
    job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"No job found with id '{job_id}'.")
    if job["status"] == "complete":
        return job["result"]
    if job["status"] == "failed":
        return JSONResponse(status_code=500, content={"status": "failed", "error": job["error"]})
    return JSONResponse(status_code=409, content={"status": job["status"]})


@app.get("/health")
async def health():
    return {"status": "ok"}


@app.get("/skill.md")
async def skill_md():
    with open("SKILL.md", "r", encoding="utf-8") as f:
        content = f.read()
    return PlainTextResponse(content, media_type="text/markdown")
