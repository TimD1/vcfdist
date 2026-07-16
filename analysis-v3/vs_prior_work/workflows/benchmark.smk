import os
from snakemake.utils import validate

configfile: os.path.join(workflow.basedir, "config", "config.yml")
validate(config, os.path.join(workflow.basedir, "config", "config_schema.yml"))

DATA = config["data_dir"]          # pre-prepared data tree (CWD-relative -> analysis-v3/data)
OUTDIR = config["outdir"]          # results tree (overridden to results-chr20 by the smoke task)
CACHE = config["cache_dir"]        # only used by the (skipped) download/prep stages
REF = config["reference"]
TRUTH = config["truth"]
DATASETS = config["datasets"]

# Constrain the {ds} wildcard to real dataset ids across all rules.
wildcard_constraints:
    ds="|".join(DATASETS),

# NOTE: download (workflows/download.smk) and prep (workflows/rules/prep.smk) are
# RETAINED for from-scratch reproducibility but SKIPPED here — the benchmark reads
# the already-prepared call sets from DATA (see config data_dir / README). To
# regenerate from scratch, run download.smk + prep.smk, then re-enable the prep
# include below.
# include: "rules/prep.smk"
include: "rules/evaluate.smk"
include: "rules/parse.smk"
include: "rules/plot.smk"


rule all:
    input:
        f"{OUTDIR}/img/fnr.pdf",
        f"{OUTDIR}/img/fdr.pdf",
        f"{OUTDIR}/img/pr_curves.pdf",
        f"{OUTDIR}/img/runtime_ram.pdf",
        expand(f"{OUTDIR}/img/confusion_{{other}}.pdf", other=OTHER_TOOLS),
