# Included by benchmark.smk; expects config, REF, TRUTH, DATASETS, DATA, OUTDIR in scope.
# Reads the already-prepared call sets directly from the DATA tree (analysis-v3/data),
# so the download/prep stages are skipped for this run.
# GNU time from the conda `time` package (macOS /usr/bin/time is BSD, no -v).
# `command` bypasses the bash `time` keyword so the external binary runs.
TIME = "command time -v"
THREADS = config["params"]["threads"]
MAX_RAM = config["params"]["max_ram_gb"]
MAX_SV = config["params"]["max_sv_len"]
CT = config["params"]["credit_threshold"]   # vcfdist -ct, set the same for v2 & v3
V3_BIN = config["params"]["vcfdist_v3_bin"]
# v3-only flags (v2.6.4 lacks --max-dist / --max-retries):
V3_EXTRA = f'-md {config["params"]["max_dist"]} -mr {config["params"]["max_retries"]}'

REF_FA = os.path.join(DATA, "refs", REF["name"])
SDF = os.path.join(DATA, "refs", REF["name"] + ".sdf")
TRUTH_DIR = os.path.join(DATA, f'{TRUTH["id"]}-{TRUTH["version"]}')
TRUTH_VCF = os.path.join(TRUTH_DIR, "split", f'{TRUTH["id"]}.most.vcf.gz')
BED = os.path.join(TRUTH_DIR, "split", config["bench_bed"])   # bench.bed | bench-chr20.bed


def _query_vcf(ds):
    """Path to a dataset's prepared query call set (<data>/<id>-<ver>/split/<id>.most.vcf.gz)."""
    return os.path.join(DATA, f'{ds}-{DATASETS[ds]["version"]}', "split", f'{ds}.most.vcf.gz')


def _vcfdist_shell(tool, binpath, extra=""):
    """Shell command for a vcfdist run; {input.*}/{wildcards.ds} are filled by Snakemake.

    `extra` carries version-specific flags (v3-only max-dist/max-retries).
    """
    return (
        f'mkdir -p {OUTDIR}/{tool} && '
        f'{TIME} {binpath} {{input.q}} {{input.t}} {{input.ref}} '
        f'--bed {{input.bed}} -t {THREADS} -r {MAX_RAM} -l {MAX_SV} -ct {CT} {extra} '
        f'-p {OUTDIR}/{tool}/{{wildcards.ds}}. '
        f'2> {OUTDIR}/{tool}/{{wildcards.ds}}.time.log'
    )


rule eval_vcfdist_v3:
    input:
        q=lambda wc: _query_vcf(wc.ds), t=TRUTH_VCF, ref=REF_FA, bed=BED,
    output:
        summary=f"{OUTDIR}/vcfdist-v3/{{ds}}.precision-recall-summary.tsv",
        curve=f"{OUTDIR}/vcfdist-v3/{{ds}}.precision-recall.tsv",
        log=f"{OUTDIR}/vcfdist-v3/{{ds}}.time.log",
    shell:
        _vcfdist_shell("vcfdist-v3", V3_BIN, V3_EXTRA)


rule eval_vcfdist_v2:
    input:
        q=lambda wc: _query_vcf(wc.ds), t=TRUTH_VCF, ref=REF_FA, bed=BED,
    output:
        summary=f"{OUTDIR}/vcfdist-v2/{{ds}}.precision-recall-summary.tsv",
        curve=f"{OUTDIR}/vcfdist-v2/{{ds}}.precision-recall.tsv",
        log=f"{OUTDIR}/vcfdist-v2/{{ds}}.time.log",
    shell:
        _vcfdist_shell("vcfdist-v2", "vcfdist")


rule eval_vcfeval:
    input:
        q=lambda wc: _query_vcf(wc.ds), t=TRUTH_VCF, sdf=SDF, bed=BED,
    output:
        tp=f"{OUTDIR}/vcfeval/{{ds}}/tp.vcf.gz",
        log=f"{OUTDIR}/vcfeval/{{ds}}.time.log",
    shell:
        r"""
        rm -rf {OUTDIR}/vcfeval/{wildcards.ds}
        mkdir -p {OUTDIR}/vcfeval
        {TIME} rtg vcfeval -b {input.t} -c {input.q} -t {input.sdf} \
            --bed-regions {input.bed} --evaluation-regions {input.bed} \
            --threads {THREADS} --ref-overlap --all-records --vcf-score-field=QUAL \
            -o {OUTDIR}/vcfeval/{wildcards.ds} 2> {OUTDIR}/vcfeval/{wildcards.ds}.time.log
        """


# --- hap.py disabled for now (benchmarking vcfdist v2/v3 + vcfeval only) -------
# Re-enable by uncommenting this rule, restoring "happy" in plot.smk ALL_TOOLS and
# the plotter TOOLS lists, and uncommenting rule parse_happy in parse.smk.
# hap.py runs in the isolated linux-64 `happy` pixi env (Python-2).
#
# rule eval_happy:
#     input:
#         q=lambda wc: _query_vcf(wc.ds), t=TRUTH_VCF, ref=REF_FA, bed=BED,
#     output:
#         vcf=f"{OUTDIR}/happy/{{ds}}.vcf.gz",
#         log=f"{OUTDIR}/happy/{{ds}}.time.log",
#     shell:
#         r"""
#         mkdir -p {OUTDIR}/happy
#         {TIME} pixi run -e happy hap.py {input.t} {input.q} -r {input.ref} -f {input.bed} \
#             --engine=xcmp --threads {THREADS} -o {OUTDIR}/happy/{wildcards.ds} \
#             2> {OUTDIR}/happy/{wildcards.ds}.time.log
#         """
