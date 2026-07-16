# Included by benchmark.smk. Produces per-(tool,dataset) shard TSVs under OUTDIR.
AGG = os.path.join(workflow.basedir, "..", "scripts", "aggregate.py")
SHARD = f"{OUTDIR}/parsed/shards"


rule parse_vcfdist_v3:
    input:
        summary=f"{OUTDIR}/vcfdist-v3/{{ds}}.precision-recall-summary.tsv",
        curve=f"{OUTDIR}/vcfdist-v3/{{ds}}.precision-recall.tsv",
        log=f"{OUTDIR}/vcfdist-v3/{{ds}}.time.log",
    output:
        counts=f"{SHARD}/vcfdist-v3.{{ds}}.counts.tsv",
        curve=f"{SHARD}/vcfdist-v3.{{ds}}.curve.tsv",
        runtime=f"{SHARD}/vcfdist-v3.{{ds}}.runtime.tsv",
    shell:
        f"mkdir -p {SHARD} && python {AGG} --tool vcfdist-v3 --dataset {{wildcards.ds}} "
        "--summary {input.summary} --curve {input.curve} --time-log {input.log} "
        "--counts-out {output.counts} --curve-out {output.curve} --runtime-out {output.runtime}"


rule parse_vcfdist_v2:
    input:
        summary=f"{OUTDIR}/vcfdist-v2/{{ds}}.precision-recall-summary.tsv",
        curve=f"{OUTDIR}/vcfdist-v2/{{ds}}.precision-recall.tsv",
        log=f"{OUTDIR}/vcfdist-v2/{{ds}}.time.log",
    output:
        counts=f"{SHARD}/vcfdist-v2.{{ds}}.counts.tsv",
        curve=f"{SHARD}/vcfdist-v2.{{ds}}.curve.tsv",
        runtime=f"{SHARD}/vcfdist-v2.{{ds}}.runtime.tsv",
    shell:
        f"mkdir -p {SHARD} && python {AGG} --tool vcfdist-v2 --dataset {{wildcards.ds}} "
        "--summary {input.summary} --curve {input.curve} --time-log {input.log} "
        "--counts-out {output.counts} --curve-out {output.curve} --runtime-out {output.runtime}"


rule parse_vcfeval:
    input:
        tp=f"{OUTDIR}/vcfeval/{{ds}}/tp.vcf.gz", log=f"{OUTDIR}/vcfeval/{{ds}}.time.log",
    output:
        counts=f"{SHARD}/vcfeval.{{ds}}.counts.tsv",
        curve=f"{SHARD}/vcfeval.{{ds}}.curve.tsv",
        runtime=f"{SHARD}/vcfeval.{{ds}}.runtime.tsv",
    shell:
        f"mkdir -p {SHARD} && python {AGG} --tool vcfeval --dataset {{wildcards.ds}} "
        f"--dir {OUTDIR}/vcfeval/{{wildcards.ds}} --time-log {{input.log}} "
        "--counts-out {output.counts} --curve-out {output.curve} --runtime-out {output.runtime}"


# --- hap.py disabled for now (re-enable alongside eval_happy in evaluate.smk) --
# rule parse_happy:
#     input:
#         vcf=f"{OUTDIR}/happy/{{ds}}.vcf.gz", log=f"{OUTDIR}/happy/{{ds}}.time.log",
#     output:
#         counts=f"{SHARD}/happy.{{ds}}.counts.tsv",
#         curve=f"{SHARD}/happy.{{ds}}.curve.tsv",
#         runtime=f"{SHARD}/happy.{{ds}}.runtime.tsv",
#     shell:
#         f"mkdir -p {SHARD} && python {AGG} --tool happy --dataset {{wildcards.ds}} "
#         f"--vcf {OUTDIR}/happy/{{wildcards.ds}}.vcf.gz --time-log {{input.log}} "
#         "--counts-out {output.counts} --curve-out {output.curve} --runtime-out {output.runtime}"
