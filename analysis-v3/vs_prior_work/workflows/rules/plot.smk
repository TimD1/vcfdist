# Included by benchmark.smk. Concatenates shards, then renders figures under OUTDIR.
SCRIPTS = os.path.join(workflow.basedir, "..", "scripts")
ALL_TOOLS = ["vcfdist-v3", "vcfdist-v2", "vcfeval"]  # "happy" disabled for now
ALL_DS = list(DATASETS)


def _shards(kind):
    return expand(f"{SHARD}/{{tool}}.{{ds}}.{kind}.tsv", tool=ALL_TOOLS, ds=ALL_DS)


rule concat_tables:
    input:
        counts=_shards("counts"), curve=_shards("curve"), runtime=_shards("runtime"),
    output:
        counts=f"{OUTDIR}/parsed/counts.tsv",
        curve=f"{OUTDIR}/parsed/pr_curve.tsv",
        runtime=f"{OUTDIR}/parsed/runtime.tsv",
    run:
        def cat(files, out):
            with open(out, "w") as o:
                for i, f in enumerate(files):
                    with open(f) as fh:
                        lines = fh.readlines()
                    o.writelines(lines if i == 0 else lines[1:])
        cat(input.counts, output.counts)
        cat(input.curve, output.curve)
        cat(input.runtime, output.runtime)


rule plot_fnr_fdr:
    input: counts=f"{OUTDIR}/parsed/counts.tsv"
    output: fnr=f"{OUTDIR}/img/fnr.pdf", fdr=f"{OUTDIR}/img/fdr.pdf"
    shell:
        f"mkdir -p {OUTDIR}/img && python {SCRIPTS}/plot_fnr_fdr.py "
        "--counts {input.counts} --fnr-out {output.fnr} --fdr-out {output.fdr}"


rule plot_pr_curves:
    input: curve=f"{OUTDIR}/parsed/pr_curve.tsv"
    output: pdf=f"{OUTDIR}/img/pr_curves.pdf"
    shell:
        f"mkdir -p {OUTDIR}/img && python {SCRIPTS}/plot_pr_curves.py "
        "--curve {input.curve} --out {output.pdf}"


rule plot_runtime_ram:
    input: runtime=f"{OUTDIR}/parsed/runtime.tsv"
    output: pdf=f"{OUTDIR}/img/runtime_ram.pdf"
    shell:
        f"mkdir -p {OUTDIR}/img && python {SCRIPTS}/plot_runtime_ram.py "
        "--runtime {input.runtime} --out {output.pdf}"


# Confusion matrices: vcfdist-v3's per-variant decision vs each other caller.
OTHER_TOOLS = [t for t in ALL_TOOLS if t != "vcfdist-v3"]


def _cm_inputs(wc):
    # Gate on the reference (v3) and the other tool's declared eval outputs; the
    # script reads summary.vcf (vcfdist) / tp.fp.fn.vcf.gz (vcfeval) co-produced there.
    ins = expand(f"{OUTDIR}/vcfdist-v3/{{ds}}.precision-recall-summary.tsv", ds=ALL_DS)
    if wc.other.startswith("vcfdist"):
        ins += expand(f"{OUTDIR}/{wc.other}/{{ds}}.precision-recall-summary.tsv", ds=ALL_DS)
    else:
        ins += expand(f"{OUTDIR}/vcfeval/{{ds}}/tp.vcf.gz", ds=ALL_DS)
    return ins


rule confusion_matrix:
    wildcard_constraints:
        other="|".join(OTHER_TOOLS),
    input:
        _cm_inputs,
    output:
        tsv=f"{OUTDIR}/parsed/confusion_{{other}}.tsv",
        pdf=f"{OUTDIR}/img/confusion_{{other}}.pdf",
    params:
        ds=" ".join(ALL_DS),
        sv=config["params"]["sv_threshold"],
    shell:
        f"mkdir -p {OUTDIR}/parsed {OUTDIR}/img && python {SCRIPTS}/confusion_matrix.py "
        f"--outdir {OUTDIR} --other-tool {{wildcards.other}} --datasets {{params.ds}} "
        "--sv-threshold {params.sv} --out-tsv {output.tsv} --out-pdf {output.pdf}"
