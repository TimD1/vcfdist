# vcfdist v3 — Benchmarking vs. Prior Work

Reproducible pixi + Snakemake harness comparing **vcfdist v3 (dev)** against
**vcfdist v2.6.4**, **vcfeval**, and **hap.py (xcmp)** on the T2T-Q100 (HG002)
truth set with HPRC / PAV / GIAB-TR query call-sets. Refactor of
`analysis-v2/vs_prior_work/`. See `docs/D7_vcfdist-v3-benchmarking-design.md`.

> **Current run:** hap.py is **disabled** (benchmarking vcfdist v2/v3 + vcfeval),
> and the benchmark reads **already-prepared data** from `../data` (config
> `data_dir` → `analysis-v3/data`), so the download and prep stages are skipped.

## Setup
```bash
pixi install
```

## Run (current: provided pre-prepared data)
```bash
pixi run test                                       # unit tests (default env)
pixi run snakemake -n -s workflows/benchmark.smk    # dry-run the DAG (default env)

# The real evaluation needs the linux-64 `bench` env (rtg-tools + vcfdist):
pixi run -e bench benchmark                          # evaluate -> parse -> plot
```
The `benchmark` task runs serially (`-j1`) so each tool's runtime/RAM is measured
without contention. Set `params.vcfdist_v3_bin` in `workflows/config/config.yml`
to your built `dev` `src/vcfdist` (default `../../src/vcfdist`); build it first.

## Regenerating data from scratch (download + prep)
The download (`workflows/download.smk`) and prep (`workflows/rules/prep.smk`)
stages are **retained for reproducibility** but not run above. To regenerate:
fill in the `TODO:` URLs in `config.yml`, run `pixi run -e bench download`, run
prep, then re-enable the `include: "rules/prep.smk"` line in `benchmark.smk`.

## Re-enabling hap.py
Uncomment `rule eval_happy` (`evaluate.smk`) and `rule parse_happy` (`parse.smk`),
and add `"happy"` back to `ALL_TOOLS` (`plot.smk`) and `TOOLS` (the three
`scripts/plot_*.py`). hap.py runs in the isolated linux-64 `happy` env.

## Figures (`results/img/`)
- `fnr.pdf` / `fdr.pdf` — FNR = FN/(FN+TP); "FPR" is reported as FDR = FP/(FP+TP) = 1−precision.
- `pr_curves.pdf` — precision-recall curves, all tools overlaid per size class × dataset.
  vcfdist curves are native (partial-credit); vcfeval/hap.py curves are QUAL-swept from classified calls.
- `runtime_ram.pdf` — wall-clock time and peak RAM per tool × dataset.

## Notes
- Size classes: SNP; INDEL (<50 bp); SV (≥50 bp, ≤1000 bp). Truvari is not included.
- `data_dir` layout: `<id>-<version>/split/<id>.most.vcf.gz`, truth
  `<id>-<version>/split/bench.bed`, and `refs/<name>` (+ `.fai`, `.sdf`).
