import os
from snakemake.utils import validate

configfile: os.path.join(workflow.basedir, "config", "config.yml")
validate(config, os.path.join(workflow.basedir, "config", "config_schema.yml"))

CACHE = config["cache_dir"]
REF = config["reference"]
TRUTH = config["truth"]
DATASETS = config["datasets"]

# Path to the pre-existing reference in the repo, moved into cache if present.
REPO_REF = os.path.join(workflow.basedir, "..", "..", "..", "data", "refs", REF["name"])


rule all:
    input:
        os.path.join(CACHE, "refs", REF["name"]),
        os.path.join(CACHE, "truth", f'{TRUTH["id"]}-{TRUTH["version"]}.vcf.gz'),
        os.path.join(CACHE, "truth", "bench.bed"),
        expand(os.path.join(CACHE, "query", "{ds}-{ver}.vcf.gz"),
               zip, ds=list(DATASETS), ver=[d["version"] for d in DATASETS.values()]),


rule reference:
    output:
        os.path.join(CACHE, "refs", REF["name"]),
    params:
        url=REF["url"], repo_ref=REPO_REF,
    shell:
        r"""
        mkdir -p $(dirname {output})
        if [ -f "{params.repo_ref}" ]; then
            echo "Moving existing repo reference into cache (no download)."
            mv "{params.repo_ref}" "{output}"
        else
            echo "Repo reference not found; downloading."
            tmp="{output}.download"
            curl -fL -o "$tmp" "{params.url}"
            case "{params.url}" in
                *.gz) gunzip -c "$tmp" > "{output}" && rm -f "$tmp" ;;
                *)    mv "$tmp" "{output}" ;;
            esac
        fi
        """


rule reference_index:
    input:
        os.path.join(CACHE, "refs", REF["name"]),
    output:
        os.path.join(CACHE, "refs", REF["name"] + ".fai"),
    shell:
        "samtools faidx {input}"


rule truth_vcf:
    output:
        os.path.join(CACHE, "truth", f'{TRUTH["id"]}-{TRUTH["version"]}.vcf.gz'),
    params:
        url=TRUTH["vcf_url"],
    shell:
        'mkdir -p $(dirname {output}) && curl -fL -o {output} "{params.url}"'


rule truth_bed:
    output:
        os.path.join(CACHE, "truth", "bench.bed"),
    params:
        url=TRUTH["bed_url"],
    shell:
        'mkdir -p $(dirname {output}) && curl -fL -o {output} "{params.url}"'


rule query_vcf:
    output:
        os.path.join(CACHE, "query", "{ds}-{ver}.vcf.gz"),
    params:
        url=lambda wc: DATASETS[wc.ds]["vcf_url"],
    shell:
        'mkdir -p $(dirname {output}) && curl -fL -o {output} "{params.url}"'
