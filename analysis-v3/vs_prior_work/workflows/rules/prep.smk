# Included by benchmark.smk; expects CACHE, TRUTH, DATASETS, REF, config in scope.

rule split_query:
    wildcard_constraints:
        ds="|".join(DATASETS),
    input:
        vcf=lambda wc: os.path.join(CACHE, "query", f'{wc.ds}-{DATASETS[wc.ds]["version"]}.vcf.gz'),
    output:
        vcf="results/prep/{ds}.most.vcf.gz",
        tbi="results/prep/{ds}.most.vcf.gz.tbi",
    params:
        maxlen=config["params"]["max_sv_len"],
    shell:
        r"""
        mkdir -p results/prep
        bcftools view -i 'TYPE="snp" || (ILEN < {params.maxlen} && ILEN > -{params.maxlen})' {input.vcf} \
            | grep -v "INV" | bgzip -f > {output.vcf}
        tabix -f -p vcf {output.vcf}
        """


rule split_truth:
    input:
        vcf=os.path.join(CACHE, "truth", f'{TRUTH["id"]}-{TRUTH["version"]}.vcf.gz'),
    output:
        vcf="results/prep/truth.most.vcf.gz",
        tbi="results/prep/truth.most.vcf.gz.tbi",
    params:
        maxlen=config["params"]["max_sv_len"],
    shell:
        r"""
        mkdir -p results/prep
        bcftools view -i 'TYPE="snp" || (ILEN < {params.maxlen} && ILEN > -{params.maxlen})' {input.vcf} \
            | grep -v "INV" | bgzip -f > {output.vcf}
        tabix -f -p vcf {output.vcf}
        """


rule build_sdf:
    input:
        ref=os.path.join(CACHE, "refs", REF["name"]),
    output:
        directory("results/prep/ref.sdf"),
    shell:
        "rm -rf {output} && rtg format -o {output} {input.ref}"
