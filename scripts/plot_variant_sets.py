import pysam
import pandas as pd
from time import perf_counter
import upsetplot
import matplotlib.pyplot as plt

names = ["truth_orig", "truth_npore"]
vcfs = [f"/home/timdunn/vcfdist/scripts/data/{name}.vcf.gz" for name in names]
cats = ["sub", "ins", "del"]
variants = {}
for cat in cats:
    variants[cat] = {}
    for name in names:
        variants[cat][name] = set()

# compute set of variants for each VCF
for name in names:
    vcf_fn = f"/home/timdunn/vcfdist/scripts/data/{name}.vcf.gz"
    vcf = pysam.VariantFile(vcf_fn, 'r')
    for record in vcf.fetch():
        if len(record.alleles) != 2:
            print("ERROR: VCFs should not contain complex variants.")
        len_diff = len(record.alleles[1]) - len(record.alleles[0])
        if len_diff > 0:
            cat = "ins"
        elif len_diff < 0:
            cat = "del"
        else:
            cat = "sub"
        variants[cat][name].add(
                f"{record.contig}:{record.start} {'|'.join(record.alleles)}")

# debug print
for cat in cats:
    for name in names:
        print(f"{name}: {len(variants[cat][name])} {cat}s")

for cat in cats:
    print("\ncat:")
    df = pd.DataFrame(columns = names, 
            index = range(len(set().union(*variants[cat].values()))))
    for idx, variant in enumerate(set().union(*variants[cat].values())):
        row = {}
        for name in names:
            row[name] = variant in variants[cat][name]
        df.loc[idx] = row
    df['count'] = 1
    counts = df.groupby(names).count().sort_values('count')
    upsetplot.plot(counts['count'], sort_by="cardinality")
    plt.savefig(f"img/{cat}_variants.png")


