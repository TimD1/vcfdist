import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

idents = sorted([ os.path.basename(f.path) for 
            f in os.scandir("../pfda-v2/benchmarking_results") ])
genomes = ['HG002', 'HG003', 'HG004']

df = pd.DataFrame()

# for i in idents[:5]:
#     for g in genomes:
#         with open(f"../pfda-v2/benchmarking_results/{i}/{i}_{g}.extended.csv") as csv:
#             for line in csv:
#                 if line[:24] == "INDEL,*,*,PASS,*,QUAL,*,":
#                     indel_recall = line.split(',')[7]
#                     indel_prec = line.split(',')[8]
#                     indel_f1 = line.split(',')[10]
#                 elif line[:22] == "SNP,*,*,PASS,*,QUAL,*,":
#                     snp_recall = line.split(',')[7]
#                     snp_prec = line.split(',')[8]
#                     snp_f1 = line.split(',')[10]
#         row = {'truth': 'bench', 'ident': i, 'genome': g, 
#                 'indel_recall': indel_recall, 'indel_prec': indel_prec, 'indel_f1': indel_f1,
#                 'snp_recall': snp_recall, 'snp_prec': snp_prec, 'snp_f1': snp_f1}
#         df = df.append(row, ignore_index=True)

for i in idents:
    for truth in ["orig", "npore"]:
        try:
            with open(f"results/{i}-all-truth_{truth}-eval421.extended.csv") as csv:
                for line in csv:
                    if line[:24] == "INDEL,*,*,PASS,*,QUAL,*,":
                        indel_recall = line.split(',')[7]
                        indel_prec = line.split(',')[8]
                        indel_f1 = line.split(',')[10]
                    elif line[:22] == "SNP,*,*,PASS,*,QUAL,*,":
                        snp_recall = line.split(',')[7]
                        snp_prec = line.split(',')[8]
                        snp_f1 = line.split(',')[10]
            row = {'truth': truth, 'ident': i, 'genome': 'HG002', 
                    'indel_recall': indel_recall, 'indel_prec': indel_prec, 'indel_f1': indel_f1,
                    'snp_recall': snp_recall, 'snp_prec': snp_prec, 'snp_f1': snp_f1}
            df = df.append(row, ignore_index=True)
        except FileNotFoundError:
            continue
 
df2 = df[df['genome'] == 'HG002']

plt.figure(figsize=(5,12))
for i in idents:
    plt.plot( list(df2[df2['ident']==i]['truth']),
        [-np.log10(1-float(x)) for x in df2[df2['ident']==i]['snp_f1']], marker='.')
plt.ylim(1.5, 3.5)
plt.yticks([1,2,3,4], labels=['90', '99', '99.9', '99.99'])
plt.title("SNP Accuracy (F1)")
plt.xlabel("Calls/Truth VCF Formats")
plt.ylabel("F1 Score (%)")
# plt.legend(idents[:5])
plt.tight_layout()
plt.savefig('img/snp_acc.png')

plt.figure(figsize=(5,12))
for i in idents:
    plt.plot( list(df2[df2['ident']==i]['truth']),
        [-np.log10(1-float(x)) for x in df2[df2['ident']==i]['indel_f1']], marker='.')
plt.ylim(1, 2)
plt.yticks([1,2], labels=['90', '99'])
plt.title("INDEL Accuracy (F1)")
plt.xlabel("Calls/Truth VCF Formats")
plt.ylabel("F1 Score (%)")
# plt.legend(idents[:5])
plt.tight_layout()
plt.savefig('img/indel_acc.png')
