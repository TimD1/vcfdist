import csv
import matplotlib.pyplot as plt

truth = "../out/truth.tsv"
calls = "../out/calls.tsv"
types = ["SUB", "INS", "DEL", "GRP"]
errtypes = ["TP", "FP", "PP", "FN"]
max_qual = 100
eps = 0.0001

# print truth summary stats
truth_data = {}
for t in types:
    for et in errtypes:
        truth_data[f"{t} {et}"] = [0]*max_qual
with open(truth) as truth_fd:
    truth_rd = csv.reader(truth_fd, delimiter="\t")
    next(truth_rd)
    for row in truth_rd:
        contig, pos, hap, ref, alt, qual, typ, errtype, \
                credit, orig_gt, cluster, loc = row
        if errtype == "PP" and float(credit) < eps: errtype = "FN"
        truth_data[f"{typ} {errtype}"][int(float(qual))] += 1
snps = sum([sum(truth_data[f'SUB {et}']) for et in errtypes])
indels = sum([sum(truth_data[f'INS {et}']) for et in errtypes]) + \
         sum([sum(truth_data[f'DEL {et}']) for et in errtypes])
print("Truth:")
print(f"  SNP:   {snps}")
print(f"  INDEL: {indels}")
print(f"  Total: {snps+indels}")

# print calls summary stats
calls_data = {}
for t in types:
    for et in errtypes:
        calls_data[f"{t} {et}"] = [0]*max_qual
with open(calls) as calls_fd:
    calls_rd = csv.reader(calls_fd, delimiter="\t")
    next(calls_rd)
    for row in calls_rd:
        contig, pos, hap, ref, alt, qual, typ, errtype, \
                credit, orig_gt, cluster, loc = row
        if errtype == "PP" and float(credit) < eps: errtype = "FP"
        calls_data[f"{typ} {errtype}"][int(float(qual))] += 1
print("\nCalls:")
print(f"  SNP:    ", end="")
for et in errtypes:
    if et == "FN":
        print(f"{et}={sum(truth_data[f'SUB {et}']):<8d}", end="")
    else:
        print(f"{et}={sum(calls_data[f'SUB {et}']):<8d}", end="")
print(" ")
print(f"  INDEL:  ", end="")
for et in errtypes:
    if et == "FN":
        print(f"{et}={sum(truth_data[f'INS {et}']) + sum(truth_data[f'DEL {et}']):<8d}", end="")
    else:
        print(f"{et}={sum(calls_data[f'INS {et}']) + sum(calls_data[f'DEL {et}']):<8d}", end="")
print(" ")

# get PP data
pp_data = {}
for t in types:
    pp_data[f"{t}"] = [0]*max_qual
with open(calls) as calls_fd:
    calls_rd = csv.reader(calls_fd, delimiter="\t")
    next(calls_rd)
    for row in calls_rd:
        contig, pos, hap, ref, alt, qual, typ, errtype, \
                credit, orig_gt, cluster, loc = row
        if errtype == "PP" and float(credit) > eps:
            pp_data[f"{typ}"][int(float(qual))] += float(credit)


fig, ax = plt.subplots(1,2)

# PESSIMISTIC: PP = FP + FN
pess_snp_p = [] # snp
pess_snp_r = []
calls_tp = 0
calls_fp = 0
calls_fn = sum(truth_data['SUB FN']) + \
        sum(calls_data['SUB TP']) + sum(calls_data['SUB PP'])
truth_tp = 0
truth_fn = 0
for i in range(max_qual-1, -1, -1):
    calls_tp += calls_data['SUB TP'][i]
    calls_fp += calls_data['SUB FP'][i] + calls_data['SUB PP'][i]
    calls_fn -= calls_data['SUB TP'][i]
    snp_p = calls_tp/max(1, calls_tp + calls_fp)
    snp_r = calls_tp/max(1, calls_tp + calls_fn)
    pess_snp_p.append(snp_p)
    pess_snp_r.append(snp_r)
pess_indel_p = [] # indel
pess_indel_r = []
calls_tp = 0
calls_fp = 0
calls_fn = sum(truth_data['INS FN']) + sum(truth_data['DEL FN']) + \
        sum(calls_data['INS TP']) + sum(calls_data['DEL TP']) + \
        sum(calls_data['INS PP']) + sum(calls_data['DEL PP'])
truth_tp = 0
truth_fn = 0
for i in range(max_qual-1, -1, -1):
    calls_tp += calls_data['INS TP'][i] + calls_data['DEL TP'][i]
    calls_fp += calls_data['INS FP'][i] + calls_data['DEL FP'][i] + \
                calls_data['INS PP'][i] + calls_data['DEL PP'][i]
    calls_fn -= calls_data['INS TP'][i] + calls_data['DEL TP'][i]
    indel_p = calls_tp/max(1, calls_tp + calls_fp)
    indel_r = calls_tp/max(1, calls_tp + calls_fn)
    pess_indel_p.append(indel_p)
    pess_indel_r.append(indel_r)

# OPTIMISTIC: PP = TP
opt_snp_p = [] # snp
opt_snp_r = []
calls_tp = 0
calls_fp = 0
calls_fn = sum(truth_data['SUB FN']) + \
           sum(calls_data['SUB TP']) + sum(calls_data['SUB PP'])
truth_tp = 0
truth_fn = 0
for i in range(max_qual-1, -1, -1):
    calls_tp += calls_data['SUB TP'][i] + calls_data['SUB PP'][i]
    calls_fp += calls_data['SUB FP'][i]
    calls_fn -= calls_data['SUB TP'][i] + calls_data['SUB PP'][i]
    snp_p = calls_tp/max(1, calls_tp + calls_fp)
    snp_r = calls_tp/max(1, calls_tp + calls_fn)
    opt_snp_p.append(snp_p)
    opt_snp_r.append(snp_r)
opt_indel_p = [] # indel
opt_indel_r = []
calls_tp = 0
calls_fp = 0
calls_fn = sum(truth_data['INS FN']) + sum(truth_data['DEL FN']) + \
        sum(calls_data['INS TP']) + sum(calls_data['DEL TP']) + \
        sum(calls_data['INS PP']) + sum(calls_data['DEL PP'])
truth_tp = 0
truth_fn = 0
for i in range(max_qual-1, -1, -1):
    calls_tp += calls_data['INS TP'][i] + calls_data['DEL TP'][i] + \
                calls_data['INS PP'][i] + calls_data['DEL PP'][i]
    calls_fp += calls_data['INS FP'][i] + calls_data['DEL FP'][i]
    calls_fn -= calls_data['INS TP'][i] + calls_data['DEL TP'][i] + \
                calls_data['INS PP'][i] + calls_data['DEL PP'][i]
    indel_p = calls_tp/max(1, calls_tp + calls_fp)
    indel_r = calls_tp/max(1, calls_tp + calls_fn)
    opt_indel_p.append(indel_p)
    opt_indel_r.append(indel_r)

# VCFDIST: PP = (credit) TP, (1-credit) FP/FN
vcfd_snp_p = [] # snp
vcfd_snp_r = []
calls_tp = 0
calls_fp = 0
calls_fn =  sum(truth_data['SUB FN']) + \
            sum(calls_data['SUB TP']) + \
            sum(calls_data['SUB PP'])
truth_tp = 0
truth_fn = 0
for i in range(max_qual-1, -1, -1):
    calls_tp += calls_data['SUB TP'][i] + pp_data['SUB'][i]
    calls_fp += calls_data['SUB FP'][i] + \
            calls_data['SUB PP'][i] - pp_data['SUB'][i]
    calls_fn -= calls_data['SUB TP'][i] + pp_data['SUB'][i]
    snp_p = calls_tp/max(1, calls_tp + calls_fp)
    snp_r = calls_tp/max(1, calls_tp + calls_fn)
    vcfd_snp_p.append(snp_p)
    vcfd_snp_r.append(snp_r)
vcfd_indel_p = [] # indel
vcfd_indel_r = []
calls_tp = 0
calls_fp = 0
calls_fn = sum(truth_data['INS FN']) + sum(truth_data['DEL FN']) + \
        sum(calls_data['INS TP']) + sum(calls_data['DEL TP']) + \
        sum(calls_data['INS PP']) + sum(calls_data['DEL PP'])
truth_tp = 0
truth_fn = 0
for i in range(max_qual-1, -1, -1):
    calls_tp += calls_data['INS TP'][i] + calls_data['DEL TP'][i] + \
                pp_data['INS'][i] + pp_data['DEL'][i]
    calls_fp += calls_data['INS FP'][i] + calls_data['DEL FP'][i] + \
            calls_data['INS PP'][i] - pp_data['INS'][i] + \
            calls_data['DEL PP'][i] - pp_data['DEL'][i]
    calls_fn -= calls_data['INS TP'][i] + calls_data['DEL TP'][i] + \
                pp_data['INS'][i] + pp_data['DEL'][i]
    indel_p = calls_tp/max(1, calls_tp + calls_fp)
    indel_r = calls_tp/max(1, calls_tp + calls_fn)
    vcfd_indel_p.append(indel_p)
    vcfd_indel_r.append(indel_r)

# plot
fig.set_size_inches(12,6)
ax[0].plot(pess_snp_r, pess_snp_p)
ax[1].plot(pess_indel_r, pess_indel_p)
ax[0].plot(opt_snp_r, opt_snp_p)
ax[1].plot(opt_indel_r, opt_indel_p)
ax[0].plot(vcfd_snp_r, vcfd_snp_p)
ax[1].plot(vcfd_indel_r, vcfd_indel_p)
ax[0].set_ylim(0.95,1)
ax[0].set_xlim(0.95,1)
ax[0].set_ylabel('SNP Precision')
ax[0].set_xlabel('SNP Recall')
ax[1].set_ylabel('INDEL Precision')
ax[1].set_xlabel('INDEL Recall')
ax[1].legend(["Pessimistic Partial Positive", "Optimistic Partial Positive", "VCFdist Partial Positive"])
plt.savefig('img/vcfdist_pr.png', dpi=200)
