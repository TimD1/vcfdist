import vcf
import numpy as np
import matplotlib.pyplot as plt

query_tsv = open("../multi_match/evals/vcfdist/pav.query.tsv")
truth_tsv = open("../multi_match/evals/vcfdist/pav.truth.tsv")

dep_gaps = []
ind_gaps = []
sg_gaps = []

# skip header
next_query = next(query_tsv, False)
next_truth = next(truth_tsv, False)
used_query = True
used_truth = True
start = 0
sc = 0
sg = 0
ctg = ""

class Variant():
    def __init__(self, line):
        if not line: # done parsing variants
            self.ctg = "None"
            self.start = 0
            self.stop = 0
            self.sc = 0
            self.sg = 0
            return

        ctg, pos, hap, ref, alt, qual, vartype, errtype, \
            credit, cluster, sc, sg, loc = line.strip().split("\t")
        self.ctg = ctg
        if len(ref) == len(alt) == 1: # SUB
            self.start = int(pos)
            self.stop = int(pos)+1
        elif len(ref) < len(alt): # INS
            self.start = int(pos)+1
            self.stop = int(pos)+1
        else: # DEL/CPX
            self.start = int(pos)+1
            self.stop = int(pos) + len(ref)-1
        self.sc = int(sc)
        self.sg = int(sg)

    def __repr__(self):
        return "Variant()"
    def __str__(self):
        return f"{self.ctg}:{self.start}-{self.stop} sc={self.sc} sg={self.sg}"

lines = 0
prev_end = 0
prev_sg_end = 0
while next_query or next_truth:

    if used_query:
        next_query = next(query_tsv, False)
        used_query = False
    if used_truth:
        next_truth = next(truth_tsv, False)
        used_truth = False
    query = Variant(next_query)
    truth = Variant(next_truth)

    # get next variant(s)
    if query.ctg == truth.ctg == ctg:
        if query.start == truth.start: # same
            used_query = used_truth = True
        elif query.start < truth.start:
            used_query = True
        else:
            used_truth = True
    elif query.ctg == ctg:
        used_query = True
    elif truth.ctg == ctg:
        used_truth = True
    else: # new contig
        ctg = truth.ctg
        print(ctg)
        prev_end = 0
        prev_sg_end = 0
    if ctg == "None": break

    # update, record gaps
    if used_query and used_truth:
        if query.sc == sc:
            if prev_end:
                dep_gaps.append(max(0, query.start - prev_end))
            prev_end = max(prev_end, max(truth.stop, query.stop))
            if query.sg == sg:
                if prev_sg_end:
                    sg_gaps.append(max(0, query.start - prev_sg_end))
                    if query.start - prev_sg_end > 10000:
                        print(query, truth, prev_sg_end)
                prev_sg_end = max(prev_sg_end, max(truth.stop, query.stop))
            else:
                prev_sg_end = 0
                sg = query.sg
        else:
            if prev_end:
                ind_gaps.append(max(0, query.start - prev_end))
            prev_end = max(prev_end, max(truth.stop, query.stop))
            prev_sg_end = 0
            sc = query.sc

    elif used_query:
        if query.sc == sc:
            if prev_end:
                dep_gaps.append(max(0, query.start - prev_end))
            prev_end = max(prev_end, query.stop)
            if query.sg == sg:
                if prev_sg_end:
                    sg_gaps.append(max(0, query.start - prev_sg_end))
                    if query.start - prev_sg_end > 10000:
                        print(query, truth, prev_sg_end)
                prev_sg_end = max(prev_sg_end, query.stop)
            else:
                prev_sg_end = 0
                sg = query.sg
        else:
            if prev_end:
                ind_gaps.append(max(0, query.start - prev_end))
            prev_end = max(prev_end, query.stop)
            prev_sg_end = 0
            sc = query.sc

    elif used_truth:
        if truth.sc == sc:
            if prev_end:
                dep_gaps.append(max(0, truth.start - prev_end))
            prev_end = max(prev_end, truth.stop)
            if truth.sg == sg:
                if prev_sg_end:
                    sg_gaps.append(max(0, truth.start - prev_sg_end))
                    if truth.start - prev_sg_end > 10000:
                        print(query, truth, prev_sg_end)
                prev_sg_end = max(prev_sg_end, truth.stop)
            else:
                prev_sg_end = 0
                sg = truth.sg
        else:
            if prev_end:
                ind_gaps.append(max(0, truth.start - prev_end))
            prev_end = max(prev_end, truth.stop)
            prev_sg_end = 0
            sc = truth.sc

    lines += 1
    # if lines > 20:
    #     break
dep_gaps = [x for x in dep_gaps if x] # filter 0
sg_gaps = [x for x in sg_gaps if x] # filter 0

plt.figure(figsize=(3.5,1.5))
bins = [1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000]
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize=5)
plt.yticks(fontsize=5)
plt.hist(sg_gaps, bins=bins, color='g', alpha=0.5, label="MULTI match: dependent", histtype='step') 
plt.hist(dep_gaps, bins=bins, color='y', alpha=0.5, label="WFA: possibly dependent", histtype='step') 
plt.hist(ind_gaps, bins=bins, color='r', alpha=0.5, label="WFA: independent", histtype='step') 
plt.xlabel("Distance between adjacent variants", fontsize=7)
plt.ylabel("Counts", fontsize=7)
plt.legend(loc='upper right', fontsize=5)
plt.tight_layout()
plt.savefig("img/gap_histogram.pdf", format="pdf")

query_tsv.close()
truth_tsv.close()
