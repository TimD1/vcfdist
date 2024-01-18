import statistics
import vcf

gaps = [10, 50, 100, 500, 1000]
POS = 0
SIZE = 1

# parse truth/query VCFs
ctgs = []
print("parse truth")
truth = {}
truth_vcf = vcf.Reader(open(f"t2t_t2t-bed.vcf", "r"))
for truth_var in truth_vcf:
    if truth_var.CHROM not in truth.keys():
        truth[truth_var.CHROM] = []
    if len(truth_var.REF) == len(truth_var.ALT): # SUB
        truth[truth_var.CHROM].append((truth_var.POS-1, len(truth_var.REF)+2))
    elif len(truth_var.REF) < len(truth_var.ALT): # INS
        truth[truth_var.CHROM].append((truth_var.POS, 1))
    else: # DEL
        truth[truth_var.CHROM].append((truth_var.POS, len(truth_var.REF)+1))
    if truth_var.CHROM not in ctgs:
        ctgs.append(truth_var.CHROM)

print("parse query")
query = {}
query_vcf = vcf.Reader(open(f"pav_t2t-bed.vcf", "r"))
for query_var in query_vcf:
    if query_var.CHROM == "chrM" or query_var.CHROM == "chrEBV": continue
    if query_var.CHROM not in query.keys():
        query[query_var.CHROM] = []
    if len(query_var.REF) == len(query_var.ALT): # SUB
        query[query_var.CHROM].append((query_var.POS-1, len(query_var.REF)+2))
    elif len(query_var.REF) < len(query_var.ALT): # INS
        query[query_var.CHROM].append((query_var.POS, 1))
    else: # DEL
        query[query_var.CHROM].append((query_var.POS, len(query_var.REF)+1))
    if query_var.CHROM not in ctgs:
        ctgs.append(query_var.CHROM)


print("clustering")
for gap in gaps:
    outfile = open(f"gap{gap}.log", "w")
    regions = []

    for ctg in ctgs:
        ti = qi = 0

        if truth[ctg][ti][POS] < query[ctg][qi][POS]:
            var = truth[ctg][ti]; ti += 1
        else:
            var = query[ctg][qi]; qi += 1
        start_pos = var[POS]
        end_pos = var[POS] + var[SIZE]


        while ti < len(truth[ctg]) or qi < len(query[ctg]): # while variants remain

            if not ti < len(truth[ctg]): # only query left
                var = query[ctg][qi]; qi += 1
                if var[POS] <= end_pos + gap: # grow cluster
                    end_pos = max(end_pos, var[POS] + var[SIZE])
                else: # new cluster
                    regions.append(end_pos - start_pos)
                    print(f"{ctg}\t{start_pos}\t{end_pos}\t{end_pos-start_pos}", file=outfile)
                    start_pos = var[POS]
                    end_pos = var[POS] + var[SIZE]

            elif not qi < len(query[ctg]): # only truth left
                var = truth[ctg][ti]; ti += 1
                if var[POS] <= end_pos + gap: # grow cluster
                    end_pos = max(end_pos, var[POS] + var[SIZE])
                else: # new cluster
                    regions.append(end_pos - start_pos)
                    print(f"{ctg}\t{start_pos}\t{end_pos}\t{end_pos-start_pos}", file=outfile)
                    start_pos = var[POS]
                    end_pos = var[POS] + var[SIZE]

            else: # both vars possible

                var = (0, 0) # choose next var
                if truth[ctg][ti][POS] < query[ctg][qi][POS]:
                    var = truth[ctg][ti]; ti += 1
                else:
                    var = query[ctg][qi]; qi += 1

                if var[POS] <= end_pos + gap: # grow cluster
                    end_pos = max(end_pos, var[POS] + var[SIZE])
                else: # new cluster
                    regions.append(end_pos - start_pos)
                    print(f"{ctg}\t{start_pos}\t{end_pos}\t{end_pos-start_pos}", file=outfile)
                    start_pos = var[POS]
                    end_pos = var[POS] + var[SIZE]

        regions.append(end_pos - start_pos)
        print(f"{ctg}\t{start_pos}\t{end_pos}\t{end_pos-start_pos}", file=outfile)

    
    print(f"\nGap {gap}")
    print(f"    regions: {len(regions)}")
    print(f"    total:   {sum(regions)}")
    print(f"    max:     {max(regions)}")
    print(f"    median:  {statistics.median(regions)}")
    print(f"    mean:    {statistics.mean(regions)}")

