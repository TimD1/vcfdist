
truth_vcf_fn = "./vcfs/t2t-q100.vcf"
query_vcf_fn = "./vcfs/pav.vcf"
out = open("diffs.txt", "w")

# filter positions with multiple variants (like WhatsHap)
print("Finding query duplicates")
query_pos = set()
query_dups = set()
with open(query_vcf_fn, "r") as query_vcf:
    for line in query_vcf:
        if line[0] == "#": continue # skip header
        fields = line.split()
        gt = fields[9].split(":")[0]
        pos = f"{fields[0]}:{fields[1]}:{gt}"
        if pos in query_pos:
            query_dups.add(pos)
        query_pos.add(pos)

print("Finding truth duplicates")
truth_pos = set()
truth_dups = set()
with open(truth_vcf_fn, "r") as truth_vcf:
    for line in truth_vcf:
        if line[0] == "#": continue # skip header
        fields = line.split()
        gt = fields[9].split(":")[0]
        pos = f"{fields[0]}:{fields[1]}:{gt}"
        print(pos)
        if pos in truth_pos:
            truth_dups.add(pos)
        truth_pos.add(pos)

query_vars = []
query_gts = {}
truth_gts = {}
count = 0
print("Parsing query VCF")
with open(query_vcf_fn, "r") as query_vcf:
    for line in query_vcf:
        if line[0] == "#": continue # skip header
        fields = line.split()
        gt = fields[9].split(":")[0]
        if gt == "1|1": continue
        if f"{fields[0]}:{fields[1]}" in query_dups: continue
        var_id = f"{fields[0]}:{fields[1]}={fields[3]}>{fields[4]}"
        query_vars.append(var_id)
        query_gts[var_id] = gt

print("Parsing truth VCF")
with open(truth_vcf_fn, "r") as truth_vcf:
    for line in truth_vcf:
        if line[0] == "#": continue # skip header
        fields = line.split()
        gt = fields[9].split(":")[0]
        if gt == "1|1": continue
        if f"{fields[0]}:{fields[1]}" in truth_dups: continue
        var_id = f"{fields[0]}:{fields[1]}={fields[3]}>{fields[4]}"
        truth_gts[var_id] = gt

print("Printing GT mismatches")
for var in query_vars:
    if var in truth_gts.keys() and query_gts[var] != truth_gts[var]:
        print(f"variant: {var}\tquery:{query_gts[var]}\ttruth:{truth_gts[var]}", file=out)
