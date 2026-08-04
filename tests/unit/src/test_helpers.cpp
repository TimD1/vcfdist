/**
 * @file test_helpers.cpp
 * @brief Shared unit-test scaffolding: global-state fixture, temporary files, in-memory builders.
 */
#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "test_helpers.h"

/* Global state fixture ***************************************************************************/

/**
 * @brief Saves the current global configuration and sets verbosity to 0.
 */
GlobalsGuard::GlobalsGuard() : saved(g) {
    g.verbosity = 0;
}

/**
 * @brief Restores the saved global configuration.
 *
 * `Globals::VERSION` and `Globals::PROGRAM` are static, so no member is const and the implicitly
 * generated copy assignment operator restores every field at once.
 */
GlobalsGuard::~GlobalsGuard() {
    g = this->saved;
}

/* Temporary files ********************************************************************************/

/**
 * @brief Creates a uniquely-named directory beneath the system temporary directory.
 * @param[in] prefix Prefix for the directory name
 * @throws ERROR if the directory cannot be created
 */
TempDir::TempDir(const std::string & prefix) {
    // mkdtemp() overwrites the six trailing X characters in place, which non-const data() allows
    std::string tmpl = (std::filesystem::temp_directory_path() / (prefix + "_XXXXXX")).string();
    if (mkdtemp(tmpl.data()) == NULL) {
        ERROR("Failed to create temporary directory '%s'", tmpl.data());
    }
    this->dir = tmpl;
}

/**
 * @brief Recursively removes the directory and everything beneath it.
 */
TempDir::~TempDir() {
    std::error_code ec;
    std::filesystem::remove_all(this->dir, ec);
}

/**
 * @brief Returns the directory path, without a trailing slash.
 * @return Absolute path of the temporary directory
 */
std::string TempDir::path() const {
    return this->dir;
}

/**
 * @brief Returns the path of a named entry within the directory.
 * @param[in] name Entry name to join onto the directory path
 * @return Absolute path of the named entry
 */
std::string TempDir::path(const std::string & name) const {
    return this->dir + "/" + name;
}

/**
 * @brief Writes a minimal single-sample VCF into a temporary directory and returns its path.
 * @param[in] dir Temporary directory that owns the written file
 * @param[in] records Record lines, written verbatim in order after the #CHROM line
 * @param[in] opts Header lines and sample name to write
 * @return Path of the written VCF
 * @throws ERROR if the VCF cannot be opened for writing
 */
std::string write_tmp_vcf(const TempDir & dir, const std::vector<std::string> & records,
        const vcf_opts & opts) {
    std::string vcf_fn = dir.path(opts.filename);
    std::ofstream out(vcf_fn);
    if (!out.is_open()) {
        ERROR("Failed to open temporary VCF '%s' for writing", vcf_fn.data());
    }
    for (const std::string & line : opts.meta) out << line << "\n";
    for (const std::string & line : opts.contigs) out << line << "\n";
    for (const std::string & line : opts.filters) out << line << "\n";
    for (const std::string & line : opts.formats) out << line << "\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << opts.sample << "\n";
    for (const std::string & rec : records) out << rec << "\n";
    out.close();
    return vcf_fn;
}

/**
 * @brief Writes a BED file into a temporary directory and returns its path.
 * @param[in] dir Temporary directory that owns the written file
 * @param[in] lines Record lines, written verbatim in order
 * @param[in] name Basename of the written file
 * @return Path of the written BED
 * @throws ERROR if the BED cannot be opened for writing
 */
std::string write_tmp_bed(const TempDir & dir, const std::vector<std::string> & lines,
        const std::string & name) {
    std::string bed_fn = dir.path(name);
    std::ofstream out(bed_fn);
    if (!out.is_open()) {
        ERROR("Failed to open temporary BED '%s' for writing", bed_fn.data());
    }
    for (const std::string & line : lines) out << line << "\n";
    out.close();
    return bed_fn;
}

/**
 * @brief Returns the path of a checked-in fixture under tests/unit/data/.
 *
 * The test binary is run both from tests/ (by pytest-workflow) and from tests/unit/build/ (by
 * hand), so each plausible location is tried in turn.
 * @param[in] name Fixture basename, such as "tiny.fasta"
 * @return Path of the fixture, relative to the current working directory
 * @throws ERROR if the fixture is not found, which is what an over-broad .gitignore looks like
 */
std::string data_path(const std::string & name) {
    std::vector<std::string> dirs;
    const char * repo = std::getenv("VCFDIST_REPO_PATH");
    if (repo != NULL) dirs.push_back(std::string(repo) + "/tests/unit/data/");
    dirs.push_back("unit/data/");       // run from tests/
    dirs.push_back("../data/");         // run from tests/unit/build/
    dirs.push_back("tests/unit/data/"); // run from the repository root
    for (const std::string & candidate : dirs) {
        if (std::filesystem::exists(candidate + name)) return candidate + name;
    }
    ERROR("Unit-test fixture '%s' not found; is tests/unit/data/ tracked by git?", name.data());
    return "";
}

/**
 * @brief Reads an entire file into a string.
 * @param[in] fn Input filename
 * @return File contents, or an empty string if the file cannot be opened
 */
std::string read_text(const std::string & fn) {
    std::ifstream in(fn);
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

/**
 * @brief Builds header options for one sample, declaring each named contig with one length.
 * @param[in] callset QUERY or TRUTH callset identifier, which names the sample
 * @param[in] contigs Contig names to declare, in order
 * @param[in] length Length declared for every contig
 * @return Header options ready for write_tmp_vcf() or parse_records()
 */
vcf_opts make_vcf_opts(int callset, const std::vector<std::string> & contigs, int length) {
    vcf_opts opts;
    opts.sample = callset_strs[callset];
    opts.contigs.clear();
    for (const std::string & ctg : contigs) {
        opts.contigs.push_back("##contig=<ID=" + ctg + ",length=" + std::to_string(length) + ">");
    }
    return opts;
}

/* Record lines ***********************************************************************************/

/**
 * @brief Joins the columns of one VCF data line with tabs.
 * @param[in] rec Column values to join
 * @return One tab-separated VCF data line, without a trailing newline
 */
std::string vcf_line(const vcf_record & rec) {
    return rec.ctg + "\t" + std::to_string(rec.pos) + "\t" + rec.id + "\t" + rec.ref + "\t" +
        rec.alt + "\t" + rec.qual + "\t" + rec.filter + "\t" + rec.info + "\t" + rec.format +
        "\t" + rec.sample;
}

/**
 * @brief Builds a data line with the given genotype, carrying QUAL 50, PASS, and phase set 1.
 * @param[in] pos 1-based VCF position
 * @param[in] ref REF allele
 * @param[in] alt ALT allele
 * @param[in] gt GT field value (e.g. "1|0", "1|.", ".|.")
 * @param[in] ctg Contig name
 * @return One tab-separated VCF data line, without a trailing newline
 */
std::string record(int pos, const std::string & ref, const std::string & alt,
        const std::string & gt, const std::string & ctg) {
    vcf_record rec;
    rec.ctg = ctg;
    rec.pos = pos;
    rec.ref = ref;
    rec.alt = alt;
    rec.sample = gt + ":1";
    return vcf_line(rec);
}

/**
 * @brief Builds a data line carrying the given FORMAT keys and sample values.
 * @param[in] pos 1-based VCF position
 * @param[in] ref REF allele
 * @param[in] alt ALT allele
 * @param[in] format FORMAT column, such as "GT" or "GT:GQ"
 * @param[in] sample Sample column, matching format field for field
 * @param[in] ctg Contig name
 * @return One tab-separated VCF data line, without a trailing newline
 */
std::string fmt_record(int pos, const std::string & ref, const std::string & alt,
        const std::string & format, const std::string & sample, const std::string & ctg) {
    vcf_record rec;
    rec.ctg = ctg;
    rec.pos = pos;
    rec.ref = ref;
    rec.alt = alt;
    rec.format = format;
    rec.sample = sample;
    return vcf_line(rec);
}

/**
 * @brief Builds a phased 1|0 SNP with the given QUAL and FILTER columns.
 * @param[in] pos 1-based VCF position
 * @param[in] qual QUAL column, such as "50" or "." for no reported quality
 * @param[in] filter FILTER column, such as "PASS", "LOWQ", or "." for no filters
 * @param[in] ctg Contig name
 * @return One tab-separated VCF data line, without a trailing newline
 */
std::string qual_filter_record(int pos, const std::string & qual, const std::string & filter,
        const std::string & ctg) {
    vcf_record rec;
    rec.ctg = ctg;
    rec.pos = pos;
    rec.qual = qual;
    rec.filter = filter;
    return vcf_line(rec);
}

/* Parse capture **********************************************************************************/

/**
 * @brief Redirects stderr to fn, truncating any existing contents.
 * @param[in] fn Path of the redirect target
 * @throws std::runtime_error if stderr cannot be duplicated or the target cannot be opened
 */
StderrToFile::StderrToFile(const std::string & fn)
        : saved_fd(dup(fileno(stderr))),
          file_fd(open(fn.data(), O_WRONLY | O_CREAT | O_TRUNC, 0644)) {
    // throw rather than redirect nowhere: a silent failure would empty the captured log and
    // fail every assertion on it, hiding the real cause behind unrelated mismatches
    if (this->saved_fd < 0 || this->file_fd < 0) {
        if (this->saved_fd >= 0) close(this->saved_fd);
        if (this->file_fd >= 0) close(this->file_fd);
        throw std::runtime_error("StderrToFile: could not redirect stderr to " + fn);
    }
    std::fflush(stderr);
    dup2(this->file_fd, fileno(stderr));
}

/**
 * @brief Flushes the redirected output and restores the original stderr.
 */
StderrToFile::~StderrToFile() {
    std::fflush(stderr);
    dup2(this->saved_fd, fileno(stderr));
    close(this->saved_fd);
    close(this->file_fd);
}

/**
 * @brief Parses VCF records with parse_variants(), capturing its stderr and output VCF.
 * @param[in] dir Temporary directory owning the fixture and captured output
 * @param[in] records VCF data lines, without trailing newlines
 * @return Surviving variants, captured log output, and the VCF written from those variants
 */
ParseResult parse_records(const TempDir & dir, const std::vector<std::string> & records) {
    return parse_records(dir, records, make_vcf_opts(), make_fasta("chr1", std::string(1000, 'A')));
}

/**
 * @brief Parses VCF records under a caller-supplied header, capturing stderr and the output VCF.
 * @param[in] dir Temporary directory owning the fixture and captured output
 * @param[in] records VCF data lines, without trailing newlines
 * @param[in] opts Header lines and sample name to write
 * @param[in] ref Reference sequence data, may be nullptr
 * @return Surviving variants, captured log output, and the VCF written from those variants
 * @note The written VCF stands in for summary.vcf: both are generated from the variants that
 *       survive parse-time filtering, so a variant absent here is absent from summary.vcf.
 */
ParseResult parse_records(const TempDir & dir, const std::vector<std::string> & records,
        const vcf_opts & opts, std::shared_ptr<fastaData> ref) {
    const std::string vcf_fn = write_tmp_vcf(dir, records, opts);
    const std::string log_fn = dir.path("parse.log");
    const std::string out_fn = dir.path("out.vcf");

    ParseResult result;
    result.vars = std::make_shared<variantData>();

    { // stderr is redirected for the parse alone, so the INFO/WARN summary can be asserted on
        StderrToFile redirect(log_fn);
        parse_variants(vcf_fn, result.vars, ref, QUERY);
    }

    result.log = read_text(log_fn);
    result.vars->write_vcf(out_fn);
    result.out_vcf = read_text(out_fn);
    return result;
}

/**
 * @brief Parses records under the given header without redirecting stderr.
 * @param[in] dir Temporary directory owning the written fixture
 * @param[in] records VCF data lines, without trailing newlines
 * @param[in] opts Header lines and sample name to write
 * @param[in] callset QUERY or TRUTH callset identifier
 */
void parse_unredirected(const TempDir & dir, const std::vector<std::string> & records,
        const vcf_opts & opts, int callset) {
    const std::string vcf_fn = write_tmp_vcf(dir, records, opts);
    std::shared_ptr<variantData> vars(new variantData());
    parse_variants(vcf_fn, vars, nullptr, callset);
}

/**
 * @brief Reports whether a captured log contains a substring.
 * @param[in] log Captured INFO/WARN output
 * @param[in] text Substring to search for
 * @return True if the log contains the substring
 */
bool logged(const std::string & log, const std::string & text) {
    return log.find(text) != std::string::npos;
}

/**
 * @brief Reports whether the captured log contains a substring.
 * @param[in] r Result of parse_records()
 * @param[in] text Substring to search for
 * @return True if the log contains the substring
 */
bool logged(const ParseResult & r, const std::string & text) {
    return logged(r.log, text);
}

/**
 * @brief Returns the variants that survived parsing on one haplotype of a contig.
 * @param[in] r Result of parse_records()
 * @param[in] hap Haplotype index (HAP1 or HAP2)
 * @param[in] ctg Contig name
 * @return Variant container for that haplotype, or nullptr if the contig is absent
 */
std::shared_ptr<ctgVariants> hap_vars(const ParseResult & r, int hap, const std::string & ctg) {
    const auto & ctg_vars = r.vars->variants[hap];
    const auto found = ctg_vars.find(ctg);
    return (found == ctg_vars.end()) ? nullptr : found->second;
}

/**
 * @brief Counts variants that survived parsing on one haplotype of a contig.
 * @param[in] r Result of parse_records()
 * @param[in] hap Haplotype index (HAP1 or HAP2)
 * @param[in] ctg Contig name
 * @return Number of surviving variants, or 0 if the contig is absent
 */
int kept_on_hap(const ParseResult & r, int hap, const std::string & ctg) {
    std::shared_ptr<ctgVariants> vars = hap_vars(r, hap, ctg);
    return (vars == nullptr) ? 0 : vars->n;
}

/**
 * @brief Counts variants that survived parsing across both haplotypes of a contig.
 * @param[in] r Result of parse_records()
 * @param[in] ctg Contig name
 * @return Total number of surviving variants
 */
int total_kept(const ParseResult & r, const std::string & ctg) {
    return kept_on_hap(r, HAP1, ctg) + kept_on_hap(r, HAP2, ctg);
}

/**
 * @brief Reports whether the written VCF holds a record at a position on a contig.
 * @param[in] r Result of parse_records()
 * @param[in] pos 1-based VCF position
 * @param[in] ctg Contig name
 * @return True if a data line at that position was written
 */
bool wrote_pos(const ParseResult & r, int pos, const std::string & ctg) {
    return count_pos(r, pos, ctg) > 0;
}

/**
 * @brief Counts the records the written VCF holds at a position on a contig.
 * @param[in] r Result of parse_records()
 * @param[in] pos 1-based VCF position
 * @param[in] ctg Contig name
 * @return Number of data lines written at that position
 */
size_t count_pos(const ParseResult & r, int pos, const std::string & ctg) {
    // the leading newline anchors the match to the start of a line, so a position never matches
    // inside another column and the header lines above the records cannot match at all
    const std::string line = "\n" + ctg + "\t" + std::to_string(pos) + "\t";
    size_t count = 0;
    for (size_t at = r.out_vcf.find(line); at != std::string::npos;
            at = r.out_vcf.find(line, at + 1)) {
        count++;
    }
    return count;
}

/**
 * @brief Returns the genotype-histogram line parse_variants() prints for a genotype and count.
 * @param[in] gt Genotype code (GT_*)
 * @param[in] count Number of records tallied under that genotype
 * @return Substring of the INFO line, derived from gt_strs so it tracks renames
 */
std::string gt_hist_line(uint8_t gt, int count) {
    std::string name = gt_strs[gt];
    while (name.size() < 3) name = " " + name; // the "%3s" in the INFO format right-justifies
    return "    " + name + ": " + std::to_string(count);
}

/**
 * @brief Returns the variant-type line parse_variants() prints for a type and count.
 * @param[in] type Variant type (TYPE_*)
 * @param[in] count Number of alleles tallied under that type across both haplotypes
 * @return Substring of the INFO line, derived from type_strs so it tracks renames
 */
std::string type_hist_line(uint8_t type, int count) {
    return "    " + type_strs[type] + ": " + std::to_string(count);
}

/* In-memory builders *****************************************************************************/

/**
 * @brief Builds a single-contig reference, bypassing the FASTA-parsing constructor.
 * @param[in] ctg Contig name
 * @param[in] seq Contig sequence
 * @param[in] uppercase Whether to uppercase the sequence, as the FASTA constructor does
 * @return Reference containing exactly this contig
 */
std::shared_ptr<fastaData> make_fasta(const std::string & ctg, const std::string & seq,
        bool uppercase) {
    return make_fasta({{ctg, seq}}, uppercase);
}

/**
 * @brief Builds a multi-contig reference from (contig, sequence) pairs, in order.
 * @param[in] seqs Contig name and sequence pairs
 * @param[in] uppercase Whether to uppercase each sequence, as the FASTA constructor does
 * @return Reference containing all given contigs
 */
std::shared_ptr<fastaData> make_fasta(
        const std::vector< std::pair<std::string, std::string> > & seqs, bool uppercase) {
    std::shared_ptr<fastaData> ref(new fastaData());
    for (const auto & [ctg, seq] : seqs) {
        std::string stored = seq;
        if (uppercase) {
            std::transform(stored.begin(), stored.end(), stored.begin(), ::toupper);
        }
        ref->fasta[ctg] = stored;
        ref->lengths[ctg] = int(stored.size());
    }
    return ref;
}

/**
 * @brief Builds a bedData holding the given [start, stop) regions on one contig, in order.
 * @param[in] ctg Contig name
 * @param[in] regions Start and stop coordinate pairs, appended in the order given
 * @return Populated BED interval container
 */
bedData make_bed(const std::string & ctg, const std::vector< std::pair<int, int> > & regions) {
    return make_bed({{ctg, regions}});
}

/**
 * @brief Builds a multi-contig bedData from (contig, regions) pairs, in order.
 *
 * add() is used rather than assigning the fields directly, so that the contig list, region
 * vectors, counts, and total size stay consistent with each other.
 * @param[in] regions Contig name and region-list pairs, appended in the order given
 * @return Populated BED interval container
 */
bedData make_bed(const std::vector< std::pair<std::string,
        std::vector< std::pair<int, int> > > > & regions) {
    bedData bed;
    for (const auto & [ctg, ctg_regions] : regions) {
        for (const auto & [start, stop] : ctg_regions) bed.add(ctg, start, stop);
    }
    return bed;
}

/**
 * @brief Builds a variantData over the given contigs, each holding an empty ctgVariants per hap.
 * @param[in] callset QUERY or TRUTH callset identifier
 * @param[in] contigs Contig names
 * @param[in] lengths Contig lengths, parallel to contigs
 * @param[in] ploidy Contig ploidies, parallel to contigs
 * @return Variant container with no variants on any contig
 * @throws ERROR if the parallel vectors have differing lengths
 */
std::shared_ptr<variantData> make_variantData(int callset,
        const std::vector<std::string> & contigs, const std::vector<int> & lengths,
        const std::vector<int> & ploidy) {
    if (contigs.size() != lengths.size() || contigs.size() != ploidy.size()) {
        ERROR("make_variantData() requires parallel vectors of equal length");
    }
    std::shared_ptr<variantData> vars(new variantData());
    vars->callset = callset;
    vars->filename = callset_strs[callset] + ".vcf";
    vars->sample = callset_strs[callset];
    vars->contigs = contigs;
    vars->lengths = lengths;
    vars->ploidy = ploidy;
    for (const std::string & ctg : contigs) {
        vars->variants[HAP1][ctg] = make_ctgVariants(ctg, {});
        vars->variants[HAP2][ctg] = make_ctgVariants(ctg, {});
    }
    return vars;
}

/**
 * @brief Builds a ctgVariants container holding the described variants, in the given order.
 *
 * Every add_var() parameter through supercluster is named explicitly, because they are all
 * defaulted past phase_set and omitting one silently binds supercluster to the wrong parameter.
 * @param[in] ctg Contig name
 * @param[in] vars Variants to append, in ascending position order
 * @return Populated variant container
 */
std::shared_ptr<ctgVariants> make_ctgVariants(const std::string & ctg,
        const std::vector<var_desc> & vars) {
    std::shared_ptr<ctgVariants> ctg_vars(new ctgVariants(ctg));
    for (const var_desc & var : vars) {
        ctg_vars->add_var(var.pos, var.rlen, var.type, var.loc, var.ref, var.alt, var.gt,
                var.qual, var.qual, var.phase_set, var.rec_idx, var.alt_idx, var.ploidy,
                var.supercluster);
    }
    return ctg_vars;
}

/**
 * @brief Builds a one-variant container with the given original and calculated genotypes.
 * @param[in] orig_gt Original genotype (GT_*) stored in orig_gts[0]
 * @param[in] calc_gt Calculated genotype (GT_*) stored in calc_gts[0]
 * @param[in] ctg Contig name
 * @param[in] pos 0-based reference start position
 * @return Container holding a single A>C substitution with the requested genotypes
 */
std::shared_ptr<ctgVariants> make_gt_var(uint8_t orig_gt, uint8_t calc_gt, const std::string & ctg,
        int pos) {
    var_desc var;
    var.pos = pos;
    var.rlen = 1;
    var.type = TYPE_SUB;
    var.ref = "A";
    var.alt = "C";
    var.gt = orig_gt;
    std::shared_ptr<ctgVariants> vars = make_ctgVariants(ctg, {var});
    vars->calc_gts[0] = calc_gt;
    return vars;
}

/**
 * @brief Builds a one-variant container of the given type with the given allele sequences.
 * @param[in] type Variant type (TYPE_*)
 * @param[in] ref Reference allele sequence, whose length becomes the variant's rlen
 * @param[in] alt Alternate allele sequence
 * @param[in] ctg Contig name
 * @param[in] pos 0-based reference start position
 * @return Container holding a single variant
 */
std::shared_ptr<ctgVariants> make_typed_var(uint8_t type, const std::string & ref,
        const std::string & alt, const std::string & ctg, int pos) {
    var_desc var;
    var.pos = pos;
    var.rlen = int(ref.size());
    var.type = type;
    var.ref = ref;
    var.alt = alt;
    return make_ctgVariants(ctg, {var});
}

/**
 * @brief Sets all six per-haplotype evaluation lanes for one variant.
 * @param[in,out] vars Variant container to modify
 * @param[in] hap Haplotype index (HAP1 or HAP2)
 * @param[in] idx Variant index
 * @param[in] errtype Error type (ERRTYPE_*)
 * @param[in] sync_group Sync group index
 * @param[in] callq Call quality
 * @param[in] ref_ed Reference edit distance
 * @param[in] query_ed Query edit distance
 * @param[in] credit Partial credit
 */
void set_hap_data(std::shared_ptr<ctgVariants> vars, int hap, int idx, uint8_t errtype,
        int sync_group, float callq, int ref_ed, int query_ed, float credit) {
    vars->errtypes[hap][idx] = errtype;
    vars->sync_group[hap][idx] = sync_group;
    vars->callq[hap][idx] = callq;
    vars->ref_ed[hap][idx] = ref_ed;
    vars->query_ed[hap][idx] = query_ed;
    vars->credit[hap][idx] = credit;
}

/**
 * @brief Sets cluster boundaries and reaches; nc defaults to clusters.size()-1.
 *
 * wf_swg_cluster() stores nc clusters as nc+1 boundaries with one reach pair per cluster, whereas
 * load_and_merge_callset_vars_across_haps() counts its trailing sentinel in nc; pass nc
 * explicitly to build the merged form.
 * @param[in,out] vars Variant container to modify
 * @param[in] clusters Variant index at the start of each cluster, plus a trailing sentinel
 * @param[in] left_reaches Leftmost reach of each cluster
 * @param[in] right_reaches Rightmost reach of each cluster
 * @param[in] nc Cluster count, or -1 to use clusters.size()-1
 */
void set_clusters(std::shared_ptr<ctgVariants> vars, const std::vector<int> & clusters,
        const std::vector<int> & left_reaches, const std::vector<int> & right_reaches, int nc) {
    vars->clusters = clusters;
    vars->left_reaches = left_reaches;
    vars->right_reaches = right_reaches;
    vars->nc = (nc < 0) ? int(clusters.size()) - 1 : nc;
}

/**
 * @brief Builds a ctgSuperclusters holding the given query and truth variant containers.
 * @param[in] qvars Query variants
 * @param[in] tvars Truth variants
 * @return Supercluster container for one contig
 */
std::shared_ptr<ctgSuperclusters> make_ctgSuperclusters(std::shared_ptr<ctgVariants> qvars,
        std::shared_ptr<ctgVariants> tvars) {
    std::shared_ptr<ctgSuperclusters> sc(new ctgSuperclusters());
    sc->callset_vars[QUERY] = qvars;
    sc->callset_vars[TRUTH] = tvars;
    return sc;
}

/**
 * @brief Builds a superclusterData over the given contigs, bypassing clustering.
 *
 * The real constructor is invoked with empty callsets so that no clustering runs, then the contig
 * and supercluster fields are overwritten. Sample and filename fields are set to QUERY and TRUTH
 * defaults.
 * @param[in] contigs Contig names
 * @param[in] lengths Contig lengths, parallel to contigs
 * @param[in] ploidy Contig ploidies, parallel to contigs
 * @param[in] superclusters Per-contig supercluster containers, parallel to contigs
 * @param[in] ref Reference sequence data, may be nullptr
 * @return Populated supercluster data
 * @throws ERROR if the parallel vectors have differing lengths
 */
std::shared_ptr<superclusterData> make_superclusterData(
        const std::vector<std::string> & contigs, const std::vector<int> & lengths,
        const std::vector<int> & ploidy,
        const std::vector< std::shared_ptr<ctgSuperclusters> > & superclusters,
        std::shared_ptr<fastaData> ref) {
    if (contigs.size() != lengths.size() || contigs.size() != ploidy.size() ||
            contigs.size() != superclusters.size()) {
        ERROR("make_superclusterData() requires parallel vectors of equal length");
    }
    std::shared_ptr<variantData> empty_query(new variantData());
    std::shared_ptr<variantData> empty_truth(new variantData());
    std::shared_ptr<superclusterData> sc_data(
            new superclusterData(empty_query, empty_truth, ref));
    sc_data->samples = {"QUERY", "TRUTH"};
    sc_data->filenames = {"query.vcf", "truth.vcf"};
    sc_data->contigs = contigs;
    sc_data->lengths = lengths;
    sc_data->ploidy = ploidy;
    for (size_t i = 0; i < contigs.size(); i++) {
        sc_data->superclusters[contigs[i]] = superclusters[i];
    }
    return sc_data;
}

/**
 * @brief Builds an alignment graph for one supercluster and truth haplotype.
 * @param[in] sc Supercluster container holding the query and truth variants
 * @param[in] ref Reference sequence data
 * @param[in] ctg Contig name
 * @param[in] truth_hap Truth haplotype index (HAP1 or HAP2)
 * @param[in] sc_idx Supercluster index within the contig
 * @return Graph ready for calc_prec_recall_aln()
 */
std::shared_ptr<Graph> make_graph(std::shared_ptr<ctgSuperclusters> sc,
        std::shared_ptr<fastaData> ref, const std::string & ctg, int truth_hap, int sc_idx) {
    return std::shared_ptr<Graph>(new Graph(sc, sc_idx, ref, ctg, truth_hap));
}

/**
 * @brief Allocates the offsets buffer that wf_swg_max_reach() requires from its caller.
 * @param[in] qlen Query sequence length
 * @param[in] tlen Truth sequence length
 * @param[in] x Substitution penalty
 * @param[in] o Gap-open penalty
 * @param[in] e Gap-extension penalty
 * @return Buffer of MATS*(max(x,o+e)+1)*(qlen+tlen-1) offsets, each initialized to -2, or an
 *   empty buffer if both sequences are empty
 */
std::vector<int> alloc_reach_offs(int qlen, int tlen, int x, int o, int e) {
    // clamp so two empty sequences give a diagonal count of 0 rather than a huge size_t
    size_t mat_len = std::max(1, qlen + tlen) - 1;
    size_t offs_size = size_t(MATS) * size_t(std::max(x, o+e) + 1) * mat_len;
    return std::vector<int>(offs_size, -2);
}
