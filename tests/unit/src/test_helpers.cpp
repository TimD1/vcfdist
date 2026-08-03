/**
 * @file test_helpers.cpp
 * @brief Shared unit-test scaffolding: global-state fixture, temporary files, in-memory builders.
 */
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>

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
    for (const std::string & record : records) out << record << "\n";
    out.close();
    return vcf_fn;
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
 * @return Buffer of MATS*(max(x,o+e)+1)*(qlen+tlen-1) offsets, each initialized to -2
 */
std::vector<int> alloc_reach_offs(int qlen, int tlen, int x, int o, int e) {
    size_t offs_size = size_t(MATS) * size_t(std::max(x, o+e) + 1) * size_t(qlen + tlen - 1);
    return std::vector<int>(offs_size, -2);
}
