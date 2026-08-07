/**
 * @file test_print.cpp
 * @brief Unit tests for print.cpp: qscore, get_ptr_repr, color wrappers, metrics, write_params.
 */
#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"

#include "../../../src/dist.h"
#include "../../../src/globals.h"
#include "../../../src/print.h"
#include "test_helpers.h"

namespace {

/**
 * @brief Reads an entire file into a string.
 * @param[in] fn Path of the file to read
 * @return File contents, or an empty string if the file cannot be opened
 */
std::string read_file(const std::string & fn) {
    std::ifstream in(fn);
    if (!in.is_open()) return "";
    std::ostringstream ostr;
    ostr << in.rdbuf();
    return ostr.str();
}

/**
 * @brief Builds a single-entry pointer map for get_ptr_repr.
 * @param[in] cell Cell to look up
 * @param[in] prev Predecessor cell the map returns
 * @return Map containing exactly the one cell-to-predecessor entry
 */
std::unordered_map<idx4,idx4> one_ptr(const idx4 & cell, const idx4 & prev) {
    std::unordered_map<idx4,idx4> ptrs;
    ptrs[cell] = prev;
    return ptrs;
}

/* qscore *****************************************************************************************/

TEST(Qscore, Prob1) {
    EXPECT_FLOAT_EQ(0.0f, qscore(1.0));
}

TEST(Qscore, Prob0p1) {
    EXPECT_FLOAT_EQ(10.0f, qscore(0.1));
}

TEST(Qscore, Prob0p01) {
    EXPECT_FLOAT_EQ(20.0f, qscore(0.01));
}

TEST(Qscore, Prob0p001) {
    EXPECT_FLOAT_EQ(30.0f, qscore(0.001));
}

TEST(Qscore, ClampHigh) {
    // -10*log10(1e-11) is 110, above the upper clamp; the cap is 100, not the Phred-60 convention
    EXPECT_FLOAT_EQ(100.0f, qscore(1e-11));
}

TEST(Qscore, Prob0) {
    // log10(0) is -inf, so the score is +inf before the upper clamp reduces it to 100
    EXPECT_FLOAT_EQ(100.0f, qscore(0.0));
}

TEST(Qscore, NegativeInput) {
    // log10 of a negative is NaN; std::max(0.0, NaN) returns its first argument, so the result is 0
    float q = qscore(-1.0);
    EXPECT_FALSE(std::isnan(q));
    EXPECT_FLOAT_EQ(0.0f, q);
}

TEST(Qscore, ProbGt1) {
    // -10*log10(2) is negative, so the lower clamp returns 0
    EXPECT_FLOAT_EQ(0.0f, qscore(2.0));
}

TEST(Qscore, Rounding) {
    // the return type is float, so fractional scores are preserved rather than rounded to integers
    float q = qscore(0.5);
    EXPECT_NEAR(3.0103f, q, 1e-4);
    EXPECT_NE(0.0f, q - std::floor(q));
}

/* get_ptr_repr ***********************************************************************************/

TEST(GetPtrRepr, NotFound) {
    idx4 cell(0, 0, 1, 1);
    std::unordered_map<idx4,idx4> ptrs;
    EXPECT_EQ("  .", get_ptr_repr(cell, ptrs));
}

TEST(GetPtrRepr, Up) {
    idx4 cell(0, 0, 3, 5);
    idx4 prev(0, 0, 2, 5);
    EXPECT_EQ("  |", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, Left) {
    idx4 cell(0, 0, 3, 5);
    idx4 prev(0, 0, 3, 4);
    EXPECT_EQ("  _", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, Diag) {
    idx4 cell(0, 0, 3, 5);
    idx4 prev(0, 0, 2, 4);
    EXPECT_EQ("  \\", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, InvalidSameMatrix) {
    // same node pair, but the predecessor is neither adjacent nor diagonal
    idx4 cell(0, 0, 5, 5);
    idx4 prev(0, 0, 1, 1);
    EXPECT_EQ(" ?1", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, QueryNode) {
    // entering query node 8 from query node 7 prints the predecessor id, zero-padded to 2 digits
    idx4 cell(8, 0, 0, 4);
    idx4 prev(7, 0, 2, 4);
    EXPECT_EQ("^07", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, TruthNode) {
    // entering truth node 8 from truth node 3 prints the predecessor id, zero-padded to 2 digits
    idx4 cell(0, 8, 3, 0);
    idx4 prev(0, 3, 3, 2);
    EXPECT_EQ("<03", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, InvalidOther) {
    // the query node changed, but not at the start of the node, so no node-id branch applies
    idx4 cell(1, 0, 5, 0);
    idx4 prev(0, 0, 4, 0);
    EXPECT_EQ(" ?2", get_ptr_repr(cell, one_ptr(cell, prev)));
}

TEST(GetPtrRepr, WidthConsistency) {
    // every branch must return exactly 3 characters for node ids below 100 so the grid stays
    // aligned, including the two invalid markers (issue #73)
    std::vector< std::pair<idx4, idx4> > cases = {
        {idx4(0, 0, 3, 5), idx4(0, 0, 2, 5)}, // up
        {idx4(0, 0, 3, 5), idx4(0, 0, 3, 4)}, // left
        {idx4(0, 0, 3, 5), idx4(0, 0, 2, 4)}, // diagonal
        {idx4(0, 0, 5, 5), idx4(0, 0, 1, 1)}, // invalid, same matrix
        {idx4(8, 0, 0, 4), idx4(7, 0, 2, 4)}, // new query node, two digits
        {idx4(1, 0, 0, 4), idx4(0, 0, 2, 4)}, // new query node, one digit
        {idx4(0, 8, 3, 0), idx4(0, 3, 3, 2)}, // new truth node, two digits
        {idx4(0, 1, 3, 0), idx4(0, 0, 3, 2)}, // new truth node, one digit
        {idx4(1, 0, 5, 0), idx4(0, 0, 4, 0)}, // invalid, different matrix
    };
    for (size_t i = 0; i < cases.size(); i++) {
        std::string repr = get_ptr_repr(cases[i].first, one_ptr(cases[i].first, cases[i].second));
        EXPECT_EQ(size_t(3), repr.size()) << "case " << i << ": '" << repr << "'";
    }

    // the not-found branch takes no pointer map entry
    std::unordered_map<idx4,idx4> empty;
    EXPECT_EQ(size_t(3), get_ptr_repr(idx4(0, 0, 1, 1), empty).size());
}

TEST(GetPtrRepr, QueryNodeWide) {
    // std::setw(2) pads but never truncates, so a node id of 100 or more widens the cell to 4
    // characters and misaligns the debug grid; this pins current output rather than correctness
    idx4 cell(124, 0, 0, 4);
    idx4 prev(123, 0, 2, 4);
    std::string repr = get_ptr_repr(cell, one_ptr(cell, prev));
    EXPECT_EQ("^123", repr);
    EXPECT_EQ(size_t(4), repr.size());
}

/* color wrappers *********************************************************************************/

TEST(Color, Green) {
    // the wrappers emit escape codes unconditionally; the isatty check lives in the COLOR_* macros
    EXPECT_EQ("\033[32m7\033[0m", GREEN(7));
}

TEST(Color, Red) {
    EXPECT_EQ("\033[31mA\033[0m", RED('A'));
}

TEST(Color, Blue) {
    EXPECT_EQ("\033[34mACGT\033[0m", BLUE(std::string("ACGT")));
}

TEST(Color, Yellow) {
    EXPECT_EQ("\033[33m-3\033[0m", YELLOW(-3));
}

TEST(Color, Purple) {
    EXPECT_EQ("\033[35mchr1\033[0m", PURPLE(std::string("chr1")));
}

TEST(Color, GreenStrEmpty) {
    // an empty payload still yields the prefix and reset suffix
    EXPECT_EQ("\033[32m\033[0m", GREEN(std::string("")));
}

/* tally_counts_by_qual ***************************************************************************/

/**
 * @brief Wraps query and truth variants of one contig in a phaseblockData.
 * @param[in] qvars Query variants
 * @param[in] tvars Truth variants
 * @return Phase block data holding the one contig, ready for tally_counts_by_qual()
 */
std::unique_ptr<phaseblockData> one_ctg(std::shared_ptr<ctgVariants> qvars,
        std::shared_ptr<ctgVariants> tvars) {
    return make_phaseblockData({"chr1"}, {1000}, {make_ctgSuperclusters(qvars, tvars)});
}

/**
 * @brief Builds a container holding one variant present on haplotype 1 only.
 * @param[in] type Variant type (TYPE_*)
 * @param[in] ref Reference allele sequence
 * @param[in] alt Alternate allele sequence
 * @param[in] ctg Contig name
 * @return Container whose single variant carries genotype 1|0 and error type ERRTYPE_UN
 */
std::shared_ptr<ctgVariants> hap1_var(edittype_t type, const std::string & ref,
        const std::string & alt, const std::string & ctg = "chr1") {
    std::shared_ptr<ctgVariants> vars = make_typed_var(type, ref, alt, ctg);
    vars->orig_gts[0] = GT_ALT1_REF;
    return vars;
}

/**
 * @brief Returns the counts of one variant type and error type across the whole quality sweep.
 * @param[in] counts Counters returned by tally_counts_by_qual()
 * @param[in] callset QUERY or TRUTH
 * @param[in] vartype Variant size class (VARTYPE_*)
 * @param[in] errtype Error type (ERRTYPE_*)
 * @return One count per quality threshold, in ascending threshold order
 */
std::vector<float> sweep(const pr_counts & counts, callset_t callset, sizeclass_t vartype,
        errtype_t errtype) {
    return callset == QUERY ? counts.query[vartype][errtype] : counts.truth[vartype][errtype];
}

TEST(TallyCountsByQual, SnpTpSingle) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = hap1_var(TYPE_SUB, "A", "G");
    set_hap_data(qvars, HAP1, 0, ERRTYPE_TP, 0, 0, 0, 0, 0);

    pr_counts counts = tally_counts_by_qual(one_ctg(qvars, make_ctgVariants("chr1", {})), 0, 1);

    EXPECT_FLOAT_EQ(1.0f, counts.query[VARTYPE_SNP][ERRTYPE_TP][0]);
    EXPECT_FLOAT_EQ(1.0f, counts.query[VARTYPE_ALL][ERRTYPE_TP][0]);
    EXPECT_FLOAT_EQ(0.0f, counts.query[VARTYPE_INDEL][ERRTYPE_TP][0]);
    EXPECT_FLOAT_EQ(0.0f, counts.query[VARTYPE_SV][ERRTYPE_TP][0]);
    EXPECT_FLOAT_EQ(0.0f, counts.query[VARTYPE_SNP][ERRTYPE_FP][0]);
    EXPECT_FLOAT_EQ(0.0f, counts.truth[VARTYPE_ALL][ERRTYPE_TP][0]);
}

TEST(TallyCountsByQual, IndelFpSingle) {
    GlobalsGuard guard;
    // a 3bp insertion is below the 50bp SV threshold, so it is an INDEL rather than an SV
    std::shared_ptr<ctgVariants> qvars = hap1_var(TYPE_INS, "A", "ACGT");
    set_hap_data(qvars, HAP1, 0, ERRTYPE_FP, 0, 0, 0, 0, 0);

    pr_counts counts = tally_counts_by_qual(one_ctg(qvars, make_ctgVariants("chr1", {})), 0, 1);

    EXPECT_FLOAT_EQ(1.0f, counts.query[VARTYPE_INDEL][ERRTYPE_FP][0]);
    EXPECT_FLOAT_EQ(1.0f, counts.query[VARTYPE_ALL][ERRTYPE_FP][0]);
    EXPECT_FLOAT_EQ(0.0f, counts.query[VARTYPE_SNP][ERRTYPE_FP][0]);
    EXPECT_FLOAT_EQ(0.0f, counts.query[VARTYPE_SV][ERRTYPE_FP][0]);
    EXPECT_FLOAT_EQ(0.0f, counts.query[VARTYPE_INDEL][ERRTYPE_TP][0]);
}

TEST(TallyCountsByQual, SvFnSingle) {
    GlobalsGuard guard;
    // a 60bp deletion is at or above the 50bp SV threshold; FN is a truth-side classification
    std::shared_ptr<ctgVariants> tvars = hap1_var(TYPE_DEL, std::string(60, 'A'), "A");
    set_hap_data(tvars, HAP1, 0, ERRTYPE_FN, 0, 1, 0, 0, 0);

    pr_counts counts = tally_counts_by_qual(one_ctg(make_ctgVariants("chr1", {}), tvars), 0, 1);

    EXPECT_FLOAT_EQ(1.0f, counts.truth[VARTYPE_SV][ERRTYPE_FN][0]);
    EXPECT_FLOAT_EQ(1.0f, counts.truth[VARTYPE_ALL][ERRTYPE_FN][0]);
    EXPECT_FLOAT_EQ(0.0f, counts.truth[VARTYPE_SNP][ERRTYPE_FN][0]);
    EXPECT_FLOAT_EQ(0.0f, counts.truth[VARTYPE_INDEL][ERRTYPE_FN][0]);
    EXPECT_FLOAT_EQ(0.0f, counts.query[VARTYPE_ALL][ERRTYPE_FN][0]);
}

TEST(TallyCountsByQual, QualSweepMonotone) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars = hap1_var(TYPE_SUB, "A", "G");
    set_hap_data(qvars, HAP1, 0, ERRTYPE_TP, 0, 3, 0, 0, 0);

    pr_counts counts = tally_counts_by_qual(one_ctg(qvars, make_ctgVariants("chr1", {})), 0, 5);

    // the variant passes every threshold up to and including its own Qscore, and none above it
    std::vector<float> expected = {1, 1, 1, 1, 0, 0};
    EXPECT_EQ(expected, sweep(counts, QUERY, VARTYPE_SNP, ERRTYPE_TP));
    EXPECT_EQ(expected, sweep(counts, QUERY, VARTYPE_ALL, ERRTYPE_TP));
}

TEST(TallyCountsByQual, TruthFnAboveQscore) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> tvars = hap1_var(TYPE_SUB, "A", "G");
    set_hap_data(tvars, HAP1, 0, ERRTYPE_TP, 0, 2, 0, 0, 0);

    pr_counts counts = tally_counts_by_qual(one_ctg(make_ctgVariants("chr1", {}), tvars), 0, 5);

    // the matching query call is filtered out above its own Qscore, turning the truth call into a
    // false negative rather than dropping it
    std::vector<float> expected_tp = {1, 1, 1, 0, 0, 0};
    std::vector<float> expected_fn = {0, 0, 0, 1, 1, 1};
    EXPECT_EQ(expected_tp, sweep(counts, TRUTH, VARTYPE_SNP, ERRTYPE_TP));
    EXPECT_EQ(expected_fn, sweep(counts, TRUTH, VARTYPE_SNP, ERRTYPE_FN));
    EXPECT_EQ(expected_fn, sweep(counts, TRUTH, VARTYPE_ALL, ERRTYPE_FN));
}

TEST(TallyCountsByQual, AcErr2To1DecrementsTruthTp) {
    GlobalsGuard guard;
    // the matched genotype was 1|1 and was forced back to the original allele count, so the
    // extra alternate allele already participated in a truth match and is compensated for here
    std::shared_ptr<ctgVariants> qvars = hap1_var(TYPE_SUB, "A", "G");
    qvars->matched_gts[0] = GT_REF_ALT1;
    qvars->ac_errtype[0] = AC_ERR_2_TO_1;
    set_hap_data(qvars, HAP1, 0, ERRTYPE_UN, 0, 2, 0, 0, 0);
    set_hap_data(qvars, HAP2, 0, ERRTYPE_TP, 0, 2, 0, 0, 0);

    pr_counts counts = tally_counts_by_qual(one_ctg(qvars, make_ctgVariants("chr1", {})), 0, 3);

    // this characterizes current behavior: with no truth call to offset it, the correction drives
    // the truth true-positive count negative
    std::vector<float> expected_tp = {-1, -1, -1, 0};
    std::vector<float> expected_fn = {1, 1, 1, 0};
    EXPECT_EQ(expected_tp, sweep(counts, TRUTH, VARTYPE_SNP, ERRTYPE_TP));
    EXPECT_EQ(expected_tp, sweep(counts, TRUTH, VARTYPE_ALL, ERRTYPE_TP));
    EXPECT_EQ(expected_fn, sweep(counts, TRUTH, VARTYPE_SNP, ERRTYPE_FN));

    // the haplotype that does carry the alternate allele is counted as a query true positive
    std::vector<float> expected_query = {1, 1, 1, 0};
    EXPECT_EQ(expected_query, sweep(counts, QUERY, VARTYPE_SNP, ERRTYPE_TP));
}

TEST(TallyCountsByQual, CalcgtSwappedHaplotype) {
    GlobalsGuard guard;
    // original 0|1 against matched 1|0 is a swap, so haplotype 1 reads the haplotype 0 lane
    std::shared_ptr<ctgVariants> qvars = hap1_var(TYPE_SUB, "A", "G");
    qvars->orig_gts[0] = GT_REF_ALT1;
    qvars->matched_gts[0] = GT_ALT1_REF;
    set_hap_data(qvars, HAP1, 0, ERRTYPE_FP, 0, 1, 0, 0, 0);
    set_hap_data(qvars, HAP2, 0, ERRTYPE_TP, 0, 4, 0, 0, 0);

    pr_counts counts = tally_counts_by_qual(one_ctg(qvars, make_ctgVariants("chr1", {})), 0, 5);

    // without the remap this would instead be a true positive counted up to threshold 4
    std::vector<float> expected_fp = {1, 1, 0, 0, 0, 0};
    std::vector<float> expected_tp = {0, 0, 0, 0, 0, 0};
    EXPECT_EQ(expected_fp, sweep(counts, QUERY, VARTYPE_SNP, ERRTYPE_FP));
    EXPECT_EQ(expected_tp, sweep(counts, QUERY, VARTYPE_SNP, ERRTYPE_TP));
}

TEST(TallyCountsByQual, ErrtypeUnknownWarnsAndSkips) {
    GlobalsGuard guard;
    TempDir dir;
    std::shared_ptr<ctgVariants> qvars = hap1_var(TYPE_SUB, "A", "G");
    std::shared_ptr<ctgVariants> tvars = hap1_var(TYPE_SUB, "A", "G");
    set_hap_data(qvars, HAP1, 0, ERRTYPE_UN, 0, 3, 0, 0, 0);
    set_hap_data(tvars, HAP1, 0, ERRTYPE_UN, 0, 3, 0, 0, 0);

    pr_counts counts;
    {
        StderrToFile capture(dir.path("warn.log"));
        counts = tally_counts_by_qual(one_ctg(qvars, tvars), 0, 5);
    }

    std::string log = read_text(dir.path("warn.log"));
    EXPECT_NE(std::string::npos, log.find("Unknown error type at QUERY chr1:100")) << log;
    EXPECT_NE(std::string::npos, log.find("Unknown error type at TRUTH chr1:100")) << log;

    // an unevaluated variant contributes nothing at all, not even the truth-side FN tail
    std::vector<float> zeros(6, 0);
    for (sizeclass_t type : EnumRange<sizeclass_t, SIZECLASS_SLOTS>{}) {
        for (errtype_t err : EnumRange<errtype_t, ERRTYPE_SLOTS>{}) {
            EXPECT_EQ(zeros, sweep(counts, QUERY, type, err))
                    << "query " << idx(type) << " " << idx(err);
            EXPECT_EQ(zeros, sweep(counts, TRUTH, type, err))
                    << "truth " << idx(type) << " " << idx(err);
        }
    }
}

TEST(TallyCountsByQual, EmptyContig) {
    GlobalsGuard guard;

    pr_counts counts = tally_counts_by_qual(
            one_ctg(make_ctgVariants("chr1", {}), make_ctgVariants("chr1", {})), 0, 2);

    std::vector<float> zeros(3, 0);
    for (sizeclass_t type : EnumRange<sizeclass_t, SIZECLASS_SLOTS>{}) {
        for (errtype_t err : EnumRange<errtype_t, ERRTYPE_SLOTS>{}) {
            EXPECT_EQ(zeros, sweep(counts, QUERY, type, err))
                    << "query " << idx(type) << " " << idx(err);
            EXPECT_EQ(zeros, sweep(counts, TRUTH, type, err))
                    << "truth " << idx(type) << " " << idx(err);
        }
    }
}

TEST(TallyCountsByQual, MultiContigSums) {
    GlobalsGuard guard;
    std::shared_ptr<ctgVariants> qvars1 = hap1_var(TYPE_SUB, "A", "G", "chr1");
    std::shared_ptr<ctgVariants> qvars2 = hap1_var(TYPE_SUB, "A", "G", "chr2");
    set_hap_data(qvars1, HAP1, 0, ERRTYPE_TP, 0, 2, 0, 0, 0);
    set_hap_data(qvars2, HAP1, 0, ERRTYPE_TP, 0, 2, 0, 0, 0);
    std::unique_ptr<phaseblockData> pb_data = make_phaseblockData({"chr1", "chr2"},
            {1000, 1000}, {make_ctgSuperclusters(qvars1, make_ctgVariants("chr1", {})),
             make_ctgSuperclusters(qvars2, make_ctgVariants("chr2", {}))});

    pr_counts counts = tally_counts_by_qual(pb_data, 0, 3);

    // both contigs accumulate into one set of counters rather than being reported separately
    std::vector<float> expected = {2, 2, 2, 0};
    EXPECT_EQ(expected, sweep(counts, QUERY, VARTYPE_SNP, ERRTYPE_TP));
    EXPECT_EQ(expected, sweep(counts, QUERY, VARTYPE_ALL, ERRTYPE_TP));
}

/* compute_pr_f1 **********************************************************************************/

TEST(ComputePrF1, PrecisionZeroQuery) {
    // no query variants at all is vacuously perfect precision, not a zero-divide
    prec_recall_f1 metrics = compute_pr_f1(0, 0, 3, 1);
    EXPECT_FLOAT_EQ(1.0f, metrics.precision);
    EXPECT_FLOAT_EQ(0.75f, metrics.recall);
    EXPECT_FLOAT_EQ(2*1.0f*0.75f / 1.75f, metrics.f1);
}

TEST(ComputePrF1, RecallZeroTruth) {
    // likewise for an empty truth set, so F1 stays defined when only one callset is empty
    prec_recall_f1 metrics = compute_pr_f1(3, 1, 0, 0);
    EXPECT_FLOAT_EQ(0.75f, metrics.precision);
    EXPECT_FLOAT_EQ(1.0f, metrics.recall);
    EXPECT_FLOAT_EQ(2*0.75f*1.0f / 1.75f, metrics.f1);
}

TEST(ComputePrF1, F1ZeroDenominator) {
    // every query call wrong and every truth call missed, so both terms are 0 and F1 is defined as 0
    prec_recall_f1 metrics = compute_pr_f1(0, 5, 0, 5);
    EXPECT_FLOAT_EQ(0.0f, metrics.precision);
    EXPECT_FLOAT_EQ(0.0f, metrics.recall);
    EXPECT_FLOAT_EQ(0.0f, metrics.f1);
}

TEST(ComputePrF1, F1Normal) {
    prec_recall_f1 metrics = compute_pr_f1(4, 0, 4, 4);
    EXPECT_FLOAT_EQ(1.0f, metrics.precision);
    EXPECT_FLOAT_EQ(0.5f, metrics.recall);
    EXPECT_NEAR(0.6667f, metrics.f1, 1e-4);
}

TEST(ComputePrF1, BothEmpty) {
    // an evaluation with no variants on either side reports perfect scores rather than NaN
    prec_recall_f1 metrics = compute_pr_f1(0, 0, 0, 0);
    EXPECT_FLOAT_EQ(1.0f, metrics.precision);
    EXPECT_FLOAT_EQ(1.0f, metrics.recall);
    EXPECT_FLOAT_EQ(1.0f, metrics.f1);
}

TEST(ComputePrF1, PrecisionKeysOffQueryRecallOffTruth) {
    // asymmetric counts pin which callset feeds which metric; a swapped denominator would still
    // pass the symmetric cases above
    prec_recall_f1 metrics = compute_pr_f1(1, 9, 1, 1);
    EXPECT_FLOAT_EQ(0.1f, metrics.precision);
    EXPECT_FLOAT_EQ(0.5f, metrics.recall);
}

TEST(ComputePrF1, NegativeTruthTp) {
    // the AC_ERR_2_TO_1 correction can drive truth TP negative, so recall goes negative and the
    // precision+recall denominator can be non-zero yet meaningless; F1 is clamped to 0 there
    prec_recall_f1 metrics = compute_pr_f1(0, 1, -2, 3);
    EXPECT_FLOAT_EQ(0.0f, metrics.precision);
    EXPECT_FLOAT_EQ(-2.0f, metrics.recall);
    EXPECT_FLOAT_EQ(0.0f, metrics.f1);
}

/* write_params ***********************************************************************************/

TEST(WriteParams, FiltersJoin) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path() + "/";
    g.filters = {"A", "B", "C"};

    write_params();

    std::string contents = read_file(g.out_prefix + "parameters.tsv");
    EXPECT_NE(std::string::npos, contents.find("\nfilters\tA,B,C\n")) << contents;
}

TEST(WriteParams, FiltersEmpty) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path() + "/";
    g.filters.clear();

    // guards the out-of-bounds g.filters[0] read that an unguarded join would perform
    write_params();

    std::string contents = read_file(g.out_prefix + "parameters.tsv");
    EXPECT_NE(std::string::npos, contents.find("\nfilters\t\n")) << contents;
}

TEST(WriteParams, FiltersSingle) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path() + "/";
    g.filters = {"PASS"};

    write_params();

    std::string contents = read_file(g.out_prefix + "parameters.tsv");
    EXPECT_NE(std::string::npos, contents.find("\nfilters\tPASS\n")) << contents;
    EXPECT_EQ(std::string::npos, contents.find("PASS,"));
}

TEST(WriteParams, FopenFail) {
    GlobalsGuard guard;
    TempDir dir;
    g.out_prefix = dir.path("missing/");

    EXPECT_EXIT(write_params(), testing::ExitedWithCode(1), "Failed to open parameters TSV file");
}

} // namespace
