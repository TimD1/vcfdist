/**
 * @file defs.h
 * @brief Global compile-time constants, type codes, and logging macros for vcfdist.
 */
#ifndef _DEF_H_
#define _DEF_H_

#include <sys/time.h>
#include <unistd.h>
#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_set>
#include <unordered_map>

class timer;
class idx4;

/* Enum-keyed containers **************************************************************************/

/** @brief Fixed-size array keyed by one scoped enum, rejecting subscripts from any other family. */
template <typename K, typename V, std::size_t N>
struct EnumArray {
    std::array<V, N> vals;

    constexpr V & operator[](K k) { return this->vals[static_cast<std::size_t>(k)]; }
    constexpr const V & operator[](K k) const { return this->vals[static_cast<std::size_t>(k)]; }
    constexpr std::size_t size() const { return N; }
    constexpr auto begin() { return this->vals.begin(); }
    constexpr auto end() { return this->vals.end(); }
    constexpr auto begin() const { return this->vals.begin(); }
    constexpr auto end() const { return this->vals.end(); }
};

/** @brief Iterable sequence of the contiguous enumerators of K, from K(0) through K(N-1). */
template <typename K, std::size_t N>
struct EnumRange {
    struct iterator {
        std::size_t i;
        constexpr K operator*() const { return static_cast<K>(this->i); }
        constexpr iterator & operator++() { ++this->i; return *this; }
        constexpr bool operator!=(const iterator & o) const { return this->i != o.i; }
    };
    constexpr iterator begin() const { return iterator{0}; }
    constexpr iterator end() const { return iterator{N}; }
};

/** @brief Converts a scoped enumerator to its underlying integer, for index arithmetic. */
template <typename K>
constexpr std::size_t idx(K k) { return static_cast<std::size_t>(k); }

// misc
#define EPSILON 1e-9 ///< Arbitrary small float value

/** @brief Which dimension of a (contig, supercluster) index pair. */
enum class idxdim_t : int8_t {
    CTG_IDX = 0, ///< Contig dimension
    SC_IDX  = 1, ///< Supercluster dimension
};
constexpr idxdim_t CTG_IDX = idxdim_t::CTG_IDX;
constexpr idxdim_t SC_IDX  = idxdim_t::SC_IDX;
constexpr std::size_t IDXDIM_SLOTS = 2; ///< Slots needed by an idxdim_t-keyed array

/** @defgroup colors ANSI terminal color macros
 *  Conditionally emit ANSI escape codes when stderr is a TTY.
 *  @{
 */
#define COLOR_RED    isatty(STDERR_FILENO) ? "\033[31m" : "" ///< ANSI red color code
#define COLOR_YELLOW isatty(STDERR_FILENO) ? "\033[33m" : "" ///< ANSI yellow color code
#define COLOR_GREEN  isatty(STDERR_FILENO) ? "\033[32m" : "" ///< ANSI green color code
#define COLOR_BLUE   isatty(STDERR_FILENO) ? "\033[34m" : "" ///< ANSI blue color code
#define COLOR_PURPLE isatty(STDERR_FILENO) ? "\033[35m" : "" ///< ANSI purple color code
#define COLOR_WHITE  isatty(STDERR_FILENO) ? "\033[0m"  : "" ///< ANSI reset/white color code
/** @} */

/** @brief Edit type of a variant relative to the reference. */
enum class edittype_t : int8_t {
    TYPE_REF = 0, ///< Reference (no variant)
    TYPE_SUB = 1, ///< Substitution (SNP)
    TYPE_INS = 2, ///< Insertion
    TYPE_DEL = 3, ///< Deletion
    TYPE_CPX = 4, ///< Complex variant
};
constexpr edittype_t TYPE_REF = edittype_t::TYPE_REF;
constexpr edittype_t TYPE_SUB = edittype_t::TYPE_SUB;
constexpr edittype_t TYPE_INS = edittype_t::TYPE_INS;
constexpr edittype_t TYPE_DEL = edittype_t::TYPE_DEL;
constexpr edittype_t TYPE_CPX = edittype_t::TYPE_CPX;
constexpr std::size_t EDITTYPE_SLOTS = 5; ///< Slots needed by an edittype_t-keyed array
constexpr int8_t TYPES = 5;               ///< Total number of variant types
static_assert(EDITTYPE_SLOTS == std::size_t(TYPES), "edittype_t slots must match TYPES");

/** @brief Variant size class. */
enum class sizeclass_t : int8_t {
    VARTYPE_SNP   = 0, ///< Single nucleotide polymorphism
    VARTYPE_INDEL = 1, ///< Small insertion or deletion
    VARTYPE_SV    = 2, ///< Structural variant (size >= sv_threshold)
    VARTYPE_ALL   = 3, ///< All variant size classes
};
constexpr sizeclass_t VARTYPE_SNP   = sizeclass_t::VARTYPE_SNP;
constexpr sizeclass_t VARTYPE_INDEL = sizeclass_t::VARTYPE_INDEL;
constexpr sizeclass_t VARTYPE_SV    = sizeclass_t::VARTYPE_SV;
constexpr sizeclass_t VARTYPE_ALL   = sizeclass_t::VARTYPE_ALL;
constexpr std::size_t SIZECLASS_SLOTS = 4; ///< Slots needed by a sizeclass_t-keyed array
constexpr int8_t VARTYPES = 4;             ///< Total number of variant size classes
static_assert(SIZECLASS_SLOTS == std::size_t(VARTYPES), "sizeclass_t slots must match VARTYPES");

/** @brief Which haplotype of a diploid sample. */
enum class hap_t : int8_t {
    HAP1 = 0, ///< First haplotype
    HAP2 = 1, ///< Second haplotype
};
constexpr hap_t HAP1 = hap_t::HAP1;
constexpr hap_t HAP2 = hap_t::HAP2;
constexpr std::size_t HAP_SLOTS = 2; ///< Slots needed by a hap_t-keyed array
constexpr int8_t HAPS = 2;           ///< Number of haplotypes
static_assert(HAP_SLOTS == std::size_t(HAPS), "hap_t slots must match HAPS");

/** @brief Returns the other haplotype of the pair. */
constexpr hap_t other_hap(hap_t h) { return h == HAP1 ? HAP2 : HAP1; }

/** @brief Location of a variant relative to the BED regions. */
enum class bedloc_t : int8_t {
    BED_OUTSIDE = 0, ///< Variant is fully outside all BED regions
    BED_INSIDE  = 1, ///< Variant is fully inside a BED region
    BED_BORDER  = 2, ///< Variant overlaps a BED region boundary
    BED_OFFCTG  = 3, ///< Variant is on a contig not present in BED file
};
constexpr bedloc_t BED_OUTSIDE = bedloc_t::BED_OUTSIDE;
constexpr bedloc_t BED_INSIDE  = bedloc_t::BED_INSIDE;
constexpr bedloc_t BED_BORDER  = bedloc_t::BED_BORDER;
constexpr bedloc_t BED_OFFCTG  = bedloc_t::BED_OFFCTG;
constexpr std::size_t BEDLOC_SLOTS = 4; ///< Subscript slots needed by a bedloc_t-keyed array

/** @brief Simplified genotype of one variant. */
enum class gt_t : int8_t {
    GT_REF       =  0, ///< Haploid reference genotype (0)
    GT_ALT1      =  1, ///< Haploid alternate genotype (1)
    GT_REF_REF   =  2, ///< Diploid reference/reference genotype (0|0)
    GT_REF_ALT1  =  3, ///< Diploid reference/alternate genotype (0|1)
    GT_ALT1_REF  =  4, ///< Diploid alternate/reference genotype (1|0)
    GT_ALT1_ALT1 =  5, ///< Diploid homozygous alternate genotype (1|1)
    GT_ALT1_ALT2 =  6, ///< Diploid compound heterozygous genotype (1|2)
    GT_ALT2_ALT1 =  7, ///< Diploid compound heterozygous genotype (2|1)
    GT_MISSING   =  8, ///< No-call genotype, every allele missing (.|. or .)
    GT_HALF      =  9, ///< Half-call genotype, exactly one allele missing (1|. or .|1)
    GT_OTHER     = 10, ///< Other/unknown genotype
};
constexpr gt_t GT_REF       = gt_t::GT_REF;
constexpr gt_t GT_ALT1      = gt_t::GT_ALT1;
constexpr gt_t GT_REF_REF   = gt_t::GT_REF_REF;
constexpr gt_t GT_REF_ALT1  = gt_t::GT_REF_ALT1;
constexpr gt_t GT_ALT1_REF  = gt_t::GT_ALT1_REF;
constexpr gt_t GT_ALT1_ALT1 = gt_t::GT_ALT1_ALT1;
constexpr gt_t GT_ALT1_ALT2 = gt_t::GT_ALT1_ALT2;
constexpr gt_t GT_ALT2_ALT1 = gt_t::GT_ALT2_ALT1;
constexpr gt_t GT_MISSING   = gt_t::GT_MISSING;
constexpr gt_t GT_HALF      = gt_t::GT_HALF;
constexpr gt_t GT_OTHER     = gt_t::GT_OTHER;
constexpr std::size_t GT_SLOTS = 11; ///< Slots needed by a gt_t-keyed array

/** @brief Benchmark error type. */
enum class errtype_t : int8_t {
    ERRTYPE_TP = 0, ///< True positive
    ERRTYPE_FP = 1, ///< False positive
    ERRTYPE_FN = 2, ///< False negative
    ERRTYPE_UN = 3, ///< Unknown (not yet evaluated)
};
constexpr errtype_t ERRTYPE_TP = errtype_t::ERRTYPE_TP;
constexpr errtype_t ERRTYPE_FP = errtype_t::ERRTYPE_FP;
constexpr errtype_t ERRTYPE_FN = errtype_t::ERRTYPE_FN;
constexpr errtype_t ERRTYPE_UN = errtype_t::ERRTYPE_UN;
constexpr std::size_t ERRTYPE_SLOTS = 4; ///< Subscript slots needed by an errtype_t-keyed array
constexpr int8_t ERRTYPES = 4;           ///< Total number of error types
static_assert(ERRTYPE_SLOTS == std::size_t(ERRTYPES), "errtype_t slots must match ERRTYPES");

/**
 * @brief A site's truth alternate allele count, then its query alternate allele count.
 *
 * The direction is absolute, so the query and truth records of a matched site carry the same value;
 * what differs is which genotype supplies which count, since a record's orig_gt is its own
 * callset's call and its matched_gt is the other callset's genotype recovered by alignment. A
 * callset's own allele count is never zero, so a query record never reaches AC_ERR_1_TO_0 or
 * AC_ERR_2_TO_0 and a truth record never reaches AC_ERR_0_TO_1 or AC_ERR_0_TO_2.
 */
enum class ac_errtype_t : int8_t {
    AC_ERR_0_TO_1 = 0, ///< 0/0 -> 0/1: 1 QUERY_FP
    AC_ERR_0_TO_2 = 1, ///< 0/0 -> 1/1: 2 QUERY_FP
    AC_ERR_1_TO_0 = 2, ///< 0/1 -> 0/0: 1 TRUTH_FN
    AC_ERR_1_TO_1 = 3, ///< 0/1 -> 0/1: 1 QUERY_TP, 1 TRUTH_TP
    AC_ERR_1_TO_2 = 4, ///< 0/1 -> 1/1: 1 QUERY_TP, 1 TRUTH_TP
    AC_ERR_2_TO_0 = 5, ///< 1/1 -> 0/0: 2 TRUTH_FN
    AC_ERR_2_TO_1 = 6, ///< 1/1 -> 0/1: 1 QUERY_TP, 1 TRUTH_TP
    AC_ERR_2_TO_2 = 7, ///< 1/1 -> 1/1: 2 QUERY_TP, 2 TRUTH_TP
    AC_UNKNOWN    = 8, ///< Unknown allele count error type
};
constexpr ac_errtype_t AC_ERR_0_TO_1 = ac_errtype_t::AC_ERR_0_TO_1;
constexpr ac_errtype_t AC_ERR_0_TO_2 = ac_errtype_t::AC_ERR_0_TO_2;
constexpr ac_errtype_t AC_ERR_1_TO_0 = ac_errtype_t::AC_ERR_1_TO_0;
constexpr ac_errtype_t AC_ERR_1_TO_1 = ac_errtype_t::AC_ERR_1_TO_1;
constexpr ac_errtype_t AC_ERR_1_TO_2 = ac_errtype_t::AC_ERR_1_TO_2;
constexpr ac_errtype_t AC_ERR_2_TO_0 = ac_errtype_t::AC_ERR_2_TO_0;
constexpr ac_errtype_t AC_ERR_2_TO_1 = ac_errtype_t::AC_ERR_2_TO_1;
constexpr ac_errtype_t AC_ERR_2_TO_2 = ac_errtype_t::AC_ERR_2_TO_2;
constexpr ac_errtype_t AC_UNKNOWN    = ac_errtype_t::AC_UNKNOWN;
constexpr std::size_t AC_ERRTYPE_SLOTS = 9; ///< Slots needed by an ac_errtype_t-keyed array
constexpr int8_t AC_ERRTYPES = 8;           ///< Total number of allele count error types
// AC_UNKNOWN is a storable sentinel rather than a real error type, so it needs a slot of its own
static_assert(AC_ERRTYPE_SLOTS == std::size_t(AC_ERRTYPES) + 1, "AC_UNKNOWN needs its own slot");

/** @brief Phase switch or flip error type. */
enum class switchtype_t : int8_t {
    SWITCHTYPE_FLIP            = 0, ///< A phase flip error
    SWITCHTYPE_SWITCH          = 1, ///< A phase switch error
    SWITCHTYPE_SWITCH_AND_FLIP = 2, ///< Combined switch and flip error
    SWITCHTYPE_SWITCH_ERR      = 3, ///< Switch error (alternate classification)
    SWITCHTYPE_FLIP_BEG        = 4, ///< Beginning boundary of a flip error region
    SWITCHTYPE_FLIP_END        = 5, ///< Ending boundary of a flip error region
    SWITCHTYPE_NONE            = 6, ///< No error
};
constexpr switchtype_t SWITCHTYPE_FLIP            = switchtype_t::SWITCHTYPE_FLIP;
constexpr switchtype_t SWITCHTYPE_SWITCH          = switchtype_t::SWITCHTYPE_SWITCH;
constexpr switchtype_t SWITCHTYPE_SWITCH_AND_FLIP = switchtype_t::SWITCHTYPE_SWITCH_AND_FLIP;
constexpr switchtype_t SWITCHTYPE_SWITCH_ERR      = switchtype_t::SWITCHTYPE_SWITCH_ERR;
constexpr switchtype_t SWITCHTYPE_FLIP_BEG        = switchtype_t::SWITCHTYPE_FLIP_BEG;
constexpr switchtype_t SWITCHTYPE_FLIP_END        = switchtype_t::SWITCHTYPE_FLIP_END;
constexpr switchtype_t SWITCHTYPE_NONE            = switchtype_t::SWITCHTYPE_NONE;
constexpr std::size_t SWITCHTYPE_SLOTS = 7; ///< Subscript slots needed by a switchtype_t-keyed array
constexpr int8_t SWITCHTYPES = 7;           ///< Total number of switch/flip error types
static_assert(SWITCHTYPE_SLOTS == std::size_t(SWITCHTYPES), "switchtype_t slots must match SWITCHTYPES");

/** @brief Pipeline stage identifying one timer. */
enum class stage_t : int8_t {
    TIME_READ       = 0, ///< Input reading stage
    TIME_CLUSTER    = 1, ///< Variant clustering stage
    TIME_ALIGN_EVAL = 2, ///< Alignment and evaluation stage
    TIME_PHASE      = 3, ///< Phasing stage
    TIME_WRITE      = 4, ///< Output writing stage
    TIME_TOTAL      = 5, ///< Total pipeline runtime
};
constexpr stage_t TIME_READ       = stage_t::TIME_READ;
constexpr stage_t TIME_CLUSTER    = stage_t::TIME_CLUSTER;
constexpr stage_t TIME_ALIGN_EVAL = stage_t::TIME_ALIGN_EVAL;
constexpr stage_t TIME_PHASE      = stage_t::TIME_PHASE;
constexpr stage_t TIME_WRITE      = stage_t::TIME_WRITE;
constexpr stage_t TIME_TOTAL      = stage_t::TIME_TOTAL;
constexpr std::size_t STAGE_SLOTS = 6; ///< Number of pipeline stage timers

/** @brief Which of the two callsets a variant or container belongs to. */
enum class callset_t : int8_t {
    QUERY = 0, ///< Query callset
    TRUTH = 1, ///< Truth callset
};
constexpr callset_t QUERY = callset_t::QUERY;
constexpr callset_t TRUTH = callset_t::TRUTH;
constexpr std::size_t CALLSET_SLOTS = 2; ///< Slots needed by a callset_t-keyed array
constexpr int8_t CALLSETS = 2;           ///< Number of callsets
static_assert(CALLSET_SLOTS == std::size_t(CALLSETS), "callset_t slots must match CALLSETS");

/** @brief Alignment backtracking pointer. Values are non-contiguous, so no EnumArray keys on it. */
enum class ptr_t : int8_t {
    PTR_INS = 1, ///< Insertion backtracking pointer
    PTR_DEL = 2, ///< Deletion backtracking pointer
    PTR_MAT = 4, ///< Match backtracking pointer
    PTR_SUB = 8, ///< Substitution backtracking pointer
};
constexpr ptr_t PTR_INS = ptr_t::PTR_INS;
constexpr ptr_t PTR_DEL = ptr_t::PTR_DEL;
constexpr ptr_t PTR_MAT = ptr_t::PTR_MAT;
constexpr ptr_t PTR_SUB = ptr_t::PTR_SUB;

/** @brief Smith-Waterman alignment matrix. */
enum class mat_t : int8_t {
    MAT_SUB = 0, ///< Substitution matrix
    MAT_INS = 1, ///< Insertion matrix
    MAT_DEL = 2, ///< Deletion matrix
};
constexpr mat_t MAT_SUB = mat_t::MAT_SUB;
constexpr mat_t MAT_INS = mat_t::MAT_INS;
constexpr mat_t MAT_DEL = mat_t::MAT_DEL;
constexpr std::size_t MAT_SLOTS = 3; ///< Subscript slots needed by a mat_t-keyed array
constexpr int8_t MATS = 3;           ///< Total number of Smith-Waterman matrices
static_assert(MAT_SLOTS == std::size_t(MATS), "mat_t slots must match MATS");

// each alignment matrix sits exactly one below the edit type it aligns, which mat_to_edittype uses
static_assert(idx(MAT_SUB) + 1 == idx(TYPE_SUB), "mat_t must sit one below edittype_t");
static_assert(idx(MAT_INS) + 1 == idx(TYPE_INS), "mat_t must sit one below edittype_t");
static_assert(idx(MAT_DEL) + 1 == idx(TYPE_DEL), "mat_t must sit one below edittype_t");

/** @brief Returns the variant edit type whose alignment matrix is the given matrix. */
constexpr edittype_t mat_to_edittype(mat_t m) {
    return static_cast<edittype_t>(static_cast<int8_t>(m) + 1);
}

/** @brief Phasing state of a variant or phase block. */
enum class phase_t : int8_t {
    PHASE_ORIG = 0, ///< Keep original haplotype assignment
    PHASE_SWAP = 1, ///< Swap haplotype assignment
    PHASE_NONE = 2, ///< No phasing information available
};
constexpr phase_t PHASE_ORIG = phase_t::PHASE_ORIG;
constexpr phase_t PHASE_SWAP = phase_t::PHASE_SWAP;
constexpr phase_t PHASE_NONE = phase_t::PHASE_NONE;
constexpr std::size_t PHASE_SLOTS = 3; ///< Slots needed by a phase_t-keyed array
constexpr int8_t PHASES = 2;           ///< Number of phasing states (ORIG and SWAP)
// PHASE_NONE is a storable sentinel rather than a phasing state, so it needs a slot of its own
static_assert(PHASE_SLOTS == std::size_t(PHASES) + 1, "PHASE_NONE needs its own slot");

/** @brief Phasing dynamic-programming backtrack pointer. */
enum class phaseptr_t : int8_t {
    PHASE_PTR_KEEP = 0, ///< Keep current phase
    PHASE_PTR_SWAP = 1, ///< Swap phase at this variant
};
constexpr phaseptr_t PHASE_PTR_KEEP = phaseptr_t::PHASE_PTR_KEEP;
constexpr phaseptr_t PHASE_PTR_SWAP = phaseptr_t::PHASE_PTR_SWAP;
constexpr std::size_t PHASEPTR_SLOTS = 2; ///< Slots needed by a phaseptr_t-keyed array

/** @brief Returns the opposite phasing state. */
constexpr phase_t other_phase(phase_t p) { return p == PHASE_ORIG ? PHASE_SWAP : PHASE_ORIG; }

/** @defgroup logging_macros Timestamped logging macros
 *  Print colored, timestamped messages to stderr. Exit on ERROR.
 *  @{
 */

/**
 * @def WARN(f_, ...)
 * @brief Prints a yellow timestamped WARNING message to stderr.
 */
#define WARN(f_, ...)                                         \
{                                                             \
    struct tm _tm123_;                                        \
    struct timeval _xxtv123_;                                 \
    gettimeofday(&_xxtv123_, NULL);                           \
    localtime_r(&_xxtv123_.tv_sec, &_tm123_);                 \
    fprintf(stderr, "%s[WARN  %s %02d:%02d:%02d]%s ",         \
            COLOR_YELLOW,  g.PROGRAM.data(), _tm123_.tm_hour, \
            _tm123_.tm_min, _tm123_.tm_sec, COLOR_WHITE);     \
    fprintf(stderr, (f_), ##__VA_ARGS__);                     \
    fprintf(stderr, "\n");                                    \
};

/**
 * @def INFO(f_, ...)
 * @brief Prints a green timestamped INFO message to stderr.
 */
#define INFO(f_, ...)                                         \
{                                                             \
    struct tm _tm123_;                                        \
    struct timeval _xxtv123_;                                 \
    gettimeofday(&_xxtv123_, NULL);                           \
    localtime_r(&_xxtv123_.tv_sec, &_tm123_);                 \
    fprintf(stderr, "%s[INFO  %s %02d:%02d:%02d]%s ",         \
            COLOR_GREEN, g.PROGRAM.data(), _tm123_.tm_hour,   \
            _tm123_.tm_min, _tm123_.tm_sec, COLOR_WHITE);     \
    fprintf(stderr, (f_), ##__VA_ARGS__);                     \
    fprintf(stderr, "\n");                                    \
};

/**
 * @def ERROR(f_, ...)
 * @brief Prints a red timestamped ERROR message to stderr and exits with code 1.
 */
#define ERROR(f_, ...)                                        \
{                                                             \
    struct tm _tm123_;                                        \
    struct timeval _xxtv123_;                                 \
    gettimeofday(&_xxtv123_, NULL);                           \
    localtime_r(&_xxtv123_.tv_sec, &_tm123_);                 \
    fprintf(stderr, "%s[ERROR %s %02d:%02d:%02d]%s ",         \
            COLOR_RED, g.PROGRAM.data(), _tm123_.tm_hour,     \
            _tm123_.tm_min, _tm123_.tm_sec, COLOR_WHITE);     \
    fprintf(stderr, (f_), ##__VA_ARGS__);                     \
    fprintf(stderr, "\n");                                    \
    std::exit(1);                                             \
};
/** @} */

#endif
