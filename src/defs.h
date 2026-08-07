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

#define CTG_IDX 0 ///< Index dimension for contig
#define SC_IDX  1 ///< Index dimension for supercluster

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

/** @defgroup type_constants Variant type constants (TYPE_*)
 *  @{
 */
#define TYPE_REF   0 ///< Reference (no variant)
#define TYPE_SUB   1 ///< Substitution (SNP)
#define TYPE_INS   2 ///< Insertion
#define TYPE_DEL   3 ///< Deletion
#define TYPE_CPX   4 ///< Complex variant
#define TYPES      5 ///< Total number of variant types
/** @} */

/** @defgroup vartype_constants Variant size-class constants (VARTYPE_*)
 *  @{
 */
#define VARTYPE_SNP   0 ///< Single nucleotide polymorphism
#define VARTYPE_INDEL 1 ///< Small insertion or deletion
#define VARTYPE_SV    2 ///< Structural variant (size >= sv_threshold)
#define VARTYPE_ALL   3 ///< All variant size classes
#define VARTYPES      4 ///< Total number of variant size classes
/** @} */

/** @defgroup hap_constants Haplotype index constants
 *  @{
 */
#define HAP1 0 ///< First haplotype index
#define HAP2 1 ///< Second haplotype index
#define HAPS 2 ///< Number of haplotypes
/** @} */

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

/** @defgroup gt_constants Genotype code constants (GT_*)
 *  @{
 */
#define GT_REF       0 ///< Haploid reference genotype (0)
#define GT_ALT1      1 ///< Haploid alternate genotype (1)
#define GT_REF_REF   2 ///< Diploid reference/reference genotype (0|0)
#define GT_REF_ALT1  3 ///< Diploid reference/alternate genotype (0|1)
#define GT_ALT1_REF  4 ///< Diploid alternate/reference genotype (1|0)
#define GT_ALT1_ALT1 5 ///< Diploid homozygous alternate genotype (1|1)
#define GT_ALT1_ALT2 6 ///< Diploid compound heterozygous genotype (1|2)
#define GT_ALT2_ALT1 7 ///< Diploid compound heterozygous genotype (2|1)
#define GT_MISSING   8 ///< No-call genotype, every allele missing (.|. or .)
#define GT_HALF      9 ///< Half-call genotype, exactly one allele missing (1|. or .|1)
#define GT_OTHER    10 ///< Other/unknown genotype
/** @} */

/** @defgroup errtype_constants Benchmark error type constants (ERRTYPE_*)
 *  @{
 */
#define ERRTYPE_TP 0 ///< True positive
#define ERRTYPE_FP 1 ///< False positive
#define ERRTYPE_FN 2 ///< False negative
#define ERRTYPE_UN 3 ///< Unknown (not yet evaluated)
#define ERRTYPES   4 ///< Total number of error types
/** @} */

/** @defgroup ac_err_constants Allele count error type constants (AC_ERR_*)
 *  A site's truth alternate allele count, then its query alternate allele count. The direction is
 *  absolute, so the query and truth records of a matched site carry the same value; what differs is
 *  which genotype supplies which count, since a record's orig_gt is its own callset's call and its
 *  calc_gt is the other callset's genotype recovered by alignment. A callset's own allele count is
 *  never zero, so a query record never reaches AC_ERR_1_TO_0 or AC_ERR_2_TO_0 and a truth record
 *  never reaches AC_ERR_0_TO_1 or AC_ERR_0_TO_2.
 *  @{
 */
#define AC_ERR_0_TO_1  0 ///< 0/0 -> 0/1: 1 QUERY_FP
#define AC_ERR_0_TO_2  1 ///< 0/0 -> 1/1: 2 QUERY_FP
#define AC_ERR_1_TO_0  2 ///< 0/1 -> 0/0: 1 TRUTH_FN
#define AC_ERR_1_TO_1  3 ///< 0/1 -> 0/1: 1 QUERY_TP, 1 TRUTH_TP
#define AC_ERR_1_TO_2  4 ///< 0/1 -> 1/1: 1 QUERY_TP, 1 TRUTH_TP 
#define AC_ERR_2_TO_0  5 ///< 1/1 -> 0/0: 2 TRUTH_FN
#define AC_ERR_2_TO_1  6 ///< 1/1 -> 0/1: 1 QUERY_TP, 1 TRUTH_TP
#define AC_ERR_2_TO_2  7 ///< 1/1 -> 1/1: 2 QUERY_TP, 2 TRUTH_TP
#define AC_UNKNOWN     8 ///< Unknown allele count error type
#define AC_ERRTYPES    8 ///< Total number of allele count error types
/** @} */

/** @defgroup switchtype_constants Phase switch/flip error type constants (SWITCHTYPE_*)
 *  @{
 */
#define SWITCHTYPE_FLIP            0 ///< A phase flip error
#define SWITCHTYPE_SWITCH          1 ///< A phase switch error
#define SWITCHTYPE_SWITCH_AND_FLIP 2 ///< Combined switch and flip error
#define SWITCHTYPE_SWITCH_ERR      3 ///< Switch error (alternate classification)
#define SWITCHTYPE_FLIP_BEG        4 ///< Beginning boundary of a flip error region
#define SWITCHTYPE_FLIP_END        5 ///< Ending boundary of a flip error region
#define SWITCHTYPE_NONE            6 ///< No error
#define SWITCHTYPES                7 ///< Total number of switch/flip error types
/** @} */

/** @brief Pipeline stage identifying one timer. */
enum class timer_t : int8_t {
    TIME_READ       = 0, ///< Input reading stage
    TIME_CLUSTER    = 1, ///< Variant clustering stage
    TIME_ALIGN_EVAL = 2, ///< Alignment and evaluation stage
    TIME_PHASE      = 3, ///< Phasing stage
    TIME_WRITE      = 4, ///< Output writing stage
    TIME_TOTAL      = 5, ///< Total pipeline runtime
};
constexpr timer_t TIME_READ       = timer_t::TIME_READ;
constexpr timer_t TIME_CLUSTER    = timer_t::TIME_CLUSTER;
constexpr timer_t TIME_ALIGN_EVAL = timer_t::TIME_ALIGN_EVAL;
constexpr timer_t TIME_PHASE      = timer_t::TIME_PHASE;
constexpr timer_t TIME_WRITE      = timer_t::TIME_WRITE;
constexpr timer_t TIME_TOTAL      = timer_t::TIME_TOTAL;
constexpr std::size_t TIMER_SLOTS = 6; ///< Number of pipeline stage timers

/** @defgroup callset_constants Callset index constants
 *  @{
 */
#define QUERY    0 ///< Query callset index
#define TRUTH    1 ///< Truth callset index
#define CALLSETS 2 ///< Number of callsets
/** @} */

/** @defgroup ptr_constants Alignment backtracking pointer constants (PTR_*)
 *  @{
 */
#define PTR_INS  1  ///< Insertion backtracking pointer
#define PTR_DEL  2  ///< Deletion backtracking pointer
#define PTR_MAT  4  ///< Match backtracking pointer
#define PTR_SUB  8  ///< Substitution backtracking pointer
/** @} */

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

/** @defgroup phase_constants Phasing state constants (PHASE_*)
 *  @{
 */
#define PHASE_ORIG 0 ///< Keep original haplotype assignment
#define PHASE_SWAP 1 ///< Swap haplotype assignment
#define PHASE_NONE 2 ///< No phasing information available
#define PHASES     2 ///< Number of phasing states (ORIG and SWAP)

#define PHASE_PTR_KEEP 0 ///< DP pointer: keep current phase
#define PHASE_PTR_SWAP 1 ///< DP pointer: swap phase at this variant
/** @} */

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
