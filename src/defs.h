#ifndef _DEF_H_
#define _DEF_H_

#include <sys/time.h>
#include <unistd.h>
#include <unordered_set>
#include <unordered_map>

class timer;
class idx3;
class idx4;

// misc
#define FALSE 0
#define TRUE  1

#define CTG_IDX 0
#define SC_IDX  1

#define COLOR_RED    isatty(STDERR_FILENO) ? "\033[31m" : ""
#define COLOR_YELLOW isatty(STDERR_FILENO) ? "\033[33m" : ""
#define COLOR_GREEN  isatty(STDERR_FILENO) ? "\033[32m" : ""
#define COLOR_BLUE   isatty(STDERR_FILENO) ? "\033[34m" : ""
#define COLOR_PURPLE isatty(STDERR_FILENO) ? "\033[35m" : ""
#define COLOR_WHITE  isatty(STDERR_FILENO) ? "\033[0m"  : ""

// variant
#define TYPE_REF 0
#define TYPE_ALL 0
#define TYPE_SUB 1
#define TYPE_INS 2
#define TYPE_DEL 3
#define TYPE_CPX 4
#define TYPE_INDEL 4
#define TYPES 5

#define VARTYPE_SNP   0
#define VARTYPE_INDEL 1
#define VARTYPE_SV    2
#define VARTYPE_ALL   3
#define VARTYPES      4

#define HAP1 0
#define HAP2 1
#define HAPS 2

#define BED_OUTSIDE 0
#define BED_INSIDE  1
#define BED_BORDER  2
#define BED_OFFCTG  3

#define GT_REF       0
#define GT_ALT1      1
#define GT_REF_REF   2
#define GT_REF_ALT1  3
#define GT_ALT1_REF  4
#define GT_ALT1_ALT1 5
#define GT_ALT1_ALT2 6
#define GT_ALT2_ALT1 7
#define GT_MISSING   8
#define GT_OTHER     9

#define ERRTYPE_TP 0 // true positive
#define ERRTYPE_FP 1 // false positive
#define ERRTYPE_FN 2 // false negative
#define ERRTYPE_UN 3 // unknown (not evaluated yet)
#define ERRTYPE_NE 3 // not evaluated (and will not be)
#define ERRTYPES   4

#define AC_ERR_0_TO_1  0 // e.g. 0/0 -> 0/1 (1 QUERY_FP)
#define AC_ERR_0_TO_2  1 // e.g. 0/0 -> 1/1 (2 QUERY_FP)
#define AC_ERR_1_TO_0  2 // e.g. 0/1 -> 0/0 (1 TRUTH_FN)
#define AC_ERR_1_TO_1  3 // e.g. 0/1 -> 0/1 (1 QUERY_TP, 1 TRUTH_TP)
#define AC_ERR_1_TO_2  4 // e.g. 0/1 -> 1/1 (1 QUERY_TP, 1 TRUTH_TP, 1 QUERY_FP)
#define AC_ERR_2_TO_0  5 // e.g. 1/1 -> 0/0 (2 TRUTH_FN)
#define AC_ERR_2_TO_1  6 // e.g. 1/1 -> 0/1 (1 QUERY_TP, 1 TRUTH_TP, 1 TRUTH_FN)
#define AC_ERR_2_TO_2  7 // e.g. 1/1 -> 1/1 (2 QUERY_TP, 2 TRUTH_TP)
#define AC_UNKNOWN     8
#define AC_ERRTYPES    8 // allele count error types

#define SWITCHTYPE_FLIP       0
#define SWITCHTYPE_SWITCH     1
#define SWITCHTYPE_SWITCH_AND_FLIP     2
#define SWITCHTYPE_SWITCH_ERR 3
#define SWITCHTYPE_FLIP_BEG   4
#define SWITCHTYPE_FLIP_END   5
#define SWITCHTYPE_NONE       6
#define SWITCHTYPES           7

// runtime timers
#define TIME_READ           0
#define TIME_EXACT_CLUSTER  1
#define TIME_ALIGN_EVAL_1   2
#define TIME_SIMPLE_CLUSTER 3
#define TIME_ALIGN_EVAL_2   4
#define TIME_PHASE          5
#define TIME_WRITE          6
#define TIME_TOTAL          7

// alignment
#define QUERY    0
#define REF      1
#define TRUTH    1
#define CALLSETS 2

#define PTR_NONE 0  // backtracking pointer flags
#define PTR_INS  1
#define PTR_DEL  2
#define PTR_MAT  4
#define PTR_SUB  8
// [14 QUERY NODE BITS] [14 TRUTH NODE BITS] [4 POINTER BITS]
// NOTE: 2^14-1 = 16383 is a fairly safe upper limit on the number of nodes in a supercluster. If we
// do hit this limit, vcfdist reports an error.
#define PTR_BITS   4  // INS, DEL, MAT, SUB
#define NODE_BITS  14 // for pointer to previous node (up to 16,383)
#define PTR_MASK   0x0000000F
#define TNODE_MASK 0x0003FFF0
#define QNODE_MASK 0xFFFC0000

#define MAT_SUB 0 // three matrices for Smith-Waterman
#define MAT_INS 1
#define MAT_DEL 2
#define MATS    3

// phasing
#define PHASE_ORIG 0
#define PHASE_SWAP 1
#define PHASE_NONE 2
#define PHASES 2

#define PHASE_PTR_KEEP 0
#define PHASE_PTR_SWAP 1

// printing
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

#endif
