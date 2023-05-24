#ifndef _DEF_H_
#define _DEF_H_

#include <sys/time.h>
#include <unordered_set>
#include <unordered_map>

class timer;
class idx1;
class idx2;

// misc
#define FAIL 0
#define PASS 1

#define FALSE 0
#define TRUE  1

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
#define VARTYPES      2

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
#define ERRTYPE_PP 3 // partial positive (reduces ED, but not TP)
#define ERRTYPE_PE 4 // phase error (0|1 -> 1|0)
#define ERRTYPE_GE 5 // genotype error (0|1 -> 1|1)
#define ERRTYPE_UN 6 // unknown
#define ERRTYPES   7

// runtime timers
#define TIME_READ            0
#define TIME_CLUST           1
#define TIME_RECLUST         2
  #define TIME_CL_GENSTR     3
  #define TIME_CL_SUBSTR     4
  #define TIME_CL_REACH      5
  #define TIME_CL_INIT       6
  #define TIME_CL_EXTEND     7
  #define TIME_CL_NEXT       8
  #define TIME_CL_MAX        9
  #define TIME_CL_BUFF_INIT  10
  #define TIME_CL_BUFF_CLEAR 11
#define TIME_REALN           12
#define TIME_SUPCLUST        13
#define TIME_ALIGN           14
  #define TIME_PR            15
    #define TIME_PR_GENPTR   16
    #define TIME_PR_ALN      17
    #define TIME_PR_INIT     18
    #define TIME_PR_EXTEND   19
    #define TIME_PR_NEXT     20
    #define TIME_PR_PATH     21
    #define TIME_PR_STAT     22
  #define TIME_SW            23
    #define TIME_SW_ALN      24
#define TIME_PHASE           25
#define TIME_WRITE           26
#define TIME_TOTAL           27

// alignment
#define QUERY    0
#define REF      1
#define TRUTH    1
#define CALLSETS 2

#define QUERY1_TRUTH1 0 // four possible alignments of truth and query haps
#define QUERY1_TRUTH2 1
#define QUERY2_TRUTH1 2
#define QUERY2_TRUTH2 3

#define PEN_SUB    0 // penalties order in boolean array *_penalties_set
#define PEN_OPEN   1
#define PEN_EXTEND 2

#define PTR_NONE    0  // backtracking pointer flags
#define PTR_INS     1
#define PTR_DEL     2
#define PTR_MAT     4
#define PTR_SUB     8
#define PTR_SWP_MAT 16
#define PTR_LPATH   32
#define PATH        32
#define PTR_RPATH   64
#define MAIN_PATH   96
#define PTR_SYNC    128

// 2 x N pointer array stores REF <-> QUERY/TRUTH
#define PTRS 0 // dimension for pointer position
#define FLAGS 1 // dimension for flags
#define PTR_DIMS 2
#define PTR_VARIANT 1
#define PTR_VAR_BEG 2
#define PTR_VAR_END 4

#define MAT_SUB 0 // three matrices for Smith-Waterman
#define MAT_INS 1
#define MAT_DEL 2
#define MATS    3

// phasing
#define PHASE_ORIG 0
#define PHASE_SWAP 1
#define PHASE_NONE 2

#define PHASE_PTR_KEEP 0
#define PHASE_PTR_SWAP 1

// printing
#define WARN(f_, ...)                                           \
{                                                               \
    struct tm _tm123_;                                          \
    struct timeval _xxtv123_;                                   \
    gettimeofday(&_xxtv123_, NULL);                             \
    localtime_r(&_xxtv123_.tv_sec, &_tm123_);                   \
    fprintf(stderr, "\033[33m[WARN  %s %02d:%02d:%02d]\033[0m ",\
            g.PROGRAM.data(), _tm123_.tm_hour,_tm123_.tm_min,   \
            _tm123_.tm_sec);                                    \
    fprintf(stderr, (f_), ##__VA_ARGS__);                       \
    fprintf(stderr, "\n");                                      \
};

#define INFO(f_, ...)                                           \
{                                                               \
    struct tm _tm123_;                                          \
    struct timeval _xxtv123_;                                   \
    gettimeofday(&_xxtv123_, NULL);                             \
    localtime_r(&_xxtv123_.tv_sec, &_tm123_);                   \
    fprintf(stderr, "\033[32m[INFO  %s %02d:%02d:%02d]\033[0m ",\
            g.PROGRAM.data(), _tm123_.tm_hour,_tm123_.tm_min,   \
            _tm123_.tm_sec);                                    \
    fprintf(stderr, (f_), ##__VA_ARGS__);                       \
    fprintf(stderr, "\n");                                      \
};

#define ERROR(f_, ...)                                           \
{                                                                \
    struct tm _tm123_;                                           \
    struct timeval _xxtv123_;                                    \
    gettimeofday(&_xxtv123_, NULL);                              \
    localtime_r(&_xxtv123_.tv_sec, &_tm123_);                    \
    fprintf(stderr, "\033[31m[ERROR %s %02d:%02d:%02d]\033[0m ", \
            g.PROGRAM.data(), _tm123_.tm_hour,_tm123_.tm_min,    \
            _tm123_.tm_sec);                                     \
    fprintf(stderr, (f_), ##__VA_ARGS__);                        \
    fprintf(stderr, "\n");                                       \
    std::exit(1);                                                \
};

#endif
