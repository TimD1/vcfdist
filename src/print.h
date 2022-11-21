#ifndef _PRINT_H_
#define _PRINT_H_

#include <string>
#include <sys/time.h>

#include "globals.h"
#include "phase.h"

#define VARTYPE_SNP   0
#define VARTYPE_INDEL 1
#define VARTYPES      2

std::string GREEN(char c);
std::string GREEN(std::string str);
std::string RED(char c);
std::string RED(std::string str);
std::string BLUE(int i);
std::string BLUE(char c);
std::string BLUE(std::string str);
std::string YELLOW(char c);
std::string YELLOW(std::string str);
std::string PURPLE(char c);
std::string PURPLE(std::string str);

void print_ptrs(std::vector< std::vector<int> > ptrs, 
        std::string alt_str, std::string ref_str);

void print_wfa_ptrs(
        std::vector<std::string> query,
        std::vector<std::string> truth,
        std::vector<int> s,
        std::vector< std::vector< std::vector<int> > > offs, 
        std::vector< std::vector< std::vector<int> > > ptrs);

void write_results(std::unique_ptr<phaseData> & phasings, int dist);

/* --------------------------------------------------------------------------- */

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
