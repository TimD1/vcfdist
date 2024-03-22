#ifndef _PRINT_H_
#define _PRINT_H_

#include <string>

#include "globals.h"
#include "phase.h"
#include "defs.h"

std::string GREEN(int i);
std::string GREEN(char c);
std::string GREEN(std::string str);
std::string RED(int i);
std::string RED(char c);
std::string RED(std::string str);
std::string BLUE(int i);
std::string BLUE(char c);
std::string BLUE(std::string str);
std::string YELLOW(int i);
std::string YELLOW(char c);
std::string YELLOW(std::string str);
std::string PURPLE(int i);
std::string PURPLE(char c);
std::string PURPLE(std::string str);

float qscore(double p_error);
void print_ref_ptrs(std::vector< std::vector<int> > ptrs);
void print_ptrs(const std::vector< std::vector<uint8_t> > & ptrs, 
        const std::string & alt_str, const std::string & ref_str);

void print_cigar(std::vector<int> cigar); 

void print_wfa_ptrs(
        const std::string & query,
        const std::string & truth,
        int s,
        const std::vector< std::vector< std::vector<uint8_t> > > & ptrs,
        const std::vector< std::vector< std::vector<int> > > & offs);

void write_results(std::unique_ptr<phaseblockData> & phasings);
void write_params();

#endif
