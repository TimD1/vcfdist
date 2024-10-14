#ifndef _PRINT_H_
#define _PRINT_H_

#include <string>

#include "globals.h"
#include "phase.h"
#include "defs.h"
#include "dist.h"

std::string GREEN(int i);
std::string GREEN(char c);
std::string GREEN(const std::string & str);
std::string RED(int i);
std::string RED(char c);
std::string RED(const std::string & str);
std::string BLUE(int i);
std::string BLUE(char c);
std::string BLUE(const std::string & str);
std::string YELLOW(int i);
std::string YELLOW(char c);
std::string YELLOW(const std::string & str);
std::string PURPLE(int i);
std::string PURPLE(char c);
std::string PURPLE(const std::string & str);

float qscore(double p_error);
void print_cigar(const std::vector<int> & cigar); 

std::string get_ptr_repr(idx4 cell, const std::unordered_map<idx4,idx4> & ptrs);
void print_graph_ptrs(const std::shared_ptr<Graph> query_graph,
        const std::unordered_map<idx4,idx4> & ptrs);

void print_wfa_ptrs(
        const std::string & query,
        const std::string & truth,
        int s,
        const std::vector< std::vector< std::vector<uint8_t> > > & ptrs,
        const std::vector< std::vector< std::vector<int> > > & offs);

void write_results(std::unique_ptr<phaseblockData> & phasings);
void write_params();

#endif
