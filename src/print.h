#ifndef _PRINT_H_
#define _PRINT_H_

#include <string>

#include "globals.h"
#include "phase.h"
#include "edit.h"
#include "defs.h"

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

void print_cigar(std::vector<int> cigar); 

void print_wfa_ptrs(
        std::vector<std::string> query,
        std::vector<std::string> truth,
        std::vector<int> s,
        std::vector< std::vector< std::vector<int> > > offs, 
        std::vector< std::vector< std::vector<int> > > ptrs);

void write_results(std::unique_ptr<phaseData> & phasings, const editData & edits);

#endif
