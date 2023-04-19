#ifndef _FASTA_H_
#define _FASTA_H_

#include <string>
#include <unordered_map>

// zlib is required for kseq
#include "zlib.h"
#include "htslib/kseq.h"
KSEQ_INIT(int, read);

class fastaData {
public:
    fastaData(FILE * ref_fasta_fp) {
        kseq_t * seq = kseq_init(fileno(ref_fasta_fp));
        while (kseq_read(seq) >= 0) {
            this->fasta[seq->name.s] = seq->seq.s;
            this->lengths[seq->name.s] = this->fasta.at(seq->name.s).size();
        }
        kseq_destroy(seq);
        fclose(ref_fasta_fp);
    }
    
    std::unordered_map<std::string,std::string> fasta;
    std::unordered_map<std::string,int> lengths;
};

#endif
