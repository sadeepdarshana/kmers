#include <iostream>
#include <fstream>
#include "vectorizer.h"
#include "kseq_cpp.h"

using namespace std;
using namespace kseq_cpp;



void p_kmer(uint8_t v) {

    int k = 3;
    uint8_t comp = 0;
    uint8_t mask = 3;
    char ch[4] = {};
    ch[k] = 0;
    for (int i = 0; i < k; i++) {
        comp = (mask & v) >> 2u * i;
        if (comp == 0)ch[k - 1 - i] = 'A';
        else if (comp == 1)ch[k - 1 - i] = 'C';
        else if (comp == 2)ch[k - 1 - i] = 'G';
        else if (comp == 3)ch[k - 1 - i] = 'T';
        mask <<= 2u;
    }
    cout<<ch<<endl;
}

int main() {


    kseq_parser kseq("seq.fasta");
    kseq.read_seq();
    std::cout << kseq.name << std::endl;
    std::cout << kseq.seq << std::endl;

    cout<<kseq.seq[41]<<kseq.seq[42]<<kseq.seq[43]<<kseq.seq[44]<<kseq.seq[45];
    vector<freq_info> list;
    //split_n_count(c,str.length(),50,20,20, &list);

    cout<<"well"<<endl;

    return 0;
}
