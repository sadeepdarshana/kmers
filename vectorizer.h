//
// Created by Sadeep on 27-Apr.
//

#ifndef PRE_PROCESSOR_VECTORIZER_H
#define PRE_PROCESSOR_VECTORIZER_H


#include <stdint.h>
#include <iostream>

struct freq_info{
    int start,end;
    int freqs[32] = {0};

    freq_info(int start, int end) : start(start), end(end) {}
};
void acgt_to_binary(const char *raw_seq, int n, uint8_t * seq);
struct freq_info kmer_counter(const uint8_t * seq,int start, int end);


#endif //PRE_PROCESSOR_VECTORIZER_H
