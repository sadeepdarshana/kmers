//
// Created by Sadeep on 27-Apr.
//

#ifndef PRE_PROCESSOR_VECTORIZER_H
#define PRE_PROCESSOR_VECTORIZER_H

#include <cmath>
#include <stdint.h>
#include <iostream>
#include <random>
#include <ctime>

struct freq_info{
    int start,end;
    int freqs[32] = {0};

    freq_info(int start, int end) : start(start), end(end) {}
};
struct freq_info4{
    int start,end;
    int freqs[136] = {0};

    freq_info4(int start, int end) : start(start), end(end) {}
};
void acgt_to_binary(const char *raw_seq, int n, uint8_t * seq);
struct freq_info kmer_counter(const uint8_t * seq,int start, int end);
void split_n_count(const char *raw_seq, int n,int length_mean,int length_stdev,int count,std::vector<freq_info> * output);
void split_n_count4(const char *raw_seq, int n,int length_mean,int length_stdev,int count,std::vector<freq_info4> * output);

#endif //PRE_PROCESSOR_VECTORIZER_H
