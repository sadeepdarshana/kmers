#include <iostream>
#include <random>
#include <time.h>
#include <fstream>
#include "vectorizer.h"

using namespace std;

vector<freq_info>fi;

void split_n_count(const char *raw_seq, int n,int length_mean,int length_stdev,int count){


    std::default_random_engine generator(time(0));;
    std::normal_distribution<double> distribution(length_mean,length_stdev);


    while(count){
        int len =  (int)distribution(generator);
        int start = rand() % n;
        if(start+len>=n||len<2)continue;

        uint8_t seq[n/4+1];
        acgt_to_binary(raw_seq, n, seq);
        freq_info segment_freqs = kmer_counter(seq, start,start+len);
        fi.push_back(segment_freqs);
        count--;
    }
}

int main() {
    std::string str;
    ifstream file ("seq.txt");
    std::getline(file,str);
    const char* c = str.c_str();

    split_n_count(c,str.length(),3,1,20);

    int n;
    n++;
    return 0;
}
