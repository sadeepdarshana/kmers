
#include "vectorizer.h"


const uint8_t kmer3_window_byte1[] = { 252,63,15,3 };
const uint8_t kmer3_window_byte2[] = { 0,0,192,240 };

const uint8_t b1_r_shifts[] = { 2,0,0,0 };
const uint8_t b1_l_shifts[] = { 0,0,2,4 };
const uint8_t b2_r_shifts[] = { 0,0,6,4 };
const uint8_t b2_l_shifts[] = { 0,0,0,0 };


const uint8_t acgt_order[] = { 'A','C','G','T' };

const uint8_t merge_map[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 7, 11, 12, 13, 3, 14, 15, 16, 13, 17, 18, 19, 10, 20, 21,
                             19, 6, 22, 23, 16, 2, 24, 25, 23, 12, 26, 27, 21, 9, 28, 27, 18, 5, 29, 25, 15, 1, 30, 29,
                             22, 11, 31, 28, 20, 8, 31, 26, 17, 4, 30, 24, 14, 0 };

void acgt_to_binary(const char* raw_seq, int n, uint8_t* seq) {
    uint8_t acgt_value[200];
    for (int i = 0; i < 4; i++) acgt_value[acgt_order[i]] = acgt_value[acgt_order[i] + 32] = i;

    for (int i = 0; i < n / 4 + 1; i++) {
        seq[i] = acgt_value[raw_seq[4 * i]] << 6u | acgt_value[raw_seq[4 * i + 1]] << 4u | acgt_value[raw_seq[4 * i + 2]] << 2u |
            acgt_value[raw_seq[4 * i + 3]];
    }

    int i = n / 4;
    seq[i] = 0;

    if (4 * i < n)   seq[i] |= (acgt_value[raw_seq[4 * i]] << 6u);
    if (4 * i + 1 < n) seq[i] |= (acgt_value[raw_seq[4 * i + 1]] << 4u);
    if (4 * i + 2 < n) seq[i] |= (acgt_value[raw_seq[4 * i + 2]] << 2u);
}

freq_info kmer_counter(const uint8_t* seq, int start, int end) { // both ends inclusive
    freq_info fi(start, end);
    uint8_t p = start % 4;

    for (int i = start; i < end - 1; i++) {
        uint8_t kmer = ((seq[i / 4] & kmer3_window_byte1[p]) >> b1_r_shifts[p]) << b1_l_shifts[p] |
            ((seq[i / 4 + 1] & kmer3_window_byte2[p]) >> b2_r_shifts[p]) << b2_l_shifts[p];
        fi.freqs[merge_map[kmer]]++;
        p = (p + 1) % 4;
    }
    return fi;
}

void split_n_count(const char* raw_seq, int n, int length_mean, int length_stdev, int count, std::vector<freq_info>* output) {
    std::default_random_engine generator(time(0));;
    std::normal_distribution<double> distribution(length_mean, length_stdev);

    while (count) {
        int len = round(distribution(generator));
        int start = rand() % n;
        if (start + len >= n || len < 2)continue;

        uint8_t* seq = new uint8_t[n / 4 + 1];
        acgt_to_binary(raw_seq, n, seq);
        freq_info segment_freqs = kmer_counter(seq, start, start + len - 1);
        output->push_back(segment_freqs);
        count--;
    }
}


/* merge_map was created using this code,

const int max_kmer = 63;
 void whatever(){
    for(int i=0;i<max_kmer+1;i++)
    {
    unsigned int rv = reverse_compliment(i);
    if(rv<i)continue;
    merge_map[i]  = last_pos;
    merge_map[rv] = last_pos++;

    }
 }
 about code uses below method
uint8_t reverse_compliment(uint8_t v)
{
    const int k = 3;

    unsigned int comp = 0;
    unsigned int mask = 3;
    for (int i=0;i<k;i++)
    {
        comp += (mask & v)>>2u*i;
        comp <<=2u;
        mask <<=2u;
    }

    return (1u<<(2u*k))-1 - (comp>>2u);
}

*/

/* code for printing a given kmer
 void p_kmer(uint8_t v){

    int k = 3;
    uint8_t comp = 0;
    uint8_t mask = 3;
    char  ch [4]={};
    ch[k] = 0;
    for (int i=0;i<k;i++)
    {
        comp = (mask & v)>>2u*i;
        if     (comp == 0)ch[k-1-i] = 'A';
        else if(comp == 1)ch[k-1-i] = 'C';
        else if(comp == 2)ch[k-1-i] = 'G';
        else if(comp == 3)ch[k-1-i] = 'T';
        mask <<=2u;
    }
}
 */





