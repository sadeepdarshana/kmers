
#include "vectorizer.h"


const uint8_t kmer3_window_byte1[] = {252,63,15,3};
const uint8_t kmer3_window_byte2[] = {0,0,192,240};

const uint8_t b1_r_shifts[] = {2,0,0,0};
const uint8_t b1_l_shifts[] = {0,0,2,4};
const uint8_t b2_r_shifts[] = {0,0,6,4};
const uint8_t b2_l_shifts[] = {0,0,0,0};

int merge_map[5000] = {-1};

const int max_kmer = 63;
uint8_t acgt['t' + 1] = {};

void acgt_to_binary(const char *raw_seq, int n, uint8_t * seq) {
    acgt['a'] = acgt['A'] = 0;
    acgt['c'] = acgt['C'] = 1;
    acgt['g'] = acgt['G'] = 2;
    acgt['t'] = acgt['T'] = 3;

    for(int i=0;i<n/4+1;i++){
        seq[i] = acgt[raw_seq[4 * i]] << 6u | acgt[raw_seq[4 * i + 1]] << 4u | acgt[raw_seq[4 * i + 2]] << 2u |
                 acgt[raw_seq[4 * i + 3]];
    }

    int i = n/4;
    seq[i] = 0;

    if(4*i<n)   seq[i] |= (acgt[raw_seq[4*i]]  <<6u);
    if(4*i+1<n) seq[i] |= (acgt[raw_seq[4*i+1]]<<4u);
    if(4*i+2<n) seq[i] |= (acgt[raw_seq[4*i+2]]<<2u);
}

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

void init(){
    int last_pos = 0;

    for(int i=0;i<max_kmer+1;i++)
    {
        unsigned int rv = reverse_compliment(i);
        if(rv<i)continue;
        merge_map[i]  = last_pos;
        merge_map[rv] = last_pos++;
    }
}



freq_info kmer_counter(const uint8_t * seq,int start, int end){
    init();
    uint8_t p = start % 4;

    freq_info fi(start,end);

    for(int i=start;i<end-1;i++){
        uint8_t kmer = ((seq[i/4] & kmer3_window_byte1[p])>>b1_r_shifts[p])<<b1_l_shifts[p] |
                ((seq[i/4+1] & kmer3_window_byte2[p])>>b2_r_shifts[p])<<b2_l_shifts[p];
        fi.freqs[merge_map[kmer]]++;
        p = (p + 1) % 4;
    }
    return fi;
}












