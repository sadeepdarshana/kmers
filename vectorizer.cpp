
#include "vectorizer.h"

const uint8_t acgt_order[] = { 'A','C','G','T' };

//3mer const
const uint8_t kmer3_window_byte1[] = { 252,63,15,3 };
const uint8_t kmer3_window_byte2[] = { 0,0,192,240 };

const uint8_t b1_r_shifts[] = { 2,0,0,0 };
const uint8_t b1_l_shifts[] = { 0,0,2,4 };
const uint8_t b2_r_shifts[] = { 0,0,6,4 };
const uint8_t b2_l_shifts[] = { 0,0,0,0 };

const uint8_t merge_map[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 7, 11, 12, 13, 3, 14, 15, 16, 13, 17, 18, 19, 10, 20, 21,
                             19, 6, 22, 23, 16, 2, 24, 25, 23, 12, 26, 27, 21, 9, 28, 27, 18, 5, 29, 25, 15, 1, 30, 29,
                             22, 11, 31, 28, 20, 8, 31, 26, 17, 4, 30, 24, 14, 0 };


//4mer const
const uint8_t kmer3_window_byte14[] = { 255,63,15,3 };
const uint8_t kmer3_window_byte24[] = { 0,192,240,255 };

const uint8_t b1_r_shifts4[] = { 0,0,0,0 };
const uint8_t b1_l_shifts4[] = { 0,2,4,6 };

const uint8_t b2_r_shifts4[] = { 0,6,4,2 };
const uint8_t b2_l_shifts4[] = { 0,0,0,0 };

const uint8_t merge_map4[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,11,31,32,33,34,35,36,37,38,39,40,41,
                              23,42,43,44,7,45,46,47,48,49,50,51,34,52,53,54,19,55,56,57,3,58,59,60,57,61,62,63,44,64,65,66,30,67,68,69,14,70,71,72,54
        ,73,74,75,41,76,77,78,26,79,80,66,10,81,82,83,51,84,85,86,37,87,88,75,22,89,90,63,6,91,92,93,47,94,95,83,33,96,97,72,18,
                              98,99,60,2,100,101,99,56,102,103,90,43,104,105,80,29,106,107,68,13,108,109,97,53,110,111,88,40,112,113,77,25,114,105,65,
                              9,115,116,95,50,117,118,85,36,119,111,74,21,120,103,62,5,121,122,92,46,123,116,82,32,124,109,71,17,125,101,59,1,126,125,
                              98,55,127,120,89,42,128,114,79,28,129,106,67,12,130,124,96,52,131,119,87,39,132,112,76,24,128,104,64,8,133,123,94,49,134
        ,117,84,35,131,110,73,20,127,102,61,4,135,121,91,45,133,115,81,31,130,108,70,16,126,100,58,0};

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
freq_info4 kmer_counter4(const uint8_t* seq, int start, int end) { // both ends inclusive
    freq_info4 fi(start, end);
    uint8_t p = start % 4;

    for (int i = start; i < end - 1; i++) {
        uint8_t kmer = ((seq[i / 4] & kmer3_window_byte14[p]) >> b1_r_shifts4[p]) << b1_l_shifts4[p] |
                       ((seq[i / 4 + 1] & kmer3_window_byte24[p]) >> b2_r_shifts4[p]) << b2_l_shifts4[p];
        fi.freqs[merge_map4[kmer]]++;
        p = (p + 1) % 4;
    }
    return fi;
}

void split_n_count(const char* raw_seq, int n, int length_mean, int length_stdev, int count, std::vector<freq_info>* output) {

    std::random_device rd;
    std::mt19937 mt(rd());
    srand(time(NULL));
    std::default_random_engine generator(rd());;
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


void split_n_count4(const char* raw_seq, int n, int length_mean, int length_stdev, int count, std::vector<freq_info4>* output) {

    std::random_device rd;
    std::mt19937 mt(rd());
    srand(time(NULL));
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(length_mean, length_stdev);

    while (count) {
        int len = round(distribution(generator));
        int start = rand() % n;
        if (start + len >= n || len < 4)continue;

        uint8_t* seq = new uint8_t[n / 4 + 1];
        acgt_to_binary(raw_seq, n, seq);
        freq_info4 segment_freqs = kmer_counter4(seq, start, start + len - 1);
        output->push_back(segment_freqs);
        count--;
    }
}


/* merge_map was created using this code,
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
uint8_t merge_map[64];
const int max_kmer = 63;
void build_merge_map(){
    int last_pos = 0;
    for(int i=0;i<max_kmer+1;i++)
    {
        unsigned int rv = reverse_compliment(i);
        if(rv<i)continue;
        merge_map[i]  = last_pos;
        merge_map[rv] = last_pos++;

    }
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





