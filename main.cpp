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



uint8_t reverse_compliment(uint8_t v)
{
    const int k = 4;
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
uint8_t merge_map[256];
const int max_kmer = 255;
void build_merge_map(){
    int last_pos = 0;
    for(int i=0;i<max_kmer+1;i++)
    {
        unsigned int rv = reverse_compliment(i);
        if(rv<i)continue;
        merge_map[i]  = last_pos;
        merge_map[rv] = last_pos++;

    }

    cout<<'{';
    for(int i=0;i<256;i++)cout<<int(merge_map[i])<<',';
    cout<<'}';
}


int main() {

    build_merge_map();
    return 0;
}
