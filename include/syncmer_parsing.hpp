#pragma once
#include <stdint.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include "utils.hpp"


typedef struct Syncmers
{
    const int16_t K;
    const int16_t S;
    const int16_t W;

    std::vector<uint64_t> buf;
    int64_t start; //start of rotating buffer
    uint64_t min; //current minimum hash for current k-mer
    uint64_t h; //last packed s-mer for current k-mer
    
    const uint64_t mask;

    const uint64_t seed;
    

    //Constructor, expects first K-mer and populates buffer
    Syncmers(int64_t K, int64_t S, uint64_t seed): 
        K(K), S(S), W(K - S + 1),  
        buf(W), mask((1ULL << (3*S)) - 1),
        seed(seed | (1ULL << 63))
        {}

    //Resets buffer and values according to given K-mer
    void reset(const uint8_t *kmer) {
        start = 0;
        h = 0;
        min = UINT64_MAX;
        int16_t j = 0;
        for(int16_t i = 0; i < K; i++) {
            h = ((h<<3) | (uint64_t)(kmer[i] & 7)) & mask;
            if(i >= S-1) {
                uint64_t hash = smer_hash(h, seed);
                buf[j++] = hash;
                if(hash < min) min = hash;
            }
        }
    }

    //Shifts the current k-mer by one character
    //Returns true if new k-mer is closed syncmer; except start & end, handled by parser
    void shift(uint8_t c) {
        bool rescan = (buf[start] == min); //If current min is shifted out, need new min
        h = ((h<<3) | (uint64_t)(c & 7)) & mask;
        uint64_t hash = smer_hash(h, seed);
        buf[start++] = hash; //Replace start with new hash and increment start position
        if(start == W) start = 0; //Loop start index

        if(rescan) { min = *std::min_element(buf.begin(), buf.end()); }
        else if(hash <= min) { min = hash; }
    } 

    bool is_syncmer() {
        const int64_t end = (start == 0 ? W : start) - 1;
        return (buf[start] == min) || (buf[end] == min); 
    }

}  Syncmers;