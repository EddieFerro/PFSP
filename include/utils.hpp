//Utils file
#pragma once
#include <string>
#include <vector>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <cstdint>
#include <numeric>

#include "MurmurHash3.h"

// splitmix64 finalizer; the s-mer hash is mix(packed value)
static inline uint64_t mix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

struct NormTable {
    uint8_t v[256];
    constexpr NormTable() : v{} {
        for (int i = 0; i < 256; ++i) v[i] = 'N';
        v['A'] = v['a'] = 'A';
        v['C'] = v['c'] = 'C';
        v['G'] = v['g'] = 'G';
        v['T'] = v['t'] = 'T';
        v['U'] = v['u'] = 'T';   // drop these two lines if you want U -> N
    }
};
inline constexpr NormTable NORM{};

struct SymbolTable {
    uint8_t v[256];
    constexpr SymbolTable() : v{} {
        for (int i = 0; i < 256; ++i) v[i] = 4;
        v[0x00] = 0;
        v[(uint8_t)'A'] = v[(uint8_t)'a'] = 1;
        v[(uint8_t)'C'] = v[(uint8_t)'c'] = 2;
        v[(uint8_t)'G'] = v[(uint8_t)'g'] = 3;
        v[(uint8_t)'N'] = v[(uint8_t)'n'] = 4;
        v[(uint8_t)'T'] = v[(uint8_t)'t'] = 5;
        v[(uint8_t)'U'] = v[(uint8_t)'u'] = 5;   // U -> T is lossless
    }
};
inline constexpr SymbolTable SYMBOL_TABLE{};
 
inline int symbol_code(uint8_t c) { return SYMBOL_TABLE.v[c]; }

static inline uint64_t murmurhash64(const void *p, size_t len, uint64_t seed) {
    uint64_t h[2];
    MurmurHash3_x64_128(p, (int)len, (uint32_t)(seed ^ (seed >> 32)), h);
    return h[0];
}

static inline uint64_t smer_hash(uint64_t h, uint64_t seed) { return h == 0 ? 0 : mix(h ^ seed); } //Make sure $^s is min s-mer


struct Writer
{
    FILE *parse = nullptr;
    FILE *gaps = nullptr;

    Writer(const std::string &gaps_path, const std::string &parse_path) {
        gaps = std::fopen(gaps_path.c_str(), "wb");
        parse = std::fopen(parse_path.c_str(), "wb");
        if(!gaps || !parse) {
            if(parse) std::fclose(parse);
            if(gaps) std::fclose(gaps);
            throw std::runtime_error("Unable to open necessary files.");
        } 
    }

    ~Writer() {
        if(parse) std::fclose(parse);
        if(gaps) std::fclose(gaps);
        parse = gaps = nullptr;
    }

    void write(uint32_t id, uint16_t gap) {
        if(std::fwrite(&id, sizeof id, 1, parse) != 1) throw std::runtime_error("Write error on temporary parse.");
        if(std::fwrite(&gap, sizeof gap, 1, gaps) != 1) throw std::runtime_error("Write error on temporary parse.");
    }

    void close() {
        bool bad = std::ferror(parse) || std::ferror(gaps);
        bad |= std::fclose(parse) != 0;  parse = nullptr;
        bad |= std::fclose(gaps)  != 0;  gaps  = nullptr;
        if (bad) throw std::runtime_error("Error writing parse/gaps output.");
    }

};

struct RecordInfo {
    std::string name;
    uint64_t length;        // original length, Ns included
    uint64_t first_entry;   // index of its first parse entry
    uint64_t n_entries;     // 0 => skipped
};

// K-mer into  PB_ packed bytes at out
// Symbols go in MSB-first; unused low bits of the last byte are zeroed, so two equal k-mers pack to equal bytes.
// Duplicate of pack function in dictionary.hpp but exposed for readers (SIFT)
inline void pack(const uint8_t *seq, uint8_t *out, int K, int bits = 3, int g_bits = std::lcm(3, 8)) {          
    const int SG = g_bits / bits; // symbols per group
    const int BG = g_bits / 8; // bytes per group
    int o = 0, i = 0;

    for (; i + SG <= K; i += SG) {
        uint64_t acc = 0;
        for (int j = 0; j < SG; ++j)
            acc = (acc << bits) | (uint64_t)symbol_code(seq[i + j]);
        for (int b = BG - 1; b >= 0; --b) out[o++] = (uint8_t)(acc >> (8 * b));
    }

    // Fewer than SG symbols remain, fall back to byte by byte.
    uint64_t acc = 0;
    int nbits = 0;
    for (; i < K; ++i) {
        acc = (acc << bits) | (uint64_t)symbol_code(seq[i]);
        nbits += bits;
        while (nbits >= 8) { nbits -= 8; out[o++] = (uint8_t)(acc >> nbits); }
    }
    if (nbits) out[o++] = (uint8_t)(acc << (8 - nbits));   // left-align, zero-pad
}


