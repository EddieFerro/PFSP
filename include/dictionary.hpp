#pragma once


// dictionary.hpp — PFSP syncmer dictionary.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "utils.hpp"

namespace pfsp {

class Dictionary {

private:
    int K_;
    int bits_;// bits per symbol
    int g_bits_; // bits per write group
    int64_t  PB_; // ceil(bits_*K_/8): fixed width of every record; the stride
    uint64_t hash_seed_;
    bool sorted_ = false; // set by sort(), also for an empty dictionary 
        
    
    // Reusable PB_-byte buffer.
    std::vector<uint8_t> scratch_;

    // Arena. one owning pointer per block
    std::vector<std::unique_ptr<uint8_t[]>> blocks_;

    // counts_[id] = occurrences in the text. counts_.size() is the next id
    std::vector<uint32_t> counts_;

    // Lookup: mix(bytes) -> construction id.
    std::unordered_map<uint64_t, uint32_t> index_;

    // Filled by sort()
    std::vector<uint32_t> order_;   // rank -> construction id
    std::vector<uint32_t> rank_;    // construction id -> rank

public:
    // Records per block. Must be a power of two so that the two-level address split below compiles to a shift and a mask.
    static constexpr uint32_t BLOCK_RECORDS = 1u << 16;

    static constexpr int64_t RADIX_CUTOFF = 32;

    // Default seed for the record fingerprint.
    static constexpr uint64_t DEFAULT_HASH_SEED = 0x51ed270b27a4d1f3ULL;

    // bits: bits per symbol. 3 for the {$,A,C,G,N,T} alphabet
    explicit Dictionary(int K, int bits = 3, uint64_t hash_seed = DEFAULT_HASH_SEED)
        : K_(K), bits_(bits), PB_((int64_t)((bits * (int64_t)K + 7) / 8)), hash_seed_(hash_seed) {
        if (K <= 0) throw std::runtime_error("Dictionary: K must be positive");
        if (bits <= 0 || bits > 8) throw std::runtime_error("Dictionary: bits must be in 1..8");
        g_bits_ = (int)std::lcm(bits_, 8);
        index_.max_load_factor(0.7f);   // below the 1.0 default: fewer probe chains
    }

    uint32_t size() const {
        return (uint32_t)counts_.size();
    }

    bool sorted() const { 
        return sorted_; 
    }

    // id -> bytes. Two-level address: block, slot
    const uint8_t *record(uint32_t id) const {
        return blocks_[id / BLOCK_RECORDS].get() + (size_t)(id & (BLOCK_RECORDS - 1)) * (size_t)PB_;
    }

    uint32_t count(uint32_t id) const { return counts_[id]; }

    // Construction id of `rec`, or -1
    int64_t find(const uint8_t *rec) const {
        auto it = index_.find(hash_record(rec));
        if (it == index_.end()) return -1;
        if (std::memcmp(record(it->second), rec, (size_t)PB_) != 0) return -1;
        return (int64_t)it->second;
    }
 
    // K-mer into  PB_ packed bytes at out
    // Symbols go in MSB-first; unused low bits of the last byte are zeroed, so two equal k-mers pack to equal bytes.
    void pack(const uint8_t *seq, uint8_t *out) const {          
        const int SG = g_bits_ / bits_; // symbols per group
        const int BG = g_bits_ / 8; // bytes per group
        int o = 0, i = 0;
 
        for (; i + SG <= K_; i += SG) {
            uint64_t acc = 0;
            for (int j = 0; j < SG; ++j)
                acc = (acc << bits_) | (uint64_t)symbol_code(seq[i + j]);
            for (int b = BG - 1; b >= 0; --b) out[o++] = (uint8_t)(acc >> (8 * b));
        }
 
        // Fewer than SG symbols remain, fall back to byte by byte.
        uint64_t acc = 0;
        int nbits = 0;
        for (; i < K_; ++i) {
            acc = (acc << bits_) | (uint64_t)symbol_code(seq[i]);
            nbits += bits_;
            while (nbits >= 8) { nbits -= 8; out[o++] = (uint8_t)(acc >> nbits); }
        }
        if (nbits) out[o++] = (uint8_t)(acc << (8 - nbits));   // left-align, zero-pad
    }


    // Insert syncmer if new, otherwise increment count. Returns id
    uint32_t add(const uint8_t *seq) {
        if (sorted()) throw std::runtime_error("Dictionary add after sort()");
        if (scratch_.size() != (size_t)PB_) scratch_.resize((size_t)PB_);
        
        pack(seq, scratch_.data());

        const uint8_t *rec = scratch_.data();

        const uint64_t h = hash_record(rec);
        auto it = index_.find(h);
        if (it != index_.end()) {
            if (std::memcmp(record(it->second), rec, (size_t)PB_) == 0) {
                ++counts_[it->second];
                return it->second;
            }
            throw std::runtime_error("Hash Collision detected!");
        }

        // New record; id is number of unique records
        if (counts_.size() >= UINT32_MAX)
            throw std::runtime_error("Dictionary: more than 2^32-2 distinct records");
        const uint32_t id = (uint32_t)counts_.size();

        // Grow the arena one block at a time
        if (id / BLOCK_RECORDS == blocks_.size())
            blocks_.emplace_back(new uint8_t[(size_t)BLOCK_RECORDS * (size_t)PB_]);

        // Put bytes into block, increment count, add id to map
        std::memcpy(slot(id), rec, (size_t)PB_);
        counts_.push_back(1);
        index_.emplace(h, id);
        return id;
    }

    void reserve(size_t distinct) { index_.reserve(distinct); }


    // MSD byte-wise radix sort
    void sort() {
        const int64_t n = (int64_t)size();
        order_.resize((size_t)n);
        std::iota(order_.begin(), order_.end(), 0u);
 
        if (n > 1) {
            std::vector<uint32_t> tmp((size_t)n);
            struct Range { int64_t lo, hi, byte; };
            std::vector<Range> work;
            work.push_back({0, n, 0});
 
            while (!work.empty()) {
                const Range r = work.back();
                work.pop_back();
 
                if (r.hi - r.lo <= RADIX_CUTOFF || r.byte == PB_) {
                    std::sort(order_.begin() + r.lo, order_.begin() + r.hi,
                              [this](uint32_t a, uint32_t b) {
                                  return std::memcmp(record(a), record(b), (size_t)PB_) < 0;
                              });
                    continue;
                }

                uint32_t count[256] = {};
                for (int64_t i = r.lo; i < r.hi; ++i)
                    ++count[record(order_[(size_t)i])[r.byte]];
 
                int64_t offs[256], acc = r.lo;
                for (int c = 0; c < 256; ++c) { offs[c] = acc; acc += count[c]; }
 
                for (int64_t i = r.lo; i < r.hi; ++i) {
                    const uint32_t id = order_[(size_t)i];
                    tmp[(size_t)offs[record(id)[r.byte]]++] = id;
                }
                std::memcpy(order_.data() + r.lo, tmp.data() + r.lo, (size_t)(r.hi - r.lo) * sizeof(uint32_t));

                int64_t lo = r.lo;
                for (int c = 0; c < 256; ++c) {
                    if (count[c] > 1)
                        work.push_back({lo, lo + (int64_t)count[c], r.byte + 1});
                    lo += (int64_t)count[c];
                }
            }
        }
 
        rank_.resize((size_t)n);
        for (int64_t r = 0; r < n; ++r) rank_[(size_t)order_[(size_t)r]] = (uint32_t)r;
        sorted_ = true;
    }


    uint32_t rank_of(uint32_t id) const {
        if (!sorted()) throw std::runtime_error("rank_of before sort()");
        return rank_[id];
    }
    uint32_t id_of_rank(uint32_t r) const {
        if (order_.empty()) throw std::runtime_error("id_of_rank before sort(), or after free_order()");
        return order_[r];
    }

    const uint8_t *record_by_rank(uint32_t r) const { return record(id_of_rank(r)); }

    // Write .dict and .counts in lex order
    void write(const std::string &dict_path, const std::string &counts_path) const {
        if (!sorted()) throw std::runtime_error("write before sort()");

        FILE *fd = std::fopen(dict_path.c_str(), "wb");
        if (!fd) throw std::runtime_error("cannot open " + dict_path);
        FILE *fc = std::fopen(counts_path.c_str(), "wb");
        if (!fc) { std::fclose(fd); throw std::runtime_error("cannot open " + counts_path); }

        // Batch into one buffer so the kernel sees a few large writes
        constexpr uint32_t BATCH = 4096;
        std::vector<uint8_t>  dbuf((size_t)BATCH * (size_t)PB_);
        std::vector<uint32_t> cbuf(BATCH);

        const uint32_t n = size();
        for (uint32_t base = 0; base < n; base += BATCH) {
            const uint32_t m = std::min(BATCH, n - base);
            for (uint32_t j = 0; j < m; ++j) {
                const uint32_t id = order_[base + j];
                std::memcpy(&dbuf[(size_t)j * (size_t)PB_], record(id), (size_t)PB_);
                cbuf[j] = counts_[id];
            }
            if (std::fwrite(dbuf.data(), (size_t)PB_, m, fd) != m ||
                std::fwrite(cbuf.data(), sizeof(uint32_t), m, fc) != m) {
                std::fclose(fd); std::fclose(fc);
                throw std::runtime_error("short write to " + dict_path);
            }
        }
        if (std::fclose(fd) != 0 || std::fclose(fc) != 0)
            throw std::runtime_error("error closing dictionary output");
    }


    // index is unneeded once a lex ordering is determined.
    void release_table() { index_ = std::unordered_map<uint64_t, uint32_t>(); }

    // order_ is only needed to write the file and to go rank to bytes
    void free_order() { order_.clear(); order_.shrink_to_fit(); }

private:

    // Mutable slot pointer
    uint8_t *slot(uint32_t id) {
        return blocks_[id / BLOCK_RECORDS].get() + (size_t)(id & (BLOCK_RECORDS - 1)) * (size_t)PB_;
    }
    uint64_t hash_record(const uint8_t *r) const {
        return murmurhash64(r, (size_t)PB_, hash_seed_);
    }


};

}  // namespace pfsp
