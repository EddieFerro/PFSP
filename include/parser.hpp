#pragma once

#include "utils.hpp"
#include "syncmer_parsing.hpp"
#include "dictionary.hpp"


namespace pfsp 
{

class Parser 
{
private:
    int64_t K_, S_;
    Syncmers syncmer_;
    Dictionary dict_;
    uint64_t n_syncmers_;
    Writer writer_;
    std::string parse_path_;
    std::string tmp_parse_path_;

    void remap_parse() {
        FILE *tmp = std::fopen(tmp_parse_path_.c_str(), "rb");
        if(!tmp) throw std::runtime_error("Unable to open temporary parse file: " + tmp_parse_path_);
        FILE *final = std::fopen(parse_path_.c_str(), "wb");
        if(!final) throw std::runtime_error("Unable to open final parse file: " + parse_path_);
        bool bad = false;
        size_t CHUNK = 1u << 20;
        std::vector<uint32_t> buf(CHUNK);
        uint64_t total = 0;

        size_t n;
        while((n = std::fread(buf.data(), sizeof(uint32_t), CHUNK, tmp)) > 0) {
            for(size_t i = 0; i < n; i++) {
                buf[i] = dict_.rank_of(buf[i]) + 1;
                
            }
            if(std::fwrite(buf.data(), sizeof(uint32_t), n, final) != n) throw std::runtime_error("Error remapping.");
            total += n;
        
        }
        bad |= std::ferror(tmp) != 0;
        std::fclose(tmp);
        bad |= std::fclose(final) != 0;
        if (bad) throw std::runtime_error("Error remapping parse.");
        if (total != n_syncmers_) throw std::runtime_error("Parse length != syncmers added.");
        std::remove(tmp_parse_path_.c_str());
    }


public:
    
    Parser(int64_t K, int64_t S, uint64_t seed, const std::string &gaps_path, const std::string &parse_path): 
            K_(K), S_(S), syncmer_(K, S, seed), dict_(K), n_syncmers_(0), writer_(gaps_path, parse_path + ".tmp"),
            parse_path_(parse_path), tmp_parse_path_(parse_path + ".tmp") {}


    void process_segment(const uint8_t *seq, size_t L) {

        if(L < (size_t) K_) {return;} //Plan for this to never happen during padding in the caller.

        syncmer_.reset(seq); // This automatically process the first K characters of the sequence.
        if(!syncmer_.is_syncmer()) throw std::runtime_error("First K-mer should always be a closed syncmer.");

        uint32_t id = dict_.add(seq);
        n_syncmers_++;
        size_t prev = 0;
        writer_.write(id, 0);

        for(size_t i = (size_t) K_; i < L; i++) {
            const size_t p = i - K_ + 1;
            syncmer_.shift(seq[i]);
            if(syncmer_.is_syncmer()) {
                id = dict_.add(seq + p);
                uint16_t gap = p - prev;
                prev = p;
                n_syncmers_++;
                writer_.write(id, gap);
            }
        }
        if (prev != L - (size_t) K_) throw std::runtime_error("Last K-mer should always be a closed syncmer.");
    }

    void finalize(const std::string &dict_path, const std::string &counts_path) {
        writer_.close();
        dict_.sort();
        dict_.release_table();
        dict_.write(dict_path, counts_path);
        dict_.free_order();
        //Free blocks here maybe
        remap_parse();
    }
    

    uint64_t num_syncmers(){ return n_syncmers_; }

    uint32_t dict_size() { return dict_.size(); }

    
};

}