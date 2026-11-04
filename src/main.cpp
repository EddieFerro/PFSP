#include <zlib.h>
#include "kseq.h"
#include "parser.hpp"
KSEQ_INIT(gzFile, gzread)

#include "CLI11.hpp"

inline void normalize(const char *src, uint8_t *dst, size_t len) {
    for (size_t i = 0; i < len; i++) dst[i] = NORM.v[(uint8_t)src[i]];
}

static double secs_since(std::chrono::steady_clock::time_point t0) {
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void write_meta(const std::string &path, int K, int S, uint64_t seed, uint64_t n_syncmers, uint32_t n_distinct, 
                const std::vector<RecordInfo> &recs) {
    FILE *f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot open " + path);
    std::fprintf(f, "#format\tpfsp-syncmer-v1\n#K\t%d\n#S\t%d\n#seed\t%llu\n"
                    "#rank_base\t1\n#gap_bytes\t2\n#id_bytes\t4\n"
                    "#n_syncmers\t%llu\n#n_distinct\t%u\n"
                    "#name\tlength\tfirst_entry\tn_entries\n",
                 K, S, (unsigned long long)seed,
                 (unsigned long long)n_syncmers, n_distinct);
    for (const auto &r : recs)
        std::fprintf(f, "%s\t%llu\t%llu\t%llu\n", r.name.c_str(),
                     (unsigned long long)r.length,
                     (unsigned long long)r.first_entry,
                     (unsigned long long)r.n_entries);
    if (std::fclose(f) != 0) throw std::runtime_error("error writing " + path);
}

static int run(int argc, char**argv) {
    CLI::App app{"pfsp: Syncmer Parsing"};

    //Default values
    std::string input = "";
    std::string output = "";
    int16_t K = 200;
    int16_t S = 20;
    uint64_t seed = 0;

    //Options
    app.add_option("input", input, "input FASTA/FASTQ file, optionally gzipped")->required();
    app.add_option("-o,--output", output, "output prefix (default: the input path)");
    app.add_option("-k,--syncmer-length", K, "syncmer length K")->capture_default_str();
    app.add_option("-s,--smer-length", S, "s-mer length S")->capture_default_str();
    app.add_option("--seed", seed, "Seed for hashing steps")->capture_default_str();
    CLI11_PARSE(app, argc, argv);

    //Param checking
    if (output.empty()) output = input;
    if (S < 1 || S >= K) throw std::runtime_error("S must be >= 1 and < K");
    if (S > 21) throw std::runtime_error("S must be <= 21 (3 bits/symbol)");

    //Open input
    gzFile fp = gzopen(input.c_str(), "r");
    if (!fp) throw std::runtime_error("Cannot open file: " + input);
    kseq_t *seq = kseq_init(fp);
    
    pfsp::Parser parser(K, S, seed, output + ".gaps", output + ".parse");


    std::vector<RecordInfo> recs;
    std::vector<uint8_t> buf;
    uint64_t n_skipped = 0;

    auto t0 = std::chrono::steady_clock::now();

    int ret;
    while((ret = kseq_read(seq)) >= 0) {
        const size_t L = seq->seq.l;
        const size_t PL = L + 2 * (size_t)S; //Padded length
        RecordInfo rec{seq->name.s, L, parser.num_syncmers(), 0};

        if(PL < (size_t)K) {
            ++n_skipped;
            recs.push_back(std::move(rec));
            continue;
        }

        if(buf.size() < PL) buf.resize(PL);
        std::memset(buf.data(), 0, S); //Pad with S 0's
        normalize(seq->seq.s, buf.data() + S, L); //Ensure ACGT
        std::memset(buf.data() + S + L, 0, S); //Pad end with S 0's

        parser.process_segment(buf.data(), PL);
        rec.n_entries = parser.num_syncmers() - rec.first_entry;
        recs.push_back(std::move(rec));
    }
    kseq_destroy(seq);
    gzclose(fp);
    if (ret == -2) throw std::runtime_error("truncated FASTQ quality string");
    if (ret < -1)  throw std::runtime_error("error reading " + input);
    std::fprintf(stderr, "[pfsp] parsed %zu records (%llu skipped), %llu syncmers in %.1fs\n",
                    recs.size(), (unsigned long long)n_skipped, 
                    (unsigned long long)parser.num_syncmers(), secs_since(t0));
    
    t0 = std::chrono::steady_clock::now();
    const uint32_t n_distinct = parser.dict_size();
    parser.finalize(output + ".dict", output + ".counts");
    std::fprintf(stderr, "[pfsp] %u distinct syncmers; finalize %.1fs\n", n_distinct, secs_since(t0));

    write_meta(output + ".meta", K, S, seed, parser.num_syncmers(), n_distinct, recs);
    return 0;
}

int main(int argc, char **argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception &e) {
        std::fprintf(stderr, "[pfsp] error: %s\n", e.what());
        return 1;
    }
}
