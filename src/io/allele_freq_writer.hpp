#ifndef ALLELE_FREQ_WRITER_H
#define ALLELE_FREQ_WRITER_H

#include <string>
#include <unordered_set>

#include "../htpv4_pos.hpp"
#include "bgz_writer.hpp"
#include "htpv4_reader.hpp"

// write allele frequencies used for gene-based tests
class AlleleFreqWriter {
  public:
    AlleleFreqWriter();
    ~AlleleFreqWriter();
    AlleleFreqWriter(const AlleleFreqWriter& other) = delete;
    AlleleFreqWriter& operator=(AlleleFreqWriter other) = delete;
    AlleleFreqWriter(AlleleFreqWriter&& other);

    void open(const std::string& filepath);

    void close();

    bool closed();

    void add_freq(const std::string& vid,
                  const std::string& chr,
                  int pos,
                  double aaf, 
                  bool is_singleton);

    void write_freqs_before(const std::string&, int before_pos);

    void clear_freqs_before(const std::string& chrom, int before_pos);

  private:
    struct AlleleFreq {
      std::string id;
      std::string chrom;
      int pos;
      double aaf;
      bool is_singleton;
      HTPv4Pos htp_pos;
    };

    std::string filepath;
    BgzWriter writer;
    std::unordered_set<std::string> variants_inserted;
    std::vector<AlleleFreq> variants;
};

#endif