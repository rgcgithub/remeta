#ifndef ALLELE_FREQUENCY_MAP_H
#define ALLELE_FREQUENCY_MAP_H

#include <string>
#include <unordered_map>

#include "io/bgz_reader.hpp"

enum is_singleton_e {
  YES,
  NO,
  UNSPECIFIED
};

struct af_info_t {
  std::string cpra;
  std::string chrom;
  int pos;
  double freq;
  is_singleton_e is_singleton;

  bool operator==(const af_info_t& other) const {
    return this->cpra == other.cpra &&
           this->chrom == other.chrom &&
           this->pos == other.pos &&
           this->freq == other.freq &&
           this->is_singleton == other.is_singleton;
  }
};

const af_info_t AF_INFO_NOT_FOUND {"", "", -1, -1, UNSPECIFIED};

/*
 * TODO: This class duplicate code from the AnnotationMap class.
 * Refactor both to eliminate duplication.
 */
class AlleleFrequencyMap {
 public:
  AlleleFrequencyMap();

  // An allele frequency file can either be in a two column format (CPRA FREQ)
  // or a five column format (CPRA FREQ IS_SINGLETON CHROM POS)
  void load(const std::string& filepath);

  void clear();

  bool loaded();

  af_info_t get_af_info(const std::string& cpra);

  static const int BUFFER_SIZE_MB = 1;

 private:
  std::string filepath;
  BgzReader af_file;
  std::unordered_map<std::string, af_info_t> af_buffer;
  std::string buffer_chrom;
  int buffer_start;
  int buffer_end;

  af_info_t next_af_info;
  bool has_next_af_info;

  af_info_t parse_af_file_line(const std::string& line);
  void load_chunk(const std::string& chrom, const int& pos);
  // for unindexed files, load one chromosome at a time
  void load_chrom(const std::string& chrom);
  bool in_current_chunk(const std::string& chrom, const int& pos);
  bool in_next_chunk(const std::string& chrom, const int& pos); 
};

#endif