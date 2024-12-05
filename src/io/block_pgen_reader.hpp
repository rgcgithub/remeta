#ifndef BLOCK_PGEN_READER_H
#define BLOCK_PGEN_READER_H

#include "../lapack_complex.hpp"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
using namespace std;

#include <Eigen/Dense>

#include "pgenlibr.h"

#include "bgz_reader.hpp"

struct pvar_variant_t {
  string chr;
  int pos;
  string id;
  int index;
};

struct psam_sample_t {
  string fid;
  string iid;
  int8_t sex;
};

struct block_pgen_region_t {
  string chr;
  int start_bp;
  int end_bp;
  int start_idx;
  int end_idx;
};

enum read_mode_e {
  REGION,
  ALL
};

class BlockPgenReader {
 public:
  BlockPgenReader() {};

  // PgenReader might not be copyable.
  // These are deleted to avoid possible bugs derived from
  // making a copy of PgenReader.
  BlockPgenReader(const BlockPgenReader&) = delete;

  BlockPgenReader& operator=(const BlockPgenReader&) = delete;

  void load(const std::string& pfile,
            const int& nthreads,
            const int& pvar_index_stride);

  void load(const string& pfile,
            const int& nthreads,
            const int& pvar_index_stride,
            const string& keep_file,
            const string& remove_file);

  void load(const string& pfile,
            const int& nthreads,
            const int& pvar_index_stride,
            const unordered_set<string>& keep_list,
            const unordered_set<string>& remove_list);

  int get_nblocks(const int& block_size);

  int get_nblocks(const int& block_size, const block_pgen_region_t& region);

  int get_nsamples();

  int get_nkeep_samples();

  unordered_set<string> get_keep_list() { return this->keep_list; }

  unordered_set<string> get_remove_list() { return this->remove_list; }

  // if a region has been set, returns number of variants in the region
  // otherwise returns the number of variants in the pvar file
  int get_nvariants();

  void set_region(const string& chr, const int& start, const int& end);

  block_pgen_region_t get_region(const string& chr, const int& start, const int& end);

  void set_extract_list(const unordered_set<string>& variant_ids);

  void set_exclude_list(const unordered_set<string>& variant_ids);

  void set_read_dosages();

  void load_block(Eigen::MatrixXf& genotype_block,
                  vector<pvar_variant_t>& variants,
                  int& nvar_read,
                  const int& block_index,
                  const int& block_size);

  void load_block(Eigen::MatrixXf& genotype_block,
                  vector<pvar_variant_t>& variants,
                  int& nvar_read,
                  const int& block_index,
                  const int& block_size,
                  const block_pgen_region_t& region);

  bool region_contains_variant(const string& vid);

 private:
  PgenReader pgen_reader;
  string pfile;
  vector<pvar_variant_t> variants; // list of variants available to load (modified by set_extract_list and set_exclude_list)
  vector<psam_sample_t> samples;
  unordered_map<string, int> variant_idx; // vid to index in pvar file
  unordered_set<string> keep_list;
  unordered_set<string> remove_list;
  unordered_set<int> extract_idx;
  read_mode_e read_mode = ALL;
  bool read_dosages = false;

  // index from chr -> pos -> index in variants vector
  // the position of the variant at variant_idx[chr][pos]
  // is guaranteed to be >= pos.
  unordered_map<string, unordered_map<int, int> > chr_pos_idx;
  // chr -> [start_bp, end_bp]
  unordered_map<string, pair<int, int> > chr_start_end;
  // store an index every idx_stride bases
  int idx_stride;
  block_pgen_region_t region;

  void load_pvar(const string& pvar);

  void load_psam(const string& psam);

  void load_keep_list(const string& keep_file);

  void load_remove_list(const string& remove_file);

  string format_sample_id(const psam_sample_t& sample);
};

#endif