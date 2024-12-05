#ifndef REF_LD_MATRIX_WRITER_H
#define REF_LD_MATRIX_WRITER_H

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include "bgz_writer.hpp"

typedef std::string variant_id;

struct ref_ld_matrix_entry_t {
  variant_id row_id;
  variant_id col_id;
  float corr;
};

struct block_block_ld_t {
  std::vector<ref_ld_matrix_entry_t> data;
};

const int32_t     BLOCK_START_MARKER = -1;
const int32_t     ROW_START_MARKER   = -2;
const int32_t     BLOCK_END_MARKER   = -3;
const float       pow_2_m7  = pow(2.f, -7.f);
const float       pow_2_m8  = pow(2.f, -8.f);
const float       pow_2_m15 = pow(2.f, -15.f);

class RefLDMatrixWriter {
 public:
  RefLDMatrixWriter(const std::string& file_prefix, const int& float_size);

  ~RefLDMatrixWriter();

  // Non-copyable
  RefLDMatrixWriter(const RefLDMatrixWriter& other) = delete;
  RefLDMatrixWriter& operator=(RefLDMatrixWriter other) = delete;

  void init_next_gene(const std::string& gene_name,
                      const float& gene_sparsity_threshold,
                      const float& buffer_sparsity_threshold,
                      const std::vector<variant_id>& gene_variant_ids,
                      const std::vector<float>& gene_variant_variances,
                      const std::vector<variant_id>& buffer_variant_ids,
                      const std::vector<float>& buffer_variant_variances,
                      const int32_t& nbuffer_blocks);

  void add_gene_gene_ld_entry(const ref_ld_matrix_entry_t& data);

  void add_gene_buffer_ld_entry(const ref_ld_matrix_entry_t& data);

  void write_gene_to_gene_ld_file();

  void write_gene_to_index(const std::string& gene_name,
                           const int64_t& first_block_block_ld_addr);

  // returns an address written
  int64_t write_block_block_ld_to_buffer_ld_file(const block_block_ld_t& block_block_ld,
                                                 const bool& within_block_ld);

  void close();

  bool closed() { return this->is_closed; }

 private:
  struct matrix_entry_data_t {
    int32_t row_idx;
    int32_t col_idx;
    float corr;
  };

  struct matrix_row_entry_data_t {
    int32_t col_idx;
    float corr;
  };

  struct gene_gene_ld_t {
    int32_t nvariants;
    int32_t nentries;
    float sparsity_threshold;
    std::vector<float> variances;
    std::vector<matrix_entry_data_t> data;
  };

  struct gene_buffer_ld_t {
    int32_t ngene_variants;
    int32_t nbuffer_variants;
    int32_t nentries;
    float sparisty_threshold;
    // nbuffer variants
    std::vector<float> variances;
    // entry buffer_idx_map[i] maps an index in buffer_variant_ids
    // to a variant index in stored in block_block_ld_t
    std::vector<int32_t> buffer_idx_map;
    int32_t nbuffer_blocks;
    std::vector<matrix_entry_data_t> data;
  };

  struct idx_entry_t {
    std::string gene;
    int64_t gene_ld_addr;
    int64_t buffer_ld_addr;
    std::vector<variant_id> gene_variants;
    std::vector<variant_id> buffer_variants;
  };

  void write_gene_gene_ld();
  void write_gene_buffer_ld();
  //void flush_block_block_ld_buffer();

  bool is_closed;
  void check_closed();

  BgzWriter gene_ld_file;
  BgzWriter buffer_ld_file;
  BgzWriter idx;
  int32_t float_size;

  std::unordered_map<variant_id, int32_t> buffer_ld_variant_idx;
  std::unordered_map<variant_id, int32_t> gene_variant_idx;
  std::unordered_map<variant_id, int32_t> buffer_variant_idx;
  gene_gene_ld_t gene_gene_ld;
  gene_buffer_ld_t gene_buffer_ld;
  std::unordered_map<std::string, idx_entry_t> idx_entries;

  //std::vector<matrix_entry_data_t> block_block_ld_buffer;
  //std::unordered_map<int32_t, std::vector<matrix_row_entry_data_t> > block_block_ld_buffer;
};

#endif