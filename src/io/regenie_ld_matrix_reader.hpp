#ifndef REGENIE_LD_MATRIX_READER_H
#define REGENIE_LD_MATRIX_READER_H

#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include <Eigen/Dense>

class RegenieLDMatrixReader {
 public:
  RegenieLDMatrixReader();
  RegenieLDMatrixReader(std::string mat_file, std::string idx_file);
  RegenieLDMatrixReader(const RegenieLDMatrixReader& other);
  RegenieLDMatrixReader& operator=(RegenieLDMatrixReader other);
  ~RegenieLDMatrixReader();

  void open(std::string mat_file, std::string idx_file);

  void close();

  bool closed() { return is_closed; }

  // Load LD matrix of all variants in gene_name and populate
  // the vector variant ids.
  void load_ld_mat(Eigen::MatrixXf& mat,
                   std::vector<std::string>& variant_ids,
                   const std::string& gene_name);

  // Load LD matrix of gene_variants from gene name. This function
  // matches the signature of the load_gene_ld_mat function of
  // RefLDMatrixReader.
  void load_gene_ld_mat(Eigen::MatrixXf& G,
                        const std::vector<std::string>& gene_variants,
                        const std::string& gene_name);

  void load_variant_ids(std::vector<std::string>& variant_ids, const string& gene_name);

  int32_t get_sample_size() { return this->sample_size; }

  std::vector<std::string> get_gene_names();

  bool contains_gene(const std::string& gene_name);

  friend void swap(RegenieLDMatrixReader& first, RegenieLDMatrixReader& second);

 private:
  void load_idx();
  void read_sample_size();
  void load_ld_mat_dense(Eigen::MatrixXf& mat);
  void load_ld_mat_sparse(Eigen::MatrixXf& mat);

  bool is_closed;

  string mat_file;
  string idx_file;
  int32_t sample_size;

  BGZF* mat_bgzf;
  BGZF* idx_bgzf;
  kstring_t buffer;

  // pointer to position of a gene in the index
  map<string, int64_t> idx_addr;
  // pointer to position of a gene in the ld matrix file.
  map<string, int64_t> mat_addr; 
};

#endif