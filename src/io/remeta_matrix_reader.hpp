/*
  Wrapper class around matrix readers for different storage formats.
  See ref_ld_matrix_reader.hpp for a description of G, C, and G_C.

  There are two supported file formats for LD matrices:
    * Reference LD matrices. These are technically covariance matrices
      of variants in a gene, and in a buffer region around a gene.
    * Regenie LD matrices. These are technically the phenotypic covariance
      adjusted covariance matrices used in the SKAT test.

  Reference matrices support load_gene_ld_mat and load_conditional_ld_mats,
  Regenie matrices only support load_gene_ld_mat (load_conditional_ld_mats
  will return empty matrices for C and G_C).
*/
#ifndef REMETA_MATRIX_READER_H
#define REMETA_MATRIX_READER_H

#include <string>
#include <vector>

#include <Eigen/Dense>

#include "ref_ld_matrix_reader.hpp"
#include "regenie_ld_matrix_reader.hpp"

class RemetaMatrixReader {
 public:
  RemetaMatrixReader();

  void open(const std::string& file_prefix);

  void close();

  bool closed();

  // only the RefLDMatrixReader class has additional variants for
  // conditional analysis.
  bool has_conditional_variants();

  void load_gene_ld_mat(Eigen::MatrixXf& G,
                        const std::vector<std::string>& gene_variants,
                        const std::string& gene_name);

  void load_conditional_ld_mats(Eigen::MatrixXf& G,
                                Eigen::MatrixXf& C,
                                Eigen::MatrixXf& G_C,
                                const std::vector<std::string>& gene_variants,
                                const std::vector<std::string>& conditional_variants,
                                const std::string& gene_name);
  
  bool contains_gene(const std::string& gene_name);

  /*
    Gives the list of variants in gene_name in the gene LD file
  */ 
  void load_gene_variant_ids(std::vector<std::string>& variant_ids,
                             const std::string& gene_name);

  /*
    Gives the list of variants in gene_name in the buffer LD file
  */
  void load_buffer_variant_ids(std::vector<std::string>& variant_ids,
                               const std::string& gene_name);

 private:
  enum reader_opened_e {
    REF_READER,
    REGENIE_READER,
    NONE
  };

  RefLDMatrixReader ref_reader;
  RegenieLDMatrixReader regenie_reader;
  reader_opened_e reader_opened;
  bool is_closed;

  void check_closed();
};

#endif