#ifndef REF_LD_MATRIX_READER_H
#define REF_LD_MATRIX_READER_H

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "bgz_reader.hpp"

/*
  Suppose we want to load an LD matrix containing variants in a gene
  and additional variants we want to condition on for conditional analysis.
  
  We can represent the covariance matrix of these variants:
     [   G     G_C  ] 
     [  G_C^T   C   ]
  where:
    * G is the covariance matrix of of variants in a gene
    * C is the covariance matrix of conditional variants
    * G_C is the cross-covariance matrix
  
  RefLDMatrixReader provides functions to load these three matrices into memory.
*/
class RefLDMatrixReader {
 public:
  RefLDMatrixReader();
  RefLDMatrixReader(const std::string& file_prefix);

  void open(const std::string& file_prefix);

  void close();

  bool closed();

  /*
    Loads the matrix G: the covariance of gene_variants
  */
  void load_gene_ld_mat(Eigen::MatrixXf& G,
                        const std::vector<string>& gene_variants,
                        const std::string& gene_name);

  void load_gene_ld_mat_sp_double(Eigen::SparseMatrix<double>& G_sp,
                                  const vector<string>& gene_variants,
                                  const string& gene_name);

  /*
    Loads the matrices for gene_name:
      * G   - the covariance matrix of variants in the gene
      * C   - the covariance of conditional_variants
      * G_C - the cross-covariance matrix of gene_variants (rows) 
              and conditional_variants (columns)

    Conditional variants can either be gene variants or buffer variants.
  */
  void load_conditional_ld_mats(Eigen::MatrixXf& G,
                                Eigen::MatrixXf& C,
                                Eigen::MatrixXf& G_C,
                                const std::vector<std::string>& gene_variants,
                                const std::vector<std::string>& conditional_variants,
                                const std::string& gene_name);

  void load_conditional_ld_mats_sp_double(Eigen::SparseMatrix<double>& G_sp,
                                          Eigen::MatrixXf& C,
                                          Eigen::MatrixXf& G_C,
                                          const std::vector<std::string>& gene_variants,
                                          const std::vector<std::string>& conditional_variants,
                                          const std::string& gene_name);

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

  bool contains_gene(const std::string& gene_name);

 private:
  BgzReader gene_ld_file;
  BgzReader buffer_ld_file;
  BgzReader index_file;
  bool has_buffer_ld_file;
  std::string file_prefix;
  bool is_closed;
  int float_size;

  struct gene_idx_t {
    int64_t gene_ld_index;
    int64_t buffer_ld_index;
    int64_t index_index;
  };

  struct gene_ld_mat_t {
    Eigen::MatrixXf* G;
    Eigen::SparseMatrix<double>* G_sp;
    // maps variant ids to the index of the stored LD matrix
    std::unordered_map<std::string, int32_t> gene_variant_idx_in_file;
    // maps variant ids to index in G (which can be smaller that what's stored)
    std::unordered_map<std::string, int32_t> gene_variant_idx_in_G;
    // maps an index in the LD file to a variant id
    std::unordered_map<int32_t, std::string> idx_in_file_gene_variant;
    // row and column names of G
    std::vector<std::string> gene_variant_ids;
  };

  struct gene_buffer_ld_mat_t {
    Eigen::MatrixXf* G_C;
    Eigen::VectorXf variances;
    // maps buffer variant ids to the index of the stored LD matrix
    std::unordered_map<std::string, int32_t> buffer_variant_idx_in_file;
    // maps buffer variant ids to column ids of G_C
    std::unordered_map<std::string, int32_t> buffer_variant_idx_in_G_C;
    // maps an index in the LD file to a variant it
    std::unordered_map<int32_t, std::string> idx_in_file_buffer_variant;
    // maps an index from the buffer ld file to a column of G_C
    std::unordered_map<int32_t, int32_t> buffer_ld_idx_to_idx_in_G_C;
    // column names
    std::vector<std::string> buffer_variant_ids;
    int32_t nbuffer_blocks;
  };

  std::unordered_map<std::string, gene_idx_t> genes;

  vector<string> last_buffer_variants_loaded;
  Eigen::MatrixXf last_buffer_mat_loaded;

  void load_idx();

  void load_gene_ld_mat(gene_ld_mat_t& gene_ld,
                        const std::vector<string>& gene_variants,
                        const std::string& gene_name,
                        bool load_sparse = false);

  void load_gene_buffer_ld_mat(gene_buffer_ld_mat_t& gene_buffer_ld,
                               const gene_ld_mat_t& gene_ld,
                               const std::vector<std::string>& buffer_variants_to_load,
                               const std::string& gene_name);

  void load_buffer_ld_mat(Eigen::MatrixXf& C,
                          const gene_buffer_ld_mat_t& gene_buffer_ld,
                          const std::vector<std::string>& buffer_variants);

  void load_conditional_ld_mats_helper(Eigen::SparseMatrix<double>& G_sp,
                                       Eigen::MatrixXf& G,
                                       Eigen::MatrixXf& C,
                                       Eigen::MatrixXf& G_C,
                                       const std::vector<std::string>& gene_variants,
                                       const std::vector<std::string>& conditional_variants,
                                       const std::string& gene_name,
                                       bool load_sparse);
};

#endif