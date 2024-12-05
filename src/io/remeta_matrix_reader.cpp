#include "remeta_matrix_reader.hpp"

#include <string>
using namespace std;

#include <Eigen/Dense>
using Eigen::MatrixXf;

#include "../util.hpp"
#include "../logging.hpp"

RemetaMatrixReader::RemetaMatrixReader()
 : reader_opened(NONE)
 , is_closed(true) {
}

void RemetaMatrixReader::open(const string& prefix) {
  this->close();
  if (util::file_exists(prefix + ".rg.ld")
      && util::file_exists(prefix + ".rg.ld.idx.gz")) {
    this->regenie_reader.open(prefix + ".rg.ld", prefix + ".rg.ld.idx.gz");
    this->reader_opened = REGENIE_READER;
    log_info("using LD matrices: " + prefix + ".{rg.ld, rg.ld.idx.gz}");
  } else if (util::file_exists(prefix + ".metamat")
    && util::file_exists(prefix + ".metamat.idx.gz")) {
    log_info("using LD matrices: " + prefix + ".{metamat, metamat.idx.gz}");
    this->regenie_reader.open(prefix + ".metamat", prefix + ".metamat.idx.gz");
    this->reader_opened = REGENIE_READER;
  } else if (util::file_exists(prefix + ".remeta.gene.ld")
    && util::file_exists(prefix + ".remeta.buffer.ld")
    && util::file_exists(prefix + ".remeta.ld.idx.gz")) {
    ref_reader.open(prefix); 
    this->reader_opened = REF_READER; 
    log_info("using LD matrices: " + prefix + ".remeta.{gene.ld, buffer.ld, ld.idx.gz}");
  } else {
    throw runtime_error("could not find LD matrices with prefix " + prefix);
  }
  this->is_closed = false;
}

void RemetaMatrixReader::close() {
  if (!this->regenie_reader.closed()) {
    this->regenie_reader.close();
  }
  if (!this->ref_reader.closed()) {
    this->ref_reader.close();
  }
  this->reader_opened = NONE;
}

bool RemetaMatrixReader::closed() {
  return this->is_closed;
}

bool RemetaMatrixReader::has_conditional_variants() {
  return this->reader_opened == REF_READER;
}

void RemetaMatrixReader::load_gene_ld_mat(MatrixXf& G,
                                          const vector<string>& gene_variants,
                                          const string& gene_name) {
  this->check_closed();
  if (this->reader_opened == REF_READER) {
    this->ref_reader.load_gene_ld_mat(G, gene_variants, gene_name);
  } else {
    this->regenie_reader.load_gene_ld_mat(G, gene_variants, gene_name);
  }
}

void RemetaMatrixReader::load_conditional_ld_mats(Eigen::MatrixXf& G,
                                                  Eigen::MatrixXf& C,
                                                  Eigen::MatrixXf& G_C,
                                                  const std::vector<std::string>& gene_variants,
                                                  const std::vector<std::string>& conditional_variants,
                                                  const std::string& gene_name) {
  this->check_closed();
  if (this->reader_opened == REF_READER) {
    this->ref_reader.load_conditional_ld_mats(G, C, G_C, gene_variants, conditional_variants, gene_name);
  } else {
    this->regenie_reader.load_gene_ld_mat(G, gene_variants, gene_name);
    C.resize(0, 0);
    G_C.resize(0, 0);
  }
}

bool RemetaMatrixReader::contains_gene(const string& gene_name) {
  this->check_closed();
  if (this->reader_opened == REF_READER) {
    return this->ref_reader.contains_gene(gene_name);
  } else {
    return this->regenie_reader.contains_gene(gene_name);
  }
}

void RemetaMatrixReader::check_closed() {
  if (this->closed()) {
    throw runtime_error("RemetaMatrixReader attempting to operate on closed file");
  }
}