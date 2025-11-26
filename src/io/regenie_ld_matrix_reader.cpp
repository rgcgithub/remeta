#include "regenie_ld_matrix_reader.hpp"

#include <algorithm>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
using namespace std;

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::MatrixXf;
#include <Eigen/Sparse>
using Eigen::SparseMatrix;
typedef Eigen::Triplet<double> Triplet;

#include "bgz_reader.hpp"

RegenieLDMatrixReader::RegenieLDMatrixReader()
  : is_closed(true)
  , mat_bgzf(nullptr)
  , idx_bgzf(nullptr)
  {}

RegenieLDMatrixReader::RegenieLDMatrixReader(string mat_file, string idx_file)
 : is_closed(false)
 , mat_file(mat_file)
 , idx_file(idx_file) {
  this->mat_bgzf = bgzf_open(mat_file.c_str(), "r");
  this->idx_bgzf = bgzf_open(idx_file.c_str(), "r");
  this->buffer = KS_INITIALIZE;

  if (this->mat_bgzf == NULL) {
    throw runtime_error("could not load: " + mat_file);
  } else if (this->idx_bgzf == NULL) {
    throw runtime_error("could not load: " + idx_file);
  }

  if (sizeof(float) != 4) {
    throw runtime_error("bad float size: sizeof(float) != 4");
  }

  this->load_idx();
  this->read_sample_size();
}

RegenieLDMatrixReader::~RegenieLDMatrixReader() {
  if (!this->closed()) {
    bgzf_close(this->idx_bgzf);
    bgzf_close(this->mat_bgzf);
    ks_free(&this->buffer);
  }
}

RegenieLDMatrixReader::RegenieLDMatrixReader(const RegenieLDMatrixReader& other) {
  this->is_closed = other.is_closed;
  this->mat_file = other.mat_file;
  this->idx_file = other.idx_file;
  this->sample_size = other.sample_size;
  this->buffer = KS_INITIALIZE;

  if (!this->is_closed) {
    this->mat_bgzf = bgzf_open(mat_file.c_str(), "r");
    this->idx_bgzf = bgzf_open(idx_file.c_str(), "r"); 
    if (this->mat_bgzf == NULL) {
      throw runtime_error("could not load: " + mat_file);
    } else if (this->idx_bgzf == NULL) {
      throw runtime_error("could not load: " + idx_file);
    }
    this->load_idx();
  }
}

RegenieLDMatrixReader& RegenieLDMatrixReader::operator=(RegenieLDMatrixReader other) {
  swap(*this, other);
  return *this;
}

void swap(RegenieLDMatrixReader& first, RegenieLDMatrixReader& second) {
  std::swap(first.is_closed, second.is_closed);
  std::swap(first.mat_file, second.mat_file);
  std::swap(first.idx_file, second.idx_file);
  std::swap(first.sample_size, second.sample_size);
  std::swap(first.mat_bgzf, second.mat_bgzf);
  std::swap(first.idx_bgzf, second.idx_bgzf);
  std::swap(first.buffer, second.buffer);
  std::swap(first.idx_addr, second.idx_addr);
  std::swap(first.mat_addr, second.mat_addr);
}

void RegenieLDMatrixReader::open(string mat_file, string idx_file) {
  RegenieLDMatrixReader reader(mat_file, idx_file);
  swap(*this, reader);
}

void RegenieLDMatrixReader::close() {
  RegenieLDMatrixReader reader;
  swap(*this, reader);
}

void RegenieLDMatrixReader::load_idx() {
  int64_t idx_addr;
  int64_t mat_addr;
  string gene;
  string line;
  stringstream ss(line);
  int ret;
  int peek = bgzf_peek(this->idx_bgzf);
  while (peek != -1) {
    if (peek == -2) {
      throw runtime_error("error loading " + this->idx_file);
    }

    idx_addr = bgzf_tell(this->idx_bgzf);
    ret = bgzf_getline(this->idx_bgzf, '\n', &this->buffer);
    if (ret < -1) {
      throw runtime_error("bgzf_getline error " + to_string(ret));
    }

    line = string(this->buffer.s);
    ss = stringstream(line);
    ss >> gene;
    ss >> mat_addr;

    this->idx_addr[gene] = idx_addr;
    this->mat_addr[gene] = mat_addr;
    peek = bgzf_peek(this->idx_bgzf);
  }
}

void RegenieLDMatrixReader::read_sample_size() {
  ssize_t bytes_read = bgzf_read(this->mat_bgzf,
                                 &(this->sample_size),
                                 sizeof(this->sample_size));
  if (bytes_read != sizeof(this->sample_size)) {
    throw runtime_error("could not read sample size from " + this->mat_file);
  }
}

void RegenieLDMatrixReader::load_ld_mat(MatrixXf& mat,
                                        vector<string>& variant_ids,
                                        const string& gene_name) {
  if (this->closed()) {
    throw runtime_error("LDMatrixReader: operating on a closed file");
  }

  variant_ids.clear();
  this->load_variant_ids(variant_ids, gene_name);

  if (this->mat_addr.count(gene_name) == 0) {
    throw runtime_error(gene_name + " is not contained in " + this->mat_file);
  }
  
  if (bgzf_seek(this->mat_bgzf, this->mat_addr[gene_name], SEEK_SET) == -1) {
    throw runtime_error("error calling seek in " + this->mat_file);
  }

  char type;
  if (bgzf_read(this->mat_bgzf, &type, 1) != 1) {
    throw runtime_error("error reading matrix type from from " + this->mat_file);
  }

  int32_t nrows;
  float sparsity_parameter;
  if (bgzf_read(this->mat_bgzf, &nrows, 4) != 4) {
    throw runtime_error("error reading nrows from " + this->mat_file);
  }
  if (bgzf_read(this->mat_bgzf, &sparsity_parameter, 4) != 4) {
    throw runtime_error("error reading sparsity_parameter from " + this->mat_file);
  }
  mat.resize(nrows, nrows);
  mat *= 0;

  if (type == 's') {
    load_ld_mat_sparse(mat);
  } else if (type == 'd') {
    load_ld_mat_dense(mat);
  } else {
    throw runtime_error("unrecognized matrix type " + to_string(type) + " in " + this->mat_file);
  }
}

void RegenieLDMatrixReader::load_gene_ld_mat(MatrixXf& G,
                                             const vector<string>& gene_variants,
                                             const string& gene_name) {
  MatrixXf full_ld_mat;
  vector<string> variant_ids;
  this->load_ld_mat(full_ld_mat, variant_ids, gene_name);
  G.resize(gene_variants.size(), gene_variants.size());
  G.setZero();

  unordered_map<string, int> variant_to_index_in_G;
  for (size_t i = 0; i < gene_variants.size(); ++i) {
    variant_to_index_in_G[gene_variants[i]] = i;
  }

  for (size_t i = 0; i < variant_ids.size(); ++i) {
    for (size_t j = 0; j < variant_ids.size(); ++j) {
      if (variant_to_index_in_G.count(variant_ids[i]) > 0
          && variant_to_index_in_G.count(variant_ids[j]) > 0) {
        G(
          variant_to_index_in_G.at(variant_ids[i]),
          variant_to_index_in_G.at(variant_ids[j])
        ) = full_ld_mat(i, j);  
      }
    }
  }
}

void RegenieLDMatrixReader::load_gene_ld_mat_sp_double(SparseMatrix<double>& G_sp,
                                                       const vector<string>& gene_variants,
                                                       const string& gene_name) {
  MatrixXf full_ld_mat;
  vector<string> variant_ids;
  this->load_ld_mat(full_ld_mat, variant_ids, gene_name);
  G_sp.resize(gene_variants.size(), gene_variants.size());

  unordered_map<string, int> variant_to_index_in_G;
  for (size_t i = 0; i < gene_variants.size(); ++i) {
    variant_to_index_in_G[gene_variants[i]] = i;
  }

  vector<Triplet> triplets;
  for (size_t i = 0; i < variant_ids.size(); ++i) {
    for (size_t j = 0; j < variant_ids.size(); ++j) {
      if (variant_to_index_in_G.count(variant_ids[i]) > 0
          && variant_to_index_in_G.count(variant_ids[j]) > 0) {
        triplets.push_back(
          Triplet(
            variant_to_index_in_G.at(variant_ids[i]),
            variant_to_index_in_G.at(variant_ids[j]),
            static_cast<double>(full_ld_mat(i, j))
          )
        );
      }
    }
  }
  G_sp.setFromTriplets(triplets.begin(), triplets.end());
}

void RegenieLDMatrixReader::load_ld_mat_dense(MatrixXf& mat) {
  float row[mat.rows()];
  for (int i = 0; i < mat.rows(); ++i) {
    if (bgzf_read(this->mat_bgzf, &row, sizeof(float)*(i+1)) != (ssize_t)sizeof(float)*(i+1)) {
      throw runtime_error("failed to read dense matrix");
    }
    for (int j = 0; j <= i; ++j) {
      mat(i, j) = row[j];
      mat(j, i) = row[j];
    }
  }
}

// helper type to load data
typedef struct {
  int32_t i;
  int32_t j;
  float data;
} spr_mat_entry_t;

void RegenieLDMatrixReader::load_ld_mat_sparse(MatrixXf& mat) {
  float variances[mat.rows()];
  if (bgzf_read(this->mat_bgzf, &variances, sizeof(float)*mat.rows()) != (ssize_t)sizeof(float)*mat.rows()) {
    throw runtime_error("failed to read sparse matrix variances");
  }

  spr_mat_entry_t entry;
  float element;
  ssize_t size = sizeof(spr_mat_entry_t);
  if (bgzf_read(this->mat_bgzf, &entry, size) != size) {
    throw runtime_error("error reading sparse matrix entry");
  }
  while (entry.i != -1 && entry.j != -1) {
    element = entry.data*sqrt(variances[entry.i])*sqrt(variances[entry.j]);
    mat(entry.i, entry.j) = element;
    mat(entry.j, entry.i) = element;
    if (bgzf_read(this->mat_bgzf, &entry, size) != size) {
      throw runtime_error("error reading sparse matrix entry");
    }
  }

  for (int i = 0; i < mat.rows(); ++i) {
    mat(i, i) = variances[i];
  }
}

void RegenieLDMatrixReader::load_variant_ids(vector<string>& variant_ids, const string& gene_name) {
  if (this->closed()) {
    throw runtime_error("LDMatrixReader: operating on a closed file");
  } 
  
  if (bgzf_seek(this->idx_bgzf, this->idx_addr[gene_name], SEEK_SET) == -1) {
    throw runtime_error("error calling seek in " + this->idx_file);
  }
  if (bgzf_getline(this->idx_bgzf, '\n', &(this->buffer)) < -1) {
    throw runtime_error("error reading " + this->idx_file);
  }

  string line(this->buffer.s);
  stringstream ss(line);
  string tmp;

  // read gene name
  ss >> tmp;
  if (tmp != gene_name) {
    throw runtime_error("bad seek in " + this->idx_file);
  }
  // read virtual pointer
  ss >> tmp;
  // read variant list
  ss >> tmp;

  ss = stringstream(tmp);
  while (getline(ss, tmp, ',')) {
    variant_ids.push_back(tmp);
  }
}

vector<string> RegenieLDMatrixReader::get_gene_names() {
  if (this->closed()) {
    throw runtime_error("LDMatrixReader: operating on a closed file");
  }

  BgzReader reader(this->idx_file);
  stringstream ss;
  string line;
  string gene;
  vector<string> gene_names;
  while (!reader.eof()) {
    line = reader.readline();
    ss = stringstream(line);
    ss >> gene;
    gene_names.push_back(gene);
  }
  return gene_names;
}

bool RegenieLDMatrixReader::contains_gene(const string& gene_name) {
  return mat_addr.count(gene_name) > 0;
}