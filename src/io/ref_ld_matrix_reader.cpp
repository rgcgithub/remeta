#include "ref_ld_matrix_reader.hpp"

#include <chrono>
#include <memory>
#include <sstream>
#include <string>
using namespace std;

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
using Eigen::MatrixXf;
using Eigen::SparseMatrix;
typedef Eigen::Triplet<double> Triplet;

#include "ref_ld_matrix_writer.hpp"
#include "../logging.hpp"
#include "../util.hpp"

RefLDMatrixReader::RefLDMatrixReader(const string& file_prefix)
 : file_prefix(file_prefix)
 , is_closed(true) {
  this->open(file_prefix);
}

RefLDMatrixReader::RefLDMatrixReader()
 : file_prefix("")
 , is_closed(true) {}

void RefLDMatrixReader::open(const string& file_prefix) {
  if (!this->is_closed) {
    this->close();
  }
  this->gene_ld_file.open(file_prefix + ".remeta.gene.ld");
  if (util::file_exists(file_prefix + ".remeta.buffer.ld")) {
    this->has_buffer_ld_file = true;
    this->buffer_ld_file.open(file_prefix + ".remeta.buffer.ld");
    string version;
    try {
      version = this->buffer_ld_file.read_bytes(13);
      this->float_size = this->buffer_ld_file.read<int32_t>();
    } catch (const runtime_error& e) {
      this->float_size = 4;
    }
  }
  this->index_file.open(file_prefix + ".remeta.ld.idx.gz");
  this->load_idx();

  if (sizeof(float) != 4) {
    throw runtime_error("bad float size: sizeof(float) != 4");
  }
  this->is_closed = false;
}

void RefLDMatrixReader::close() {
  this->gene_ld_file.close();
  this->buffer_ld_file.close();
  this->index_file.close();
  this->is_closed = true;
}

bool RefLDMatrixReader::closed() { 
  return this->is_closed; 
}

void RefLDMatrixReader::load_gene_ld_mat(gene_ld_mat_t& gene_ld,
                                         const vector<string>& gene_variants,
                                         const std::string& gene_name,
                                         bool load_sparse) {
  if (this->genes.count(gene_name) == 0) {
    throw runtime_error(
      "gene " + gene_name + " is not found in index " 
      + this->file_prefix + ".remeta.ld.idx.gz"
    );
  }
  gene_ld.gene_variant_ids = gene_variants;

  vector<string> variant_ids;
  this->load_gene_variant_ids(variant_ids, gene_name);
  for (size_t i = 0; i < variant_ids.size(); ++i) {
    gene_ld.gene_variant_idx_in_file[variant_ids[i]] = static_cast<int32_t>(i);
    gene_ld.idx_in_file_gene_variant[static_cast<int32_t>(i)] = variant_ids[i];
  }

  int nmissing = 0;
  for (size_t i = 0; i < gene_variants.size(); ++i) {
    if (gene_ld.gene_variant_idx_in_file.count(gene_variants[i]) == 0) {
      ++nmissing;
    }
    if (gene_ld.gene_variant_idx_in_G.count(gene_variants[i]) > 0) {
      throw runtime_error("duplicate variant in load_gene_ld_mat: this is a bug");
    }

    gene_ld.gene_variant_idx_in_G[gene_variants[i]] = 
      static_cast<int32_t>(gene_ld.gene_variant_idx_in_G.size());
  }
  // if (nmissing > 0.5*gene_variants.size()) {
  //   log_warning("LD matrix file has fewer than 50% of requested variants for " + gene_name);
  // }

  if (!load_sparse) {
    gene_ld.G->resize(
      gene_ld.gene_variant_idx_in_G.size(),
      gene_ld.gene_variant_idx_in_G.size()
    );
    gene_ld.G->setZero();
  }

  this->gene_ld_file.seek_addr(genes[gene_name].gene_ld_index);

  int32_t nvariants = this->gene_ld_file.read<int32_t>();
  int32_t nentries = this->gene_ld_file.read<int32_t>();
  // sparsity threshold
  this->gene_ld_file.read<float>();

  vector<Triplet> triplets;
  vector<float> var(gene_ld.gene_variant_idx_in_G.size(), 0.0f);
  for (int32_t i = 0; i < nvariants; ++i) {
    string variant = gene_ld.idx_in_file_gene_variant[i];
    if (gene_ld.gene_variant_idx_in_G.count(variant) > 0) {
      int32_t k = gene_ld.gene_variant_idx_in_G[variant];
      if (load_sparse) {
        var[k] = this->gene_ld_file.read<float>();
        triplets.push_back(
          Triplet(k, k, static_cast<double>(var[k]))
        );
      } else {
        (*gene_ld.G)(k, k) = this->gene_ld_file.read<float>();
      }
    } else {
      this->gene_ld_file.read<float>();
    }
  }

  int i1;
  int i2;
  string vid1;
  string vid2;
  float corr;
  float s;
  double cov;
  for (int32_t i = 0; i < nentries; ++i) {
    vid1 = gene_ld.idx_in_file_gene_variant[this->gene_ld_file.read<int32_t>()];
    vid2 = gene_ld.idx_in_file_gene_variant[this->gene_ld_file.read<int32_t>()];
    corr = this->gene_ld_file.read<float>();
    if (gene_ld.gene_variant_idx_in_G.count(vid1) > 0 
        && gene_ld.gene_variant_idx_in_G.count(vid2) > 0) {
      i1 = gene_ld.gene_variant_idx_in_G[vid1];
      i2 = gene_ld.gene_variant_idx_in_G[vid2];
      if (load_sparse) {
        cov = static_cast<double>(sqrt(var[i1])*sqrt(var[i2])*corr);
        triplets.push_back(Triplet(i1, i2, cov));
        triplets.push_back(Triplet(i2, i1, cov));
      } else {
        s = sqrt((*gene_ld.G)(i1, i1))*sqrt((*gene_ld.G)(i2, i2));
        (*gene_ld.G)(i1, i2) = corr*s;
        (*gene_ld.G)(i2, i1) = corr*s;
      }
    }
  }

  if (load_sparse) {
    gene_ld.G_sp->resize(gene_ld.gene_variant_idx_in_G.size(), gene_ld.gene_variant_idx_in_G.size());
    gene_ld.G_sp->setFromTriplets(triplets.begin(), triplets.end());
  }
}

void RefLDMatrixReader::load_gene_ld_mat(Eigen::MatrixXf& G,
                                         const std::vector<string>& gene_variants,
                                         const std::string& gene_name) {
  std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

  gene_ld_mat_t gene_ld;
  gene_ld.G = &G;
  gene_ld.G_sp = nullptr;
  this->load_gene_ld_mat(gene_ld,
                         gene_variants,
                         gene_name,
                         false);

  std::chrono::duration<double> d = std::chrono::steady_clock::now() - start;
  log_debug("load_gene_ld_mat: " + to_string(d.count()));
}

void RefLDMatrixReader::load_gene_ld_mat_sp_double(Eigen::SparseMatrix<double>& G_sp,
                                                   const vector<string>& gene_variants,
                                                   const string& gene_name) {
  std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

  gene_ld_mat_t gene_ld;
  gene_ld.G = nullptr;
  gene_ld.G_sp = &G_sp;
  this->load_gene_ld_mat(gene_ld,
                         gene_variants,
                         gene_name,
                         true);

  std::chrono::duration<double> d = std::chrono::steady_clock::now() - start;
  log_debug("load_gene_ld_mat_sp_double: " + to_string(d.count()));
}

void RefLDMatrixReader::load_conditional_ld_mats(MatrixXf& G,
                                                 MatrixXf& C,
                                                 MatrixXf& G_C,
                                                 const vector<string>& gene_variants,
                                                 const vector<string>& conditional_variants,
                                                 const string& gene_name) {
  SparseMatrix<double> G_sp;
  this->load_conditional_ld_mats_helper(G_sp, G, C, G_C, gene_variants, conditional_variants, gene_name, false);
}


void RefLDMatrixReader::load_conditional_ld_mats_sp_double(SparseMatrix<double>& G_sp,
                                                           MatrixXf& C,
                                                           MatrixXf& G_C,
                                                           const vector<string>& gene_variants,
                                                           const vector<string>& conditional_variants,
                                                           const string& gene_name) {
  MatrixXf G;
  this->load_conditional_ld_mats_helper(G_sp, G, C, G_C, gene_variants, conditional_variants, gene_name, true);
}


void RefLDMatrixReader::load_conditional_ld_mats_helper(SparseMatrix<double>& G_sp,
                                                        MatrixXf& G,
                                                        MatrixXf& C,
                                                        MatrixXf& G_C,
                                                        const vector<string>& gene_variants,
                                                        const vector<string>& conditional_variants,
                                                        const string& gene_name,
                                                        bool load_sparse) {
  std::chrono::time_point<std::chrono::steady_clock> start;
  std::chrono::duration<double> d;

  // first we need to determine which conditional variants fall in a gene region or a buffer region
  vector<string> gv;
  this->load_gene_variant_ids(gv, gene_name);
  unordered_set<string> all_gene_variants;
  for ( const string& v : gv ) {
    all_gene_variants.insert(v);
  }
  unordered_set<string> all_buffer_variants;
  this->load_buffer_variant_ids(gv, gene_name);
  for ( const string& v : gv) {
    all_buffer_variants.insert(v);
  }

  vector<string> gene_variants_to_load;
  vector<string> buffer_variants_to_load;
  unordered_set<string> gene_variants_seen;
  for ( const string& v : gene_variants ) {
    gene_variants_to_load.push_back(v);
    gene_variants_seen.insert(v);
  }
  for ( const string& v : conditional_variants ) {
    if (all_gene_variants.count(v) > 0 && gene_variants_seen.count(v) == 0) {
      gene_variants_to_load.push_back(v);
    } else if (all_buffer_variants.count(v) > 0) {
      buffer_variants_to_load.push_back(v);
    }
  }
  sort(buffer_variants_to_load.begin(), buffer_variants_to_load.end());

  // then we need to load intermediate LD matrices G1, C1, G_C1 with
  // gene variances in G1, and buffer variants in C1.
  gene_ld_mat_t gene_ld;
  MatrixXf G1;
  SparseMatrix<double> G1_sp;
  gene_ld.G = &G1;
  gene_ld.G_sp = &G1_sp;

  start = std::chrono::steady_clock::now();
  this->load_gene_ld_mat(gene_ld,
                         gene_variants_to_load,
                         gene_name,
                         load_sparse);
  d = std::chrono::steady_clock::now() - start;
  log_debug("load_gene_ld_mat: " + to_string(d.count()));

  gene_buffer_ld_mat_t gene_buffer_ld;
  MatrixXf G_C1;
  gene_buffer_ld.G_C = &G_C1;

  start = std::chrono::steady_clock::now();
  this->load_gene_buffer_ld_mat(gene_buffer_ld,
                                gene_ld,
                                buffer_variants_to_load,
                                gene_name);
  d = std::chrono::steady_clock::now() - start;
  log_debug("load_gene_buffer_ld_mat: " + to_string(d.count()));

  this->buffer_ld_file.seek_addr(this->genes[gene_name].buffer_ld_index);

  start = std::chrono::steady_clock::now();
  MatrixXf C1;
  this->load_buffer_ld_mat(C1, gene_buffer_ld, buffer_variants_to_load);
  d = std::chrono::steady_clock::now() - start;
  log_debug("load_buffer_ld_mat: " + to_string(d.count()));

  // finally we need to move gene variants in the set of conditional variants
  // from G1 to G, and remove them from C and G_C.
  vector<Triplet> triplets;
  if (!load_sparse) {
    G = G1(
      Eigen::seqN(0, gene_variants.size()),
      Eigen::seqN(0, gene_variants.size())
    );
  } else {
    G_sp.resize(gene_variants.size(), gene_variants.size());
    vector<Triplet> triplets;
    for (int k = 0; k < gene_ld.G_sp->outerSize(); ++k) {
      for (SparseMatrix<double>::InnerIterator it(*(gene_ld.G_sp), k); it; ++it) {
        if (it.row() < (int)gene_variants.size() && it.col() < (int)gene_variants.size()) {
          triplets.push_back(Triplet(it.row(), it.col(), it.value()));
        }
      }
    }
    G_sp.setFromTriplets(triplets.begin(), triplets.end());
  }
  G_C.resize(gene_variants.size(), conditional_variants.size());
  C.resize(conditional_variants.size(), conditional_variants.size());
  G_C.setZero();
  C.setZero();

  string v1;
  string v2;
  for (size_t j = 0; j < conditional_variants.size(); ++j) {
    v2 = conditional_variants[j];

    // construct the covariance matrix of the conditional variants
    for (size_t i = 0; i <= j; ++i) {
      v1 = conditional_variants[i];

      if (i == j && gene_buffer_ld.buffer_variant_idx_in_G_C.count(v1) > 0) {
        C(i, i) = gene_buffer_ld.variances(gene_buffer_ld.buffer_variant_idx_in_G_C.at(v1));
      } else if (gene_ld.gene_variant_idx_in_G.count(v1) > 0 
                 && gene_ld.gene_variant_idx_in_G.count(v2) > 0) {
        if (load_sparse) {
          C(i, j) = (*gene_ld.G_sp).coeff(
            gene_ld.gene_variant_idx_in_G[v1],
            gene_ld.gene_variant_idx_in_G[v2]
          );
        } else {
          C(i, j) = (*gene_ld.G)(gene_ld.gene_variant_idx_in_G[v1], gene_ld.gene_variant_idx_in_G[v2]);
        }
        C(j, i) = C(i, j);
      } else if (gene_buffer_ld.buffer_variant_idx_in_G_C.count(v1) > 0
                 && gene_buffer_ld.buffer_variant_idx_in_G_C.count(v2) > 0) {
        C(i, j) = C1(
          gene_buffer_ld.buffer_variant_idx_in_G_C.at(v1),
          gene_buffer_ld.buffer_variant_idx_in_G_C.at(v2)
        );
        C(j, i) = C (i, j);
      } else if (gene_buffer_ld.buffer_variant_idx_in_G_C.count(v1) > 0
                 && gene_ld.gene_variant_idx_in_G.count(v2) > 0) {
        C(i, j) = (*gene_buffer_ld.G_C)(
          gene_ld.gene_variant_idx_in_G.at(v2),
          gene_buffer_ld.buffer_variant_idx_in_G_C.at(v1)
        );
        C(j, i) = C(i, j);
      } else if (gene_ld.gene_variant_idx_in_G.count(v1) > 0
                 && gene_buffer_ld.buffer_variant_idx_in_G_C.count(v2) > 0) {
        C(i, j) = (*gene_buffer_ld.G_C)(
          gene_ld.gene_variant_idx_in_G.at(v1),
          gene_buffer_ld.buffer_variant_idx_in_G_C.at(v2)
        );
        C(j, i) = C(i, j);
      } else {
        C(i, j) = 0;
        C(j, i) = 0;
        //throw runtime_error("failed to find covariance for "  + v1 + " " + v2 + ": this is a bug");
      }
    }

    // construct the cross covariance matrix of gene-conditional variants
    for (size_t i = 0; i < gene_variants.size(); ++i) {
      v1 = gene_variants[i];
      if (gene_ld.gene_variant_idx_in_G.count(v2)) {
        if (load_sparse) {
          G_C(i, j) = (*gene_ld.G_sp).coeff(
            gene_ld.gene_variant_idx_in_G[v1],
            gene_ld.gene_variant_idx_in_G[v2]
          );
        } else {
          G_C(i, j) = (*gene_ld.G)(
            gene_ld.gene_variant_idx_in_G[v1],
            gene_ld.gene_variant_idx_in_G[v2]
          );
        }

      } else if (gene_buffer_ld.buffer_variant_idx_in_G_C.count(v2)) {
        G_C(i, j) = (*gene_buffer_ld.G_C)(
          gene_ld.gene_variant_idx_in_G[v1],
          gene_buffer_ld.buffer_variant_idx_in_G_C[v2]
        );
      }
    }
  }
}

void RefLDMatrixReader::load_gene_buffer_ld_mat(gene_buffer_ld_mat_t& gene_buffer_ld,
                                                const gene_ld_mat_t& gene_ld,
                                                const std::vector<std::string>& buffer_variants_to_load,
                                                const std::string& gene_name) {
  vector<string> variant_ids;
  this->load_buffer_variant_ids(variant_ids, gene_name);
  for (size_t i = 0; i < variant_ids.size(); ++i) {
    gene_buffer_ld.buffer_variant_idx_in_file[variant_ids[i]] = (int32_t)i;
    gene_buffer_ld.idx_in_file_buffer_variant[(int32_t)i] = variant_ids[i];
  }
  gene_buffer_ld.buffer_variant_ids = buffer_variants_to_load;

  for (size_t i = 0; i < buffer_variants_to_load.size(); ++i) {
    if (gene_buffer_ld.buffer_variant_idx_in_file.count(buffer_variants_to_load[i]) == 0) {
      log_error(
        "variant " + buffer_variants_to_load[i] + " is missing from gene " 
        + gene_name + " in buffer LD file: this is a bug",
        1
      );
    }
    if (gene_ld.gene_variant_idx_in_G.count(buffer_variants_to_load[i]) > 0) {
      throw runtime_error("duplicate variant in load_gene_buffer_ld_mat: this is a bug");
    }

    gene_buffer_ld.buffer_variant_idx_in_G_C[buffer_variants_to_load[i]] = 
      (int32_t)gene_buffer_ld.buffer_variant_idx_in_G_C.size();
  }
  gene_buffer_ld.G_C->resize(
    gene_ld.gene_variant_ids.size(),
    gene_buffer_ld.buffer_variant_ids.size()
  );
  gene_buffer_ld.G_C->setZero();

  // ngene_variants
  this->gene_ld_file.read<int32_t>();
  int32_t nbuffer_variants = this->gene_ld_file.read<int32_t>();
  int32_t nentries = this->gene_ld_file.read<int32_t>();
  // sparsity threshold
  this->gene_ld_file.read<float>();

  gene_buffer_ld.variances.resize(buffer_variants_to_load.size());
  for (int32_t i = 0; i < nbuffer_variants; ++i) {
    string variant = gene_buffer_ld.idx_in_file_buffer_variant[i];
    if (gene_buffer_ld.buffer_variant_idx_in_G_C.count(variant) > 0) {
      gene_buffer_ld.variances(gene_buffer_ld.buffer_variant_idx_in_G_C[variant])
        = this->gene_ld_file.read<float>();
    } else {
      this->gene_ld_file.read<float>();
    }
  }

  int32_t w;
  for (int32_t i = 0; i < nbuffer_variants; ++i) {
    w = this->gene_ld_file.read<int32_t>();
    if (gene_buffer_ld.buffer_variant_idx_in_G_C.count(variant_ids[i]) > 0) {
      gene_buffer_ld.buffer_ld_idx_to_idx_in_G_C[w] =
        gene_buffer_ld.buffer_variant_idx_in_G_C[variant_ids[i]];
    }
  }
  gene_buffer_ld.nbuffer_blocks = this->gene_ld_file.read<int32_t>();

  int i1;
  int i2;
  string vid1;
  string vid2;
  float corr;;
  for (int32_t i = 0; i < nentries; ++i) {
    i1 = this->gene_ld_file.read<int32_t>();
    i2 = this->gene_ld_file.read<int32_t>();
    vid1 = gene_ld.idx_in_file_gene_variant.at(i1);
    vid2 = gene_buffer_ld.idx_in_file_buffer_variant[i2];
    if (this->float_size == 1) {
      corr = this->gene_ld_file.read<int8_t>()*pow_2_m7;
      corr = corr == 0 ? 1 : corr;
    } else if (this->float_size == 2) {
      corr = this->gene_ld_file.read<int16_t>()*pow_2_m15;
      corr = corr == 0 ? 1 : corr;
    } else if (this->float_size == 4) {
      corr = this->gene_ld_file.read<float>();
    } else {
      throw runtime_error("bad float size");
    }

    if (gene_ld.gene_variant_idx_in_G.count(vid1) > 0 
        && gene_buffer_ld.buffer_variant_idx_in_G_C.count(vid2) > 0) {
      i1 = gene_ld.gene_variant_idx_in_G.at(vid1);
      i2 = gene_buffer_ld.buffer_variant_idx_in_G_C[vid2];
      if (gene_ld.G->size() > 0) {
        (*gene_buffer_ld.G_C)(i1, i2) = corr*sqrt((*gene_ld.G)(i1, i1))
                                          *sqrt(gene_buffer_ld.variances(i2));
      } else {
        (*gene_buffer_ld.G_C)(i1, i2) = corr*sqrt((*gene_ld.G_sp).coeff(i1, i1))
                                          *sqrt(gene_buffer_ld.variances(i2));
      }
    }
  }
}

void RefLDMatrixReader::load_buffer_ld_mat(MatrixXf& C,
                                           const gene_buffer_ld_mat_t& gene_buffer_ld,
                                           const std::vector<std::string>& buffer_variants) {
  C.resize(buffer_variants.size(), buffer_variants.size());
  C.setZero();

  if (buffer_variants.size() == 0) {
    return;
  }

  if (buffer_variants.size() == this->last_buffer_variants_loaded.size()) {
    bool is_same_ld_mat = true;
    for (size_t i = 0; i < buffer_variants.size(); ++i) {
      is_same_ld_mat = is_same_ld_mat
        && buffer_variants[i] == this->last_buffer_variants_loaded[i];
    }

    if (is_same_ld_mat) {
      C = this->last_buffer_mat_loaded;
      return;
    }
  }

  int32_t blocks_seen = 0;
  int32_t i;
  int32_t j;
  int32_t i1;
  int32_t j1;
  int32_t n_entries;
  float corr;
  while (blocks_seen <= gene_buffer_ld.nbuffer_blocks && !this->buffer_ld_file.eof()) {
    i = this->buffer_ld_file.read<int32_t>();
    if (i == BLOCK_START_MARKER) {
      ++blocks_seen;
    }

    while (!this->buffer_ld_file.eof()) {
      i = this->buffer_ld_file.read<int32_t>();
      if (i == BLOCK_END_MARKER) {
        break;
      } else if (i != ROW_START_MARKER) {
        throw runtime_error("corrupted buffer LD file");
      }
      n_entries = this->buffer_ld_file.read<int32_t>();
      i = this->buffer_ld_file.read<int32_t>();

      if (gene_buffer_ld.buffer_ld_idx_to_idx_in_G_C.count(i) == 0) {
        for (int32_t l = 0; l < n_entries; ++l) {
          j = this->buffer_ld_file.read<int32_t>();
          if (this->float_size == 1) {
            corr = this->buffer_ld_file.read<int8_t>()*pow_2_m7;
            corr = corr == 0 ? 1 : corr;
          } else if (this->float_size == 2) {
            corr = this->buffer_ld_file.read<int16_t>()*pow_2_m15;
            corr = corr == 0 ? 1 : corr;
          } else if (this->float_size == 4) {
            corr = this->buffer_ld_file.read<float>();
          } else {
            throw runtime_error("bad float size");
          }
        }
      } else {
        i1 = gene_buffer_ld.buffer_ld_idx_to_idx_in_G_C.at(i);
        for (int32_t l = 0; l < n_entries; ++l) {
          j = this->buffer_ld_file.read<int32_t>();
          if (this->float_size == 1) {
            corr = this->buffer_ld_file.read<int8_t>()*pow_2_m7;
            corr = corr == 0 ? 1 : corr;
          } else if (this->float_size == 2) {
            corr = this->buffer_ld_file.read<int16_t>()*pow_2_m15;
            corr = corr == 0 ? 1 : corr;
          } else if (this->float_size == 4) {
            corr = this->buffer_ld_file.read<float>();
          } else {
            throw runtime_error("bad float size");
          }

          if (gene_buffer_ld.buffer_ld_idx_to_idx_in_G_C.count(j) > 0) {
            j1 = gene_buffer_ld.buffer_ld_idx_to_idx_in_G_C.at(j);
            C(i1, j1) = corr
                        * sqrt(gene_buffer_ld.variances(i1))
                        * sqrt(gene_buffer_ld.variances(j1));
            C(j1, i1) = C(i1, j1);
          }
        }
      }
    }
  }
  this->last_buffer_mat_loaded = C;
  this->last_buffer_variants_loaded = buffer_variants;
}

void RefLDMatrixReader::load_gene_variant_ids(vector<string>& variant_ids,
                                              const string& gene_name) {
  this->index_file.seek_addr(this->genes[gene_name].index_index);
  string line = this->index_file.readline();
  vector<string> cols = util::str_split(line, "\t");
  variant_ids = util::str_split(cols[3], ",");
}

void RefLDMatrixReader::load_buffer_variant_ids(vector<string>& variant_ids,
                                                const string& gene_name) {
  this->index_file.seek_addr(this->genes[gene_name].index_index);
  string line = this->index_file.readline();
  vector<string> cols = util::str_split(line, "\t");
  variant_ids = util::str_split(cols[4], ",");
}

void RefLDMatrixReader::load_idx() {
  string line;
  string gene;
  stringstream ss;
  int64_t ld_addr;
  int64_t buffer_addr;
  int64_t idx_addr = this->index_file.tell();
  while (!this->index_file.eof()) {
    line = this->index_file.readline();
    ss = stringstream(line);
    ss >> gene;
    ss >> ld_addr;
    ss >> buffer_addr;
    genes[gene] = gene_idx_t {
      ld_addr,
      buffer_addr,
      idx_addr
    };
    idx_addr = this->index_file.tell();
  }
}

bool RefLDMatrixReader::contains_gene(const string& gene_name) {
  return this->genes.count(gene_name) > 0;
}