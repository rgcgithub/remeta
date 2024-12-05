#include "run_compute_ref_ld.hpp"

#if defined(WITH_MKL)
#include <mkl.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <chrono>
#include <numeric>
#include <algorithm>

#include <Eigen/Dense>
using Eigen::MatrixXf;
using Eigen::VectorXf;
#include <Eigen/Sparse>
using Eigen::SparseMatrix;

#include "pgenlibr.h"

#include "htpv4_pos.hpp"
#include "io/block_pgen_reader.hpp"
#include "io/bgz_reader.hpp"
#include "io/ref_ld_block_mapper.hpp"
#include "io/ref_ld_matrix_writer.hpp"
#include "logging.hpp"
#include "parameter_checks.hpp"
namespace pc = parameter_checks;
#include "util.hpp"
#include "variant_filter.hpp"

struct ld_variant_t {
  variant_id id;
  HTPv4Pos pos;
  double variance;
  double mean;

  inline bool operator<(const ld_variant_t& other) {
    if (this->pos == other.pos) {
      return this->id < other.id;
    } else {
      return this->pos < other.pos;
    }
  }
};

void parse_gene_file_entry(string& gene,
                           string& chr,
                           int& start,
                           int& end,
                           const string& line,
                           const string& gene_file) {
  vector<string> row = util::str_split(line, "\t ");
  if (row.size() < 4) {
    throw runtime_error("gene_file " + gene_file + " has fewer than 4 columns");
  }
  gene = row[0];
  chr = row[1];
  start = stoi(row[2]);
  end = stoi(row[3]);
}

string format_region(string chrom, int start, int end) {
  return chrom + ":" + to_string(start) + "-" + to_string(end);
}

inline void extract_block_corr(vector<ref_ld_matrix_entry_t>& matrix_entries,
                               unordered_set<variant_id>& block1_variants_inserted,
                               unordered_set<variant_id>& block2_variants_inserted,
                               const MatrixXf& block1,
                               const MatrixXf& block2,
                               const vector<pvar_variant_t>& block1_variants,
                               const vector<pvar_variant_t>& block2_variants,
                               const int& nvar_read1,
                               const int& nvar_read2,
                               const unordered_map<variant_id, ld_variant_t>& block1_var_info,
                               const unordered_map<variant_id, ld_variant_t>& block2_var_info,
                               const float& r2_thr,
                               const bool& is_same_block) {
  if (nvar_read1 == 0 || nvar_read2 == 0) {
    return ;
  }

  variant_id id1;
  variant_id id2;
  ld_variant_t v1;
  ld_variant_t v2;
  float cov;
  float corr;

  vector<ld_variant_t> block1_ld_variants;
  vector<ld_variant_t> block2_ld_variants;
  vector<int> block1_keep_idx;
  vector<int> block2_keep_idx;
  for (int i = 0; i < nvar_read1; ++i) {
    id1 = block1_variants[i].id;
    if (block1_var_info.count(id1) > 0) {
      block1_ld_variants.push_back(block1_var_info.at(id1));
      block1_keep_idx.push_back(i);
    }
  }

  for (int j = 0; j < nvar_read2; ++j) {
    id2 = block2_variants[j].id;
    if (block2_var_info.count(id2) > 0) {
      block2_ld_variants.push_back(block2_var_info.at(id2));
      block2_keep_idx.push_back(j);
    }
  }
  MatrixXf GtG = block1(Eigen::all, block1_keep_idx).transpose() * block2(Eigen::all, block2_keep_idx) / block1.rows();

  for (size_t i = 0; i < block1_keep_idx.size(); ++i) {
    v1 = block1_ld_variants[i];
    id1 = block1_variants[block1_keep_idx[i]].id;
    if (v1.variance == 0) {
      continue;
    }

    size_t end;
    if (is_same_block) {
      end = i;
    } else {
      end = block2_keep_idx.size();
    }

    for (size_t j = 0; j < end; ++j) {
      v2 = block2_ld_variants[j];
      id2 = block2_variants[block2_keep_idx[j]].id;
      if (id1 == id2 || v2.variance == 0) {
        continue;
      }

      cov = (GtG(i,j) - v1.mean*v2.mean);
      corr = cov / sqrt(v1.variance*v2.variance);

      if (corr*corr >= r2_thr) {
        block1_variants_inserted.insert(id1);
        block2_variants_inserted.insert(id2);
        matrix_entries.push_back({id1, id2, corr});
      }
    }
  }
}

inline VectorXf compute_block_mean(const MatrixXf& Gblock) {
  VectorXf c = Eigen::VectorXf::Ones(Gblock.rows()) / Gblock.rows();
  VectorXf block_means = c.transpose() * Gblock;
  return block_means;
}

inline VectorXf compute_block_variance(const MatrixXf& Gblock, const VectorXf& block_means) {
  VectorXf block_variances = Eigen::VectorXf::Zero(Gblock.cols());
  for (ssize_t p = 0; p < Gblock.cols(); ++p) {
    block_variances(p) = Gblock.col(p).transpose() * Gblock.col(p);
    block_variances(p) /= Gblock.rows();
    block_variances(p) -= block_means(p)*block_means(p);
  }
  return block_variances;
}

inline VectorXf compute_block_mean(const MatrixXf& Gblock, const vector<int>& keep_idx) {
  VectorXf c = Eigen::VectorXf::Ones(Gblock.rows()) / Gblock.rows();
  VectorXf block_means = c.transpose() * Gblock;
  return block_means(keep_idx);
}

inline VectorXf compute_block_variance(const MatrixXf& Gblock, const VectorXf& block_means, const vector<int>& keep_idx) {
  VectorXf block_variances = Eigen::VectorXf::Zero(keep_idx.size());
  for (size_t i = 0; i < keep_idx.size(); ++i) {
    block_variances(i) = Gblock.col(keep_idx[i]).transpose() * Gblock.col(keep_idx[i]);
    block_variances(i) /= Gblock.rows();
    block_variances(i) -= block_means(i)*block_means(i);
  }
  return block_variances;
}

void run_compute_ref_ld(const string& target_pfile,
                        const string& gene_file,
                        const string& chrom,
                        const string& out_prefix,
                        const std::string& buffer_pfile,
                        const double& buffer_mb,
                        const double& buffer_cm,
                        const std::string& genetic_map_file,
                        const double& target_r2,
                        const double& buffer_r2,
                        const int& float_size,
                        const int& block_size,
                        const int& nthreads,
                        const string& target_extract,
                        const string& target_exclude,
                        const string& target_keep,
                        const string& target_remove,
                        const string& buffer_extract,
                        const string& buffer_exclude,
                        const bool& skip_buffer) {
  pc::check_file_exists(target_pfile + ".pgen");
  pc::check_file_exists(target_pfile + ".pvar");
  pc::check_file_exists(target_pfile + ".psam");
  pc::check_file_exists(gene_file);
  if (buffer_pfile != "") {
    pc::check_file_exists(buffer_pfile + ".pgen");
    pc::check_file_exists(buffer_pfile + ".pvar");
    pc::check_file_exists(buffer_pfile + ".psam");
  }
  pc::check_float_size(float_size, buffer_r2);
  if (target_extract != "") {
    pc::check_file_exists(target_extract);
  }
  if (target_exclude != "") {
    pc::check_file_exists(target_exclude);
  }
  if (target_keep != "") {
    pc::check_file_exists(target_keep);
  }
  if (target_remove != "") {
    pc::check_file_exists(target_remove);
  }
  if (buffer_extract != "") {
    pc::check_file_exists(buffer_extract);
  }
  if (buffer_exclude != "") {
    pc::check_file_exists(buffer_exclude);
  }

  pc::check_block_mapper_params(buffer_mb, buffer_cm, genetic_map_file);
  RefLDBlockMapper block_mapper;
  ld_block_map_t ld_block_map;
  if (genetic_map_file != "") {
    ld_block_map = block_mapper.map_ld_blocks_cm(gene_file,
                                                 genetic_map_file,
                                                 chrom,
                                                 buffer_cm);
  } else {
    ld_block_map = block_mapper.map_ld_blocks_mb(gene_file, chrom, buffer_mb);
  }
  unordered_map<gene_id, int64_t> gene_first_block_addr;
  unordered_set<gene_id> genes_skipped;

#if defined(WITH_MKL)
  mkl_set_num_threads(nthreads);
#endif

#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
  Eigen::setNbThreads(nthreads);
#endif

  log_info("loading target pgen/pvar/psam files: " + target_pfile);
  BlockPgenReader target_reader;
  target_reader.load(target_pfile, nthreads, 10000, target_keep, target_remove);
  int nsamples = target_reader.get_nkeep_samples();

  BlockPgenReader buffer_reader;
  if (buffer_pfile != "") {
    log_info("loading buffer pgen/pvar/psam files: " + buffer_pfile);
    buffer_reader.load(buffer_pfile, 
                       nthreads,
                       10000,
                       target_reader.get_keep_list(), 
                       target_reader.get_remove_list());
  } else {
    log_info("loading buffer pgen/pvar/psam files: " + target_pfile);
    buffer_reader.load(target_pfile,
                       nthreads,
                       10000,
                       target_reader.get_keep_list(),
                       target_reader.get_remove_list());
  }

  if (buffer_reader.get_nkeep_samples() != target_reader.get_nkeep_samples()) {
    log_error("buffer file is missing samples in the target file (kept " 
      + to_string(target_reader.get_nkeep_samples()) 
      + " in target but " + to_string(buffer_reader.get_nkeep_samples()) 
      + " in buffer)",
      1
    );
  }

  if (target_extract != "") {
    BgzReader extract_reader(target_extract);
    unordered_set<variant_id> extract_list;
    while (!extract_reader.eof()) {
      string line = extract_reader.readline();
      extract_list.insert(line);
    }
    target_reader.set_extract_list(extract_list);
  }
  if (target_exclude != "") {
    BgzReader exclude_reader(target_exclude);
    unordered_set<variant_id> exclude_list;
    while (!exclude_reader.eof()) {
      string line = exclude_reader.readline();
      exclude_list.insert(line);
    }
    target_reader.set_exclude_list(exclude_list);
  }
  if (buffer_extract != "") {
    BgzReader extract_reader(buffer_extract);
    unordered_set<variant_id> extract_list;
    while (!extract_reader.eof()) {
      string line = extract_reader.readline();
      extract_list.insert(line);
    }
    buffer_reader.set_extract_list(extract_list);
  }
  if (buffer_exclude != "") {
    BgzReader exclude_reader(buffer_exclude);
    unordered_set<variant_id> exclude_list;
    while (!exclude_reader.eof()) {
      string line = exclude_reader.readline();
      exclude_list.insert(line);
    }
    buffer_reader.set_exclude_list(exclude_list);
  }

  BgzReader genes(gene_file);
  RefLDMatrixWriter out(out_prefix, float_size);

  MatrixXf Gblock1(nsamples, block_size);
  MatrixXf Gblock2(nsamples, block_size);
  vector<pvar_variant_t> block1_variants;
  vector<pvar_variant_t> block2_variants;
  int nvar_read1;
  int nvar_read2;
  unordered_map<variant_id, ld_variant_t> all_buffer_variants;
  int ld_block_idx = 0;
  size_t gene_variants_saved = 0;
  size_t gene_records_written = 0;
  size_t gene_buffer_records_written = 0;
  size_t buffer_records_written = 0;

  std::chrono::time_point<std::chrono::steady_clock> start, gene_start, t;
  std::chrono::duration<double> d;
  start = std::chrono::steady_clock::now();

  while (!genes.eof()) {
    string line = genes.readline();

    string gene;
    string chr;
    int start_bp;
    int end_bp;
    try {
      parse_gene_file_entry(gene, chr, start_bp, end_bp, line, gene_file);
    } catch(std::exception& e) {
      log_error("failed to parse entry in " + gene_file 
                + " (line: " + line + ")", 1);
    }

    if (chr != chrom) {
      continue;
    }

    int buffer_start_bp = -1;
    int buffer_end_bp = -1;
    for ( const ld_block_t& ld_block : ld_block_map.gene_ld_blocks[gene] ) {
      if (buffer_start_bp == -1) {
        buffer_start_bp = ld_block.start_bp;
        buffer_end_bp = ld_block.end_bp;
      } else {
        buffer_start_bp = min(buffer_start_bp, ld_block.start_bp);
        buffer_end_bp = max(buffer_end_bp, ld_block.end_bp);
      }
    }

    while (!skip_buffer && buffer_start_bp > ld_block_map.sorted_ld_block_pairs[ld_block_idx].max_pos()
           && ld_block_idx < (int)ld_block_map.sorted_ld_block_pairs.size()) {
      t = std::chrono::steady_clock::now();
      LDBlockPair block_pair = ld_block_map.sorted_ld_block_pairs[ld_block_idx];

      block_pgen_region_t region1 = buffer_reader.get_region(
        chrom,
        block_pair.ld_block1.start_bp,
        block_pair.ld_block1.end_bp
      );
      block_pgen_region_t region2 = buffer_reader.get_region(
        chrom,
        block_pair.ld_block2.start_bp,
        block_pair.ld_block2.end_bp
      );

      log_info("computing LD matrix for block pairs "
        + format_region(chrom, block_pair.ld_block1.start_bp, block_pair.ld_block1.end_bp)
        + " "
        + format_region(chrom, block_pair.ld_block2.start_bp, block_pair.ld_block2.end_bp)
      );

      unordered_set<variant_id> block1_variants_inserted;
      unordered_set<variant_id> block2_variants_inserted;
      for (int i = 0; i < buffer_reader.get_nblocks(block_size, region1); ++i) {
        block_block_ld_t block_block_ld;
        buffer_reader.load_block(Gblock1, block1_variants, nvar_read1, i, block_size, region1);
        int j_end;
        if (block_pair.ld_block1 == block_pair.ld_block2) {
          j_end = i+1;
        } else {
          j_end = buffer_reader.get_nblocks(block_size, region2);
        }
        for (int j = 0; j < j_end; ++j) {
          vector<ref_ld_matrix_entry_t> block_matrix_entries;
          buffer_reader.load_block(Gblock2, block2_variants, nvar_read2, j, block_size, region2);
          extract_block_corr(block_matrix_entries,
                             block1_variants_inserted,
                             block2_variants_inserted,
                             Gblock1,
                             Gblock2,
                             block1_variants,
                             block2_variants,
                             nvar_read1,
                             nvar_read2,
                             all_buffer_variants,
                             all_buffer_variants,
                             buffer_r2,
                             block_pair.ld_block1 == block_pair.ld_block2 && i == j);

          for ( const ref_ld_matrix_entry_t& e :  block_matrix_entries ) {
            if (all_buffer_variants.count(e.row_id) > 0 && all_buffer_variants.count(e.col_id) > 0) {
              ++buffer_records_written;
              block_block_ld.data.push_back(e);
            }
          }
        }

        int64_t addr = out.write_block_block_ld_to_buffer_ld_file(
          block_block_ld, block_pair.ld_block1 == block_pair.ld_block2
        );

        for ( const string& g : ld_block_map.ld_block_genes[block_pair.ld_block1] ) {
          if (gene_first_block_addr.count(g) == 0 && genes_skipped.count(g) == 0) {
            if (addr == -1) {
              log_error("bad address to LD block: this is a bug", 1);
            }
            gene_first_block_addr[g] = addr;
            out.write_gene_to_index(g, addr);
          }
        }
      }
      ++ld_block_idx;
      d = std::chrono::steady_clock::now() - t;
      log_info("finished in " + to_string(d.count()) + " seconds");
    }

    log_info("computing LD matrix for gene " + gene);
    log_info(" * target region: " + chr + ":" + to_string(start_bp) + "-" + to_string(end_bp));

    if (!skip_buffer) {
      log_info(" * buffer region: " + chr + ":" 
        + to_string(buffer_start_bp) + "-" + to_string(buffer_end_bp));
    }

    gene_start = std::chrono::steady_clock::now();
    target_reader.set_region(chr, start_bp, end_bp);
    buffer_reader.set_region(chr, buffer_start_bp, buffer_end_bp);
    unordered_map<variant_id, ld_variant_t> buffer_variants;
    unordered_map<variant_id, ld_variant_t> target_variants;
    vector<variant_id> target_variants_seen;
    vector<variant_id> buffer_variants_seen;

    log_info(" * computing target variances");
    t = std::chrono::steady_clock::now();
    for (int i = 0; i < target_reader.get_nblocks(block_size); ++i) {
      target_reader.load_block(Gblock1, block1_variants, nvar_read1, i, block_size);
      VectorXf block_mean = compute_block_mean(Gblock1);
      VectorXf block_variance = compute_block_variance(Gblock1, block_mean);
      for (int j = 0; j < nvar_read1; ++j) {
        if (block_mean(j) > 0) {
          target_variants[block1_variants[j].id] = {
            block1_variants[j].id,
            HTPv4Pos(block1_variants[j].chr, block1_variants[j].pos),
            block_variance(j),
            block_mean(j)
          };
          target_variants_seen.push_back(block1_variants[j].id);
        }
      }
    }
    d = std::chrono::steady_clock::now() - t;
    log_info("     done (" + to_string(d.count()) + " seconds)");
    log_info(" * found " + to_string(target_variants.size()) + " target variants");

    if (target_variants.size() == 0) {
      log_warning("gene " + gene + " does not have any variants in target pfile, skipping...");
      genes_skipped.insert(gene);
      continue;
    }

    std::chrono::time_point<std::chrono::steady_clock> t1;
    std::chrono::duration<double> d1 = t1 - t1;

    if (!skip_buffer) {
      log_info(" * computing buffer variances");
      t = std::chrono::steady_clock::now();
      for (int i = 0; i < buffer_reader.get_nblocks(block_size); ++i) {
        buffer_reader.load_block(Gblock1, block1_variants, nvar_read1, i, block_size);
        vector<int> new_buffer_idx;
        for (int j = 0; j < nvar_read1; ++j) {
          if ((target_variants.count(block1_variants[j].id) == 0) && (buffer_variants.count(block1_variants[j].id) == 0)) {
            new_buffer_idx.push_back(j);
            buffer_variants_seen.push_back(block1_variants[j].id);
          } else if (buffer_variants.count(block1_variants[j].id) > 0) {
            buffer_variants[block1_variants[j].id] 
              = all_buffer_variants.at(block1_variants[j].id);
            buffer_variants_seen.push_back(block1_variants[j].id);
          }
        }

        VectorXf block_mean = compute_block_mean(Gblock1, new_buffer_idx);
        VectorXf block_variance = compute_block_variance(Gblock1, block_mean, new_buffer_idx);
        for (size_t k = 0; k < new_buffer_idx.size(); ++k) {
          int j = new_buffer_idx[k];
          buffer_variants[block1_variants[j].id] = {
            block1_variants[j].id,
            HTPv4Pos(block1_variants[j].chr, block1_variants[j].pos),
            block_variance(k),
            block_mean(k)
          };
        }
      }
      d = std::chrono::steady_clock::now() - t;
      log_info("     done (" + to_string(d.count()) + " seconds)");
      log_info(" * found " + to_string(buffer_variants.size()) + " candidate variants in buffer");
    }

    unordered_set<variant_id> target_variants_inserted;
    vector<ref_ld_matrix_entry_t> target_matrix_entries;
    log_info(" * computing target-target covariance");
    t = std::chrono::steady_clock::now();
    for (int i = 0; i < target_reader.get_nblocks(block_size); ++i) {
      target_reader.load_block(Gblock1, block1_variants, nvar_read1, i, block_size);
      for (int j = 0; j <= i; ++j) {
        target_reader.load_block(Gblock2, block2_variants, nvar_read2, j, block_size);
        extract_block_corr(target_matrix_entries,
                           target_variants_inserted,
                           target_variants_inserted,
                           Gblock1,
                           Gblock2,
                           block1_variants,
                           block2_variants,
                           nvar_read1,
                           nvar_read2,
                           target_variants,
                           target_variants,
                           target_r2,
                           i == j);
      }
    }
    d = std::chrono::steady_clock::now() - t;
    log_info("     done (" + to_string(d.count()) + " seconds)");

    unordered_set<variant_id> buffer_variants_inserted;
    vector<ref_ld_matrix_entry_t> buffer_target_matrix_entries;
    if (!skip_buffer) {
      log_info(" * computing buffer-target covariance");
      t = std::chrono::steady_clock::now();
      for (int i = 0; i < buffer_reader.get_nblocks(block_size); ++i) {
        buffer_reader.load_block(Gblock1, block1_variants, nvar_read1, i, block_size);
        for (int j = 0; j < target_reader.get_nblocks(block_size); ++j) {;
          target_reader.load_block(Gblock2, block2_variants, nvar_read2, j, block_size);
          extract_block_corr(buffer_target_matrix_entries,
                             buffer_variants_inserted,
                             target_variants_inserted,
                             Gblock1,
                             Gblock2,
                             block1_variants,
                             block2_variants,
                             nvar_read1,
                             nvar_read2,
                             buffer_variants,
                             target_variants,
                             buffer_r2,
                             false);
        }
      }
      d = std::chrono::steady_clock::now() - t;
      log_info("     done (" + to_string(d.count()) + " seconds)");
      log_info(" * keeping " + to_string(buffer_variants_inserted.size()) + " buffer variants in LD with target variants");
    }

    log_info(" * writing gene LD output file");
    t = std::chrono::steady_clock::now();
    vector<variant_id> gene_variants;
    vector<float> gene_variant_variances;
    vector<variant_id> window_variants;
    vector<float> window_variant_variances;
    for ( const variant_id& vid : target_variants_seen ) {
      gene_variants.push_back(vid);
      gene_variant_variances.push_back((float)target_variants[vid].variance);
    }
    gene_variants_saved += target_variants_seen.size();
    for ( const variant_id& vid : buffer_variants_seen ) {
      if (target_variants.count(vid) == 0 && buffer_variants_inserted.count(vid) > 0) {
        window_variants.push_back(vid);
        window_variant_variances.push_back((float)buffer_variants[vid].variance);
        all_buffer_variants[vid] = buffer_variants[vid];
      }
    }
    out.init_next_gene(gene,
                       target_r2,
                       buffer_r2,
                       gene_variants,
                       gene_variant_variances,
                       window_variants,
                       window_variant_variances,
                       (int32_t)ld_block_map.gene_ld_blocks[gene].size());

    for (const ref_ld_matrix_entry_t& e : target_matrix_entries) {
      if (target_variants.count(e.row_id) > 0 && target_variants.count(e.col_id) > 0) {
        ++gene_records_written;
        out.add_gene_gene_ld_entry(e);
      }
    }

    for ( const ref_ld_matrix_entry_t& e : buffer_target_matrix_entries ) {
      if (buffer_variants_inserted.count(e.row_id) > 0) {
        ++gene_buffer_records_written;
        out.add_gene_buffer_ld_entry(e);
      }
    }

    out.write_gene_to_gene_ld_file();
    d = std::chrono::steady_clock::now() - t;

    if (skip_buffer) {
      int64_t addr = out.write_block_block_ld_to_buffer_ld_file(
        block_block_ld_t {}, true
      );
      out.write_gene_to_index(gene, addr);
    }

    if (!skip_buffer && buffer_variants_inserted.size() == 0) {
      int64_t addr = out.write_block_block_ld_to_buffer_ld_file(
        block_block_ld_t {}, true
      );
      gene_first_block_addr[gene] = addr;
      out.write_gene_to_index(gene, addr);
    }
    log_info("     done (" + to_string(d.count()) + " seconds)");


    d = std::chrono::steady_clock::now() - gene_start;
    log_info("finished in " + to_string(d.count()) + " seconds");
  }

  log_info("computing remaining buffer LD blocks");
  while (!skip_buffer && ld_block_idx < (int)ld_block_map.sorted_ld_block_pairs.size()) {
    t = std::chrono::steady_clock::now();
    LDBlockPair block_pair = ld_block_map.sorted_ld_block_pairs[ld_block_idx];

    block_pgen_region_t region1 = buffer_reader.get_region(
      chrom,
      block_pair.ld_block1.start_bp,
      block_pair.ld_block1.end_bp
    );
    block_pgen_region_t region2 = buffer_reader.get_region(
      chrom,
      block_pair.ld_block2.start_bp,
      block_pair.ld_block2.end_bp
    );

    log_info("computing LD matrix for block pairs "
      + format_region(chrom, block_pair.ld_block1.start_bp, block_pair.ld_block1.end_bp)
      + " "
      + format_region(chrom, block_pair.ld_block2.start_bp, block_pair.ld_block2.end_bp)
    );

    unordered_set<variant_id> block1_variants_inserted;
    unordered_set<variant_id> block2_variants_inserted;
    for (int i = 0; i < buffer_reader.get_nblocks(block_size, region1); ++i) {
      block_block_ld_t block_block_ld;
      
      buffer_reader.load_block(Gblock1, block1_variants, nvar_read1, i, block_size, region1);
      int j_end;
      if (block_pair.ld_block1 == block_pair.ld_block2) {
        j_end = i+1;
      } else {
        j_end = buffer_reader.get_nblocks(block_size, region2);
      }
      for (int j = 0; j < j_end; ++j) {
        vector<ref_ld_matrix_entry_t> block_matrix_entries;
        buffer_reader.load_block(Gblock2, block2_variants, nvar_read2, j, block_size, region2);
        extract_block_corr(block_matrix_entries,
                            block1_variants_inserted,
                            block2_variants_inserted,
                            Gblock1,
                            Gblock2,
                            block1_variants,
                            block2_variants,
                            nvar_read1,
                            nvar_read2,
                            all_buffer_variants,
                            all_buffer_variants,
                            buffer_r2,
                            block_pair.ld_block1 == block_pair.ld_block2 && i == j);

        for ( const ref_ld_matrix_entry_t& e :  block_matrix_entries ) {
          if (all_buffer_variants.count(e.row_id) > 0 && all_buffer_variants.count(e.col_id) > 0) {
            ++buffer_records_written;
            block_block_ld.data.push_back(e);
          }
        }
      }

      int64_t addr = out.write_block_block_ld_to_buffer_ld_file(
        block_block_ld, block_pair.ld_block1 == block_pair.ld_block2
      );

      for ( const string& g : ld_block_map.ld_block_genes[block_pair.ld_block1] ) {
        if (gene_first_block_addr.count(g) == 0 && genes_skipped.count(g) == 0) {
          if (addr == -1) {
            log_error("bad address to LD block: this is a bug", 1);
          }
          gene_first_block_addr[g] = addr;
          out.write_gene_to_index(g, addr);
        }
      }
    }
    ++ld_block_idx;
    d = std::chrono::steady_clock::now() - t;
    log_info("finished in " + to_string(d.count()) + " seconds");
  }

  out.close();
  d = std::chrono::steady_clock::now() - start;
  log_info("total runtime: " + to_string(d.count() / 60.0) + " minutes");
  log_info("gene variants saved: " + to_string(gene_variants_saved));
  log_info("buffer variants saved: " + to_string(all_buffer_variants.size()));
  log_info("gene-gene entries written: " + to_string(gene_records_written));
  log_info("gene-buffer entries written: " + to_string(gene_buffer_records_written));
  log_info("buffer-buffer entries written: " + to_string(buffer_records_written));
}