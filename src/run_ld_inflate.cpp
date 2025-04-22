#include "run_ld_inflate.hpp"

#include <sstream>
using namespace std;

#include <boost/format.hpp>

#include <Eigen/Dense>
using Eigen::MatrixXf;

#include "io/bgz_writer.hpp"
#include "io/remeta_matrix_reader.hpp"
#include "parameter_checks.hpp"
namespace pc = parameter_checks;

void run_ld_inflate(const string& mat_prefix,
                    const string& out_prefix,
                    const string& gene_name,
                    const string& extract) {
  pc::check_ld_files(vector<string>{mat_prefix}, 1);

  RemetaMatrixReader reader;
  reader.open(mat_prefix);
  if (!reader.contains_gene(gene_name)) {
    log_error(gene_name + " not found", 1);
  }

  vector<string> variants;
  unordered_map<string, int> variant_idx;
  BgzReader extract_file(extract);
  string line;
  while (!extract_file.eof()) {
    line = extract_file.readline();
    variants.push_back(line);
    variant_idx[line] = (int)variants.size() - 1;
  }
  extract_file.close();
  if (variants.size() == 0) {
    log_error("no variants found in " + extract, 1);
  }

  vector<string> gene_variants;
  reader.load_gene_variant_ids(gene_variants, gene_name);

  vector<string> buffer_variants;
  reader.load_buffer_variant_ids(buffer_variants, gene_name);

  vector<string> gene_variants_to_load;
  vector<string> buffer_variants_to_load;
  for (const string& vid: variants) {
    if (find(gene_variants.begin(), gene_variants.end(), vid) != gene_variants.end()) {
      gene_variants_to_load.push_back(vid);
    } else if (find(buffer_variants.begin(), buffer_variants.end(), vid) != buffer_variants.end()) {
      buffer_variants_to_load.push_back(vid);
    }
  }
  size_t variants_to_load = gene_variants_to_load.size() + buffer_variants_to_load.size();
  log_info("found " + to_string(variants_to_load) + " of " + to_string(variants.size()) + " variants in the LD matrix");
  if (variants_to_load == 0) {
    log_error("none of the variants in " + extract + " are in the LD matrix", 1);
  }

  log_info("loading LD matrix...");
  MatrixXf G, C, G_C;
  reader.load_conditional_ld_mats(G,
                                  C,
                                  G_C,
                                  gene_variants_to_load,
                                  buffer_variants_to_load,
                                  gene_name);

  MatrixXf out_matrix = MatrixXf::Zero(variants_to_load, variants_to_load);
  int u1, u2;
  for (size_t v1 = 0; v1 < gene_variants_to_load.size(); ++v1) {
    u1 = variant_idx[gene_variants_to_load[v1]];
    for (size_t v2 = 0; v2 < gene_variants_to_load.size(); ++v2) {
      u2 = variant_idx[gene_variants_to_load[v2]];
      out_matrix(u1, u2) = G(v1, v2);
      out_matrix(u2, u1) = out_matrix(u1, u2);
    }

    for (size_t v2 = 0; v2 < buffer_variants_to_load.size(); ++v2) {
      u2 = variant_idx[buffer_variants_to_load[v2]];
      out_matrix(u1, u2) = G_C(v1, v2);
      out_matrix(u2, u1) = out_matrix(u1, u2);
    }
  }

  for (size_t v1 = 0; v1 < buffer_variants_to_load.size(); ++v1) {
    u1 = variant_idx[buffer_variants_to_load[v1]];
    for (size_t v2 = 0; v2 < buffer_variants_to_load.size(); ++v2) {
      u2 = variant_idx[buffer_variants_to_load[v2]];
      out_matrix(u1, u2) = C(v1, v2);
      out_matrix(u2, u1) = out_matrix(u1, u2);
    }
  }

  log_info("writing output...");
  BgzWriter out(out_prefix + ".raw.ld.gz", "w");
  for (ssize_t v1 = 0; v1 < out_matrix.rows(); ++v1) {
    for (ssize_t v2 = 0; v2 < out_matrix.cols()-1; ++v2) {
      out.write((boost::format("%1$g") % out_matrix(v1, v2)).str() + "\t");
    }
    out.write((boost::format("%1$g") % out_matrix(v1, out_matrix.cols()-1)).str() + "\n");
  }
  out.close();

  out.open(out_prefix + ".snplist.gz", "w");
  for (const string& vid: variants) {
    out.write(vid + "\n");
  }
  out.close();

  log_info("done!");
}