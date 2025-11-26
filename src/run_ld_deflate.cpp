#include "run_ld_deflate.hpp"

#include <string>
using std::string;

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "io/bgz_reader.hpp"
#include "io/ref_ld_matrix_writer.hpp"
#include "logging.hpp"
#include "util.hpp"

void run_ld_deflate(const string& ld_file,
                    const int& sample_size,
                    const double& sparsity_threshold,
                    const string& out) {
  BgzReader ld_reader(ld_file);
  RefLDMatrixWriter ld_writer(out, 1);

  while (!ld_reader.eof()) {
    string gene = ld_reader.readline();
    log_info("processing " + gene);

    string chrom = ld_reader.readline();
    vector<string> variant_ids = util::str_split(ld_reader.readline(), " ");
    vector<string> tmp = util::str_split(ld_reader.readline(), " ");
    if (variant_ids.size() == 0) {
      log_warning("gene " + gene + " does not have any variants");
      ld_reader.readline();
      continue;
    }
    if (tmp.size() != variant_ids.size()) {
      log_error("number of variants does not match the number of variances", 1);
    }

    vector<float> variances(variant_ids.size());
    for (size_t i = 0; i < variant_ids.size(); ++i) {
      try {
        variances[i] = stof(tmp[i]);
      } catch (const std::invalid_argument& e) {
        log_error("invalid variance value " + tmp[i] + " for gene " + gene, 1);
      }

      if (variances[i] <= 0) {
        log_error("found non-positive variance for variant " + variant_ids[i] + " in gene " + gene, 1);
      }
    }
    ld_writer.init_next_gene(gene,
                             sparsity_threshold,
                             sparsity_threshold,
                             variant_ids,
                             variances,
                             vector<variant_id>(),
                             vector<float>(),
                             0);

    string n_ld_entries_str = ld_reader.readline();
    int n_ld_entries = 0;
    try {
      n_ld_entries = stoi(n_ld_entries_str);
    } catch (const std::invalid_argument& e) {
      log_error("could not read number of ld entries for gene " + gene, 1);
    }

    for (int i = 0; i < n_ld_entries; ++i) {
      string ld_entry = ld_reader.readline();
      vector<string> ld_entries = util::str_split(ld_entry, " ");
      if (ld_entries.size() != 3) {
        log_error("invalid ld entry for gene " + gene + "(" + ld_entry + ")", 1);
      }
      try {
        int32_t i = stoi(ld_entries[0]);
        int32_t j = stoi(ld_entries[1]);
        float corr = stof(ld_entries[2]);
        if (abs(corr*corr) >= sparsity_threshold) {
          ld_writer.add_gene_gene_ld_entry(ref_ld_matrix_entry_t{variant_ids[i], variant_ids[j], corr});
        }
      } catch (const std::invalid_argument& e) {
        log_error("could not parse ld entries " + to_string(i + 1) + " for gene " + gene, 1);
      }
    }
    ld_writer.write_gene_to_gene_ld_file();

    int64_t addr = ld_writer.write_block_block_ld_to_buffer_ld_file(
        block_block_ld_t {}, true
    );
    ld_writer.write_gene_to_index(gene, addr);
  }
  ld_writer.close();
  log_info("done!");
}