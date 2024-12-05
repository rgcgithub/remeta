#ifndef RUN_HTP_H
#define RUN_HTP_H

#include "lapack_complex.hpp"

#include <string>
#include <vector>

#include "io/htpv4_reader.hpp"
#include "enums.hpp"

void run_htp(
  const std::vector<std::string>& htp_files,
  const std::vector<std::string>& ld_prefixes,
  const std::vector<std::string>& cohorts,
  const std::string& anno_file,
  const std::string& set_list_file,
  const std::string& mask_def_file,
  const std::string& trait_name,
  const std::string& trait_type,
  const std::string& out_prefix,
  const std::vector<double>& burden_aaf_bins,
  const std::string& burden_singleton_def,
  const std::string& burden_weight_strategy,
  const bool& skip_burden,
  const double& skato_max_aaf,
  const std::vector<double> skato_rho_values,
  const int& skato_min_aac,
  const std::string& skato_weight_strategy,
  const bool& skip_skato,
  const double& acatv_max_aaf,
  const int& acatv_min_aac,
  const std::string& acatv_weight_strategy,
  const bool& skip_acatv,
  const std::string& condition_list_file,
  const std::vector<std::string>& condition_htp_files,
  const std::string& af_strategy,
  const std::string& af_file,
  const double& spa_pval,
  const double& spa_ccr,
  const std::string& chr,
  const std::string& gene,
  const std::string& extract_file,
  const std::string& exclude_file,
  const std::vector<std::string>& sources,
  const int& threads,
  const bool& write_cohort_burden_tests,
  const bool& write_mask_snplist,
  const bool& ignore_mask_ld,
  const bool& recompute_score,
  const bool& keep_variants_not_in_ld_mat
);

#endif
