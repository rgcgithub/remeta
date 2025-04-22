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
  double burden_mask_spa_pval,
  double burden_mask_spa_ccr,
  double burden_sv_spa_pval,
  double burden_sv_spa_ccr,
  bool skip_burden,
  double skato_max_aaf,
  const std::vector<double>& skato_rho_values,
  int skato_min_aac,
  const std::string& skato_weight_strategy,
  double skato_mask_spa_pval,
  double skato_mask_spa_ccr,
  double skato_sv_spa_pval,
  double skato_sv_spa_ccr,
  bool skip_skato,
  double acatv_max_aaf,
  int acatv_min_aac,
  const std::string& acatv_weight_strategy,
  double acatv_sv_spa_pval,
  double acatv_sv_spa_ccr,
  bool skip_acatv,
  const std::string& condition_list_file,
  const std::vector<std::string>& condition_htp_files,
  int max_cond_var_per_gene,
  const std::string& af_strategy,
  const std::string& af_file,
  const std::string& chr,
  const std::string& gene,
  const std::string& extract_file,
  const std::string& exclude_file,
  const std::vector<std::string>& sources,
  int threads,
  bool write_cohort_burden_tests,
  bool write_mask_snplist,
  bool ignore_mask_ld,
  bool recompute_score,
  bool keep_variants_not_in_ld_mat,
  bool write_freqs
);

#endif
