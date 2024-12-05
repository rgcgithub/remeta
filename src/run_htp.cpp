#include "run_htp.hpp"

#if defined(WITH_MKL)
#include <mkl.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <string>
#include <vector>
#include <unordered_set>

#include "io/bgz_reader.hpp"
#include "io/gene_set_reader.hpp"
#include "io/htpv4_reader.hpp"
#include "io/htpv4_writer.hpp"
#include "io/regenie_ld_matrix_reader.hpp"
#include "meta/htp_meta_analyzer.hpp"
#include "meta/set_meta_analyzer.hpp"
#include "enums.hpp"
#include "htpv4_pos.hpp"
#include "logging.hpp"
#include "parameter_checks.hpp"
namespace pc = parameter_checks;

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
) {
  pc::check_htp_files(htp_files);
  pc::check_file_exists(anno_file);
  pc::check_file_exists(set_list_file);
  pc::check_file_exists(mask_def_file);
  pc::check_anno_file_indexed(anno_file);
  pc::check_pval(spa_pval);
  pc::check_case_control_ratio(spa_ccr);
  pc::check_cohorts(cohorts, htp_files);
  pc::check_aac(skato_min_aac);
  pc::check_af(skato_max_aaf);
  pc::check_aac(acatv_min_aac);
  pc::check_af(acatv_max_aaf);
  pc::check_aaf_bins(burden_aaf_bins);

  bool ld_required = !ignore_mask_ld || !keep_variants_not_in_ld_mat || condition_list_file != "";
  if (ld_required) {
    pc::check_ld_files(ld_prefixes, htp_files.size());
  }

  trait_type_e trait = pc::check_trait_type(trait_type);

  af_strategy_e af;
  if (af_file != "") {
    pc::check_file_exists(af_file);
    af = USE_EXTERNAL_AF;
  } else {
    af = pc::check_af_strategy(af_strategy);
  }

  if (extract_file != "") {
    pc::check_file_exists(extract_file);
  }
  if (exclude_file != "") {
    pc::check_file_exists(exclude_file);
  }

  VariantFilter vf;
  if (!recompute_score) {
    vf.set_info_has_score();
  }
  if (extract_file != "") {
    vf.set_extract_file(extract_file);
  }
  if (exclude_file != "") {
    vf.set_exclude_file(exclude_file);
  }
  if (sources.size() > 0) {
    vf.set_info_source_is_one_of(sources);
  }

#if defined(WITH_MKL)
  mkl_set_num_threads(threads);
#endif

#if defined(_OPENMP)
  omp_set_num_threads(threads);
  Eigen::setNbThreads(threads);
#endif

  HTPMetaAnalyzer meta(
    trait_name,
    trait,
    cohorts,
    af,
    mask_def_file,
    anno_file,
    ld_prefixes,
    spa_pval,
    spa_ccr
  );

  if (!skip_burden) {
    weight_strategy_e burden_weights = pc::check_weight_strategy(burden_weight_strategy);
    singleton_def_e singleton = pc::check_singleton_def(burden_singleton_def);
    meta.set_run_burden(
      burden_aaf_bins,
      singleton,
      burden_weights
    );
  }
  if (!skip_skato) {
    weight_strategy_e skato_weights = pc::check_weight_strategy(skato_weight_strategy);
    meta.set_run_skato(
      vector<double>{skato_max_aaf},
      skato_rho_values,
      skato_weights,
      skato_min_aac
    );
  }
  if (!skip_acatv) {
    weight_strategy_e acatv_weights = pc::check_weight_strategy(acatv_weight_strategy);
    meta.set_run_acatv(
      vector<double>{acatv_max_aaf},
      acatv_weights,
      acatv_min_aac
    );
  }

  if (af_file != "") {
    meta.set_af_file(af_file);
  }

  if (condition_list_file != "") {
    log_info("loading conditional variants...");
    pc::check_conditional_analysis_files(condition_list_file, condition_htp_files, htp_files.size());
    vector<vector<htpv4_record_t> > conditional_variants =
      util::get_htp_variants(condition_list_file, condition_htp_files, chr);

    bool found_one_conditional_variant = false;
    for (size_t i = 0; i < conditional_variants.size(); ++ i) {
      log_info("found " + to_string(conditional_variants[i].size()) + " variants to condition on in study " + to_string(i+1));
      if (conditional_variants[i].size() > 0) {
        found_one_conditional_variant = true;
      }
    }
    if (!found_one_conditional_variant) {
      log_warning("no variants found for conditional analysis");
    }

    meta.set_conditional_variants(conditional_variants);
  }

  if (recompute_score) {
    meta.set_recompute_score();
  }
  if (keep_variants_not_in_ld_mat) {
    meta.set_keep_variants_not_in_ld_mat();
  }
  if (ignore_mask_ld) {
    meta.set_ignore_mask_ld();
  }

  HTPv4Writer burden_out; 
  if (write_cohort_burden_tests) {
    burden_out.open(out_prefix + ".burden_per_cohort.remeta.gz", "w");
    meta.set_write_cohort_burden_tests(burden_out);
  }

  BgzWriter snplist_out;
  if (write_mask_snplist) {
    snplist_out.open(out_prefix + ".snplist.remeta.gz", "w");
    meta.set_write_mask_snplist(snplist_out);
  }

  HTPv4Writer out(out_prefix + ".remeta.gz", "w");
  out.writeheader();
  vector<HTPv4Reader> readers;
  for (const string& htp_file : htp_files) {
    readers.push_back(HTPv4Reader(htp_file));
  }

  run_set_meta_analysis(meta,
                        vf,
                        readers,
                        set_list_file,
                        out,
                        chr,
                        gene);

  burden_out.close();
  snplist_out.close();
}
