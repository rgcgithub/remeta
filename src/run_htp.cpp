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

#include "io/allele_freq_writer.hpp"
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
  const vector<string>& htp_files,
  const vector<string>& ld_prefixes,
  const vector<string>& cohorts,
  const string& anno_file,
  const string& set_list_file,
  const string& mask_def_file,
  const string& trait_name,
  const string& trait_type,
  const string& out_prefix,
  const vector<double>& burden_aaf_bins,
  const string& burden_singleton_def,
  const string& burden_weight_strategy,
  double burden_mask_spa_pval,
  double burden_mask_spa_ccr,
  double burden_sv_spa_pval,
  double burden_sv_spa_ccr,
  bool skip_burden,
  double skato_max_aaf,
  const vector<double>& skato_rho_values,
  int skato_min_aac,
  const string& skato_weight_strategy,
  double skato_mask_spa_pval,
  double skato_mask_spa_ccr,
  double skato_sv_spa_pval,
  double skato_sv_spa_ccr,
  int skato_collapse_below_aac,
  bool skip_skato,
  double acatv_max_aaf,
  int acatv_min_aac,
  const string& acatv_weight_strategy,
  double acatv_sv_spa_pval,
  double acatv_sv_spa_ccr,
  bool skip_acatv,
  const string& condition_list_file,
  const vector<string>& condition_htp_files,
  int max_cond_var_per_gene,
  const string& af_strategy,
  const string& af_file,
  const string& chr,
  const string& gene,
  const string& extract_file,
  const vector<string>& cohort_extract_files,
  const string& exclude_file,
  const vector<string>& sources,
  int threads,
  bool write_cohort_burden_tests,
  bool write_mask_snplist,
  bool ignore_mask_ld,
  bool recompute_score,
  bool keep_variants_not_in_ld_mat,
  bool write_freqs,
  const string& collapse_anno_name
) {
  pc::check_htp_files(htp_files);
  pc::check_file_exists(anno_file);
  pc::check_file_exists(set_list_file);
  pc::check_file_exists(mask_def_file);
  pc::check_anno_file_indexed(anno_file);
  pc::check_pval(burden_mask_spa_pval);
  pc::check_pval(burden_sv_spa_pval);
  pc::check_case_control_ratio(burden_mask_spa_ccr);
  pc::check_case_control_ratio(burden_sv_spa_ccr);
  pc::check_pval(skato_mask_spa_pval);
  pc::check_pval(skato_sv_spa_pval);
  pc::check_case_control_ratio(skato_mask_spa_ccr);
  pc::check_case_control_ratio(skato_sv_spa_ccr);
  pc::check_pval(acatv_sv_spa_pval);
  pc::check_case_control_ratio(acatv_sv_spa_ccr);
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
  vector<string> parsed_cohort_extract_files(htp_files.size(), "");
  if (cohort_extract_files.size() > 0) {
    parsed_cohort_extract_files = pc::check_cohort_extract_files(cohort_extract_files, cohorts);
  }
  if (exclude_file != "") {
    pc::check_file_exists(exclude_file);
  }

  VariantFilter vf;
  vf.set_effect_not_na();
  if (!recompute_score) {
    vf.set_info_has_score();
  }
  if (extract_file != "") {
    vf.set_extract_file(extract_file);
  }
  for (size_t i = 0; i < parsed_cohort_extract_files.size(); ++ i) {
    if (parsed_cohort_extract_files[i] != "") {
      vf.set_cohort_extract_file(parsed_cohort_extract_files[i], i);
    }
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
    ld_prefixes
  );

  if (!skip_burden) {
    weight_strategy_e burden_weights = pc::check_weight_strategy(burden_weight_strategy);
    singleton_def_e singleton = pc::check_singleton_def(burden_singleton_def);
    meta.set_run_burden(
      burden_aaf_bins,
      singleton,
      burden_weights,
      burden_mask_spa_pval,
      burden_mask_spa_ccr,
      burden_sv_spa_pval,
      burden_sv_spa_ccr
    );
  }
  if (!skip_skato) {
    weight_strategy_e skato_weights = pc::check_weight_strategy(skato_weight_strategy);
    meta.set_run_skato(
      vector<double>{skato_max_aaf},
      skato_rho_values,
      skato_weights,
      skato_min_aac,
      skato_collapse_below_aac,
      skato_mask_spa_pval,
      skato_mask_spa_ccr,
      skato_sv_spa_pval,
      skato_sv_spa_ccr
    );
  }
  if (!skip_acatv) {
    weight_strategy_e acatv_weights = pc::check_weight_strategy(acatv_weight_strategy);
    meta.set_run_acatv(
      vector<double>{acatv_max_aaf},
      acatv_weights,
      acatv_min_aac,
      acatv_sv_spa_pval,
      acatv_sv_spa_ccr
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

    meta.set_conditional_variants(conditional_variants, max_cond_var_per_gene);
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
  if (collapse_anno_name != "") {
    meta.set_collapse_anno(collapse_anno_name);
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

  AlleleFreqWriter freq_out;
  if (write_freqs) {
    freq_out.open(out_prefix + ".freqs.remeta.gz");
    meta.set_write_freqs(freq_out);
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
