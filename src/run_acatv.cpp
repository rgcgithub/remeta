#include "run_acatv.hpp"

#include <string>
#include <vector>
using std::vector;
using std::string;

#if defined(WITH_MKL)
#include <mkl.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "io/gene_set_reader.hpp"
#include "io/htpv4_reader.hpp"
#include "io/htpv4_writer.hpp"
#include "meta/htp_meta_analyzer.hpp"
#include "logging.hpp"
#include "htpv4_pos.hpp"
#include "parameter_checks.hpp"
namespace pc = parameter_checks;
#include "variant_filter.hpp"

void run_acatv(const vector<string>& htp_files,
               const vector<string>& cohorts,
               const string& anno_file,
               const string& set_list_file,
               const string& mask_def_file,
               const string& trait_name,
               const string& trait_type,
               const string& out_prefix,
               const string& af_strategy,
               const string& af_file,
               const string& weight_strategy,
               const int& min_aac,
               const double& max_aaf,
               const string& chr,
               const string& gene,
               const string& extract_file,
               const string& exclude_file,
               const vector<string>& sources,
               const vector<string>& ld_prefixes,
               const string& condition_list_file,
               const vector<string>& condition_htp_files,
               const int& threads,
	       const bool& write_mask_snplist) {
  pc::check_htp_files(htp_files);
  pc::check_file_exists(anno_file);
  pc::check_file_exists(set_list_file);
  pc::check_file_exists(mask_def_file);
  pc::check_anno_file_indexed(anno_file);
  pc::check_cohorts(cohorts, htp_files);
  pc::check_aac(min_aac);
  pc::check_af(max_aaf);

  trait_type_e trait = pc::check_trait_type(trait_type);
  weight_strategy_e weight = pc::check_weight_strategy(weight_strategy);

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
  vf.set_effect_not_na();
  vf.set_info_has_se();
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
    0,
    0
  );
  meta.set_run_acatv({max_aaf}, weight, min_aac);

  if (af_file != "") {
    meta.set_af_file(af_file);
  }

  if (condition_htp_files.size() > 0 || condition_list_file != "") {
    log_info("loading conditional variants...");
    pc::check_ld_files(ld_prefixes, htp_files.size());
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
      log_error("no variants found for conditional analysis", 1);
    }

    meta.set_conditional_variants(conditional_variants);
  }

  BgzWriter snplist_out;
  if (write_mask_snplist) {
    snplist_out.open(out_prefix + ".acatv.snplist.remeta.gz", "w");
    meta.set_write_mask_snplist(snplist_out);
  }

  HTPv4Writer out(out_prefix + ".acatv.remeta.gz", "w");
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

  snplist_out.close();
}
