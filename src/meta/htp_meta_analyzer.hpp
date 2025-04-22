#ifndef HTP_META_ANALYZER_H
#define HTP_META_ANALYZER_H

#include "../lapack_complex.hpp"

#include <string>
#include <vector>

#include <Eigen/Dense>

#include "set_meta_analyzer.hpp"
#include "../io/allele_freq_writer.hpp"
#include "../io/gene_set_reader.hpp"
#include "../io/bgz_reader.hpp"
#include "../io/bgz_writer.hpp"
#include "../io/htpv4_reader.hpp"
#include "../io/remeta_matrix_reader.hpp"
#include "../stat/misc.hpp"
#include "../stat/spa.hpp"
#include "../stat/tests.hpp"
#include "../allele_frequency_map.hpp"
#include "../annotation_map.hpp"
#include "../enums.hpp"
#include "../mask_set.hpp"
#include "../logging.hpp"
#include "../util.hpp"

typedef std::string variant_id;

enum htp_test_e {
  SKATO,
  BURDEN,
  ACATV
};

const double SPA_NEVER_APPLIED = std::numeric_limits<double>::lowest();
const double SPA_ALWAYS_APPLIED = std::numeric_limits<double>::max();

struct variant_metadata {
  double aaf;
  double weight;
  bool   is_singleton;
};

class HTPMetaVariant {
  public:
    HTPMetaVariant() {};

    HTPMetaVariant(const variant_id& id,
                   const std::string& chrom,
                   const int& pos,
                   const int& nstudies);
    // these get set when HTPMetaVariant is added to the class
    variant_id id;
    std::string chrom;
    int pos;
    int nstudies;
    Eigen::VectorXd scores;
    Eigen::VectorXd betas;
    Eigen::VectorXd ses;
    Eigen::VectorXd aafs; // alt allele frequencies
    Eigen::VectorXd cases_het;
    Eigen::VectorXd cases_alt;
    Eigen::VectorXd controls_het;
    Eigen::VectorXd controls_alt;
    Eigen::VectorXd sample_sizes;
    Eigen::VectorXd scorev;
    Eigen::VectorXi macs;
    Eigen::VectorXi ncases;
    Eigen::VectorXi ncontrols;
    std::vector<int> cohort_idx;
    int max_mac;
    double aac;
    double an;
    double max_aaf;

    // these get set when meta_analyze_gene is called
    double aaf;
    double burden_weight;
    double skato_weight;
    double acatv_weight;

    // this gets update when meta_analyze_gene is called
    bool is_singleton;
};

struct burden_params_t {
  std::vector<double> af_bins;
  singleton_def_e singleton_def;
  weight_strategy_e weight_strategy;
  double mask_spa_z_score;
  double mask_spa_case_control_ratio;
  double sv_spa_z_score;
  double sv_spa_case_control_ratio;

  bool params_set() { return af_bins.size() > 0; }
};

struct skato_params_t {
  std::vector<double> af_bins;
  std::vector<double> rho_values;
  weight_strategy_e weight_strategy;
  int min_aac;
  double mask_spa_z_score;
  double mask_spa_case_control_ratio;
  double sv_spa_z_score;
  double sv_spa_case_control_ratio;

  bool params_set() { return af_bins.size() > 0; }
};

struct acatv_params_t {
  std::vector<double> af_bins;
  weight_strategy_e weight_strategy;
  int min_aac;
  double sv_spa_z_score;
  double sv_spa_case_control_ratio;

  bool params_set() { return af_bins.size() > 0; }
};

class HTPMetaAnalyzer : public SetMetaAnalyzer {
 public:
  HTPMetaAnalyzer(const std::string& trait_name,
                  const trait_type_e& trait_type,
                  const std::vector<std::string>& cohorts,
                  const af_strategy_e& af_strategy,
                  const std::string& mask_def_file,
                  const std::string& annotation_file,
                  const std::vector<std::string>& ld_file_prefixes);

  void set_af_file(const string& af_file);

  void set_conditional_variants(const vector<vector<htpv4_record_t> >& study_conditional_variants, int max_cond_var_per_gene);

  void set_run_burden(const std::vector<double>& af_bins,
                      const singleton_def_e& singleton_def,
                      const weight_strategy_e& weight_strategy,
                      double mask_spa_pval,
                      double mask_spa_case_control_ratio,
                      double sv_spa_pval,
                      double sv_spa_case_control_ratio);

  void set_run_skato(const std::vector<double>& af_bins,
                     const std::vector<double>& rho_values,
                     const weight_strategy_e& weight_strategy,
                     int min_aac,
                     double mask_spa_pval,
                     double mask_spa_case_control_ratio,
                     double sv_spa_pval,
                     double sv_spa_case_control_ratio);

  void set_run_acatv(const std::vector<double>& af_bins,
                     const weight_strategy_e& weight_strategy,
                     int min_aac,
                     double sv_spa_pval,
                     double sv_spa_case_control_ratio);

  void set_write_cohort_burden_tests(HTPv4Writer& burden_writer);

  void set_write_mask_snplist(BgzWriter& snplist_writer);

  void set_ignore_mask_ld();

  void set_recompute_score();

  void set_keep_variants_not_in_ld_mat();

  void set_write_freqs(AlleleFreqWriter& freq_writer);

  void add_line(const htpv4_record_t& rec, const int& study_index);

  vector<htpv4_record_t> meta_analyze_gene(Gene g);

  void clear_before(const string& chrom, const int& before_pos);

 private:
  std::string trait_name;
  trait_type_e trait_type;
  double max_aaf;
  std::vector<std::string> cohorts;
  af_strategy_e af_strategy;
  weight_strategy_e weight_strategy;
  singleton_def_e singleton_def;
  std::string mask_def_file;
  std::string annotation_file;
  std::vector<std::string> ld_file_prefixes;
  double mask_spa_z_score;
  double mask_spa_case_control_ratio;
  double sv_spa_z_score;
  double sv_spa_case_control_ratio;
  std::string af_file;
  AlleleFrequencyMap af_map;
  int nstudies;
  vector<double> af_bins; // all bins across burden, SKATO, and ACATV
  double mask_max_aaf;

  unordered_map<variant_id, HTPMetaVariant> variants;
  // mask_set must appear before anno_map so that it is initialized first
  MaskSet mask_set;
  AnnotationMap anno_map;
  vector<RemetaMatrixReader> ld_mat_readers;

  // set by set_run_$TEST
  burden_params_t burden_params;
  skato_params_t skato_params;
  acatv_params_t acatv_params;

  // set by set_conditional_variants
  vector<variant_id> conditional_variants;
  unordered_set<variant_id> conditional_variants_set;
  Eigen::MatrixXd conditional_scores;
  Eigen::MatrixXd conditional_scorev;
  Eigen::MatrixXd conditional_betas;
  Eigen::MatrixXd conditional_ses;
  Eigen::VectorXd conditional_scores_sum;
  Eigen::VectorXd conditional_scorev_sum;
  Eigen::VectorXd conditional_chisq_stats;
  vector<int> conditional_sorted_variants_idx;
  int max_cond_var_per_gene;

  // set by set_write_cohort_burden_tests
  bool write_cohort_burden_tests;
  HTPv4Writer* burden_writer;

  // set by set_write_mask_snplist
  bool write_mask_snplist;
  BgzWriter* snplist_writer;

  bool ignore_mask_ld;
  bool recompute_score;
  bool keep_variants_not_in_ld_mat;

  // set by set_write_freqs
  bool write_freqs;
  AlleleFreqWriter* freq_writer;

  void check_info(const htpv4_record_t& rec);
  double get_score(const htpv4_record_t& rec);
  double get_scorev(const htpv4_record_t& rec);
  bool ld_mat_required();
};

#endif
