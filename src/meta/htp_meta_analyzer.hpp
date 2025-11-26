#ifndef HTP_META_ANALYZER_H
#define HTP_META_ANALYZER_H

#include "../lapack_complex.hpp"

#include <string>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

#include "set_meta_analyzer.hpp"
#include "../io/allele_freq_writer.hpp"
#include "../io/anno_reader.hpp"
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

    // this gets set when get_variants is called
    double aaf;

    // this gets set when get_variants is called
    bool is_singleton;
};

class GeneSumStats {
 public:
  GeneSumStats(int nvariants, int nmasks, int nstudies);

  int nvariants;
  int nmasks;
  int nstudies;

  // Indicator matrix. Entry i, j is 1 if variants_found[j]
  // is contained in mask_set.get_mask(i).
  Eigen::MatrixXd masks;
  // Matrix types hold sum stats per study
  // Vector types hold sum stats aggregated across studies
  Eigen::MatrixXi mask_ncases_per_study;
  Eigen::MatrixXi mask_ncontrols_per_study;
  Eigen::MatrixXd mask_betas;
  Eigen::MatrixXd mask_ses;
  Eigen::MatrixXd mask_aafs;
  Eigen::MatrixXd mask_cases_het;
  Eigen::MatrixXd mask_cases_alt;
  Eigen::MatrixXd mask_controls_het;
  Eigen::MatrixXd mask_controls_alt;
  Eigen::VectorXd mask_cf;
  Eigen::VectorXd mask_pval;
  Eigen::VectorXd mask_nhet;
  Eigen::VectorXd mask_nhom_alt;
  Eigen::VectorXd mask_nhom_ref;
  Eigen::VectorXd mask_case_control_ratio;
  Eigen::VectorXd mask_aaf;
  Eigen::MatrixXi mask_max_mac;
  std::vector<std::vector<int> > mask_indices;
  std::vector<std::unordered_set<int> > mask_cohort_idx;
  std::vector<std::vector<std::string> > mask_variants;
  vector<int> mask_ncases;
  vector<int> mask_ncontrols;
  std::vector<Eigen::SparseMatrix<double> > mask_selectors;
  std::vector<Eigen::SparseMatrix<double> > mask_collapsors;

  // Single variant stats
  Eigen::MatrixXd scores;
  Eigen::MatrixXd scorev;
  Eigen::MatrixXd betas;
  Eigen::MatrixXd ses;
  Eigen::MatrixXd aafs; // per cohort aafs
  Eigen::MatrixXd ncases;
  Eigen::MatrixXd cases_het;
  Eigen::MatrixXd cases_alt;
  Eigen::MatrixXd ncontrols;
  Eigen::MatrixXd controls_het;
  Eigen::MatrixXd controls_alt;
  Eigen::VectorXd meta_aafs; // aafs for masks
  Eigen::VectorXd meta_aacs;
  Eigen::VectorXi is_singleton;
  Eigen::VectorXd burden_weights;
  Eigen::VectorXd skato_weights;
  Eigen::VectorXd acatv_weights;
  Eigen::VectorXd scores_sum;
  Eigen::VectorXd scorev_sum;

  std::vector<double> study_sample_size;
  Eigen::VectorXd skato_var_to_collapse;
  Eigen::VectorXd collapse_anno_indicators;
  double max_ccr;
  std::vector<std::string> gene_conditional_variants;
  Eigen::MatrixXd gene_cond_scores;
  Eigen::MatrixXd gene_cond_scorev;
  std::vector<std::unordered_set<variant_id> > ld_variants_found;

  Eigen::SparseMatrix<double> gene_ld_mat;
  Eigen::MatrixXd gene_buffer_ld_mat;
  Eigen::MatrixXd buffer_ld_mat;

  bool found_at_least_one_conditional_variant;
  bool found_at_least_one_ld_mat;
  string conditional_variant_info;

  Eigen::VectorXd variant_cf_sqrt;
  Eigen::VectorXd variant_z_scores;
  Eigen::VectorXd burden_scores;
  Eigen::VectorXd skato_scores;

  std::vector<Eigen::SparseMatrix<double> > cohort_ld_mats;
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
  int collapse_below_aac;
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
                     int collapse_below_aac,
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

  void set_collapse_anno(const string& anno_name);

  void add_line(const htpv4_record_t& rec, const int& study_index);

  vector<htpv4_record_t> meta_analyze_gene(Gene g);

  void clear_before(const string& chrom, const int& before_pos);

  // populate variants_found and variant_annotations based on the htp entries collected with add_line
  void get_variants(vector<variant_id>& variants_found,
                    vector<int>& variant_annotations,
                    Gene g);

  GeneSumStats collect_sum_stats(const vector<variant_id>& variants_found,
                                 const vector<int>& variant_annotations,
                                 const MaskSet& mask_set,
                                 const Gene& g);

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

  // set by set_collapse_anno
  int collapse_anno_int;
  bool collapse_anno;

  void check_info(const htpv4_record_t& rec);
  double get_score(const htpv4_record_t& rec);
  double get_scorev(const htpv4_record_t& rec);
  bool ld_mat_required();

  void write_snplist(const GeneSumStats& sum_stats, const Gene& g);

  void compute_conditional_variant_stats(GeneSumStats& sum_stats);

  void compute_variant_spa(GeneSumStats& sum_stats);

  void compute_burden(vector<htpv4_record_t>& burden_results,
                      GeneSumStats& sum_stats,
                      const Eigen::SparseMatrix<double>& ld_mat_meta,
                      const htpv4_record_t& result_template,
                      int mask_idx);

  void compute_skato(vector<htpv4_record_t>& skato_results,
                     const GeneSumStats& sum_stats,
                     const Eigen::SparseMatrix<double>& ld_mat_meta,
                     const htpv4_record_t& result_template,
                     int mask_idx);

  void compute_acatv(vector<htpv4_record_t>& acatv_results,
                     const GeneSumStats& sum_stats,
                     const htpv4_record_t& result_template,
                     int mask_idx);
};

#endif
