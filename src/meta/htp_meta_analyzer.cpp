// Run some or all gene-based tests in a single pass.
#include "htp_meta_analyzer.hpp"

#include <algorithm>
#include <chrono>
#include <exception>
#include <set>
#include <string>
#include <unordered_map>
using std::set;
using std::string;
using std::unordered_map;

#include <boost/format.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <Eigen/Dense>
using Eigen::DiagonalMatrix;
using Eigen::MatrixXd;
using Eigen::MatrixXf;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix<float> SpMatF;
typedef Eigen::Triplet<double> Triplet;

bool vector_contains(const vector<double>& vec, const double& val) {
  return std::find(vec.begin(), vec.end(), val) != vec.end();
}

bool unordered_set_contains(const std::unordered_set<std::string>& vec, const std::string& val) {
  return vec.find(val) != vec.end();
}

HTPMetaVariant::HTPMetaVariant(const variant_id& id,
                               const string& chrom,
                               const int& pos,
                               const int& nstudies)
 : id(id)
 , chrom(chrom)
 , pos(pos)
 , nstudies(nstudies)  {
  this->scores         = VectorXd::Zero(nstudies);
  this->betas          = VectorXd::Zero(nstudies);
  this->ses            = VectorXd::Zero(nstudies);
  this->aafs           = VectorXd::Zero(nstudies);
  this->cases_het      = VectorXd::Zero(nstudies);
  this->cases_alt      = VectorXd::Zero(nstudies);
  this->controls_het   = VectorXd::Zero(nstudies);
  this->controls_alt   = VectorXd::Zero(nstudies);
  this->sample_sizes   = VectorXd::Zero(nstudies);
  this->scorev         = VectorXd::Zero(nstudies);
  this->macs           = VectorXi::Zero(nstudies);
  this->ncases         = VectorXi::Zero(nstudies);
  this->ncontrols      = VectorXi::Zero(nstudies);
  this->aac            = 0;
  this->an             = 0;
  this->max_aaf        = 0;
  this->max_mac        = 0;
  this->aaf            = 0;
  this->is_singleton   = false;
}

GeneSumStats::GeneSumStats(int nvariants, int nmasks, int nstudies)
: nvariants(nvariants)
, nmasks(nmasks)
, nstudies(nstudies)
, max_ccr(0)
, found_at_least_one_conditional_variant(false)
, found_at_least_one_ld_mat(false) {
  this->masks                    = MatrixXd(nmasks, nvariants);
  this->mask_ncases_per_study    = MatrixXi::Zero(nmasks, nstudies);
  this->mask_ncontrols_per_study = MatrixXi::Zero(nmasks, nstudies);
  this->mask_betas               = MatrixXd::Zero(nmasks, nstudies);
  this->mask_ses                 = MatrixXd::Zero(nmasks, nstudies);
  this->mask_aafs                = MatrixXd::Zero(nmasks, nstudies);
  this->mask_cases_het           = MatrixXd::Zero(nmasks, nstudies);
  this->mask_cases_alt           = MatrixXd::Zero(nmasks, nstudies);
  this->mask_controls_het        = MatrixXd::Zero(nmasks, nstudies);
  this->mask_controls_alt        = MatrixXd::Zero(nmasks, nstudies);
  this->mask_cf                  = VectorXd::Ones(nmasks);
  this->mask_pval                = VectorXd::Constant(nmasks, -1);  // -1 indicates no result
  this->mask_nhet                = VectorXd::Zero(nmasks);
  this->mask_nhom_alt            = VectorXd::Zero(nmasks);
  this->mask_nhom_ref            = VectorXd::Zero(nmasks);
  this->mask_case_control_ratio  = VectorXd::Zero(nmasks);
  this->mask_aaf                 = VectorXd::Zero(nmasks);
  this->mask_max_mac             = VectorXi::Zero(nmasks);
  this->mask_indices             = vector<vector<int> >(nmasks);
  this->mask_cohort_idx          = vector<unordered_set<int> >(nmasks);
  this->mask_variants            = vector<vector<string> >(nmasks);
  this->mask_ncases              = vector<int>(nmasks, 0);
  this->mask_ncontrols           = vector<int>(nmasks, 0);
  this->mask_selectors           = vector<SpMat>(nmasks);
  this->mask_collapsors          = vector<SpMat>(nmasks);
  this->scores                   = MatrixXd::Zero(nvariants, nstudies);
  this->scorev                   = MatrixXd::Zero(nvariants, nstudies);
  this->betas                    = MatrixXd::Zero(nvariants, nstudies);
  this->ses                      = MatrixXd::Zero(nvariants, nstudies);
  this->aafs                     = MatrixXd::Zero(nvariants, nstudies);
  this->ncases                   = MatrixXd::Zero(nvariants, nstudies);
  this->cases_het                = MatrixXd::Zero(nvariants, nstudies);
  this->cases_alt                = MatrixXd::Zero(nvariants, nstudies);
  this->ncontrols                = MatrixXd::Zero(nvariants, nstudies);
  this->controls_het             = MatrixXd::Zero(nvariants, nstudies);
  this->controls_alt             = MatrixXd::Zero(nvariants, nstudies);
  this->meta_aafs                = VectorXd::Zero(nvariants);
  this->meta_aacs                = VectorXd::Zero(nvariants);
  this->is_singleton             = VectorXi::Zero(nvariants);
  this->burden_weights           = VectorXd::Zero(nvariants);
  this->skato_weights            = VectorXd::Zero(nvariants);
  this->acatv_weights            = VectorXd::Zero(nvariants);
  this->study_sample_size        = vector<double>(nstudies, 0);
  this->skato_var_to_collapse    = VectorXd::Zero(nvariants);
  this->collapse_anno_indicators = VectorXd::Zero(nvariants);
  this->max_ccr                  = 0;
  this->ld_variants_found        = vector<unordered_set<variant_id> >(nstudies);
}

HTPMetaAnalyzer::HTPMetaAnalyzer(const string& trait_name,
                                 const trait_type_e& trait_type,
                                 const vector<string>& cohorts,
                                 const af_strategy_e& af_strategy,
                                 const string& mask_def_file,
                                 const string& annotation_file,
                                 const vector<string>& ld_file_prefixes)
 : trait_name(trait_name)
 , trait_type(trait_type)
 , cohorts(cohorts)
 , af_strategy(af_strategy)
 , mask_def_file(mask_def_file)
 , annotation_file(annotation_file)
 , ld_file_prefixes(ld_file_prefixes)
 , mask_spa_z_score(SPA_NEVER_APPLIED)
 , mask_spa_case_control_ratio(0)
 , sv_spa_z_score(SPA_NEVER_APPLIED)
 , sv_spa_case_control_ratio(0)
 , af_file("")
 , nstudies(cohorts.size())
 , mask_max_aaf(0)
 , variants(unordered_map<variant_id, HTPMetaVariant>())
 , mask_set(mask_def_file, vector<double>{1})
 , anno_map(annotation_file, mask_set)
 , ld_mat_readers(vector<RemetaMatrixReader>())
 , conditional_variants(vector<variant_id>())
 , conditional_scores(Eigen::MatrixXd(0, 0))
 , conditional_scorev(Eigen::MatrixXd(0, 0))
 , conditional_betas(Eigen::MatrixXd(0, 0))
 , conditional_ses(Eigen::MatrixXd(0, 0))
 , write_cohort_burden_tests(false)
 , burden_writer(nullptr)
 , write_mask_snplist(false)
 , snplist_writer(nullptr)
 , ignore_mask_ld(false)
 , recompute_score(false)
 , keep_variants_not_in_ld_mat(false)
 , write_freqs(false)
 , freq_writer(nullptr)
 , collapse_anno_int(-1)
 , collapse_anno(false) {
  for (const string& prefix : ld_file_prefixes) {
    ld_mat_readers.push_back(
      RemetaMatrixReader()
    );
    ld_mat_readers[ld_mat_readers.size()-1].open(prefix);
  }
}

void HTPMetaAnalyzer::set_af_file(const string& af_file) {
  log_info("using external allele frequencies");
  this->af_strategy = USE_EXTERNAL_AF;
  this->af_file = af_file;
  this->af_map.load(af_file);
}

void HTPMetaAnalyzer::set_conditional_variants(const vector<vector<htpv4_record_t> >& study_conditional_variants, int max_cond_var_per_gene) {
  this->conditional_variants.clear();
  unordered_map<variant_id, int> variant_idx;
  for (size_t i = 0; i < study_conditional_variants.size(); ++i) {
    for ( const htpv4_record_t& rec : study_conditional_variants[i] ) {
      if (variant_idx.count(rec.name) == 0) {
        variant_idx[rec.name] = variant_idx.size();
        this->conditional_variants.push_back(rec.name);
        this->conditional_variants_set.insert(rec.name);
      }
    }
  }

  this->conditional_scores = MatrixXd::Zero(variant_idx.size(), study_conditional_variants.size());
  this->conditional_scorev = MatrixXd::Zero(variant_idx.size(), study_conditional_variants.size());
  this->conditional_betas = MatrixXd::Zero(variant_idx.size(), study_conditional_variants.size());
  this->conditional_ses = MatrixXd::Zero(variant_idx.size(), study_conditional_variants.size());
  for (size_t i = 0; i < study_conditional_variants.size(); ++i) {
    for ( const htpv4_record_t& rec : study_conditional_variants[i] ) {
      this->check_info(rec);
      this->conditional_scores(variant_idx[rec.name], i) = this->get_score(rec);
      this->conditional_scorev(variant_idx[rec.name], i) = this->get_scorev(rec);
      this->conditional_betas(variant_idx[rec.name], i) = HTPv4Reader::get_beta(rec);
      this->conditional_ses(variant_idx[rec.name], i) = HTPv4Reader::get_se(rec);
    }
  }
  this->conditional_scores_sum = this->conditional_scores.rowwise().sum();
  this->conditional_scorev_sum = this->conditional_scorev.rowwise().sum();
  this->conditional_chisq_stats = this->conditional_scores_sum.array().square() / this->conditional_scorev_sum.array();
  this->conditional_sorted_variants_idx.resize(this->conditional_variants.size());
  iota(this->conditional_sorted_variants_idx.begin(), this->conditional_sorted_variants_idx.end(), 0);
  std::stable_sort(this->conditional_sorted_variants_idx.begin(), this->conditional_sorted_variants_idx.end(),
       [&](const int& a, const int& b) {
         return this->conditional_chisq_stats(a) > this->conditional_chisq_stats(b);
       });
  this->max_cond_var_per_gene = max_cond_var_per_gene;
}

void HTPMetaAnalyzer::set_run_burden(const vector<double>& af_bins,
                                     const singleton_def_e& singleton_def,
                                     const weight_strategy_e& weight_strategy,
                                     double mask_spa_pval,
                                     double mask_spa_case_control_ratio,
                                     double sv_spa_pval,
                                     double sv_spa_case_control_ratio) {
  this->burden_params.af_bins.clear();
  if (singleton_def != OMIT) {
    this->burden_params.af_bins.push_back(SINGLETON);

    if (std::find(this->af_bins.begin(), this->af_bins.end(), SINGLETON) == this->af_bins.end()) {
      this->af_bins.push_back(SINGLETON);
    }
  }
  for ( const double& af : af_bins ) {
    this->burden_params.af_bins.push_back(af);

    if (std::find(this->af_bins.begin(), this->af_bins.end(), af) == this->af_bins.end()) {
      this->af_bins.push_back(af);
    }
    this->mask_max_aaf = max(this->mask_max_aaf, af);
  }
  this->burden_params.singleton_def = singleton_def;
  this->burden_params.weight_strategy = weight_strategy;
  std::sort(this->af_bins.begin(), this->af_bins.end());
  this->mask_set = MaskSet(this->mask_def_file, this->af_bins);

  this->burden_params.mask_spa_case_control_ratio = mask_spa_case_control_ratio;
  this->mask_spa_case_control_ratio = max(
    this->mask_spa_case_control_ratio,
    this->burden_params.mask_spa_case_control_ratio
  );
  if (mask_spa_pval <= 0) {
    this->burden_params.mask_spa_z_score = SPA_NEVER_APPLIED;
  } else if (mask_spa_pval >= 1) {
    this->burden_params.mask_spa_z_score = SPA_ALWAYS_APPLIED;
  } else {
    boost::math::normal s;
    this->burden_params.mask_spa_z_score = quantile(s, mask_spa_pval / 2); // two-sided pvalue
  }
  this->mask_spa_z_score = max(
    this->mask_spa_z_score,
    this->burden_params.sv_spa_z_score
  );

  this->burden_params.sv_spa_case_control_ratio = sv_spa_case_control_ratio;
  if (sv_spa_pval <= 0) {
    this->burden_params.sv_spa_z_score = SPA_NEVER_APPLIED;
  } else if (sv_spa_pval >= 1) {
    this->burden_params.sv_spa_z_score = SPA_ALWAYS_APPLIED;
  } else {
    boost::math::normal s;
    this->burden_params.sv_spa_z_score = quantile(s, sv_spa_pval / 2);
  }
  this->sv_spa_z_score = max(
    this->sv_spa_z_score,
    this->burden_params.sv_spa_z_score
  );
}

void HTPMetaAnalyzer::set_run_skato(const vector<double>& af_bins,
                                    const vector<double>& rho_values,
                                    const weight_strategy_e& weight_strategy,
                                    int min_aac,
                                    int collapse_below_aac,
                                    double mask_spa_pval,
                                    double mask_spa_case_control_ratio,
                                    double sv_spa_pval,
                                    double sv_spa_case_control_ratio) {
  this->skato_params.af_bins.clear();
  for ( const double& af : af_bins ) {
    this->skato_params.af_bins.push_back(af);

    if (std::find(this->af_bins.begin(), this->af_bins.end(), af) == this->af_bins.end()) {
      this->af_bins.push_back(af);
    }
    this->mask_max_aaf = max(this->mask_max_aaf, af);
  }
  this->skato_params.rho_values = rho_values;
  this->skato_params.weight_strategy = weight_strategy;
  this->skato_params.min_aac = min_aac;
  this->skato_params.collapse_below_aac = collapse_below_aac;
  std::sort(this->af_bins.begin(), this->af_bins.end());
  this->mask_set = MaskSet(this->mask_def_file, this->af_bins);

  this->skato_params.mask_spa_case_control_ratio = mask_spa_case_control_ratio;
  this->mask_spa_case_control_ratio = max(
    this->mask_spa_case_control_ratio,
    this->skato_params.mask_spa_case_control_ratio
  );
  if (mask_spa_pval <= 0) {
    this->skato_params.mask_spa_z_score = SPA_NEVER_APPLIED;
  } else if (mask_spa_pval >= 1) {
    this->skato_params.mask_spa_z_score = SPA_ALWAYS_APPLIED;
  } else {
    boost::math::normal s;
    this->skato_params.mask_spa_z_score = quantile(s, mask_spa_pval / 2);
  }
  this->mask_spa_z_score = max(
    this->mask_spa_z_score,
    this->skato_params.mask_spa_z_score
  );

  this->skato_params.sv_spa_case_control_ratio = sv_spa_case_control_ratio;
  this->sv_spa_case_control_ratio = max(
    this->sv_spa_case_control_ratio,
    this->skato_params.sv_spa_case_control_ratio
  );
  if (sv_spa_pval <= 0) {
    this->skato_params.sv_spa_z_score = SPA_NEVER_APPLIED;
  } else if (sv_spa_pval >= 1) {
    this->skato_params.sv_spa_z_score = SPA_ALWAYS_APPLIED;
  } else {
    boost::math::normal s;
    this->skato_params.sv_spa_z_score = quantile(s, sv_spa_pval / 2);
  }
  this->sv_spa_z_score = max(
    this->sv_spa_z_score,
    this->skato_params.sv_spa_z_score
  );
}

void HTPMetaAnalyzer::set_run_acatv(const vector<double>& af_bins,
                                    const weight_strategy_e& weight_strategy,
                                    int min_aac,
                                    double sv_spa_pval,
                                    double sv_spa_case_control_ratio) {
  this->acatv_params.af_bins.clear();
  for ( const double& af : af_bins ) {
    this->acatv_params.af_bins.push_back(af);

    if (std::find(this->af_bins.begin(), this->af_bins.end(), af) == this->af_bins.end()) {
      this->af_bins.push_back(af);
    }
    this->mask_max_aaf = max(this->mask_max_aaf, af);
  }
  this->acatv_params.weight_strategy = weight_strategy;
  this->acatv_params.min_aac = min_aac;
  std::sort(this->af_bins.begin(), this->af_bins.end());
  this->mask_set = MaskSet(this->mask_def_file, this->af_bins);

  this->acatv_params.sv_spa_case_control_ratio = sv_spa_case_control_ratio;
  this->sv_spa_case_control_ratio = max(
    this->sv_spa_case_control_ratio,
    this->acatv_params.sv_spa_case_control_ratio
  );
  if (sv_spa_pval <= 0) {
    this->acatv_params.sv_spa_z_score = SPA_NEVER_APPLIED;
  } else if (sv_spa_pval >= 1) {
    this->acatv_params.sv_spa_z_score = SPA_ALWAYS_APPLIED;
  } else {
    boost::math::normal s;
    this->acatv_params.sv_spa_z_score = quantile(s, sv_spa_pval / 2);
  }
  this->sv_spa_z_score = max(
    this->sv_spa_z_score,
    this->acatv_params.sv_spa_z_score
  );
}

void HTPMetaAnalyzer::add_line(const htpv4_record_t& rec, const int& study_index) {
  if (rec.pval == HTPv4_NA) { // typically happens when a variant fails Firth
    log_debug(rec.name + " does not have a p-value, skipping...");
    return ;
  }

  // occurs when trying to update older htpv4 files with the SKATV info field
  if (rec.info.count("SKATV") > 0 && rec.info.at("SKATV") == "inf") {
    log_debug(rec.name + " has bad skatv value, skipping...");
    return ;
  }
  double scorev = this->get_scorev(rec);
  if (scorev <= 0 || boost::math::isnan(scorev) || boost::math::isinf(scorev)) {
    log_debug(rec.name + " in study " + this->cohorts[study_index] + " has bad skatv value, skipping...");
    return ;
  }
  this->check_info(rec);

  double aac;
  double an;
  if (this->trait_type == QT) {
    aac = rec.cases_het + 2*rec.cases_alt;
    an = 2*rec.num_cases;
  } else {
    aac = rec.cases_het + rec.controls_het + 2*(rec.cases_alt + rec.controls_alt);
    an = 2*(rec.num_cases + rec.num_controls);
  }

  if (this->variants.count(rec.name) == 0) {
    this->variants[rec.name] = HTPMetaVariant(rec.name, rec.chr, rec.pos, this->nstudies);
    this->variants[rec.name].is_singleton = (aac == 1);
  }

  this->variants[rec.name].scores[study_index] = this->get_score(rec);
  this->variants[rec.name].scorev[study_index] = this->get_scorev(rec);
  this->variants[rec.name].aafs[study_index] = rec.aaf;
  this->variants[rec.name].cases_het[study_index] = rec.cases_het;
  this->variants[rec.name].cases_alt[study_index] = rec.cases_alt;
  this->variants[rec.name].controls_het[study_index] = rec.controls_het;
  this->variants[rec.name].controls_alt[study_index] = rec.controls_alt;
  this->variants[rec.name].sample_sizes[study_index] = this->trait_type == BT ? rec.num_cases + rec.num_controls : rec.num_cases;
  this->variants[rec.name].ncases[study_index] = rec.num_cases;
  if (this->trait_type == BT) {
    this->variants[rec.name].ncontrols[study_index] = rec.num_controls;
  }
  this->variants[rec.name].aac += aac;
  this->variants[rec.name].an += an;
  this->variants[rec.name].max_aaf = max(variants[rec.name].max_aaf, rec.aaf);

  this->variants[rec.name].cohort_idx.push_back(study_index);

  if (rec.aaf <= 0.5) {
    this->variants[rec.name].macs[study_index] = aac;
  } else {
    this->variants[rec.name].macs[study_index] = an - aac;
  }
  this->variants[rec.name].max_mac = max(this->variants[rec.name].max_mac, this->variants[rec.name].macs[study_index]);

  this->variants[rec.name].betas[study_index] = HTPv4Reader::get_beta(rec);
  this->variants[rec.name].ses[study_index] = HTPv4Reader::get_se(rec);

  if (this->burden_params.singleton_def == WITHIN) {
    this->variants[rec.name].is_singleton =
      this->variants[rec.name].is_singleton && (aac == 1);
  } else if (this->burden_params.singleton_def == ACROSS) {
    this->variants[rec.name].is_singleton = (this->variants[rec.name].aac == 1);
  } else {
    this->variants[rec.name].is_singleton = false;
  }
}

void HTPMetaAnalyzer::clear_before(const string& chrom, const int& before_pos) {
  if (this->write_freqs) {
    log_debug("writing allele frequencies...");
    this->freq_writer->write_freqs_before(chrom, before_pos);
    this->freq_writer->clear_freqs_before(chrom, before_pos);
  }
  log_debug("clearing records chr " + chrom + " before position " + to_string(before_pos));
  auto it = this->variants.begin();
  while (it != this->variants.end()) {
    if (it->second.chrom == chrom && it->second.pos <= before_pos) {
      it = this->variants.erase(it);
    } else {
      ++it;
    }
  }
}

void HTPMetaAnalyzer::get_variants(vector<variant_id>& variants_found,
                                   vector<int>& variant_annotations,
                                   Gene g) {
  int variants_without_af = 0;
  for (const string& vid : *g.get_variants()) {
    if (this->variants.count(vid) > 0) {
      HTPMetaVariant *variant = &this->variants[vid];

      if (this->af_strategy == USE_MAX_AF) {
        variant->aaf = variant->max_aaf;
      } else if (this->af_strategy == USE_OVERALL_AF) {
        variant->aaf = variant->aac / variant->an;
      } else if (this->af_strategy == USE_EXTERNAL_AF) {
        af_info_t af_info = this->af_map.get_af_info(vid);
        if (af_info == AF_INFO_NOT_FOUND) {
          ++variants_without_af;
          continue;
        }
        variant->aaf = af_info.freq;
        variant->is_singleton = af_info.is_singleton == UNSPECIFIED ? variant->is_singleton : af_info.is_singleton == YES;
      }

      if (this->write_freqs) {
        this->freq_writer->add_freq(vid, variant->chrom, variant->pos, variant->aaf, variant->is_singleton);
      }

      int anno = anno_map.get_annotation(vid, g.get_name(), variant->chrom, variant->pos);
      if (anno == -1) {
        continue;
      }

      if (variant->aaf <= this->mask_max_aaf) {
        variants_found.push_back(vid);
        variant_annotations.push_back(anno);
      }
    }
  }
  if (this->af_strategy == USE_EXTERNAL_AF && variants_without_af > 0) {
    log_warning(
      "found " + to_string(variants_without_af) +
      " variants without allele frequencies, "
      "these will be omitted from masks"
    );
  }
}

GeneSumStats HTPMetaAnalyzer::collect_sum_stats(const vector<variant_id>& variants_found,
                                                const vector<int>& variant_annotations,
                                                const MaskSet& mask_set,
                                                const Gene& g) {
  std::chrono::time_point<std::chrono::steady_clock> start;
  std::chrono::duration<double> d;
  string gene_name = g.get_name();

  int nmasks = mask_set.n_masks();
  int nstudies = this->cohorts.size();
  GeneSumStats sum_stats(variants_found.size(), nmasks, nstudies);

  boost::math::beta_distribution<> beta(1, 25);
  double weight_aaf = 0;
  double weight_maf = 0;
  log_debug(" * scanning through variants...");
  // Indicator matrix. Entry i, j is 1 if variants_found[j]
  // is contained in mask_set.get_mask(i).
  bool collapse_anno_found = false;
  for (size_t j = 0; j < variants_found.size(); ++j) {
    HTPMetaVariant variant = this->variants[variants_found[j]];
    sum_stats.scores.row(j) = variant.scores;
    sum_stats.scorev.row(j) = variant.scorev;
    sum_stats.betas.row(j) = variant.betas;
    sum_stats.ses.row(j) = variant.ses;
    sum_stats.aafs.row(j) = variant.aafs;
    sum_stats.cases_het.row(j) = variant.cases_het;
    sum_stats.cases_alt.row(j) = variant.cases_alt;
    sum_stats.ncases.row(j) = variant.ncases.cast<double>();
    sum_stats.controls_het.row(j) = variant.controls_het;
    sum_stats.controls_alt.row(j) = variant.controls_alt;
    sum_stats.ncontrols.row(j) = variant.ncontrols.cast<double>();
    sum_stats.meta_aafs(j) = variant.aaf;
    sum_stats.meta_aacs(j) = variant.aac;
    sum_stats.is_singleton(j) = variant.is_singleton ? 1 : 0;

    weight_aaf = variant.aac / variant.an;
    weight_maf = weight_aaf > 0.5 ? 1 - weight_aaf : weight_aaf;
    if (this->burden_params.weight_strategy == USE_BETA_WEIGHTS) {
      sum_stats.burden_weights[j] = pdf(beta, weight_maf);
    } else {
      sum_stats.burden_weights[j] = 1;
    }
    if (this->skato_params.weight_strategy == USE_BETA_WEIGHTS && variant.aac >= this->skato_params.min_aac) {
      sum_stats.skato_weights[j] = pdf(beta, weight_maf);
    } else if (variant.aac >= this->skato_params.min_aac) {
      sum_stats.skato_weights[j] = 1;
    } else {
      sum_stats.skato_weights[j] = 0;
    }
    if (this->acatv_params.weight_strategy == USE_BETA_WEIGHTS && variant.aac >= this->acatv_params.min_aac) {
      sum_stats.acatv_weights[j] = pow(pdf(beta, weight_maf), 2) * weight_maf * (1 - weight_maf);
    } else if (variant.aac >= this->acatv_params.min_aac) {
      sum_stats.acatv_weights[j] = 1;
    } else {
      sum_stats.acatv_weights[j] = 0;
    }

    if (this->conditional_variants_set.count(variant.id) > 0) {
      sum_stats.scores.row(j).setZero();
      sum_stats.scorev.row(j).setZero();
    }

    if (variant.aac < this->skato_params.collapse_below_aac) {
      sum_stats.skato_var_to_collapse[j] = 1;
    } else {
      sum_stats.skato_var_to_collapse[j] = 0;
    }

    for (int i = 0; i < nmasks; ++i) {
      if (mask_set.mask_contains(i, variant_annotations[j], variant.aaf, variant.is_singleton)) {
        sum_stats.masks(i, j) = 1;
        sum_stats.mask_max_mac(i) = max(sum_stats.mask_max_mac(i), variant.max_mac);
        sum_stats.mask_indices[i].push_back(j);
        sum_stats.mask_cohort_idx[i].insert(variant.cohort_idx.begin(), variant.cohort_idx.end());

        for (int k = 0; k < nstudies; ++k) {
          sum_stats.mask_ncases_per_study(i, k) = max(sum_stats.mask_ncases_per_study(i, k), variant.ncases(k));
          sum_stats.mask_ncontrols_per_study(i, k) = max(sum_stats.mask_ncontrols_per_study(i, k), variant.ncontrols(k));
        }

        if (this->write_mask_snplist) {
          sum_stats.mask_variants[i].push_back(variants_found[j]);
        }
      } else {
        sum_stats.masks(i, j) = 0;
      }
    }

    for (int k = 0; k < nstudies; ++k) {
      sum_stats.study_sample_size[k] = max(sum_stats.study_sample_size[k], variant.sample_sizes[k]);
    }

    if (this->collapse_anno && variant_annotations[j] == this->collapse_anno_int) {
      sum_stats.collapse_anno_indicators[j] = 1;
      collapse_anno_found = true;
    } else {
      sum_stats.collapse_anno_indicators[j] = 0;
    }
  }

  log_debug(" * construcing selectors and collapsors...");
  for (int m = 0; m < nmasks; ++m) {
    sum_stats.mask_ncases[m] = sum_stats.mask_ncases_per_study.row(m).sum();
    sum_stats.mask_ncontrols[m] = sum_stats.mask_ncontrols_per_study.row(m).sum();
    sum_stats.max_ccr = max(sum_stats.max_ccr, (double)sum_stats.mask_ncases[m] / (double)sum_stats.mask_ncontrols[m]);

    vector<Triplet> selector_triplet_list(sum_stats.mask_indices[m].size());

    int nvar_to_collapse = sum_stats.masks.row(m) * sum_stats.skato_var_to_collapse;
    int nvar_to_keep = sum_stats.masks.row(m).sum() - nvar_to_collapse;
    int collapsor_idx = 0;
    vector<Triplet> collapsor_triplet_list(sum_stats.mask_indices[m].size());

    for (size_t k = 0; k < sum_stats.mask_indices[m].size(); ++k) {
      selector_triplet_list.push_back(Triplet(k, sum_stats.mask_indices[m][k], 1));

      if (sum_stats.skato_var_to_collapse(sum_stats.mask_indices[m][k]) > 0 || sum_stats.skato_weights(sum_stats.mask_indices[m][k]) == 0) {
        collapsor_triplet_list[k] = Triplet(nvar_to_keep, k, 1);
      } else {
        collapsor_triplet_list[k] = Triplet(collapsor_idx, k, 1);
        ++collapsor_idx;
      }
    }

    SpMat selector(sum_stats.mask_indices[m].size(), variants_found.size());
    selector.setFromTriplets(selector_triplet_list.begin(), selector_triplet_list.end());
    sum_stats.mask_selectors[m] = selector;

    SpMat collapsor(nvar_to_keep + 1, sum_stats.mask_indices[m].size());
    collapsor.setFromTriplets(collapsor_triplet_list.begin(), collapsor_triplet_list.end());
    sum_stats.mask_collapsors[m] = collapsor;
  }

  if (this->conditional_variants.size() > 0) {
    log_debug(" * checking conditional variants");
    unordered_set<string> available_conditional_variants;
    vector<int> gene_conditional_variants_idx;
    vector<string> variant_ids;
    for (int s = 0; s < nstudies; ++s) {
      if (ld_mat_readers[s].contains_gene(gene_name)) {
        ld_mat_readers[s].load_gene_variant_ids(variant_ids, gene_name);
        available_conditional_variants.insert(variant_ids.begin(), variant_ids.end());

        ld_mat_readers[s].load_buffer_variant_ids(variant_ids, gene_name);
        available_conditional_variants.insert(variant_ids.begin(), variant_ids.end());
      }
    }

    for (const int& idx: this->conditional_sorted_variants_idx) {
      if ((int)sum_stats.gene_conditional_variants.size() < this->max_cond_var_per_gene &&
          available_conditional_variants.count(this->conditional_variants[idx]) > 0) {
        sum_stats.gene_conditional_variants.push_back(this->conditional_variants[idx]);
        gene_conditional_variants_idx.push_back(idx);
      }
    }
    sum_stats.gene_cond_scores = this->conditional_scores(gene_conditional_variants_idx, Eigen::all);
    sum_stats.gene_cond_scorev = this->conditional_scorev(gene_conditional_variants_idx, Eigen::all);
  }

  start = std::chrono::steady_clock::now();
  sum_stats.found_at_least_one_ld_mat = false;
  sum_stats.found_at_least_one_conditional_variant = false;
  sum_stats.gene_ld_mat = SpMat(variants_found.size(), variants_found.size());
  sum_stats.gene_buffer_ld_mat = MatrixXd(variants_found.size(), sum_stats.gene_conditional_variants.size());
  sum_stats.buffer_ld_mat = MatrixXd(sum_stats.gene_conditional_variants.size(), sum_stats.gene_conditional_variants.size());
  for (int s = 0; s < nstudies; ++s) {
    if (this->ld_mat_required() && !ld_mat_readers[s].contains_gene(gene_name)) {
      log_warning("LD matrix for " + gene_name + " is not found in cohort " + this->cohorts[s] + ", skipping...");
      sum_stats.scores.col(s).setZero();
      continue;
    }

    SpMat Gsp;
    // MatrixXf Gf;
    MatrixXf G_Cf;
    MatrixXf Cf;
    //MatrixXd G;       // gene_ld_mat
    MatrixXd G_C;     // gene_buffer_ld_mat
    MatrixXd C;       // inverse of buffer_ld_mat

    if (sum_stats.gene_conditional_variants.size() > 0) {
      log_debug("loading conditional LD matrices");
      ld_mat_readers[s].load_conditional_ld_mats_sp_double(
        Gsp,
        Cf,
        G_Cf,
        variants_found,
        sum_stats.gene_conditional_variants,
        gene_name
      );
      C = Cf.cast<double>();
      Cf.resize(0, 0);
      G_C = G_Cf.cast<double>();
      G_Cf.resize(0, 0);
      sum_stats.found_at_least_one_conditional_variant = sum_stats.found_at_least_one_conditional_variant || (C.array() > 0).any();
    } else if (this->ld_mat_required()) {
      log_debug("loading marginal LD matrices");
      ld_mat_readers[s].load_gene_ld_mat_sp_double(Gsp, variants_found, gene_name);
    } else {
      log_debug("ignoring LD matrices");
      Gsp = SpMat(variants_found.size(), variants_found.size());
    }

    VectorXd G_diag = VectorXd(Gsp.diagonal());
    if (this->ignore_mask_ld) {
      log_debug("assuming diagonal LD matrices");
      Gsp = SpMat(G_diag.asDiagonal());
    }

    VectorXd G_diag_update = VectorXd::Zero(G_diag.size());
    for (size_t i = 0; i < variants_found.size(); ++i) {
      if (G_diag(i) > 0 && sum_stats.scores(i, s) != 0) {
        sum_stats.ld_variants_found[s].insert(variants_found[i]);
        sum_stats.found_at_least_one_ld_mat = true;
      } else if (G_diag(i) == 0 && sum_stats.scores(i, s) != 0 && !this->keep_variants_not_in_ld_mat) {
        log_warning(variants_found[i] + " does not have an entry in LD matrix for cohort " + this->cohorts[s]);
        sum_stats.scores(i, s) = 0;
      } else if (G_diag(i) == 0 && sum_stats.scores(i, s) != 0 && this->keep_variants_not_in_ld_mat) {
        sum_stats.ld_variants_found[s].insert(variants_found[i]);
        G_diag_update(i) = sum_stats.scorev(i, s);
        sum_stats.found_at_least_one_ld_mat = true;
      }
    }
    if (this->keep_variants_not_in_ld_mat) {
      G_diag = G_diag + G_diag_update;
      Gsp = Gsp + SpMat(G_diag_update.asDiagonal());
    }

    // The LD matrix exists but contains none of the requested variants
    if (Gsp.nonZeros() == 0) {
      log_warning("LD matrix for cohort " + this->cohorts[s] + " does not contain any of the requested variants");
      continue;
    }

    if (this->burden_params.params_set() || this->skato_params.params_set()) {
      log_debug("computing AAFs and genotype counts");
#if defined(_OPENMP)
      int nthreads = Eigen::nbThreads();
      Eigen::setNbThreads(1);
      #pragma omp parallel for schedule(dynamic)
#endif
      for (int m = 0; m < sum_stats.nmasks; ++m) {
        if (sum_stats.mask_indices[m].size() == 0 || (sum_stats.mask_ncases_per_study(m, s) + sum_stats.mask_ncontrols_per_study(m, s) == 0)) continue;
        SpMat Gsp_mask = sum_stats.mask_selectors[m] * Gsp * sum_stats.mask_selectors[m].transpose();

        sum_stats.mask_cases_het(m, s) = stat::misc::estimate_nhets(
          sum_stats.cases_het(sum_stats.mask_indices[m], s),
          sum_stats.aafs(sum_stats.mask_indices[m], s),
          Gsp_mask
        );
        sum_stats.mask_cases_alt(m, s) = stat::misc::estimate_nalts(
          sum_stats.cases_alt(sum_stats.mask_indices[m], s),
          sum_stats.aafs(sum_stats.mask_indices[m], s),
          Gsp_mask
        );
        if (this->trait_type == BT) {
          sum_stats.mask_controls_het(m, s) = stat::misc::estimate_nhets(
            sum_stats.controls_het(sum_stats.mask_indices[m], s),
            sum_stats.aafs(sum_stats.mask_indices[m], s),
            Gsp_mask
          );
          sum_stats.mask_controls_alt(m, s) = stat::misc::estimate_nalts(
            sum_stats.controls_alt(sum_stats.mask_indices[m], s),
            sum_stats.aafs(sum_stats.mask_indices[m], s),
            Gsp_mask
          );
        }
        if (this->trait_type == QT) {
          sum_stats.mask_aafs(m, s) = max(
            1.0 / (2.0*sum_stats.mask_ncases_per_study(m, s)),
            (sum_stats.mask_cases_het(m, s) + 2*sum_stats.mask_cases_alt(m, s)) / (2.0*sum_stats.mask_ncases_per_study(m, s))
          );
        } else {
          sum_stats.mask_aafs(m, s) = max(
            1.0 / (2.0*sum_stats.mask_ncases_per_study(m, s) + 2.0*sum_stats.mask_ncontrols_per_study(m, s)),
            (sum_stats.mask_cases_het(m, s) + sum_stats.mask_controls_het(m, s) + 2*sum_stats.mask_cases_alt(m, s) + 2*sum_stats.mask_controls_alt(m, s))
              / (2.0*sum_stats.mask_ncases_per_study(m, s) + 2.0*sum_stats.mask_ncontrols_per_study(m, s))
          );
        }
      }
#if defined(_OPENMP)
      Eigen::setNbThreads(nthreads);
#endif
    }

    if (!sum_stats.found_at_least_one_conditional_variant) {
      log_debug("rescaling matrices");
      stat::misc::rescale_cov(Gsp, G_diag, sum_stats.scorev.col(s));
    } else {
      log_debug("rescaling matrices");
      stat::misc::rescale_cov_block(
        Gsp,
        C,
        G_C,
        G_diag,
        sum_stats.scorev.col(s).array(),
        sum_stats.gene_cond_scorev.col(s).array()
      );
    }

    if (this->collapse_anno && collapse_anno_found) {
      stat::misc::collapse_anno(Gsp, sum_stats.scores.col(s), sum_stats.scorev.col(s), sum_stats.ses.col(s), sum_stats.collapse_anno_indicators);
      G_diag = Gsp.diagonal();
    }

    sum_stats.gene_ld_mat += Gsp;
    if (sum_stats.found_at_least_one_conditional_variant) {
      sum_stats.gene_buffer_ld_mat += G_C;
      sum_stats.buffer_ld_mat += C;
    }

    if (this->burden_params.params_set() || sum_stats.found_at_least_one_conditional_variant) {
      log_debug("computing effect sizes");
      VectorXd study_scores = sum_stats.scores.col(s);

      G_diag = Gsp.diagonal();
      SpMat beta_cov;
      if (sum_stats.found_at_least_one_conditional_variant) {
        Eigen::LDLT<MatrixXd> C_dcmp = C.ldlt();
        Gsp = Gsp - G_C * C_dcmp.solve(G_C.transpose());
        study_scores = study_scores - G_C * C_dcmp.solve(sum_stats.gene_cond_scores.col(s));
      }
      VectorXd D = (G_diag.array() > 0).select(
        sum_stats.ses.col(s).array() / G_diag.array().sqrt(),
        VectorXd::Zero(sum_stats.ses.rows())
      ).matrix();
      beta_cov = D.asDiagonal() * Gsp * D.asDiagonal();
      VectorXd beta_weights = (beta_cov.diagonal().array() > 0).select(
        beta_cov.diagonal().array().inverse(),
        0
      );

      for (int m = 0; m < sum_stats.nmasks; ++m) {
        double burden_score = sum_stats.masks.row(m) * study_scores;
        if (burden_score == 0) continue;
        //SpMat Gsp_mask = mask_selectors[m] * Gsp * mask_selectors[m].transpose();

        double burden_score_std = sqrt((sum_stats.masks.row(m) * Gsp * sum_stats.masks.row(m).transpose()).sum());
        double burden_chi = pow(burden_score / burden_score_std, 2);
        double cf = 1.0;
        double case_control_ratio = (double)sum_stats.mask_ncases[m] / (double)sum_stats.mask_ncontrols[m];
        if (this->trait_type == BT
            && case_control_ratio < this->burden_params.mask_spa_case_control_ratio
            && -abs(burden_score) / burden_score_std < this->burden_params.mask_spa_z_score) {

          double nhet = sum_stats.mask_cases_het(m, s) + sum_stats.mask_controls_het(m, s);
          double nhom_alt = sum_stats.mask_cases_alt(m, s) + sum_stats.mask_controls_alt(m, s);
          double nhom_ref = max(0.0, sum_stats.mask_ncases_per_study(m ,s) + sum_stats.mask_ncontrols_per_study(m ,s) - nhet - nhom_alt);
          double spa_chi = stat::spa::compute_spa_chival_from_geno_counts(burden_score,
                                                                          burden_score_std,
                                                                          nhom_ref,
                                                                          nhet,
                                                                          nhom_alt,
                                                                          sum_stats.mask_ncases_per_study(m ,s),
                                                                          sum_stats.mask_ncontrols_per_study(m ,s),
                                                                          false);
          if (spa_chi != stat::spa::SPA_FAILED) {
            cf = max(1.0, burden_chi / spa_chi);
          }
        }

        double weight_sum = (sum_stats.masks.row(m) * beta_weights).value();
        double beta_var = (
          (sum_stats.masks.row(m).array()*beta_weights.transpose().array()).matrix()/weight_sum 
          * beta_cov
          * (sum_stats.masks.row(m).transpose().array()*beta_weights.array()).matrix()/weight_sum
        ).value();
        double mask_se = sqrt( beta_var );
        sum_stats.mask_betas(m, s) = mask_se * burden_score / burden_score_std;

        stat::tests::test_result_t burden_pval = stat::tests::wst_burden(
          burden_score*burden_score, cf*burden_score_std*burden_score_std
        );
        if (burden_pval == stat::tests::TEST_FAILED) {
          log_warning("burden test failed for " + g.get_name() + "." + mask_set.get_mask_alt(m) + " in study " + to_string(s+1));
          continue;
        }
        double se = stat::misc::get_se_from_beta_log10p(sum_stats.mask_betas(m, s), burden_pval.log10p);
        sum_stats.mask_ses(m, s) = se;

        if (this->write_cohort_burden_tests) {
          string bin_name = mask_set.get_mask_alt(m);
          htpv4_record_t result = {
            g.get_name() + "." + bin_name,
            g.get_chrom(),
            g.get_start(),
            "ref",
            bin_name,
            this->trait_name,
            this->cohorts[s],
            "REMETA-BURDEN-META",
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            HTPv4_NA,
            map<string, string>()
          };
          result.num_cases = sum_stats.mask_ncases_per_study(m, s);
          result.cases_het = sum_stats.mask_cases_het(m, s);
          result.cases_alt = sum_stats.mask_cases_alt(m, s);
          result.cases_ref = max(0.0, sum_stats.mask_ncases_per_study(m, s) - result.cases_het - result.cases_alt);
          if (this->trait_type == BT) {
            result.num_controls = sum_stats.mask_ncontrols_per_study(m, s);
            result.controls_het = sum_stats.mask_controls_het(m, s);
            result.controls_alt = sum_stats.mask_controls_alt(m, s);
            result.controls_ref = max(0.0, sum_stats.mask_ncontrols_per_study(m, s) - result.controls_het - result.controls_alt);
          }
          result.aaf = sum_stats.mask_aafs(m, s);
          result.info[HTPv4_ESTIMATED_GENOTYPE_COUNTS_FLAG] = "";

          result.pval = burden_pval.pval;
          result.info["LOG10P"] = to_string(-burden_pval.log10p);

          if (this->trait_type == QT) {
            result.effect = sum_stats.mask_betas(m, s);
            result.lci_effect = sum_stats.mask_betas(m, s) - 1.96 * sum_stats.mask_ses(m, s);
            result.uci_effect = sum_stats.mask_betas(m, s) + 1.96 * sum_stats.mask_ses(m, s);
          } else {
            result.effect = exp(sum_stats.mask_betas(m, s));
            result.lci_effect = exp(sum_stats.mask_betas(m, s) - 1.96 * sum_stats.mask_ses(m, s));
            result.uci_effect = exp(sum_stats.mask_betas(m, s) + 1.96 * sum_stats.mask_ses(m, s));
            result.info["BETA"] = to_string(sum_stats.mask_betas(m, s));
          }
          result.info["SE"] = to_string(se);
          this->burden_writer->writerec(result);
        }
      }
    }
  }
  d = std::chrono::steady_clock::now() - start;
  log_debug(" * loading ld matrices: " + to_string(d.count()) + " seconds");

  for (int m = 0; m < sum_stats.nmasks; ++m) {
    double mask_aac = 0;
    double mask_an = 0;
    for (int s = 0; s < nstudies; ++s) {
      mask_aac += sum_stats.mask_cases_het(m, s) + 2*sum_stats.mask_cases_alt(m, s);
      mask_an += 2*sum_stats.mask_ncases_per_study(m, s);

      if (this->trait_type == BT) {
        mask_aac += sum_stats.mask_controls_het(m, s) + 2*sum_stats.mask_controls_alt(m, s);
        mask_an += 2*sum_stats.mask_ncontrols_per_study(m, s);
      }
    }
    sum_stats.mask_aaf[m] = mask_aac / mask_an;

    sum_stats.mask_nhet[m] = sum_stats.mask_cases_het.row(m).sum() + sum_stats.mask_controls_het.row(m).sum();
    sum_stats.mask_nhom_alt[m] = sum_stats.mask_cases_alt.row(m).sum() + sum_stats.mask_controls_alt.row(m).sum();
    sum_stats.mask_nhom_ref[m] = max(0.0, sum_stats.mask_ncases[m] + sum_stats.mask_ncontrols[m] - sum_stats.mask_nhet[m] - sum_stats.mask_nhom_alt[m]);
    sum_stats.mask_case_control_ratio[m] = (double)sum_stats.mask_ncases[m] / (double)sum_stats.mask_ncontrols[m];
  }

  sum_stats.scores_sum = sum_stats.scores.rowwise().sum();
  sum_stats.scorev_sum = sum_stats.scorev.rowwise().sum();

  return sum_stats;
}

void HTPMetaAnalyzer::write_snplist(const GeneSumStats& sum_stats, const Gene& g) {
  log_debug("writing mask snplist");
  for (int m = 0; m < sum_stats.nmasks; ++m) {
    if (sum_stats.mask_variants[m].size() == 0) continue;

    this->snplist_writer->write(g.get_name() + "." + mask_set.get_mask_alt(m) + "\t");
    for (size_t j = 0; j < sum_stats.mask_variants[m].size() - 1; ++j) {
      this->snplist_writer->write(sum_stats.mask_variants[m][j] + ",");
    }
    if (sum_stats.mask_variants[m].size() > 0) {
      this->snplist_writer->write(sum_stats.mask_variants[m][sum_stats.mask_variants[m].size() - 1] + "\n");
    }
  }
}

void HTPMetaAnalyzer::compute_conditional_variant_stats(GeneSumStats& sum_stats) {
  sum_stats.conditional_variant_info = "";
  if (sum_stats.found_at_least_one_conditional_variant) {
    log_debug("computing conditional statistics");

    Eigen::LDLT<MatrixXd> buffer_ld_mat_dcmp = MatrixXd(sum_stats.buffer_ld_mat).ldlt();
    MatrixXd gene_buffer_ld_mat_dense = MatrixXd(sum_stats.gene_buffer_ld_mat);
    VectorXd cond_scores_sum = sum_stats.gene_cond_scores.rowwise().sum();

    sum_stats.scores_sum = sum_stats.scores_sum - gene_buffer_ld_mat_dense * buffer_ld_mat_dcmp.solve(cond_scores_sum);
    sum_stats.gene_ld_mat = sum_stats.gene_ld_mat - gene_buffer_ld_mat_dense * buffer_ld_mat_dcmp.solve(gene_buffer_ld_mat_dense.transpose());
    sum_stats.scorev_sum = MatrixXd(sum_stats.gene_ld_mat).diagonal();

    for ( const variant_id& vid : sum_stats.gene_conditional_variants ) {
      if (sum_stats.conditional_variant_info.size() > 0) {
        sum_stats.conditional_variant_info += "," + vid; 
      } else {
        sum_stats.conditional_variant_info = vid;
      }
    }
  }
}

void HTPMetaAnalyzer::compute_variant_spa(GeneSumStats& sum_stats) {
  sum_stats.variant_cf_sqrt = VectorXd::Ones(sum_stats.nvariants);
  sum_stats.variant_z_scores = VectorXd::Zero(sum_stats.nvariants);
  if (this->sv_spa_z_score != SPA_NEVER_APPLIED && sum_stats.max_ccr < this->sv_spa_case_control_ratio) {
    log_debug("applying per variant spa");
    double z;
    double spa_chi;
    double cf;
    for (int j = 0; j < sum_stats.nvariants; ++j) {
      z = -abs(sum_stats.scores_sum(j) / sqrt(sum_stats.scorev_sum(j)));
      sum_stats.variant_z_scores(j) = z;
      if (z < this->sv_spa_z_score) {
        spa_chi = stat::spa::compute_spa_chival_from_geno_counts(
          sum_stats.scores_sum(j),
          sqrt(sum_stats.scorev_sum(j)),
          (sum_stats.ncontrols.row(j) + sum_stats.ncases.row(j) - (sum_stats.controls_het.row(j) + sum_stats.controls_alt.row(j) + sum_stats.cases_het.row(j) + sum_stats.cases_alt.row(j))).sum(),
          (sum_stats.controls_het.row(j) + sum_stats.cases_het.row(j)).sum(),
          (sum_stats.controls_alt.row(j) + sum_stats.cases_alt.row(j)).sum(),
          (sum_stats.ncases.row(j)).sum(),
          (sum_stats.ncontrols.row(j)).sum(),
          false
        );
        cf = max(1.0, z*z / spa_chi);
        if (spa_chi != stat::spa::SPA_FAILED && cf > 1) {
          sum_stats.variant_cf_sqrt(j) = sqrt(cf);
        }
      }
    }
  }
}

void HTPMetaAnalyzer::compute_burden(vector<htpv4_record_t>& burden_results,
                                     GeneSumStats& sum_stats, // modifies mask_cf
                                     const SpMat& ld_mat_meta,
                                     const htpv4_record_t& result_template,
                                     int mask_idx) {
  log_debug("computing burden");
  int m = mask_idx;
  double Qburden = sum_stats.masks.row(m) * sum_stats.burden_scores;
  Qburden *= Qburden;

  SpMat D(sum_stats.mask_indices[m].size(), sum_stats.mask_indices[m].size());
  D.reserve(VectorXi::Constant(sum_stats.mask_indices[m].size(), 1));
  for (size_t v = 0; v < sum_stats.mask_indices[m].size(); ++v) {
    double c = 1;
    if (sum_stats.variant_z_scores(sum_stats.mask_indices[m][v]) < this->burden_params.sv_spa_z_score
        && sum_stats.max_ccr < this->burden_params.sv_spa_case_control_ratio) {
      c *= sum_stats.variant_cf_sqrt(sum_stats.mask_indices[m][v]);
    }
    D.insert(v, v) = c*sum_stats.burden_weights(sum_stats.mask_indices[m])[v];
  }
  SpMat burden_ld_mat = D * ld_mat_meta * D;

  VectorXd variant_cf_sqrt_mask = VectorXd::Ones(sum_stats.mask_indices[m].size());
  if (sum_stats.max_ccr < this->burden_params.sv_spa_case_control_ratio) {
    variant_cf_sqrt_mask = (sum_stats.variant_z_scores(sum_stats.mask_indices[m]).array() < this->burden_params.sv_spa_z_score).select(
      sum_stats.variant_cf_sqrt(sum_stats.mask_indices[m]).array(),
      1
    );
  }
  double burden_score = sum_stats.masks.row(m) * sum_stats.scores_sum;
  double burden_score_std = sqrt(
    variant_cf_sqrt_mask.transpose() * ld_mat_meta * variant_cf_sqrt_mask
  );
  double burden_z_score = -abs(burden_score / burden_score_std);
  double burden_chi = pow(burden_score / burden_score_std, 2);
  double cf = 1.0;
  if (this->trait_type == BT
      && sum_stats.mask_case_control_ratio[m] < this->burden_params.mask_spa_case_control_ratio
      && burden_z_score < this->burden_params.mask_spa_z_score) {
    log_debug("computing burden mask correction factor");
    double spa_chi = stat::spa::compute_spa_chival_from_geno_counts(burden_score,
                                                                    burden_score_std,
                                                                    sum_stats.mask_nhom_ref[m],
                                                                    sum_stats.mask_nhet[m],
                                                                    sum_stats.mask_nhom_alt[m],
                                                                    sum_stats.mask_ncases[m],
                                                                    sum_stats.mask_ncontrols[m],
                                                                    false);
    if (spa_chi != stat::spa::SPA_FAILED) {
      cf = max(1.0, burden_chi / spa_chi);
    }
    burden_ld_mat = burden_ld_mat * cf;
  }
  sum_stats.mask_cf[m] = cf;

  htpv4_record_t burden_result = result_template;
  burden_result.model = "REMETA-BURDEN-META";
  burden_result.cases_het = sum_stats.mask_cases_het.row(m).sum();
  burden_result.cases_alt = sum_stats.mask_cases_alt.row(m).sum();
  burden_result.cases_ref = max(0.0, sum_stats.mask_ncases[m] - burden_result.cases_het - burden_result.cases_alt);
  if (this->trait_type == BT) {
    burden_result.controls_het = sum_stats.mask_controls_het.row(m).sum();
    burden_result.controls_alt = sum_stats.mask_controls_alt.row(m).sum();
    burden_result.controls_ref = max(0.0, sum_stats.mask_ncontrols[m] - burden_result.controls_het - burden_result.controls_alt);
  }

  burden_result.aaf = sum_stats.mask_aaf[m];
  burden_result.info["Direction"] = "";
  burden_result.info[HTPv4_ESTIMATED_GENOTYPE_COUNTS_FLAG] = "";
  for (int k = 0; k < this->nstudies; ++k) {
    if (sum_stats.mask_betas(m, k) > 0) {
      burden_result.info["Direction"] += "+";
    } else if (sum_stats.mask_betas(m, k) < 0) {
      burden_result.info["Direction"] += "-";
    } else {
      burden_result.info["Direction"] += "?";
    }
  }
  pair<double, double> effects = stat::misc::meta_analyze_effects(
    sum_stats.mask_betas.row(m),
    sum_stats.mask_ses.row(m)
  );
  stat::tests::test_result_t burden_pval = stat::tests::wst_burden(Qburden, burden_ld_mat.sum());
  burden_result.pval = burden_pval.pval;
  sum_stats.mask_pval[m] = burden_pval.pval;
  double se = stat::misc::get_se_from_beta_log10p(effects.first, burden_pval.log10p);
  burden_result.info["LOG10P"] = to_string(-burden_pval.log10p);
  burden_result.info["CF"] = to_string(cf);

  if (this->trait_type == QT) {
    burden_result.effect = effects.first;
    burden_result.lci_effect = effects.first - 1.96 * se;
    burden_result.uci_effect = effects.first + 1.96 * se;
  } else {
    burden_result.effect = exp(effects.first);
    burden_result.lci_effect = exp(effects.first - 1.96 * se);
    burden_result.uci_effect = exp(effects.first + 1.96 * se);
    burden_result.info["BETA"] = to_string(effects.first);
  }
  burden_result.info["SE"] = to_string(se);
  burden_result.info["NVAR"] = to_string(sum_stats.mask_indices[m].size());

  if (sum_stats.conditional_variant_info != "") {
    burden_result.info["CONDITIONAL"] = sum_stats.conditional_variant_info;
  }
  if (burden_pval == stat::tests::TEST_FAILED) {
    log_warning("burden test failed for " + burden_result.name);
  } else {
    burden_results.push_back(burden_result);
  }
}

void HTPMetaAnalyzer::compute_skato(vector<htpv4_record_t>& skato_results,
                                    const GeneSumStats& sum_stats,
                                    const SpMat& ld_mat_meta,
                                    const htpv4_record_t& result_template,
                                    int mask_idx) {
  log_debug("computing SKATO");
  int m = mask_idx;

  SpMat D(sum_stats.mask_indices[m].size(), sum_stats.mask_indices[m].size());
  D.reserve(VectorXi::Constant(sum_stats.mask_indices[m].size(), 1));
  for (size_t v = 0; v < sum_stats.mask_indices[m].size(); ++v) {
    double c = 1;
    if (sum_stats.variant_z_scores(sum_stats.mask_indices[m][v]) < this->skato_params.sv_spa_z_score
        && sum_stats.max_ccr < this->skato_params.sv_spa_case_control_ratio) {
      c *= sum_stats.variant_cf_sqrt(sum_stats.mask_indices[m][v]);
    }
    D.insert(v, v) = c*sum_stats.skato_weights(sum_stats.mask_indices[m])[v];
  }
  //D.makeCompressed();
  MatrixXd skato_ld_mat = sum_stats.mask_collapsors[m] * D * ld_mat_meta * D * sum_stats.mask_collapsors[m].transpose();

  VectorXd variant_cf_sqrt_mask = VectorXd::Ones(sum_stats.mask_indices[m].size());
  if (sum_stats.max_ccr < this->skato_params.sv_spa_case_control_ratio) {
    variant_cf_sqrt_mask = (sum_stats.variant_z_scores(sum_stats.mask_indices[m]).array() < this->skato_params.sv_spa_z_score).select(
      sum_stats.variant_cf_sqrt(sum_stats.mask_indices[m]).array(),
      1
    );
  }
  double burden_score = sum_stats.masks.row(m) * sum_stats.scores_sum;
  double burden_score_std = sqrt(
    variant_cf_sqrt_mask.transpose() * ld_mat_meta * variant_cf_sqrt_mask
  );
  double burden_z_score = -abs(burden_score / burden_score_std);
  double burden_chi = pow(burden_score / burden_score_std, 2);
  double cf = 1.0;
  if (this->trait_type == BT
      && sum_stats.mask_case_control_ratio[m] < this->skato_params.mask_spa_case_control_ratio
      && burden_z_score < this->skato_params.mask_spa_z_score) {
    log_debug("computing SKATO mask correction factor");
    double spa_chi = stat::spa::compute_spa_chival_from_geno_counts(burden_score,
                                                                    burden_score_std,
                                                                    sum_stats.mask_nhom_ref[m],
                                                                    sum_stats.mask_nhet[m],
                                                                    sum_stats.mask_nhom_alt[m],
                                                                    sum_stats.mask_ncases[m],
                                                                    sum_stats.mask_ncontrols[m],
                                                                    false);
    if (spa_chi != stat::spa::SPA_FAILED) {
      cf = max(1.0, burden_chi / spa_chi);
    }
    skato_ld_mat = skato_ld_mat * cf;
  }

  htpv4_record_t skat_result = result_template;
  skat_result.model = "REMETA-SKAT-META";
  htpv4_record_t skato_result = result_template;
  skato_result.model = "REMETA-SKATO-ACAT-META";

  vector<double> log10_pvals;
  double Qskat = (sum_stats.mask_collapsors[m] * sum_stats.mask_selectors[m] * sum_stats.skato_scores).array().square().sum();
  double Qburden = (sum_stats.mask_selectors[m] * sum_stats.skato_scores).array().sum();
  Qburden *= Qburden;
  for (const double& rho : this->skato_params.rho_values) {
    stat::tests::test_result_t test_result = stat::tests::skato(
      Qskat,
      Qburden,
      rho,
      skato_ld_mat
    );
    if (test_result == stat::tests::TEST_FAILED) {
      continue;
    }

    log10_pvals.push_back(test_result.log10p);
    if (rho == 0) {
      skat_result.pval = test_result.pval;
      skat_result.info["LOG10P"] = to_string(-test_result.log10p);
      skat_result.info["CF"] = to_string(cf);
    }
  }

  int nvar = 0;
  for (ssize_t v = 0; v < skato_ld_mat.rows(); ++v) {
    if (skato_ld_mat(v, v) > 0) {
      ++nvar;
    }
  }
  skato_result.info["NVAR"] = to_string(nvar);

  stat::tests::test_result_t skato_acat = stat::tests::acat(log10_pvals);
  skato_result.pval = skato_acat.pval;
  skato_result.info["LOG10P"] = to_string(-skato_acat.log10p);
  skato_result.info["CF"] = to_string(cf);
  if (sum_stats.conditional_variant_info != "") {
    skato_result.info["CONDITIONAL"] = sum_stats.conditional_variant_info;
  }
  if (skato_acat == stat::tests::TEST_FAILED || log10_pvals.size() == 0) {
    log_warning("SKATO test failed for " + skato_result.name);
  } else {
    skato_results.push_back(skato_result);
  }

  if (skat_result.pval == HTPv4_NA) {
    stat::tests::test_result_t test_result = stat::tests::skato(Qskat, Qburden, 0, ld_mat_meta);
    skat_result.pval = test_result.pval;
    skat_result.info["LOG10P"] = to_string(-test_result.log10p);
    skat_result.info["CF"] = to_string(cf);
  }

  if (sum_stats.conditional_variant_info != "") {
    skat_result.info["CONDITIONAL"] = sum_stats.conditional_variant_info;
  }
  if (skat_result.pval < 0) {
    log_warning("SKAT test failed for " + skat_result.name);
  } else {
    skato_results.push_back(skat_result);
  }
}

void HTPMetaAnalyzer::compute_acatv(vector<htpv4_record_t>& acatv_results,
                                    const GeneSumStats& sum_stats,
                                    const htpv4_record_t& result_template,
                                    int mask_idx) {
  log_debug("computing ACATV");
  htpv4_record_t acatv_result = result_template;
  acatv_result.model = "REMETA-ACATV-META";

  vector<double> mask_log10_pvals;
  vector<double> mask_weights;
  int nvar = 0;
  double z = 0;

  VectorXd Qacatv = sum_stats.masks.row(mask_idx).transpose().array() * sum_stats.scores_sum.array();
  VectorXd weights = sum_stats.masks.row(mask_idx).transpose().array() * sum_stats.acatv_weights.array();

  for (ssize_t v = 0; v < Qacatv.size(); ++v) {
    double c = 1;
    if (sum_stats.variant_z_scores(v) < this->acatv_params.sv_spa_z_score
        && sum_stats.mask_case_control_ratio(mask_idx) < this->acatv_params.sv_spa_case_control_ratio) {
      c *= sum_stats.variant_cf_sqrt(v);
    }
    if (weights[v] > 0 && sum_stats.scorev_sum(v) > 0) {
      z = -abs(Qacatv[v]) / (c*sqrt(sum_stats.scorev_sum(v)));
      mask_weights.push_back(weights[v]);
      mask_log10_pvals.push_back(
        min(
          log10(0.99),
          stat::tests::normal(z, true).log10p
        )
      );
      ++nvar;
    }
  }

  if (nvar > 0) {
    stat::tests::test_result_t test_result = stat::tests::weighted_acat(mask_log10_pvals, mask_weights);
    acatv_result.pval = test_result.pval;
    acatv_result.info["LOG10P"] = to_string(-test_result.log10p);
    acatv_result.info["NVAR"] = to_string(nvar);

    if (sum_stats.conditional_variant_info != "") {
      acatv_result.info["CONDITIONAL"] = sum_stats.conditional_variant_info;
    }
    acatv_results.push_back(acatv_result);
  }
}

void HTPMetaAnalyzer::set_write_cohort_burden_tests(HTPv4Writer& burden_writer) {
  this->write_cohort_burden_tests = true;
  this->burden_writer = &burden_writer;
  this->burden_writer->writeheader();
}

void HTPMetaAnalyzer::set_write_mask_snplist(BgzWriter& snplist_writer) {
  this->write_mask_snplist = true;
  this->snplist_writer = &snplist_writer;
}

void HTPMetaAnalyzer::set_ignore_mask_ld() {
  this->ignore_mask_ld = true;
}

void HTPMetaAnalyzer::set_recompute_score() {
  this->recompute_score = true;
}

void HTPMetaAnalyzer::set_keep_variants_not_in_ld_mat() {
  this->keep_variants_not_in_ld_mat = true;
}

void HTPMetaAnalyzer::set_write_freqs(AlleleFreqWriter& freq_writer) {
  this->write_freqs = true;
  this->freq_writer = &freq_writer;
}

void HTPMetaAnalyzer::check_info(const htpv4_record_t& rec) {
  if (rec.info.find("SCORE") == rec.info.end() && !this->recompute_score) {
    throw runtime_error("missing SCORE from INFO field of: " + rec.name);
  }
}

void HTPMetaAnalyzer::set_collapse_anno(const string& annotation) {
  this->collapse_anno_int = this->mask_set.anno_to_int(annotation);
  this->collapse_anno = true;
}

double HTPMetaAnalyzer::get_score(const htpv4_record_t& rec) {
  double beta = HTPv4Reader::get_beta(rec);
  if (this->recompute_score && rec.info.count("SCORE") == 0) {
    double se = HTPv4Reader::get_se(rec);
    return beta / (se * se);
  } else {
    // Some versions of regenie report the score of the minor allele
    // instead of the alternate allele. In those cases the score is flipped.
    // Since the beta is reported with respect to the alt allele, we can use
    // it to check the sign.
    double score = stod(rec.info.at("SCORE"));
    return std::copysign(score, beta);
  }
}

double HTPMetaAnalyzer::get_scorev(const htpv4_record_t& rec) {
  double score = this->get_score(rec);
  double pval = rec.pval;
  double scorev = 0;
  boost::math::normal dist(0, 1);
  if (rec.info.count("SKATV") == 0) {
    double z = quantile(dist, max(pval/2, std::numeric_limits<double>::min()));
    scorev = pow(score / z, 2);
  } else {
    // The scorev info field only has 6 decimal places. It will be truncated to 0
    // if scorev<0.000001. To future-proof, we need to clip the value to some small
    // number above 0. If this triggers a massive p-value from the score test,
    // it will get re-adjusted by regenie's pvalue.
    scorev = max(stod(rec.info.at("SKATV")), 1e-10);

    // check that p-value from score test matches reported p-value
    double p = 2*cdf(dist, -abs(score)/sqrt(scorev));
    if (log10(p) < log10(pval) - 2 && pval > 1e-300) {
      double z = quantile(dist, max(pval/2, std::numeric_limits<double>::min()));
      scorev = pow(score / z, 2);
    }
  }

  // Check if variance of score statistic is unusually high. This signals
  // that something went wrong in the Firth step.
  if (this->trait_type == BT) {
    double aac = rec.cases_het + rec.controls_het + 2*(rec.cases_alt + rec.controls_alt);
    double mu = (double)rec.num_cases / (double)(rec.num_cases + rec.num_controls);

    // this corresponds to a variance ~2 orders of magnitude larger than expected
    if (aac == 1 && scorev > 100*aac*mu*(1-mu)) {
      scorev = -1;
    }
  }

  return scorev;
}

bool HTPMetaAnalyzer::ld_mat_required() {
  return !this->ignore_mask_ld
    || !this->keep_variants_not_in_ld_mat
    || this->conditional_variants.size() > 0;
}

vector<htpv4_record_t> HTPMetaAnalyzer::meta_analyze_gene(Gene g) {
  log_debug("entering HTPMetaAnalyzer::meta_analyze_gene");
  log_debug("meta-analyzing " + g.get_name());

  std::chrono::time_point<std::chrono::steady_clock> start;
  std::chrono::duration<double> d;

  start = std::chrono::steady_clock::now();
  vector<variant_id> variants_found;
  vector<int> variant_annotations;
  this->get_variants(variants_found, variant_annotations, g);
  d = std::chrono::steady_clock::now() - start;
  log_debug("set variants: " + to_string(d.count()) + " seconds");
  log_debug(to_string(variants_found.size()) + " variants found");

  start = std::chrono::steady_clock::now();
  log_debug("collecting summary statistics...");
  GeneSumStats sum_stats = this->collect_sum_stats(variants_found, variant_annotations, this->mask_set, g);
  d = std::chrono::steady_clock::now() - start;
  log_debug("collect sum stats: " + to_string(d.count()) + " seconds");

  if (sum_stats.masks.sum() == 0) {
    log_warning("no variants found in " + g.get_name() + " masks");
    return vector<htpv4_record_t>();
  }

  if (this->ld_mat_required() && !sum_stats.found_at_least_one_ld_mat) {
    log_warning("no study contains an LD matrix for " + g.get_name() + ", skipping...");
    return vector<htpv4_record_t>();
  }

  if (this->write_mask_snplist) {
    this->write_snplist(sum_stats, g);
  }

  if (sum_stats.found_at_least_one_conditional_variant) {
    this->compute_conditional_variant_stats(sum_stats);
  }

  this->compute_variant_spa(sum_stats);

  log_debug("computing gene tests...");
  start = std::chrono::steady_clock::now();
  sum_stats.burden_scores = sum_stats.scores_sum.array() * sum_stats.burden_weights.array();
  sum_stats.skato_scores =  sum_stats.scores_sum.array() * sum_stats.skato_weights.array();
  vector<htpv4_record_t> burden_results;
  vector<htpv4_record_t> skato_results;
  vector<htpv4_record_t> acatv_results;

  htpv4_record_t result_template_gene {
    g.get_name(),
    g.get_chrom(),
    g.get_start(),
    "ref",
    "alt",
    this->trait_name,
    "",
    "",
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    map<string, string>()
  };
  result_template_gene.info["SOURCE"] = "REMETA-GB";

  for (int m = 0; m < sum_stats.nmasks; ++m) {
    string mask_cohort_meta = util::format_cohort_meta(
      this->cohorts,
      vector<int>(sum_stats.mask_cohort_idx[m].begin(),
      sum_stats.mask_cohort_idx[m].end())
    );
    string bin_name = mask_set.get_mask_alt(m);

    htpv4_record_t result_template = result_template_gene;
    result_template.name = g.get_name() + "." + bin_name;
    result_template.alt = bin_name;
    result_template.cohort = mask_cohort_meta;
    result_template.info["MAX_MAC"] = to_string((int)sum_stats.mask_max_mac(m));

    result_template.num_cases = sum_stats.mask_ncases[m];
    if (this->trait_type == BT) {
      result_template.num_controls = sum_stats.mask_ncontrols[m];
    }

    log_debug("running mask " + bin_name);
    if (sum_stats.mask_indices[m].size() == 0) {
      log_warning("mask " + g.get_name() + "." + bin_name + " does not have any variants");
      continue;
    }

    SpMat ld_mat_meta = sum_stats.mask_selectors[m] * sum_stats.gene_ld_mat * sum_stats.mask_selectors[m].transpose();
    if ((ld_mat_meta.diagonal().array() == 0).all()) {
      // likely including imputed variants in mask, or as mismatch between variants in the LD matrix
      // and variants in the summary statistics file
      log_warning("found zero LD matrix for " + g.get_name() + "." + bin_name + ", skipping...");
      continue;
    }

    log_debug("computing p-values");
    if (vector_contains(this->burden_params.af_bins, this->mask_set.get_mask_bin(m))) {
      this->compute_burden(burden_results, sum_stats, ld_mat_meta, result_template, m);
    }

    if (vector_contains(this->skato_params.af_bins, this->mask_set.get_mask_bin(m))) {
      this->compute_skato(skato_results, sum_stats, ld_mat_meta, result_template, m);
    }

    if (vector_contains(this->acatv_params.af_bins, this->mask_set.get_mask_bin(m))) {
      this->compute_acatv(acatv_results, sum_stats, result_template, m);
    }
  }

  vector<htpv4_record_t> results;
  for (const htpv4_record_t& rec : burden_results) {
    results.push_back(rec);
  }
  for (const htpv4_record_t& rec : skato_results) {
    results.push_back(rec);
  }
  for (const htpv4_record_t& rec : acatv_results) {
    results.push_back(rec);
  }
  d = std::chrono::steady_clock::now() - start;
  log_debug("computing gene tests: " + to_string(d.count()) + " seconds");
  return results;
}
