// Run some or all gene-based tests in a single pass.
#include "htp_meta_analyzer.hpp"

#include <algorithm>
#include <chrono>
#include <exception>
#include <string>
#include <unordered_map>
using std::string;
using std::unordered_map;

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
typedef Eigen::Triplet<double> Triplet;

bool vector_contains(const vector<double>& vec, const double& val) {
  return std::find(vec.begin(), vec.end(), val) != vec.end();
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
  this->burden_weight  = 0;
  this->skato_weight   = 0;
  this->acatv_weight   = 0;
  this->is_singleton   = false;
}

HTPMetaAnalyzer::HTPMetaAnalyzer(const string& trait_name,
                                 const trait_type_e& trait_type,
                                 const vector<string>& cohorts,
                                 const af_strategy_e& af_strategy,
                                 const string& mask_def_file,
                                 const string& annotation_file,
                                 const vector<string>& ld_file_prefixes,
                                 const double& spa_pval,
                                 const double& spa_case_control_ratio)
 : trait_name(trait_name)
 , trait_type(trait_type)
 , cohorts(cohorts)
 , af_strategy(af_strategy)
 , mask_def_file(mask_def_file)
 , annotation_file(annotation_file)
 , ld_file_prefixes(ld_file_prefixes)
 , spa_case_control_ratio(spa_case_control_ratio)
 , af_file("")
 , nstudies(cohorts.size())
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
 , keep_variants_not_in_ld_mat(false) {
  for (const string& prefix : ld_file_prefixes) {
    ld_mat_readers.push_back(
      RemetaMatrixReader()
    );
    ld_mat_readers[ld_mat_readers.size()-1].open(prefix);
  }

  if (spa_pval <= 0) {
    this->spa_z_score = SPA_NEVER_APPLIED;
  } else if (spa_pval >= 1) {
    this->spa_z_score = SPA_ALWAYS_APPLIED;
  } else {
    boost::math::normal s;
    this->spa_z_score = quantile(s, spa_pval / 2); // two-sided pvalue
  }
}

void HTPMetaAnalyzer::set_af_file(const string& af_file) {
  log_info("using external allele frequencies");
  this->af_strategy = USE_EXTERNAL_AF;
  this->af_file = af_file;
  this->af_map.load(af_file);
}

void HTPMetaAnalyzer::set_conditional_variants(vector<vector<htpv4_record_t> > study_conditional_variants) {
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
}

void HTPMetaAnalyzer::set_run_burden(const vector<double>& af_bins,
                                     const singleton_def_e& singleton_def,
                                     const weight_strategy_e& weight_strategy) {
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
}

void HTPMetaAnalyzer::set_run_skato(const vector<double>& af_bins,
                                    const vector<double>& rho_values,
                                    const weight_strategy_e& weight_strategy,
                                    const int& min_aac) {
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
  std::sort(this->af_bins.begin(), this->af_bins.end());
  this->mask_set = MaskSet(this->mask_def_file, this->af_bins);
}

void HTPMetaAnalyzer::set_run_acatv(const vector<double>& af_bins,
                                    const weight_strategy_e& weight_strategy,
                                    const int& min_aac) {
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
    log_debug(rec.name + " in study " + to_string(study_index + 1) + " has bad skatv value, skipping...");
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

vector<htpv4_record_t> HTPMetaAnalyzer::meta_analyze_gene(Gene g) {
  log_debug("entering HTPMetaAnalyzer::meta_analyze_gene");
  log_debug("meta-analyzing " + g.get_name());

  std::chrono::time_point<std::chrono::steady_clock> start;
  std::chrono::duration<double> d;

  vector<variant_id>                    variants_found;
  unordered_map<variant_id, int>        variant_idx;
  vector<int>                           variant_annotations;

  boost::math::beta_distribution<> beta(1, 25);
  start = std::chrono::steady_clock::now();
  int variants_without_af = 0;
  for (const string& vid : *g.get_variants()) {
    if (this->variants.count(vid) > 0) {
      HTPMetaVariant *variant = &this->variants[vid];

      int anno = anno_map.get_annotation(vid, g.get_name());
      if (anno == -1) {
        continue;
      }

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

      if (variant->aaf <= this->mask_max_aaf) {
        variant_idx[vid] = variants_found.size();
        variants_found.push_back(vid);
        variant_annotations.push_back(anno);

        double weight_aaf = variant->aac / variant->an;
        double weight_maf = weight_aaf > 0.5 ? 1 - weight_aaf : weight_aaf;
        if (this->burden_params.weight_strategy == USE_BETA_WEIGHTS) {
          variant->burden_weight = pdf(beta, weight_maf);
        } else {
          variant->burden_weight = 1;
        }
        if (this->skato_params.weight_strategy == USE_BETA_WEIGHTS && variant->aac >= this->skato_params.min_aac) {
          variant->skato_weight = pdf(beta, weight_maf);
        } else if (variant->aac >= this->skato_params.min_aac) {
          variant->skato_weight = 1;
        } else {
          variant->skato_weight = 0;
        }
        if (this->acatv_params.weight_strategy == USE_BETA_WEIGHTS && variant->aac >= this->acatv_params.min_aac) {
          variant->acatv_weight = pow(pdf(beta, weight_maf), 2) * weight_maf * (1 - weight_maf);
        } else if (variant->aac >= this->acatv_params.min_aac) {
          variant->acatv_weight = 1;
        } else {
          variant->acatv_weight = 0;
        }
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
  d = std::chrono::steady_clock::now() - start;
  log_debug("collecting variants: " + to_string(d.count()) + " seconds");
  log_debug(to_string(variants_found.size()) + " variants found");

  start = std::chrono::steady_clock::now();
  int nmasks = this->mask_set.n_masks();
  int nstudies = this->cohorts.size();
  // Indicator matrix. Entry i, j is 1 if variants_found[j]
  // is contained in mask_set.get_mask(i).
  MatrixXd masks(nmasks, variants_found.size());
  MatrixXi mask_ncases_per_study = MatrixXi::Zero(nmasks, nstudies);
  MatrixXi mask_ncontrols_per_study = MatrixXi::Zero(nmasks, nstudies);
  MatrixXd scores(variants_found.size(), nstudies);
  MatrixXd scorev(variants_found.size(), nstudies);
  MatrixXd betas(variants_found.size(), nstudies);
  MatrixXd ses(variants_found.size(), nstudies);
  MatrixXd aafs(variants_found.size(), nstudies);
  MatrixXd cases_het(variants_found.size(), nstudies);
  MatrixXd cases_alt(variants_found.size(), nstudies);
  MatrixXd controls_het(variants_found.size(), nstudies);
  MatrixXd controls_alt(variants_found.size(), nstudies);
  VectorXd burden_weights(variants_found.size());
  VectorXd skato_weights(variants_found.size());
  VectorXd acatv_weights(variants_found.size());
  MatrixXi masks_max_mac = VectorXi::Zero(nmasks);
  vector<vector<int> > mask_indices(nmasks);
  vector<unordered_set<int> > mask_cohort_idx(nmasks);
  vector<double> study_sample_size(nstudies, 0);
  vector<vector<string> > mask_variants(nmasks);
  for (size_t j = 0; j < variants_found.size(); ++j) {
    HTPMetaVariant variant = this->variants[variants_found[j]];
    burden_weights[j] = variant.burden_weight;
    skato_weights[j] = variant.skato_weight;
    acatv_weights[j] = variant.acatv_weight;
    scores.row(j) = variant.scores;
    scorev.row(j) = variant.scorev;
    betas.row(j) = variant.betas;
    ses.row(j) = variant.ses;
    aafs.row(j) = variant.aafs;
    cases_het.row(j) = variant.cases_het;
    cases_alt.row(j) = variant.cases_alt;
    controls_het.row(j) = variant.controls_het;
    controls_alt.row(j) = variant.controls_alt;

    if (this->conditional_variants_set.count(variant.id) > 0) {
      scores.row(j).setZero();
      scorev.row(j).setZero();
    }

    for (int i = 0; i < nmasks; ++i) {
      if (mask_set.mask_contains(i, variant_annotations[j], variant.aaf, variant.is_singleton)) {
        masks(i, j) = 1;
        masks_max_mac(i) = max(masks_max_mac(i), variant.max_mac);
        mask_indices[i].push_back(j);
        mask_cohort_idx[i].insert(variant.cohort_idx.begin(), variant.cohort_idx.end());

        for (int k = 0; k < nstudies; ++k) {
          mask_ncases_per_study(i, k) = max(mask_ncases_per_study(i, k), variant.ncases(k));
          mask_ncontrols_per_study(i, k) = max(mask_ncontrols_per_study(i, k), variant.ncontrols(k));
        }

        if (this->write_mask_snplist) {
          mask_variants[i].push_back(variants_found[j]);
        }
      } else {
        masks(i, j) = 0;
      }
    }

    for (int k = 0; k < nstudies; ++k) {
      study_sample_size[k] = max(study_sample_size[k], variant.sample_sizes[k]);
    }
  }

  vector<int> mask_ncases(nmasks, 0);
  vector<int> mask_ncontrols(nmasks, 0);
  for (int i = 0; i < nmasks; ++i) {
    mask_ncases[i] = mask_ncases_per_study.row(i).sum();
    mask_ncontrols[i] = mask_ncontrols_per_study.row(i).sum();
  }

  d = std::chrono::steady_clock::now() - start;
  log_debug("building masks: " + to_string(d.count()) + " seconds");

  if (masks.sum() == 0) {
    log_warning("no variants found in " + g.get_name() + " masks");
    return vector<htpv4_record_t>();
  }

  if (this->write_mask_snplist) {
    log_debug("writing mask snplist");
    for (int m = 0; m < nmasks; ++m) {
      if (mask_variants[m].size() == 0) continue;

      this->snplist_writer->write(g.get_name() + "." + mask_set.get_mask_alt(m) + "\t");
      for (size_t j = 0; j < mask_variants[m].size() - 1; ++j) {
        this->snplist_writer->write(mask_variants[m][j] + ",");
      }
      if (mask_variants[m].size() > 0) {
        this->snplist_writer->write(mask_variants[m][mask_variants[m].size() - 1] + "\n");
      }
    }
  }

  start = std::chrono::steady_clock::now();
  vector<double> sample_sizes;
  vector<unordered_set<variant_id> > ld_variants_found(nstudies);
  MatrixXd mask_betas = MatrixXd::Zero(nmasks, nstudies);
  MatrixXd mask_ses = MatrixXd::Zero(nmasks, nstudies);
  MatrixXd mask_aafs = MatrixXd::Zero(nmasks, nstudies);
  MatrixXd mask_cases_het = MatrixXd::Zero(nmasks, nstudies);
  MatrixXd mask_cases_alt = MatrixXd::Zero(nmasks, nstudies);
  MatrixXd mask_controls_het = MatrixXd::Zero(nmasks, nstudies);
  MatrixXd mask_controls_alt = MatrixXd::Zero(nmasks, nstudies);
  bool found_at_least_one_ld_mat = false;
  bool found_at_least_one_conditional_variant = false;
  SpMat gene_ld_mat(variants_found.size(), variants_found.size());
  SpMat gene_buffer_ld_mat(variants_found.size(), this->conditional_variants.size());
  SpMat buffer_ld_mat(this->conditional_variants.size(), this->conditional_variants.size());
  unordered_set<string> gene_conditional_variants;
  for (int s = 0; s < nstudies; ++s) {
    if (this->ld_mat_required() && !ld_mat_readers[s].contains_gene(g.get_name())) {
      log_warning("LD matrix for " + g.get_name() + " is not found in study " + to_string(s+1));
      continue;
    }

    MatrixXf Gf;
    MatrixXf G_Cf;
    MatrixXf Cf;
    MatrixXd G;       // gene_ld_mat
    MatrixXd G_C;     // gene_buffer_ld_mat
    MatrixXd C;       // inverse of buffer_ld_mat

    if (this->conditional_variants.size() > 0) {
      log_debug("loading conditional LD matrices");
      ld_mat_readers[s].load_conditional_ld_mats(
        Gf,
        Cf,
        G_Cf,
        variants_found,
        this->conditional_variants,
        g.get_name()
      );
      G = Gf.cast<double>();
      Gf.resize(0, 0);
      C = Cf.cast<double>();
      Cf.resize(0, 0);
      G_C = G_Cf.cast<double>();
      G_Cf.resize(0, 0);
      found_at_least_one_conditional_variant = found_at_least_one_conditional_variant || (C.array() > 0).any();
    } else if (this->ld_mat_required()) {
      log_debug("loading marginal LD matrices");
      ld_mat_readers[s].load_gene_ld_mat(Gf, variants_found, g.get_name());
      G = Gf.cast<double>();
      Gf.resize(0, 0);
    } else {
      log_debug("assuming diagonal LD matrices");
      G = MatrixXd::Zero(variants_found.size(), variants_found.size());
    }

    if (this->ignore_mask_ld) {
      MatrixXd G_diag = G.diagonal();
      G = G_diag.asDiagonal();
    }

    for (size_t i = 0; i < variants_found.size(); ++i) {
      if (G(i, i) > 0 && scores(i, s) != 0) {
        ld_variants_found[s].insert(variants_found[i]);
        found_at_least_one_ld_mat = true;
      } else if (G(i, i) == 0 && scores(i, s) != 0 && !this->keep_variants_not_in_ld_mat) {
        //log_warning(variants_found[i] + " does not have an entry in LD matrix for study " + to_string(s+1));
        scores(i, s) = 0;
      } else if (G(i, i) == 0 && scores(i, s) != 0 && this->keep_variants_not_in_ld_mat) {
        ld_variants_found[s].insert(variants_found[i]);
        G(i, i) = scorev(i, s);
        found_at_least_one_ld_mat = true;
      }
    }

    // The LD matrix exists but contains none of the requested variants
    if (G.isZero(0)) {
      log_warning("LD matrix for study " + to_string(s+1) + " does not contain any of the requested variants");
      continue;
    }

    if (this->burden_params.params_set() || this->skato_params.params_set()) {
      log_debug("computing AAFs and genotype counts");
#if defined(_OPENMP)
      int nthreads = Eigen::nbThreads();
      Eigen::setNbThreads(1);
      #pragma omp parallel for schedule(dynamic)
#endif
      for (int m = 0; m < nmasks; ++m) {
        if (mask_indices[m].size() == 0) continue;
        mask_aafs(m, s) = max(
          1.0 / (2.0*study_sample_size[s]),
          stat::misc::estimate_burden_aaf(aafs(mask_indices[m], s), G, mask_indices[m])
        );
        mask_cases_het(m, s) = stat::misc::estimate_nhets(
          cases_het(mask_indices[m], s),
          aafs(mask_indices[m], s),
          G,
          mask_indices[m]
        );
        mask_cases_alt(m, s) = stat::misc::estimate_nalts(
          cases_alt(mask_indices[m], s),
          aafs(mask_indices[m], s),
          G,
          mask_indices[m]
        );
        if (this->trait_type == BT) {
          mask_controls_het(m, s) = stat::misc::estimate_nhets(
            controls_het(mask_indices[m], s),
            aafs(mask_indices[m], s),
            G,
            mask_indices[m]
          );
          mask_controls_alt(m, s) = stat::misc::estimate_nalts(
            controls_alt(mask_indices[m], s),
            aafs(mask_indices[m], s),
            G,
            mask_indices[m]
          );
        }
      }
#if defined(_OPENMP)
      Eigen::setNbThreads(nthreads);
#endif
    }

    if (this->burden_params.params_set()) {
      log_debug("computing effect sizes");
      MatrixXd beta_cov;
      MatrixXd C_inv;
      if (found_at_least_one_conditional_variant) {
        stat::misc::rescale_cov_block(
          G,
          C,
          G_C,
          ses.col(s),
          this->conditional_ses.col(s)
        );
        C_inv = MatrixXd(C).completeOrthogonalDecomposition().pseudoInverse().sparseView();
        betas.col(s) = betas.col(s) - G_C * C_inv * this->conditional_betas.col(s);
        beta_cov = G - G_C * C_inv * G_C.transpose();
      }
      for (int m = 0; m < nmasks; ++m) {
        if (!found_at_least_one_conditional_variant) {
          pair<double, double> effects = stat::misc::wst_get_effect_size_unconditional(
            betas(mask_indices[m], s),
            ses(mask_indices[m], s),
            G,
            mask_indices[m]
          );
          mask_betas(m, s) = effects.first;
          mask_ses(m, s) = effects.second;
        } else {
          pair<double, double> effects = stat::misc::wst_get_effect_size_unconditional(
            betas(mask_indices[m], s),
            beta_cov(mask_indices[m], mask_indices[m]).diagonal(),
            beta_cov,
            mask_indices[m]
          );
          mask_betas(m, s) = effects.first;
          mask_ses(m, s) = effects.second;
        }
      }
    }

    if (!found_at_least_one_conditional_variant) {
      log_debug("rescaling matrices");
      stat::misc::rescale_cov(G, scorev.col(s));
    } else {
      log_debug("rescaling matrices");
      stat::misc::rescale_cov_block(
        G,
        C,
        G_C,
        scorev.col(s).array(),
        this->conditional_scorev.col(s).array()
      );
      for (size_t i = 0; i < this->conditional_variants.size(); ++i) {
        if (C(i, i) > 0) {
          gene_conditional_variants.insert(this->conditional_variants[i]);
        }
      }
    }

    gene_ld_mat += G.sparseView();
    if (found_at_least_one_conditional_variant) {
      gene_buffer_ld_mat += G_C.sparseView();
      buffer_ld_mat += C.sparseView();
    }

    if (this->write_cohort_burden_tests) {
      VectorXd study_scores = scores.col(s);
      if (found_at_least_one_conditional_variant) {
        C = C.completeOrthogonalDecomposition().pseudoInverse().sparseView();
        G = G - G_C * C * G_C.transpose();
        study_scores = study_scores - G_C * C * this->conditional_scores.col(s);
      }

      for (int m = 0; m < nmasks; ++m) {
        double burden_score = masks.row(m) * study_scores;
        if (burden_score == 0) continue;
        double burden_score_std = sqrt(G(mask_indices[m], mask_indices[m]).sum());
        double burden_chi = pow(burden_score / burden_score_std, 2);
        double cf = 1.0;
        double case_control_ratio = (double)mask_ncases[m] / (double)mask_ncontrols[m];
        if (this->trait_type == BT 
            && case_control_ratio < this->spa_case_control_ratio
            && -abs(burden_score) / burden_score_std < this->spa_z_score) {

          double nhet = mask_cases_het(m, s) + mask_controls_het(m, s);
          double nhom_alt = mask_cases_alt(m, s) + mask_controls_alt(m, s);
          double nhom_ref = max(0.0, mask_ncases_per_study(m ,s) + mask_ncontrols_per_study(m ,s) - nhet - nhom_alt);
          double spa_chi = stat::spa::compute_spa_chival_from_geno_counts(burden_score,
                                                                          burden_score_std,
                                                                          nhom_ref,
                                                                          nhet,
                                                                          nhom_alt,
                                                                          mask_ncases_per_study(m ,s),
                                                                          mask_ncontrols_per_study(m ,s),
                                                                          false);
          if (spa_chi != stat::spa::SPA_FAILED) {
            cf = max(1.0, burden_chi / spa_chi);
          }
        }
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
        result.num_cases = mask_ncases_per_study(m, s);
        result.cases_het = mask_cases_het(m, s);
        result.cases_alt = mask_cases_alt(m, s);
        result.cases_ref = max(0.0, mask_ncases_per_study(m, s) - result.cases_het - result.cases_alt);
        if (this->trait_type == BT) {
          result.num_controls = mask_ncontrols_per_study(m, s);
          result.controls_het = mask_controls_het(m, s);
          result.controls_alt = mask_controls_alt(m, s);
          result.controls_ref = max(0.0, mask_ncontrols_per_study(m, s) - result.controls_het - result.controls_alt);
        }
        result.aaf = mask_aafs(m, s);
        result.info[HTPv4_ESTIMATED_GENOTYPE_COUNTS_FLAG] = "";

        stat::tests::test_result_t burden_pval = stat::tests::chisq_df1(burden_chi/cf);
        result.pval = burden_pval.pval;
        result.info["LOG10P"] = to_string(-burden_pval.log10p);

        double se = stat::misc::get_se_from_beta_pval(mask_betas(m, s), result.pval);
        if (this->trait_type == QT) {
          result.effect = mask_betas(m, s);
          result.lci_effect = mask_betas(m, s) - 1.96 * se;
          result.uci_effect = mask_betas(m, s) + 1.96 * se;
        } else {
          result.effect = exp(mask_betas(m, s));
          result.lci_effect = exp(mask_betas(m, s) - 1.96 * se);
          result.uci_effect = exp(mask_betas(m, s) + 1.96 * se);
          result.info["BETA"] = to_string(mask_betas(m, s));
        }
        result.info["SE"] = to_string(se);
        this->burden_writer->writerec(result);
      }
    }
  }
  d = std::chrono::steady_clock::now() - start;
  log_debug("loading ld matrices: " + to_string(d.count()) + " seconds");

  if (this->ld_mat_required() && !found_at_least_one_ld_mat) {
    log_warning("no study contains an LD matrix for " + g.get_name() + ", skipping...");
    return vector<htpv4_record_t>();
  }
  VectorXd scores_sum = scores.rowwise().sum();
  VectorXd scorev_sum = scorev.rowwise().sum();
  string conditional_variant_info = "";
  if (found_at_least_one_conditional_variant) {
    log_debug("computing conditional statistics");

    buffer_ld_mat = MatrixXd(buffer_ld_mat).completeOrthogonalDecomposition().pseudoInverse().sparseView();
    scores_sum = scores_sum - gene_buffer_ld_mat * buffer_ld_mat * this->conditional_scores.rowwise().sum();
    SpMat gene_buffer_ld_mat_T = gene_buffer_ld_mat.transpose();
    gene_ld_mat = gene_ld_mat - gene_buffer_ld_mat * buffer_ld_mat * gene_buffer_ld_mat_T;
    scorev_sum = MatrixXd(gene_ld_mat).diagonal();

    for ( const variant_id& vid : this->conditional_variants ) {
      if (gene_conditional_variants.count(vid) > 0 && conditional_variant_info.size() > 0) {
        conditional_variant_info += "," + vid; 
      } else if (gene_conditional_variants.count(vid) > 0) {
        conditional_variant_info = vid;
      }
    }
  }

  log_debug("computing test statistics");
  start = std::chrono::steady_clock::now();
  VectorXd burden_scores = scores_sum.array() * burden_weights.array();
  VectorXd Qburden = (masks * burden_scores).array().square().matrix();
  VectorXd skato_scores =  scores_sum.array() * skato_weights.array();
  VectorXd Qskato_burden = (masks * skato_scores).array().square().matrix();
  VectorXd Qskato_skat = masks * skato_scores.array().square().matrix();
  vector<htpv4_record_t> burden_results;
  vector<htpv4_record_t> skato_results;
  vector<htpv4_record_t> acatv_results;
  for (int i = 0; i < nmasks; ++i) {
    string mask_cohort_meta = util::format_cohort_meta(this->cohorts, vector<int>(mask_cohort_idx[i].begin(), mask_cohort_idx[i].end()));
    string bin_name = mask_set.get_mask_alt(i);

    htpv4_record_t result_template {
      g.get_name() + "." + bin_name,
      g.get_chrom(),
      g.get_start(),
      "ref",
      bin_name,
      this->trait_name,
      mask_cohort_meta,
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
    result_template.info["MAX_MAC"] = to_string((int)masks_max_mac(i));
    result_template.info["SOURCE"] = "REMETA-GB";

    result_template.num_cases = mask_ncases[i];
    if (this->trait_type == BT) {
      result_template.num_controls = mask_ncontrols[i];
    }

    log_debug("running mask " + bin_name);
    if (mask_indices[i].size() == 0) {
      log_warning("mask " + g.get_name() + "." + bin_name + " does not have any variants");
      continue;
    }

    MatrixXd ld_mat_meta;
    double cf = 1;
    double mask_aaf = 0;
    log_debug("building meta ld matrix of size " + to_string(mask_indices[i].size()));
    vector<Triplet> triplet_list(mask_indices[i].size());
    for (size_t k = 0; k < mask_indices[i].size(); ++k) {
      triplet_list.push_back(Triplet(k, mask_indices[i][k], 1));
    }
    SpMat selector(mask_indices[i].size(), variants_found.size());
    selector.setFromTriplets(triplet_list.begin(), triplet_list.end());
    ld_mat_meta = selector * gene_ld_mat * selector.transpose();
    if ((ld_mat_meta.diagonal().array() == 0).all()) {
      // likely including imputed variants in mask, or as mismatch between variants in the LD matrix
      // and variants in the summary statistics file
      log_warning("found zero LD matrix for " + g.get_name() + "." + bin_name + ", skipping...");
      continue;
    }

    log_debug("computing burden aaf");
    double mask_aac = 0;
    double mask_an = 0;
    for (int s = 0; s < nstudies; ++s) {
      mask_aac += 2*mask_aafs(i, s)*study_sample_size[s];
      mask_an += 2*study_sample_size[s];
    }
    mask_aaf = mask_aac / mask_an;

    // this has to happen before the LD matrix is reweighted
    log_debug("computing correction factor");
    double burden_score = masks.row(i) * scores_sum;
    double burden_score_std = sqrt(ld_mat_meta.sum());
    double burden_chi = pow(burden_score / burden_score_std, 2);
    double case_control_ratio = (double)mask_ncases[i] / (double)mask_ncontrols[i];
    if (this->trait_type == BT
        && case_control_ratio < this->spa_case_control_ratio
        && -abs(burden_score) / burden_score_std < this->spa_z_score) {

      double nhet = mask_cases_het.row(i).sum() + mask_controls_het.row(i).sum();
      double nhom_alt = mask_cases_alt.row(i).sum() + mask_controls_alt.row(i).sum();
      double nhom_ref = max(0.0, mask_ncases[i] + mask_ncontrols[i] - nhet - nhom_alt);
      double spa_chi = stat::spa::compute_spa_chival_from_geno_counts(burden_score,
                                                                      burden_score_std,
                                                                      nhom_ref,
                                                                      nhet,
                                                                      nhom_alt,
                                                                      mask_ncases[i],
                                                                      mask_ncontrols[i],
                                                                      false);
      // double spa_chi = stat::spa::compute_spa_chival_from_sum_stats(burden_score,
      //                                                               burden_score_std,
      //                                                               mask_aaf,
      //                                                               mask_ncases[i],
      //                                                               mask_ncontrols[i],
      //                                                               false);
      if (spa_chi != stat::spa::SPA_FAILED) {
        cf = max(1.0, burden_chi / spa_chi);
      }
    }
    ld_mat_meta *= cf;

    log_debug("computing p-values");
    if (vector_contains(this->burden_params.af_bins, this->mask_set.get_mask_bin(i))) {
      log_debug("computing burden");
      if (Qburden[i] == 0) {
        continue;
      }
      MatrixXd burden_ld_mat = burden_weights(mask_indices[i]).asDiagonal()
                              * ld_mat_meta
                              * burden_weights(mask_indices[i]).asDiagonal();

      htpv4_record_t burden_result = result_template;
      burden_result.model = "REMETA-BURDEN-META";
      burden_result.cases_het = mask_cases_het.row(i).sum();
      burden_result.cases_alt = mask_cases_alt.row(i).sum();
      burden_result.cases_ref = max(0.0, mask_ncases[i] - burden_result.cases_het - burden_result.cases_alt);
      if (this->trait_type == BT) {
        burden_result.controls_het = mask_controls_het.row(i).sum();
        burden_result.controls_alt = mask_controls_alt.row(i).sum();
        burden_result.controls_ref = max(0.0, mask_ncontrols[i] - burden_result.controls_het - burden_result.controls_alt);
      }
      burden_result.aaf = mask_aaf;
      burden_result.info["Direction"] = "";
      burden_result.info[HTPv4_ESTIMATED_GENOTYPE_COUNTS_FLAG] = "";
      for (int k = 0; k < this->nstudies; ++k) {
        if (mask_betas(i, k) > 0) {
          burden_result.info["Direction"] += "+";
        } else if (mask_betas(i, k) < 0) {
          burden_result.info["Direction"] += "-";
        } else {
          burden_result.info["Direction"] += "?";
        }
      }
      pair<double, double> effects = stat::misc::meta_analyze_effects(
        mask_betas.row(i),
        mask_ses.row(i)
      );
      stat::tests::test_result_t burden_pval = stat::tests::wst_burden(Qburden[i], burden_ld_mat);
      burden_result.pval = burden_pval.pval;
      double se = stat::misc::get_se_from_beta_pval(effects.first, burden_result.pval);
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
      burden_result.info["NVAR"] = to_string(mask_indices[i].size());

      if (conditional_variant_info != "") {
        burden_result.info["CONDITIONAL"] = conditional_variant_info;
      }
      if (burden_pval == stat::tests::TEST_FAILED) {
        log_warning(burden_result.name + " test failed");
      } else {
        burden_results.push_back(burden_result);
      }
    }

    if (vector_contains(this->skato_params.af_bins, this->mask_set.get_mask_bin(i))) {
      log_debug("computing SKATO");
      MatrixXd skato_ld_mat = skato_weights(mask_indices[i]).asDiagonal()
                              * ld_mat_meta
                              * skato_weights(mask_indices[i]).asDiagonal();

      htpv4_record_t skat_result = result_template;
      skat_result.model = "REMETA-SKAT-META";
      htpv4_record_t skato_result = result_template;
      skato_result.model = "REMETA-SKATO-ACAT-META";

      vector<double> log10_pvals;
      for (const double& rho : this->skato_params.rho_values) {
        stat::tests::test_result_t test_result = stat::tests::skato(
          Qskato_skat[i],
          Qskato_burden[i],
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
      if (conditional_variant_info != "") {
        skato_result.info["CONDITIONAL"] = conditional_variant_info;
      }
      if (skato_acat == stat::tests::TEST_FAILED || log10_pvals.size() == 0) {
        log_warning("SKATO test failed for " + skato_result.name);
      } else {
        skato_results.push_back(skato_result);
      }

      if (skat_result.pval == HTPv4_NA) {
        stat::tests::test_result_t test_result = stat::tests::skato(Qskato_skat[i], Qskato_burden[i], 0, ld_mat_meta);
        skat_result.pval = test_result.pval;
        skat_result.info["LOG10P"] = to_string(-test_result.log10p);
        skat_result.info["CF"] = to_string(cf);
      }

      if (conditional_variant_info != "") {
        skat_result.info["CONDITIONAL"] = conditional_variant_info;
      }
      if (skat_result.pval < 0) {
        log_warning("SKAT test failed for " + skat_result.name);
      } else {
        skato_results.push_back(skat_result);
      }
    }

    if (vector_contains(this->acatv_params.af_bins, this->mask_set.get_mask_bin(i))) {
      log_debug("computing ACATV");
      htpv4_record_t acatv_result = result_template;
      acatv_result.model = "REMETA-ACATV-META";

      vector<double> mask_log10_pvals;
      vector<double> mask_weights;
      int nvar = 0;
      double z = 0;

      VectorXd Qacatv = masks.row(i).transpose().array() * scores_sum.array();
      VectorXd weights = masks.row(i).transpose().array() * acatv_weights.array();
      for (ssize_t v = 0; v < Qacatv.size(); ++v) {
        if (weights[v] > 0 && scorev_sum(v) > 0) {
          z = -abs(Qacatv[v]) / sqrt(scorev_sum(v));
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

        if (conditional_variant_info != "") {
          acatv_result.info["CONDITIONAL"] = conditional_variant_info;
        }
        acatv_results.push_back(acatv_result);
      }
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
  return results;
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

void HTPMetaAnalyzer::check_info(const htpv4_record_t& rec) {
  if (rec.info.find("SCORE") == rec.info.end() && !this->recompute_score) {
    throw runtime_error("missing SCORE from INFO field of: " + rec.name);
  }
}

double HTPMetaAnalyzer::get_score(const htpv4_record_t& rec) {
  if (this->recompute_score && rec.info.count("SCORE") == 0) {
    double beta = HTPv4Reader::get_beta(rec);
    double se = HTPv4Reader::get_se(rec);
    return beta / (se * se);
  } else {
    return stod(rec.info.at("SCORE"));
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
