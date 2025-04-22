#include "es_meta_analyzer.hpp"

#include <algorithm>
#include <boost/format.hpp>
#include <cmath>
#include <set>
#include <stdexcept>
#include <iostream>
using namespace std;

#include "../io/htpv4_reader.hpp"
#include "../stat/tests.hpp"
#include "../logging.hpp"
#include "../util.hpp"

ESMetaAnalyzer::ESMetaAnalyzer(const trait_type_e& trait_type,
                               const string& trait_name,
                               const vector<string>& cohorts,
                               const int& nstudies)
 : trait_type(trait_type)
 , trait_name(trait_name)
 , cohorts(cohorts)
 , nstudies(nstudies)
 , include_study_aafs(false) {}

void ESMetaAnalyzer::add_line(const htpv4_record_t& rec, 
                              const int& study_index) {
  log_debug("entering ESMetaAnalyzer::add_line");
  log_debug("processing " + rec.name + " from study " + to_string(study_index));
  if (!HTPv4Reader::has_beta(rec)) {
    throw runtime_error("REGENIE_SE, META_SE or, EXTERNAL_GWAS_SE, or SE info field is required for ESMA: " + rec.name);
  }

  htpv4_key_t key;
  if (util::is_cpra(rec.name)) {
    key = rec.name;
  } else {
    key = rec.name + "." + rec.model;
  }
  if (this->keys.find(key) == this->keys.end()) {
    this->keys.insert(key);
    this->key_order.push_back(key);

    es_variant_t variant;
    variant.key = key;
    variant.betas = vector<double>(this->nstudies, 0);
    variant.standard_errors = vector<double>(this->nstudies, 0);
    variant.effect_directions = string(this->nstudies, '?');
    variant.info_scores = vector<double>(this->nstudies, HTPv4_NA);
    variant.aafs = vector<double>(this->nstudies, 0);
    variant.sample_sizes = vector<double>(this->nstudies, 0);
    variant.sources = vector<string>(this->nstudies, "");

    variant.meta_result = {
      rec.name,
      rec.chr,
      rec.pos,
      rec.ref,
      rec.alt,
      this->trait_name,
      "",
      "",
      HTPv4_NA,
      HTPv4_NA,
      HTPv4_NA,
      HTPv4_NA,
      HTPv4_NA,
      0,
      0,
      0,
      0,
      HTPv4_NA,
      HTPv4_NA,
      HTPv4_NA,
      HTPv4_NA,
      map<string, string>()
    };
  
    if (this->is_bt()) {
      variant.meta_result.num_controls = 0;
      variant.meta_result.controls_ref = 0;
      variant.meta_result.controls_het = 0;
      variant.meta_result.controls_alt = 0;
    } else {
      variant.meta_result.num_controls = HTPv4_NA;
      variant.meta_result.controls_ref = HTPv4_NA;
      variant.meta_result.controls_het = HTPv4_NA;
      variant.meta_result.controls_alt = HTPv4_NA;    
    }

    this->variants[key] = variant;
  }

  es_variant_t *variant = &this->variants[key];
  variant->betas[study_index] = HTPv4Reader::get_beta(rec);
  variant->standard_errors[study_index] = HTPv4Reader::get_se(rec);
  variant->effect_directions[study_index] = variant->betas[study_index] > 0 ? '+' : '-';
  variant->aafs[study_index] = rec.aaf;
  if (rec.info.count("SOURCE") != 0) {
    variant->sources[study_index] = rec.info.at("SOURCE");
  }
  if (rec.info.count("INFO") != 0) {
    variant->info_scores[study_index] = stod(rec.info.at("INFO"));
  }
  variant->cohort_idx.push_back(study_index);
  variant->models.insert(rec.model);

  variant->sample_sizes[study_index] = rec.num_cases;
  variant->meta_result.num_cases += rec.num_cases;
  variant->meta_result.cases_ref += rec.cases_ref;
  variant->meta_result.cases_het += rec.cases_het;
  variant->meta_result.cases_alt += rec.cases_alt;
  if (this->is_bt()) {
    variant->sample_sizes[study_index] += rec.num_controls;
    variant->meta_result.num_controls += rec.num_controls;
    variant->meta_result.controls_ref += rec.controls_ref;
    variant->meta_result.controls_het += rec.controls_het;
    variant->meta_result.controls_alt += rec.controls_alt;
  }
  log_debug("leaving ESMetaAnalyzer::add_line");
}

vector<htpv4_record_t> ESMetaAnalyzer::meta_analyze_before(const string& chr,
                                                           const int& pos) {
  log_debug("entering ESMetaAnalyzer::meta_analyzer_before");
  log_debug("meta-analyzing records on chr " + chr + " before position " + to_string(pos));
  vector<htpv4_record_t> results;
  for (htpv4_key_t& key : this->key_order) {

    string rec_chr = this->variants[key].meta_result.chr;
    int rec_pos = this->variants[key].meta_result.pos;
    if (rec_chr != chr || rec_pos >= pos) {
      continue;
    }

    es_variant_t v = this->variants[key];

    double aac = 0;
    double an = 0;
    for (int i = 0; i < this->nstudies; ++i) {
      aac += 2*v.sample_sizes[i]*v.aafs[i];
      an += 2*v.sample_sizes[i];
    }
    v.meta_result.aaf = aac / an;

    vector<double> weights;
    double weight_sum = 0;
    double beta = 0;
    double weight = 0;
    double sse = 0;
    for (int i = 0; i < this->nstudies; ++i) {
      if (v.standard_errors[i] != 0) {
        sse = v.standard_errors[i];
        weight = 1 / pow(sse, 2);
        weights.push_back(weight);
        weight_sum += weight;
        beta += weight * v.betas[i];      
      }
    }

    beta /= weight_sum;
    double se = sqrt(1 / weight_sum);
    double z = beta / se;

    if (boost::math::isnan(z)) {
      log_warning("meta-analysis failed for: " + v.meta_result.name);
      v.meta_result.pval = HTPv4_NA;
    } else {
      stat::tests::test_result_t test_result = stat::tests::normal(-abs(z), true);
      test_result.pval = min(test_result.pval, 1 - 1e-6);
      v.meta_result.pval = test_result.pval;
      v.meta_result.info["LOG10P"] = to_string(-test_result.log10p);
    }

    if (this->is_bt()) {
      v.meta_result.lci_effect = exp(beta - 1.96 * se);
      v.meta_result.uci_effect = exp(beta + 1.96 * se);
      v.meta_result.effect = exp(beta);
      v.meta_result.info["BETA"] = to_string(beta);
    } else {
      v.meta_result.lci_effect = beta - 1.96 * se;
      v.meta_result.uci_effect = beta + 1.96 * se;
      v.meta_result.effect = beta;
    }
    v.meta_result.cohort = util::format_cohort_meta(this->cohorts, v.cohort_idx);
    v.meta_result.info["META"] = "ESMA";
    v.meta_result.info["META_SE"] = to_string(se);
    v.meta_result.info["Direction"] = v.effect_directions;

    v.meta_result.model = "";
    for (const string& m : v.models) {
      if (v.meta_result.model != "") {
        v.meta_result.model  += "__";
      }
      v.meta_result.model += m;
    }
    v.meta_result.model += "-META";

    if (this->include_study_aafs) {
      string aafs = "";
      for (size_t i = 0; i < v.aafs.size(); ++i) {
        if (aafs != "") {
          aafs += ",";
        }
        if (v.aafs[i] != 0) {
          aafs += to_string(v.aafs[i]);
        }
      }
      v.meta_result.info["AAF"] = aafs;
    }

    string source_list = "";
    string source = "";
    string category = "";
    bool at_least_one_variant_has_source = false;
    for (int i = 0; i < this->nstudies; ++i) {
      if (v.sources[i] != "") {
        at_least_one_variant_has_source = true;
        if (source_list != "") {
          source_list += "," ;
        }
        source_list += source_map.get_source_short(v.sources[i]);

        if (source == "") {
          source = v.sources[i];
        } else if (source != v.sources[i]) {
          source = "MULTI";
        }

        if (category == "") {
          category = source_map.get_category(v.sources[i]);
        } else if (category != source_map.get_category(v.sources[i])) {
          category = "MULTI";
        }
      }
    }
    if (source == "MULTI" && category != "MULTI" && category != "") {
      source = "MULTI-" + category;
    }

    if (at_least_one_variant_has_source) {
      v.meta_result.info["SOURCE"] = source;
      v.meta_result.info["SOURCE_LIST"] = source_list;
    }

    string info = "";
    bool at_least_one_variant_has_info = false;
    double min_info = 2;
    boost::format fmter("%.3f");
    for (int i = 0; i < this->nstudies; ++i) {
      if (v.info_scores[i] != HTPv4_NA && v.info_scores[i] < min_info) {
        min_info = v.info_scores[i];
      }

      if (v.info_scores[i] != HTPv4_NA) {
        at_least_one_variant_has_info = true;
        if (info != "") {
          info += ",";
        }
        info += str(fmter % v.info_scores[i]);
      }
    }
    if (at_least_one_variant_has_info) {
      v.meta_result.info["INFO_LIST"] = info;
      v.meta_result.info["INFO"] = str(fmter % min_info);
    }

    results.push_back(v.meta_result);
  }
  log_debug("leaving ESMetaAnalyzer::meta_analyzer_before");
  return results;
}

void ESMetaAnalyzer::clear_before(const string& chr, const int& pos) {
  log_debug("entering ESMetaAnalyzer::clear_before");
  log_debug("clearing records on chr " + chr + " before position " + to_string(pos));
  vector<htpv4_key_t> new_key_order;
  for (htpv4_key_t& key : this->key_order) {
    string rec_chr = this->variants[key].meta_result.chr;
    int rec_pos = this->variants[key].meta_result.pos;
    if (rec_chr == chr && rec_pos < pos) {
      this->keys.erase(key);
      this->variants.erase(key);
    } else {
      new_key_order.push_back(key);
    }
  }
  std::swap(new_key_order, this->key_order);
  log_debug("leaving ESMetaAnalyzer::clear_before");
}