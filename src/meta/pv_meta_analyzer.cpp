#include "pv_meta_analyzer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
using namespace std;

#include "../stat/tests.hpp"
#include "../logging.hpp"

PVMetaAnalyzer::PVMetaAnalyzer(pvma_method_e method,
                               trait_type_e trait_type,
                               string trait_name,
                               vector<string> cohorts,
                               bool unweighted) 
  : method(method)
  , trait_type(trait_type)
  , trait_name(trait_name)
  , cohorts(cohorts)
  , unweighted(unweighted) { }

void PVMetaAnalyzer::add_line(const htpv4_record_t& rec, const int& study_index) {
  log_debug("entering PVMetaAnalyzer::add_line");
  log_debug("processing " + rec.name + " from study " + to_string(study_index));

  htpv4_key_t key = rec.name + "." + rec.model;
  if (this->keys.find(key) == this->keys.end()) {
    this->keys.insert(key);
    this->key_order.push_back(key);
  }
  this->htpv4_recs[key].push_back(rec);
  this->htpv4_cohort_idx[key].push_back(study_index);
  log_debug("leaving PVMetaAnalyzer::add_line");
}

string pvma_method_to_string(const pvma_method_e& method) {
  if (method == STOUFFERS) {
    return "stouffers";
  } else if (method == FISHERS) {
    return "fishers";
  } else {
    return "";
  }
}

vector<htpv4_record_t> PVMetaAnalyzer::meta_analyze_before(const string& chr,
                                                           const int& pos) {
  log_debug("entering PVMetaAnalyzer::meta_analyze_before");
  log_debug("meta-analyzing records on chr " + chr + " before position " + to_string(pos));
  vector<htpv4_record_t> results;
  for (htpv4_key_t& key : this->key_order) {
    string rec_chr = this->htpv4_recs[key][0].chr;
    int rec_pos = this->htpv4_recs[key][0].pos;
    if (rec_chr != chr || rec_pos >= pos) {
      continue;
    }

    vector<double> log10_pvals;
    vector<double> sample_sizes;
    vector<int> cohort_idx;
    double sample_sizes_sum = 0.0;
    int num_cases = 0;
    int cases_ref = 0;
    int cases_het = 0;
    int cases_alt = 0;
    int num_controls = 0;
    int controls_ref = 0;
    int controls_het = 0;
    int controls_alt = 0;
    rec_chr = this->htpv4_recs[key][0].chr;
    rec_pos = this->htpv4_recs[key][0].pos;
    string source_list = "";
    string source = "";
    for (htpv4_record_t rec : this->htpv4_recs[key]) {
      if (rec.chr != rec_chr || rec.pos != rec_pos) {
        string expected = rec_chr + ":" + to_string(rec_pos);
        string found = rec.chr + ":" + to_string(rec.pos);
        log_warning("found inconsistent position for " + key + " (found " + found + " but expected " + expected + ")");
      }

      if (rec.pval >= std::numeric_limits<double>::min()) {
        log10_pvals.push_back(log10(rec.pval));
      } else if (rec.info.count("LOG10P")) {
        log10_pvals.push_back(-stod(rec.info.at("LOG10P")));
      } else if (rec.pval != HTPv4_NA) {
        log10_pvals.push_back(log10(std::numeric_limits<double>::min()));
      } else {
        log_warning("found invalid p-value for " + key);
        continue;
      }

      num_cases += rec.num_cases;
      cases_ref += rec.cases_ref;
      cases_het += rec.cases_het;
      cases_alt += rec.cases_alt;
      num_controls += rec.num_controls;
      if (this->is_bt()) {
        controls_ref += rec.controls_ref;
        controls_het += rec.controls_het;
        controls_alt += rec.controls_alt;
      }

      if (this->unweighted) {
        sample_sizes.push_back(1.0);
        sample_sizes_sum += 1.0;
      } else if (this->is_qt()) {
        sample_sizes.push_back(rec.num_cases);
        sample_sizes_sum += rec.num_cases;
      } else if (this->is_bt()) {
        sample_sizes.push_back(
          4.0 / ((1.0/rec.num_cases + 1.0/rec.num_controls))
        );
        sample_sizes_sum += (
          4.0 / ((1.0/rec.num_cases + 1.0/rec.num_controls))
        );
      }

      if (rec.info.count("SOURCE") > 0) {
        if (source == "") {
          source = rec.info.at("SOURCE");
        } else if (source != "" && source != rec.info.at("SOURCE")) {
          source = "MULTI";
        }

        if (source_list != "") {
          source_list += ",";
        }
        source_list += source_map.get_source_short(rec.info.at("SOURCE"));
      }
    }

    stat::tests::test_result_t test_result;
    if (this->method == STOUFFERS) {
      vector<double> weights;
      for (const double& s : sample_sizes) {
        weights.push_back(sqrt(s));
      }
      test_result = stat::tests::stouffers(log10_pvals, weights);
    } else if (this->method == FISHERS) {
      test_result = stat::tests::fishers(log10_pvals);
    } else {
      throw runtime_error("invalid pvma method");
    }

    htpv4_record_t rec0 = this->htpv4_recs[key][0];
    results.push_back(
      htpv4_record_t {
        rec0.name,
        rec0.chr,
        rec0.pos,
        rec0.ref,
        rec0.alt,
        this->trait_name,
        util::format_cohort_meta(this->cohorts, this->htpv4_cohort_idx[key]),
        rec0.model + "-META",
        HTPv4_NA,
        HTPv4_NA,
        HTPv4_NA,
        test_result.pval,
        HTPv4_NA,
        num_cases,
        HTPv4_NA,
        HTPv4_NA,
        HTPv4_NA,
        (this->is_bt() ? num_controls : HTPv4_NA),
        HTPv4_NA,
        HTPv4_NA,
        HTPv4_NA,
        map<string, string>()
      }
    );
    int last_idx = results.size() - 1;
    results[last_idx].info["LOG10P"] = to_string(-test_result.log10p);
    results[last_idx].info["META"] = "PVMA_" + pvma_method_to_string(this->method);

    if (source != "") {
      results[last_idx].info["SOURCE"] = source;
      results[last_idx].info["SOURCE_LIST"] = source_list;
    }
  }
  log_debug("leaving PVMetaAnalyzer::meta_analyze_before");
  return results;
}

void PVMetaAnalyzer::clear_before(const string& chr, const int& pos) {
  log_debug("entering PVMetaAnalyzer::clear_before");
  log_debug("clearing records chr " + chr + " before position " + to_string(pos));
  vector<htpv4_key_t> new_key_order;
  for (htpv4_key_t& key : this->key_order) {
    string rec_chr = this->htpv4_recs[key][0].chr;
    int rec_pos = this->htpv4_recs[key][0].pos;
    if (rec_chr == chr && rec_pos < pos) {
      this->keys.erase(key);
      this->htpv4_recs.erase(key);
    } else {
      new_key_order.push_back(key);
    }
  }
  std::swap(new_key_order, this->key_order);
  log_debug("leaving PVMetaAnalyzer::clear_before");
}