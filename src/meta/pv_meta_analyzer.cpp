#include "pv_meta_analyzer.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
using namespace std;

#include "../stat/tests.hpp"
#include "../logging.hpp"

PVMetaAnalyzer::PVMetaAnalyzer(pvma_method_e method,
                               trait_type_e trait_type,
                               string trait_name,
                               vector<string> cohorts,
                               bool unweighted,
                               bool two_sided) 
  : method(method)
  , trait_type(trait_type)
  , trait_name(trait_name)
  , cohorts(cohorts)
  , unweighted(unweighted)
  , two_sided(two_sided) { }

void PVMetaAnalyzer::add_line(const htpv4_record_t& rec, const int& study_index) {
  log_debug("entering PVMetaAnalyzer::add_line");
  log_debug("processing " + rec.name + " from study " + to_string(study_index));

  htpv4_key_t key;
  if (util::is_cpra(rec.name)) {
    key = rec.name;
  } else {
    key = rec.name + "." + rec.model;
  }
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
    vector<double> signs;
    vector<int> cohort_idx;
    set<string> models;
    double sample_sizes_sum = 0.0;
    int num_cases = 0;
    int cases_ref = 0;
    int cases_het = 0;
    int cases_alt = 0;
    int num_controls = 0;
    int controls_ref = 0;
    int controls_het = 0;
    int controls_alt = 0;
    double aac = 0;
    double an = 0;
    bool has_genotype_counts = true;
    rec_chr = this->htpv4_recs[key][0].chr;
    rec_pos = this->htpv4_recs[key][0].pos;
    string source_list = "";
    string source = "";
    string category = "";
    for (htpv4_record_t rec : this->htpv4_recs[key]) {
      if (rec.chr != rec_chr || rec.pos != rec_pos) {
        string expected = rec_chr + ":" + to_string(rec_pos);
        string found = rec.chr + ":" + to_string(rec.pos);
        log_warning("found inconsistent position for " + key + " (found " + found + " but expected " + expected + ")");
      }

      if (rec.pval >= std::numeric_limits<double>::min()) {
        log10_pvals.push_back(log10(rec.pval));
      } else if (rec.pval != HTPv4_NA && rec.info.count("LOG10P") && rec.info.at("LOG10P") != "NA") {
        log10_pvals.push_back(-stod(rec.info.at("LOG10P")));
      } else if (rec.pval != HTPv4_NA) {
        log10_pvals.push_back(log10(std::numeric_limits<double>::min()));
      } else {
        log_warning("found invalid p-value for " + key);
        continue;
      }

      if (this->two_sided && HTPv4Reader::has_beta(rec)) {
        signs.push_back(HTPv4Reader::get_beta(rec) > 0 ? 1.0 : -1.0);
      }
      models.insert(rec.model);

      bool has_genotype_counts = rec.cases_ref != HTPv4_NA && has_genotype_counts;
      num_cases += rec.num_cases;
      cases_ref += rec.cases_ref != HTPv4_NA ? rec.cases_ref : 0;
      cases_het += rec.cases_het != HTPv4_NA ? rec.cases_het : 0;
      cases_alt += rec.cases_alt != HTPv4_NA ? rec.cases_alt : 0;
      num_controls += rec.num_controls;
      if (this->is_bt()) {
        controls_ref += rec.controls_ref != HTPv4_NA ? rec.controls_ref : 0;
        controls_het += rec.controls_het != HTPv4_NA ? rec.controls_het : 0;
        controls_alt += rec.controls_alt != HTPv4_NA ? rec.controls_alt : 0;
      }
      if (rec.aaf != HTPv4_NA) {
        double ss = this->is_bt() ? rec.num_cases + rec.num_controls : rec.num_cases;
        aac += 2*ss*rec.aaf;
        an += 2*ss;
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

        if (category == "") {
          category = source_map.get_category(rec.info.at("SOURCE"));
        } else if (category != source_map.get_category(rec.info.at("SOURCE"))) {
          category = "MULTI";
        }
      }
    }

    if (source == "MULTI" && category != "MULTI" && category != "") {
      source = "MULTI-" + category;
    }

    if (log10_pvals.size() == 0) {
      log_warning("no valid p-values found for " + key);
      continue;
    }

    stat::tests::test_result_t test_result;
    if (this->method == STOUFFERS) {
      vector<double> weights;
      double weight_sum = 0;
      for (const double& s : sample_sizes) {
        weights.push_back(sqrt(s));
        weight_sum += s;
      }

      bool all_var_have_sign = log10_pvals.size() == signs.size();
      if (all_var_have_sign) {
        test_result = stat::tests::stouffers_two_sided(log10_pvals, weights, signs);
      } else {
        test_result = stat::tests::stouffers(log10_pvals, weights);
      }

    } else if (this->method == FISHERS) {
      test_result = stat::tests::fishers(log10_pvals);
    } else {
      throw runtime_error("invalid pvma method");
    }

    htpv4_record_t rec0 = this->htpv4_recs[key][0];
    string model = "";
    for (const string& m : models) {
      if (model != "") {
        model  += "__";
      }
      model += m;
    }
    model += "-META";

    results.push_back(
      htpv4_record_t {
        rec0.name,
        rec0.chr,
        rec0.pos,
        rec0.ref,
        rec0.alt,
        this->trait_name,
        util::format_cohort_meta(this->cohorts, this->htpv4_cohort_idx[key]),
        model,
        HTPv4_NA,
        HTPv4_NA,
        HTPv4_NA,
        test_result.pval,
        has_genotype_counts ? aac / an : HTPv4_NA,
        num_cases,
        has_genotype_counts ? static_cast<double>(cases_ref) : HTPv4_NA,
        has_genotype_counts ? static_cast<double>(cases_het) : HTPv4_NA,
        has_genotype_counts ? static_cast<double>(cases_alt) : HTPv4_NA,
        (this->is_bt() ? num_controls : HTPv4_NA),
        this->is_bt() && has_genotype_counts ? static_cast<double>(controls_ref) : HTPv4_NA,
        this->is_bt() && has_genotype_counts ? static_cast<double>(controls_het) : HTPv4_NA,
        this->is_bt() && has_genotype_counts ? static_cast<double>(controls_alt) : HTPv4_NA,
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