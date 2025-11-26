#ifndef ES_META_ANALYZER
#define ES_META_ANALYZER

#include <string>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

#include <boost/math/distributions/normal.hpp>

#include "../io/htpv4_reader.hpp"
#include "../enums.hpp"
#include "../source_map.hpp"
#include "variant_meta_analyzer.hpp"

// The Name+Model columns in an HTPv4 file should be unique, so
// let's use these as a unique key.
typedef string htpv4_key_t;

struct es_variant_t {
  htpv4_key_t key;
  vector<double> betas;
  vector<double> standard_errors;
  vector<double> aafs;
  vector<double> sample_sizes;
  vector<double> info_scores;
  vector<int> cohort_idx;
  set<string> models;
  string effect_directions;
  vector<string> sources;
  htpv4_record_t meta_result;
};

class ESMetaAnalyzer : public VariantMetaAnalyzer {
 public:
  ESMetaAnalyzer(const trait_type_e& trait_type,
                 const string& trait_name,
                 const vector<string>& cohorts,
                 int nstudies);

  void add_line(const htpv4_record_t& htpv4_rec, const int& study_index);

  // meta-analyze all records that occur on chromosome chr and < before_pos
  vector<htpv4_record_t> meta_analyze_before(const string& chr, const int& before_pos);

  // clear all records that occur on chromosome chr and < before_pos
  void clear_before(const string& chr, const int& before_pos);

  bool is_bt() { return this->trait_type == BT; }

  bool is_qt() { return this->trait_type == QT; }

  // include per study aafs in info field
  void set_include_study_aafs() { this->include_study_aafs = true; }

  void set_source_definitions(string source_def) { this->source_map.load(source_def); }

 private:
  trait_type_e trait_type;
  string trait_name;
  vector<string> cohorts;
  int nstudies;
  bool include_study_aafs;
  SourceMap source_map;

  boost::math::normal s;

  unordered_set<htpv4_key_t> keys;
  vector<htpv4_key_t> key_order;
  unordered_map<htpv4_key_t, es_variant_t> variants;
};

#endif