#ifndef GENEP_META_ANALYZER_H
#define GENEP_META_ANALYZER_H

#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

#include "../enums.hpp"
#include "../io/htpv4_reader.hpp"
#include "variant_meta_analyzer.hpp"

typedef string gene_id;

struct genep_record_t {
  string name;
  string chr;
  int pos;
  string trait;
  string cohort;

  string acatv_source;
  string burden_source;
  string skato_source;
  string sbat_source;

  // one vector per group in genep file
  vector<vector<double> > acatv_log10_pvals;
  vector<vector<double> > burden_log10_pvals;
  vector<vector<double> > skato_log10_pvals;
  vector<vector<double> > sbat_log10_pvals;
};

enum model_e {
  BURDEN,
  ACATV,
  SKATO,
  SBAT,
  NONE
};

class GenePMetaAnalyzer : public VariantMetaAnalyzer {
 public:
  GenePMetaAnalyzer(const string& genep_file,
                    const string& burden_model,
                    const string& acatv_model,
                    const string& skato_model,
                    const string& sbat_model);

  void add_line(const htpv4_record_t& rec, const int& study_index);

  vector<htpv4_record_t> meta_analyze_before(const string& chr, const int& before_pos);

  void clear_before(const string& chr, const int& before_pos);

 private:
  string burden_model;
  string acatv_model;
  string skato_model;
  string sbat_model;

  // mask -> index in genep_groups
  map<string, vector<int> > mask_map;
  vector<string> genep_groups;
  unordered_map<gene_id, genep_record_t> genep_records;
  vector<htpv4_record_t> htpv4_records;

  void load_genep_file(const string& genep_file);

  model_e parse_model(const htpv4_record_t& rec);

  pair<htpv4_record_t, double> meta_analyze_genep_group(const genep_record_t& gene,
                                                        const string& model,
                                                        const string& genep_group,
                                                        const string& source,
                                                        const vector<double>& log10_pvals);
};

#endif