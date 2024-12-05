#ifndef REMETA_PV_META_ANALYZER_H
#define REMETA_PV_META_ANALYZER_H

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

#include "variant_meta_analyzer.hpp"
#include "../io/htpv4_reader.hpp"
#include "../util.hpp"
#include "../enums.hpp"
#include "../source_map.hpp"

// The Name+Model columns in an HTPv4 file should be unique, so
// let's use these as a unique key.
typedef string htpv4_key_t;

class PVMetaAnalyzer : public VariantMetaAnalyzer {
  public:
    PVMetaAnalyzer(pvma_method_e method,
                   trait_type_e trait_type,
                   string trait_name,
                   vector<string> cohorts,
                   bool unweighted);

    void add_line(const htpv4_record_t& htpv4_rec, const int& study_index);

    // meta-analyze all records that occur on chromosome chr and <= before_pos
    vector<htpv4_record_t> meta_analyze_before(const string& chr, 
                                               const int& before_pos);

    // clear all records that occur on chromosome chr and <= before_pos
    void clear_before(const string& chr, const int& before_pos);

    bool is_bt() { return this->trait_type == BT; }

    bool is_qt() { return this->trait_type == QT; }

    void set_source_definitions(string source_def) { this->source_map.load(source_def); }

  private:
    pvma_method_e method;
    trait_type_e trait_type;
    string trait_name;
    vector<string> cohorts;
    bool unweighted;
    SourceMap source_map;

    unordered_set<htpv4_key_t> keys;
    vector<htpv4_key_t> key_order;
    unordered_map<htpv4_key_t, vector<htpv4_record_t> > htpv4_recs;
    unordered_map<htpv4_key_t, vector<int> > htpv4_cohort_idx;
};

#endif