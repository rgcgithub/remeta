#ifndef VARIANT_META_ANALYZER
#define VARIANT_META_ANALYZER

#include <string>
#include <vector>
using namespace std;

#include "../io/htpv4_reader.hpp"
#include "../io/htpv4_writer.hpp"
#include "../variant_filter.hpp"

class VariantMetaAnalyzer {
 public:
  virtual void add_line(const htpv4_record_t& htpv4_rec,
                        const int& study_index) = 0;
  
  // meta-analyze all records that occur on chromosome chr and < before_pos
  virtual vector<htpv4_record_t> meta_analyze_before(const string& chr,
                                                     const int& before_pos) = 0;
  
  // clear all records that occur on chromosome chr and < before_pos
  virtual void clear_before(const string& chr, const int& before_pos) = 0;
};

void run_variant_meta_analysis(VariantMetaAnalyzer& meta,
                               VariantFilter& vf,
                               vector<HTPv4Reader>& in_files,
                               HTPv4Writer& out,
                               const string& chr);

#endif