#ifndef SET_META_ANALYZER_H
#define SET_META_ANALYZER_H

#include <string>
#include <vector>
using namespace std;

#include "../io/gene_set_reader.hpp"
#include "../io/htpv4_reader.hpp"
#include "../io/htpv4_writer.hpp"
#include "../variant_filter.hpp"

class SetMetaAnalyzer {
 public:
  virtual void add_line(const htpv4_record_t& htpv4_rec,
                        const int& study_index) = 0;

  virtual vector<htpv4_record_t> meta_analyze_gene(Gene g) = 0;

  virtual void clear_before(const string& chrom,
                            const int& before_pos) = 0;
};

void run_set_meta_analysis(SetMetaAnalyzer& meta,
                           VariantFilter& vf,
                           vector<HTPv4Reader>& in_files,
                           const string& set_list_file,
                           HTPv4Writer& out,
                           const string& chr,
                           const string& gene);

#endif