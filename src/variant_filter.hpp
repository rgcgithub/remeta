#ifndef VARIANT_FILTER_H
#define VARIANT_FILTER_H

#include <memory>
#include <string>
#include <unordered_set>
using namespace std;

#include "io/htpv4_reader.hpp"

typedef string variant_id;

class VariantFilter {
 public:
  VariantFilter();

  void set_extract_file(const string& extract_file);

  void set_exclude_file(const string& exclude_file);

  void set_info_source_is_one_of(const vector<string>& sources);

  void set_effect_not_na();

  void set_effect_is_na();

  void set_info_has_se();

  void set_info_has_score();

  void set_skip_genep();

  void set_keep_remeta_gene_only();

  void set_keep_remeta_gene_and_sbat();

  bool include_htp_record(const htpv4_record_t& rec);

  bool exclude_htp_record(const htpv4_record_t& rec);

  bool include_variant(const variant_id& vid);

  bool exclude_variant(const variant_id& vid);

  bool extract_file_is_set() { return this->extract_file != ""; }

  bool exclude_file_is_set() { return this->exclude_file != ""; }

 private:
  string extract_file;
  string exclude_file;
  unordered_set<string> sources;
  bool effect_not_na;
  bool effect_is_na;
  bool info_has_se;
  bool info_has_score;
  bool skip_genep;
  bool keep_remeta_gene_only;
  bool keep_remeta_gene_and_sbat;
  // avoid an expensive copy
  shared_ptr<unordered_set<variant_id> > extract_ids_ptr;
  shared_ptr<unordered_set<variant_id> > exclude_ids_ptr;
};

#endif