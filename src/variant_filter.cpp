#include "variant_filter.hpp"

#include "io/bgz_reader.hpp"

VariantFilter::VariantFilter()
 : extract_file("")
 , exclude_file("")
 , sources({})
 , effect_not_na(false)
 , effect_is_na(false)
 , info_has_se(false)
 , info_has_score(false)
 , skip_genep(false)
 , keep_remeta_gene_only(false)
 , keep_remeta_gene_and_sbat(false)
 , extract_ids_ptr(new unordered_set<string>)
 , exclude_ids_ptr(new unordered_set<string>) {}

void VariantFilter::set_extract_file(const string& extract_file) {
  this->extract_file = extract_file;
  BgzReader reader(extract_file);
  while (!reader.eof()) {
    extract_ids_ptr->insert(reader.readline());
  }
}

void VariantFilter::set_exclude_file(const string& exclude_file) {
  this->exclude_file = exclude_file;
  BgzReader reader(exclude_file);
  while (!reader.eof()) {
    exclude_ids_ptr->insert(reader.readline());
  }
}

void VariantFilter::set_info_source_is_one_of(const vector<string>& sources) {
  this->sources.clear();
  for (const string& s : sources) {
    this->sources.insert(s);
  }
}

void VariantFilter::set_effect_not_na() {
  if (this->effect_is_na) {
    throw runtime_error("VariantFilter: cannot set both effect not NA and effect is NA");
  }
  this->effect_not_na = true;
}

void VariantFilter::set_effect_is_na() {
  if (this->effect_not_na) {
    throw runtime_error("VariantFilter: cannot set both effect not NA and effect is NA");
  }
  this->effect_is_na = true;
}

void VariantFilter::set_info_has_se() {
  this->info_has_se = true;
}

void VariantFilter::set_info_has_score() {
  this->info_has_score = true;
}

void VariantFilter::set_skip_genep() {
  this->skip_genep = true;
}

void VariantFilter::set_keep_remeta_gene_only() {
  this->keep_remeta_gene_only = true;
}

void VariantFilter::set_keep_remeta_gene_and_sbat() {
  this->keep_remeta_gene_and_sbat = true;
}

bool VariantFilter::include_htp_record(const htpv4_record_t& rec) {
  if (this->sources.size() > 0 && (rec.info.count("SOURCE") == 0  || this->sources.count(rec.info.at("SOURCE")) == 0)) {
    return false;
  } else if (this->effect_not_na && !HTPv4Reader::has_beta(rec)) {
    return false;
  } else if (this->effect_is_na && rec.effect != HTPv4_NA) {
    return false;
  } else if (this->info_has_se && !HTPv4Reader::has_se(rec)) {
    return false;
  } else if (this->info_has_score && rec.info.count("SCORE") == 0) {
    return false;
  } else if (this->skip_genep && rec.name.find(".GENE") == rec.name.length() - 5 && rec.model.find("SBAT") == string::npos) {
    return false;
  } else if (keep_remeta_gene_only) {
    return rec.name.find(".GENE") == string::npos || rec.model.find("REMETA") != string::npos;
  } else if (keep_remeta_gene_and_sbat) {
    return rec.name.find(".GENE") == string::npos 
           || rec.model.find("REMETA") != string::npos
           || (
                (rec.model.rfind("ADD-WGR-BURDEN-SBAT_POS", 0) == 0 || rec.model.rfind("ADD-WGR-BURDEN-SBAT_NEG", 0) == 0) && rec.alt == "set"
              ); 
  } else if (this->extract_file_is_set()) {
    return this->extract_ids_ptr->count(rec.name) > 0;
  } else if (this->exclude_file_is_set()) {
    return this->exclude_ids_ptr->count(rec.name) == 0;
  } else {
    return true;
  }
}

bool VariantFilter::exclude_htp_record(const htpv4_record_t& rec) {
  return !this->include_htp_record(rec);
}

bool VariantFilter::include_variant(const variant_id& vid) {
  if (this->extract_file_is_set()) {
    return this->extract_ids_ptr->count(vid) > 0;
  } else if (this->exclude_file_is_set()) {
    return this->exclude_ids_ptr->count(vid) == 0;
  } else {
    return true;
  }
}

bool VariantFilter::exclude_variant(const variant_id& vid) {
  return !this->include_variant(vid);
}