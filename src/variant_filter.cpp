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
 , exclude_ids_ptr(new unordered_set<string>)
 , cohorts_with_extract_file_ptr(new unordered_set<int>)
 , variant_cohort_idx_ptr(new unordered_map<string, unordered_set<int> >) {}

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

void VariantFilter::set_cohort_extract_file(const string& extract_file, int cohort_idx) {
  this->cohorts_with_extract_file_ptr->insert(cohort_idx);
  BgzReader reader(extract_file);
  while (!reader.eof()) {
    (*this->variant_cohort_idx_ptr)[reader.readline()].insert(cohort_idx);
  }
}

bool VariantFilter::include_htp_record(const htpv4_record_t& rec, int cohort_idx) {
  // Handle cohort-specific extract files first
  if (this->cohorts_with_extract_file_ptr->count(cohort_idx) > 0) {
    return this->exclude_ids_ptr->count(rec.name) == 0 &&
           this->variant_cohort_idx_ptr->count(rec.name) > 0 &&
           this->variant_cohort_idx_ptr->at(rec.name).count(cohort_idx) > 0;
  }

  // Check source filter
  if (!this->sources.empty()) {
    auto source_it = rec.info.find("SOURCE");
    if (source_it == rec.info.end() || this->sources.count(source_it->second) == 0) {
      return false;
    }
  }

  // Check effect filters
  if (this->effect_not_na && !HTPv4Reader::has_beta(rec)) {
    return false;
  }
  if (this->effect_is_na && rec.effect != HTPv4_NA) {
    return false;
  }

  // Check info field requirements
  if (this->info_has_se && !HTPv4Reader::has_se(rec)) {
    return false;
  }
  if (this->info_has_score && rec.info.count("SCORE") == 0) {
    return false;
  }

  // Check gene-specific filters
  const bool is_genep_or_sbat_record = rec.name.size() >= 5 && rec.name.rfind(".GENE") == rec.name.size() - 5;

  if (this->skip_genep && is_genep_or_sbat_record) {
    const bool is_sbat = rec.model.find("SBAT_NEG") != string::npos ||
                          rec.model.find("SBAT-META_NEG") != string::npos ||
                          rec.model.find("SBAT_POS") != string::npos ||
                          rec.model.find("SBAT-META_POS") != string::npos;
    if (!is_sbat) {
      return false;
    }
  }

  // Keep REMETA gene-based tests but drop REGENIE's
  // Also keep single variant results
  const bool is_genetest = rec.name.find(".GENE") != string::npos;
  if (this->keep_remeta_gene_only) {
    return !is_genetest || rec.model.find("REMETA") != string::npos;
  }

  if (this->keep_remeta_gene_and_sbat) {
    if (!is_genetest) {
      return true;
    }

    // Allow REMETA models
    if (rec.model.find("REMETA") != string::npos) {
      return true;
    }

    // Allow specific SBAT models
    const bool is_sbat_burden = (rec.model.rfind("ADD-WGR-BURDEN-SBAT_POS", 0) == 0 ||
                                 rec.model.rfind("ADD-WGR-BURDEN-SBAT_NEG", 0) == 0) &&
                                rec.alt == "set";
    return is_sbat_burden;
  }

  // Check extract/exclude files
  if (this->extract_file_is_set()) {
    return this->extract_ids_ptr->count(rec.name) > 0;
  }
  if (this->exclude_file_is_set()) {
    return this->exclude_ids_ptr->count(rec.name) == 0;
  }

  return true;
}

bool VariantFilter::exclude_htp_record(const htpv4_record_t& rec, int cohort_idx) {
  return !this->include_htp_record(rec, cohort_idx);
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