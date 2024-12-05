#include "genep_meta_analyzer.hpp"

#include <sstream>

#include "../io/bgz_reader.hpp"
#include "../stat/tests.hpp"
#include "../htpv4_pos.hpp"
#include "../logging.hpp"

void htpv4_name_to_gene_mask(const string& name, string& gene, string& mask) {
  size_t start = name.find(".GENE") + 5;
  if (start == string::npos) {
    gene = "";
    mask = "";
    return;
  }

  size_t end = name.find(".", start+1);
  if (end == string::npos) {
    gene = name.substr(0, start);
    mask = "";
    return;
  }

  gene = name.substr(0, start);
  mask = name.substr(start+1, end-start-1);
}

void htpv4_rec_to_gene_mask(const htpv4_record_t& rec, string& gene, string& mask) {
  if (rec.model.find("SBAT") != string::npos) {
    gene = rec.name;
    mask = "";
    return;
  } else if (rec.ref != "ref" || rec.alt == "set") {
    gene = "";
    mask = "";
    return;
  } else {
    size_t name_end = rec.name.find("." + rec.alt);
    if (name_end == string::npos) {
      log_error("expected Name field to end with mask name for gene-based tests", 1);
    }
    size_t mask_end = rec.name.find(".", name_end + 1);
    gene = rec.name.substr(0, name_end);
    mask = rec.name.substr(name_end+1, mask_end - name_end - 1);
    return;
  }
}

string update_source(const string& current_source, const htpv4_record_t& rec) {
  if (rec.info.count("SOURCE") == 0 && current_source != "") {
    return "MULTI";
  } else if (rec.info.count("SOURCE") == 0) {
    return "";
  } else if (current_source == "") {
    return rec.info.at("SOURCE");
  } else if (rec.info.at("SOURCE") != current_source) {
    return "MULTI";
  } else {
    return current_source;
  }
}

string get_genep_source(const vector<string>& sources) {
  string source = sources[0];
  for (const string& s : sources) {
    if (s != "" && source == "") {
      source = s;
    } else if (s != "" && s != source) {
      source = "MULTI";
    }
  }
  return source;
}

GenePMetaAnalyzer::GenePMetaAnalyzer(const string& genep_file,
                                     const string& burden_model,
                                     const string& acatv_model,
                                     const string& skato_model,
                                     const bool& include_sbat)
 : burden_model(burden_model)
 , acatv_model(acatv_model)
 , skato_model(skato_model)
 , include_sbat(include_sbat) {
  if (genep_file != "") {
    this->load_genep_file(genep_file);
  }
 }

void GenePMetaAnalyzer::add_line(const htpv4_record_t& rec, const int& study_index) {
  log_debug("entering GenePMetaAnalyzer::add_line");
  log_debug("processing " + rec.name);
  this->htpv4_records.push_back(rec);
  model_e model = this->parse_model(rec);
  if (model == NONE) {
    return;
  }

  string gene;
  string mask;
  // htpv4_name_to_gene_mask(rec.name, gene, mask);
  htpv4_rec_to_gene_mask(rec, gene, mask);
  // if (mask == "" || this->mask_map.at(mask).size() == 0) {
  //   return;
  // }

  if (this->genep_records.count(gene) == 0) {
    genep_record_t genep_rec;
    genep_rec.name = gene;
    genep_rec.chr = rec.chr,
    genep_rec.pos = rec.pos,
    genep_rec.trait =  rec.trait,
    genep_rec.cohort = rec.cohort;
    genep_rec.acatv_log10_pvals.resize(genep_groups.size(), vector<double>());
    genep_rec.burden_log10_pvals.resize(genep_groups.size(), vector<double>());
    genep_rec.skato_log10_pvals.resize(genep_groups.size(), vector<double>());
    genep_rec.sbat_log10_pvals.resize(genep_groups.size(), vector<double>());
    this->genep_records[gene] = genep_rec;
  }

  double log10p;
  if (rec.info.count("LOG10P") && rec.pval <= std::numeric_limits<double>::min()) {
    log10p = -stod(rec.info.at("LOG10P"));
  } else {
    log10p = log10(max(rec.pval, std::numeric_limits<double>::min()));
  }

  if (mask != "" && this->mask_map.count(mask) > 0) {
    for (const int& genep_idx : this->mask_map.at(mask)) {
      if (model == ACATV) {
        log_debug("adding to ACATV pvalues " + this->genep_groups[genep_idx]);
        this->genep_records[gene].acatv_log10_pvals[genep_idx].push_back(log10p);
        this->genep_records[gene].acatv_source = update_source(this->genep_records[gene].acatv_source, rec);
      } else if (model == BURDEN) {
        log_debug("adding to BURDEN pvalues " + this->genep_groups[genep_idx]);
        this->genep_records[gene].burden_log10_pvals[genep_idx].push_back(log10p);
        this->genep_records[gene].burden_source = update_source(this->genep_records[gene].burden_source, rec);
      } else if (model == SKATO) {
        log_debug("adding to SKATO pvalues " + this->genep_groups[genep_idx]);
        this->genep_records[gene].skato_log10_pvals[genep_idx].push_back(log10p);
        this->genep_records[gene].skato_source = update_source(this->genep_records[gene].skato_source, rec);
      }
    }
  } else if (this->include_sbat && model == SBAT) {
    // check for gene_p group in rec_model
    for (size_t i = 0; i < this->genep_groups.size(); ++i) {
      if (rec.model == "ADD-WGR-BURDEN-SBAT_POS_" + this->genep_groups[i]
        || rec.model == "ADD-WGR-BURDEN-SBAT_POS_" + this->genep_groups[i] + "-META"
        || rec.model == "ADD-WGR-BURDEN-SBAT_NEG_" + this->genep_groups[i]
        || rec.model == "ADD-WGR-BURDEN-SBAT_NEG_" + this->genep_groups[i] + "-META") {
        log_debug("adding to SBAT pvalues " + this->genep_groups[i]);
        this->genep_records[gene].sbat_log10_pvals[i].push_back(log10p);
        this->genep_records[gene].sbat_source = update_source(this->genep_records[gene].sbat_source, rec);
      }
    }
  }
  log_debug("leaving GenePMetaAnalyzer::add_line");
}

vector<htpv4_record_t> GenePMetaAnalyzer::meta_analyze_before(const string& chr, const int& pos) {
  log_debug("entering GenePMetaAnalyzer::meta_analyze_before");
  log_debug("computing gene_p of genes on chr " + chr + " before position " + to_string(pos));
  vector<htpv4_record_t> results;
  HTPv4Pos end_pos(chr, pos);
  for (const htpv4_record_t& rec : this->htpv4_records) {
    if (HTPv4Pos(rec.chr, rec.pos) < end_pos) {
      results.push_back(rec);
    }
  }

  for (const auto& rec : this->genep_records) {
    if (HTPv4Pos(rec.second.chr, rec.second.pos) < end_pos) {
      for (size_t i = 0; i < this->genep_groups.size(); ++i) {
        vector<double> genep_log10_pvals;
        if (rec.second.acatv_log10_pvals[i].size() > 0) {
          pair<htpv4_record_t, double> acatv_genep =
            this->meta_analyze_genep_group(rec.second, this->acatv_model, this->genep_groups[i], rec.second.acatv_source, rec.second.acatv_log10_pvals[i]);
            results.push_back(acatv_genep.first);
            genep_log10_pvals.push_back(acatv_genep.second);
        }
        if (rec.second.burden_log10_pvals[i].size() > 0) {
          pair<htpv4_record_t, double> burden_genep =
            this->meta_analyze_genep_group(rec.second, this->burden_model, this->genep_groups[i], rec.second.burden_source, rec.second.burden_log10_pvals[i]);
          results.push_back(burden_genep.first);
          genep_log10_pvals.push_back(burden_genep.second);
        }
        if (rec.second.skato_log10_pvals[i].size() > 0) {
          pair<htpv4_record_t, double> skato_genep =
            this->meta_analyze_genep_group(rec.second, this->skato_model, this->genep_groups[i], rec.second.skato_source, rec.second.skato_log10_pvals[i]);
          results.push_back(skato_genep.first);
          genep_log10_pvals.push_back(skato_genep.second);
        }
        if (rec.second.sbat_log10_pvals[i].size() > 0){
          pair<htpv4_record_t, double> sbat_genep =
            this->meta_analyze_genep_group(rec.second, "BURDEN-SBAT-META", this->genep_groups[i], rec.second.sbat_source, rec.second.sbat_log10_pvals[i]);
          results.push_back(sbat_genep.first);
          genep_log10_pvals.push_back(sbat_genep.second);
        }
        if (genep_log10_pvals.size() > 0) {
          pair<htpv4_record_t, double> genep =
            this->meta_analyze_genep_group(
              rec.second,
              "GENE_P",
              this->genep_groups[i],
              get_genep_source(vector<string>{rec.second.acatv_source, rec.second.burden_source, rec.second.skato_source, rec.second.sbat_source}),
              genep_log10_pvals
            );
          results.push_back(genep.first);
        }
      }
    }
  }
  log_debug("leaving GenePMetaAnalyzer::meta_analyze_before");
  return results;
}

void GenePMetaAnalyzer::clear_before(const string& chr, const int& pos) {
  log_debug("entering GenePMetaAnalyzer::clear_before");
  log_debug("clearing records on chr " + chr + " before position " + to_string(pos));
  HTPv4Pos end_pos(chr, pos);

  for (auto it = this->htpv4_records.begin(); it != this->htpv4_records.end();) {
    if (HTPv4Pos(it->chr, it->pos) < end_pos) {
      it = this->htpv4_records.erase(it);
    } else {
      ++it;
    }
  }

  for (auto it = this->genep_records.begin(); it != this->genep_records.end();) {
    if (HTPv4Pos(it->second.chr, it->second.pos) < end_pos) {
      it = this->genep_records.erase(it);
    } else {
      ++it;
    }
  }
  log_debug("leaving GenePMetaAnalyzer::clear_before");
}

pair<htpv4_record_t, double> GenePMetaAnalyzer::meta_analyze_genep_group(const genep_record_t& genep_record,
                                                                         const string& model,
                                                                         const string& genep_group,
                                                                         const string& source,
                                                                         const vector<double>& log10_pvals) {
  if (log10_pvals.size() == 0) {
    log_error("bad call to GenePMetaAnalyzer::meta_analyze_genep_group - this is a bug", 1);
  }
  htpv4_record_t result {
    genep_record.name,
    genep_record.chr,
    genep_record.pos,
    "ref",
    "set",
    genep_record.trait,
    genep_record.cohort,
    "",
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    HTPv4_NA,
    map<string, string>()
  };
  if (source != "") {
    result.info["SOURCE"] = source;
  }
  result.model = model + "_" + genep_group;
  stat::tests::test_result_t test_result = stat::tests::acat(log10_pvals);
  result.pval = test_result.pval;
  result.info["LOG10P"] = to_string(-test_result.log10p);
  return make_pair(result, test_result.log10p);
}

void GenePMetaAnalyzer::load_genep_file(const string& genep_file) {
  BgzReader reader(genep_file);

  string line;
  string genep_label;
  string mask_label;
  string tmp;
  stringstream ss;
  while (!reader.eof()) {
    line = reader.readline();
    ss = stringstream(line);

    ss >> genep_label;
    this->genep_groups.push_back(genep_label);
    ss >> tmp;
    ss = stringstream(tmp);
    while (getline(ss, mask_label, ',')) {
      this->mask_map[mask_label].push_back(this->genep_groups.size() - 1);
    }
  }
}

model_e GenePMetaAnalyzer::parse_model(const htpv4_record_t& rec) {
  if (rec.model == this->skato_model && rec.alt != "set") {
    return SKATO;
  } else if (rec.model == this->acatv_model && rec.alt != "set") {
    return ACATV;
  } else if (rec.model == this->burden_model && rec.alt != "set") {
    return BURDEN;
  } else if ( (rec.model.rfind("ADD-WGR-BURDEN-SBAT_POS", 0) == 0 || rec.model.rfind("ADD-WGR-BURDEN-SBAT_NEG", 0) == 0) && rec.alt == "set") {
    return SBAT;
  } else {
    return NONE;
  }
}