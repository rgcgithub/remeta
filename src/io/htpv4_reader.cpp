#include "htpv4_reader.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "../logging.hpp"
#include "../util.hpp"

string htpv4_to_string(int i) {
  if (i == HTPv4_NA || i == -9) {
    return "NA";
  } else {
    return to_string(i);
  }
}

string htpv4_to_string(double d) {
  if (d == HTPv4_NA || d == -9) {
    return "NA";
  } else {
    return (boost::format("%1$g") % d).str();
  }
}

string htpv4_rounded_double_to_string(double d, bool keep_as_float) {
  if (d == HTPv4_NA || d == -9) {
    return "NA";
  } else if (keep_as_float) {
    return (boost::format("%1.2f") % d).str();
  } else {
    return to_string((int)d);
  }
}

string htpv4_make_info_string(const map<string, string>& info) {
  string s = "";
  if (info.size() == 0) {
    return "NA";
  }
  for (auto it = info.begin(); it != info.end(); ++it) {
    if (it->first == HTPv4_ESTIMATED_GENOTYPE_COUNTS_FLAG || it->first == HTPv4_PVAL_STRING) {
      continue;
    } else if (it->second != "") {
      s += it->first + "=" + it->second;
    } else {
      s += it->first;
    }
    if (next(it) != info.end()) {
      s += ";";
    }
  }
  return s;
}

string htpv4_pval_to_string(double pval, const string& mlog10p, const string& pval_string) {
  if (pval_string != "") {
    return pval_string;
  } else if ( mlog10p != "" && (pval == HTPv4_NA || pval == -9 || pval <= std::numeric_limits<double>::min()) ) {
    if (mlog10p.find("-") != string::npos) {
      log_error("htpv4_pval_to_string expected -log10p but got " + mlog10p, 1);
    }

    // we want to write 10^{-log10p} as 10^{a}*10^{-b} where 0 < a < 1 and b is an integer
    // if log10p = m.n where m and n are the parts before and after the decimal, we can write
    // 10^{-log10p} = 10^{-m - 0.n} = 10^{-m} * 10^{-0.n} = 10^{-m} * a
    vector<string> pval_parts = util::str_split(mlog10p, ".");
    double a = pow(10, -stod("0." + pval_parts[1]));
    double m = -stod(pval_parts[0]);
    if (a < 1) {
      a *= 10;
      m -= 1;
    }

    string base = (boost::format("%1.5f") % a).str();
    while (base.back() == '0' || base.back() == '.') {
      base.pop_back();
    }
    string exp = (boost::format("%1i") % m).str();
    while(exp[0] == '0') {
      exp = exp.substr(1);
    }
    return base + "e" + exp;
  } else {
    return htpv4_to_string(pval);
  }
}

double htpv4_pval_string_to_log10p(string pval) {
  if (pval == "NA" || pval == "-9") {
    return HTPv4_NA;
  }

  double log10p = HTPv4_NA;
  boost::algorithm::to_lower(pval);
  vector<string> pval_parts = util::str_split(pval, "e");
  if (pval_parts.size() == 2) {
    log10p = log10(stod(pval_parts[0])) + stod(pval_parts[1]);
  }
  return log10p;
}

HTPv4Reader::HTPv4Reader(string filepath)
  : BgzReader(filepath) 
  , first_line("") {
  // skip header lines
  while (first_line == "" || first_line[0] == '#' || first_line.find("Name") != string::npos) {
    if (BgzReader::eof()) {
      break;
    }
    first_line = this->readline();
  }
  if (first_line.find("Name") != string::npos) {
    first_line = "";
  }
}

void HTPv4Reader::open(string filepath) {
  BgzReader::open(filepath);
  this->first_line = "";
  while (first_line == "" || first_line[0] == '#' || first_line.find("Name") != string::npos) {
    if (BgzReader::eof()) {
      break;
    }
    first_line = this->readline();
  }
  if (first_line.find("Name") != string::npos) {
    first_line = "";
  }
}

htpv4_record_t HTPv4Reader::readrec() {
  string line;
  if (this->first_line != "") {
    line = this->first_line;
    this->first_line = "";
  } else if (this->eof()) {
    throw out_of_range("read past end of file");
  } else if (this->closed()) {
    throw runtime_error("reading a closed file");
  } else {
    line = this->readline();
  }
  htpv4_record_t rec;
  rec.info = map<string, string>();
  stringstream ss(line);
  string tmp;

  ss >> rec.name >> rec.chr >> rec.pos >> rec.ref >> rec.alt >> rec.trait >>
      rec.cohort >> rec.model;

  ss >> tmp;
  try {
    rec.effect = (tmp == "NA") ? HTPv4_NA : stod(tmp);
  } catch (const std::out_of_range& e) {
    log_warning("effect column out of range for " + rec.name + " (effect: " + tmp + ")");
    rec.effect = HTPv4_NA;
  }

  ss >> tmp;
  try {
    rec.lci_effect = (tmp == "NA") ? HTPv4_NA : stod(tmp);
  } catch (const std::out_of_range& e) {
    log_warning("lci_effect column out of range for " + rec.name + " (lci_effect: " + tmp + ")");
    rec.lci_effect = HTPv4_NA;
  }

  ss >> tmp;
  try {
    rec.uci_effect = (tmp == "NA") ? HTPv4_NA : stod(tmp);
  } catch (const std::out_of_range& e) {
    log_warning("uci_effect column out of range for  " + rec.name + " (uci_effect: " + tmp + ")");
    rec.uci_effect = HTPv4_NA;
  }

  ss >> tmp;
  double log10p = HTPv4_NA;
  try {
    rec.pval = (tmp == "NA") ? HTPv4_NA : stod(tmp);
    log10p = -log10(rec.pval);
  } catch (const std::out_of_range& e) {
    rec.pval = boost::lexical_cast<double>(tmp);
    log_warning("pvalue underflow for " + rec.name + " (pval: " + tmp + ")");

    // try to compute a log p-value
    log10p = -htpv4_pval_string_to_log10p(tmp);
  }
  rec.info[HTPv4_PVAL_STRING] = tmp;

  bool geno_are_float = false;
  ss >> tmp;
  rec.aaf = (tmp == "NA") ? HTPv4_NA : stod(tmp);

  ss >> tmp;
  rec.num_cases = (tmp == "NA") ? HTPv4_NA : stoi(tmp);

  ss >> tmp;
  rec.cases_ref = (tmp == "NA") ? HTPv4_NA : stod(tmp);
  geno_are_float = tmp.find(".") != string::npos;

  ss >> tmp;
  rec.cases_het = (tmp == "NA") ? HTPv4_NA : stod(tmp);

  ss >> tmp;
  rec.cases_alt = (tmp == "NA") ? HTPv4_NA : stod(tmp);

  ss >> tmp;
  rec.num_controls = (tmp == "NA") ? HTPv4_NA : stoi(tmp);

  ss >> tmp;
  rec.controls_ref = (tmp == "NA") ? HTPv4_NA : stod(tmp);

  ss >> tmp;
  rec.controls_het = (tmp == "NA") ? HTPv4_NA : stod(tmp);

  ss >> tmp;
  rec.controls_alt = (tmp == "NA") ? HTPv4_NA : stod(tmp);

  if (ss.fail()) {
    throw runtime_error("failed to parse htpv4 record at " + rec.name);
  }

  
  if (ss.eof()) {
    return rec;
  }

  ss >> tmp;
  ss = stringstream(tmp);
  string info_field;
  while (getline(ss, info_field, ';')) {
    size_t pos = info_field.find("=");
    if (pos != string::npos) {
      rec.info[info_field.substr(0, pos)] =
          info_field.substr(pos + 1, info_field.size());
    } else {
      rec.info[info_field] = string("");
    }
    if (ss.fail()) {
      throw runtime_error("failed to parse htpv4 info at at " + rec.name + " info: " + tmp);
    }
  }

  if (log10p != HTPv4_NA && rec.info.count("LOG10P") == 0) {
    rec.info["LOG10P"] = to_string(log10p);
  }

  if (geno_are_float) {
    rec.info[HTPv4_ESTIMATED_GENOTYPE_COUNTS_FLAG] = "";
  }
  return rec;
}

void HTPv4Reader::seek(string chrom, hts_pos_t position) {
  this->first_line = "";
  BgzReader::seek(chrom, position);
}

bool HTPv4Reader::eof() {
  return this->first_line == "" && BgzReader::eof();
}

double HTPv4Reader::get_beta(const htpv4_record_t& rec) {
  if (rec.info.count("REGENIE_BETA") > 0) {
    return stod(rec.info.at("REGENIE_BETA"));
  } else if (rec.info.count("BETA") > 0) {
    return stod(rec.info.at("BETA"));
  } else if (rec.info.count("EXTERNAL_GWAS_BETA") > 0) {
    return stod(rec.info.at("EXTERNAL_GWAS_BETA"));
  } else if (rec.info.count("INGENIE_BETA")) {
    return stod(rec.info.at("INGENIE_BETA"));
  } else if (rec.num_controls == HTPv4_NA) { // is a qt
    return rec.effect;
  } else if (rec.num_controls > 0 && rec.effect != 0) { // is a bt
    return log(rec.effect);
  } else {
    throw runtime_error("could not find beta for: " + rec.name);
  }
}

bool HTPv4Reader::has_beta(const htpv4_record_t& rec) {
  try {
    HTPv4Reader::get_beta(rec);
    return true;
  } catch (const runtime_error& e) {
    return false;
  }
}

double HTPv4Reader::get_se(const htpv4_record_t& rec) {
  if (rec.info.count("REGENIE_SE") > 0) {
    return stod(rec.info.at("REGENIE_SE"));
  } else if (rec.info.count("SE") > 0) {
    return stod(rec.info.at("SE"));
  } else if (rec.info.count("META_SE") > 0) {
    return stod(rec.info.at("META_SE"));
  } else if (rec.info.count("EXTERNAL_GWAS_SE") > 0 ) {
    return stod(rec.info.at("EXTERNAL_GWAS_SE"));
  } else if (rec.info.count("INGENIE_SE")) {
    return stod(rec.info.at("INGENIE_SE"));
  } else {
    throw runtime_error("could not find se for: " + rec.name);
  }
}

bool HTPv4Reader::has_se(const htpv4_record_t& rec) {
  try {
    HTPv4Reader::get_se(rec);
    return true;
  } catch (const runtime_error& e) {
    return false;
  }
}