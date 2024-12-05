#include "htpv4_reader.hpp"

#include <boost/lexical_cast.hpp>

#include "../logging.hpp"

HTPv4Reader::HTPv4Reader(string filepath)
  : BgzReader(filepath) 
  , first_line("") {
  // skip header lines
  while (first_line == "" || first_line[0] == '#' || first_line.find("Name") != string::npos) {
    if (this->eof()) {
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
    if (this->eof()) {
      break;
    }
    first_line = this->readline();
  }
  if (first_line.find("Name") != string::npos) {
    first_line = "";
  }
}

htpv4_record_t HTPv4Reader::readrec() {
  string line;;
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
  try {
    rec.pval = (tmp == "NA") ? HTPv4_NA : stod(tmp);
  } catch (const std::out_of_range& e) {
    rec.pval = boost::lexical_cast<double>(tmp);
    log_warning("pvalue underflow for " + rec.name + " (pval: " + tmp + ")");
  }

  ss >> tmp;
  rec.aaf = (tmp == "NA") ? HTPv4_NA : stod(tmp);

  ss >> tmp;
  rec.num_cases = (tmp == "NA") ? HTPv4_NA : stoi(tmp);

  ss >> tmp;
  rec.cases_ref = (tmp == "NA") ? HTPv4_NA : stod(tmp);

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

  rec.info = map<string, string>();
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
  return rec;
}

void HTPv4Reader::seek(string chrom, hts_pos_t position) {
  this->first_line = "";
  BgzReader::seek(chrom, position);
}

double HTPv4Reader::get_beta(const htpv4_record_t& rec) {
  if (rec.info.count("REGENIE_BETA") > 0) {
    return stod(rec.info.at("REGENIE_BETA"));
  } else if (rec.info.count("BETA") > 0) {
    return stod(rec.info.at("BETA"));
  } else if (rec.info.count("EXTERNAL_GWAS_BETA") > 0) {
    return stod(rec.info.at("EXTERNAL_GWAS_BETA"));
  } else if (rec.num_controls == HTPv4_NA) { // is a qt
    return rec.effect;
  } else if (rec.num_controls > 0) { // is a bt
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
  } else if (rec.info.count("EXTERNAL_GWAS_SE") > 0 ) {
    return stod(rec.info.at("EXTERNAL_GWAS_SE"));
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