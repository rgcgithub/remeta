#include "htpv4_writer.hpp"

#include <boost/format.hpp>

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

string htpv4_make_info_string(map<string, string> info) {
  string s = "";
  if (info.size() == 0) {
    return "NA";
  }
  for (auto it = info.begin(); it != info.end(); ++it) {
    if (it->first == HTPv4_ESTIMATED_GENOTYPE_COUNTS_FLAG) {
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

void HTPv4Writer::writerec(htpv4_record_t rec) {
  bool keep_as_float = rec.info.count(HTPv4_ESTIMATED_GENOTYPE_COUNTS_FLAG) > 0;
  string line = 
    rec.name+ "\t"
    + rec.chr + "\t"
    + htpv4_to_string(rec.pos) + "\t"
    + rec.ref + "\t"
    + rec.alt + "\t"
    + rec.trait + "\t"
    + rec.cohort + "\t"
    + rec.model + "\t"
    + htpv4_to_string(rec.effect) + "\t"
    + htpv4_to_string(rec.lci_effect) + "\t"
    + htpv4_to_string(rec.uci_effect) + "\t"
    + htpv4_to_string(rec.pval) + "\t"
    + htpv4_to_string(rec.aaf) + "\t"
    + htpv4_to_string(rec.num_cases) + "\t"
    + htpv4_rounded_double_to_string(rec.cases_ref, keep_as_float) + "\t"
    + htpv4_rounded_double_to_string(rec.cases_het, keep_as_float) + "\t"
    + htpv4_rounded_double_to_string(rec.cases_alt, keep_as_float) + "\t"
    + htpv4_to_string(rec.num_controls) + "\t"
    + htpv4_rounded_double_to_string(rec.controls_ref, keep_as_float) + "\t"
    + htpv4_rounded_double_to_string(rec.controls_het, keep_as_float) + "\t"
    + htpv4_rounded_double_to_string(rec.controls_alt, keep_as_float) + "\t"
    + htpv4_make_info_string(rec.info) + "\n";

  BgzWriter::write(line);
}

void HTPv4Writer::writeheader() {
  string header = 
    "Name\tChr\tPos\tRef\tAlt\tTrait\tCohort\tModel\t"
    "Effect\tLCI_Effect\tUCI_Effect\tPval\tAAF\t"
    "Num_Cases\tCases_Ref\tCases_Het\tCases_Alt\t"
    "Num_Controls\tControls_Ref\tControls_Het\tControls_Alt\t"
    "Info\n";
  BgzWriter::write(header);
}