#include "htpv4_writer.hpp"

#include <cmath>

#include "../util.hpp"

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
    + htpv4_pval_to_string(
      rec.pval,
      rec.info.count("LOG10P") > 0 ? rec.info.at("LOG10P") : "",
      rec.info.count(HTPv4_PVAL_STRING) > 0 ? rec.info.at(HTPv4_PVAL_STRING) : ""
    ) + "\t"
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