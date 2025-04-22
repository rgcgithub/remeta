#include "tabixed_anno_reader.hpp"

#include <string>
using std::string;
#include <vector>
using std::vector;

#include "../util.hpp"

void TabixedAnnoReader::open(const string& filepath) {
  this->filepath = filepath;
  this->reader.open(filepath);
}

void TabixedAnnoReader::close() {
  this->reader.close();
}

bool TabixedAnnoReader::closed() {
  return this->reader.closed();
}

bool TabixedAnnoReader::eof() {
  return this->reader.eof();
}

bool TabixedAnnoReader::indexed() {
  return this->reader.indexed();
}

void TabixedAnnoReader::seek(const string& chrom, int pos) {
  this->reader.seek(chrom, pos);
}

annorec_t TabixedAnnoReader::readrec() {
  string line = this->reader.readline();
  vector<string> fields = util::str_split(line, "\t ");
  if (fields.size() < 5) {
    throw std::runtime_error("TabixedAnnoReader expected at least 5 fields in tabixed annotation file");
  }
  annorec_t rec {
    .name = fields[0],
    .chrom = fields[3],
    .gene = fields[1],
    .annotation = fields[2],
    .pos = std::stoi(fields[4])
  };
  return rec;
}