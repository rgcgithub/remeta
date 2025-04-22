#include "anno_reader.hpp"

#include <string>
using std::string;
#include <vector>
using std::vector;

#include "../logging.hpp"
#include "../util.hpp"

void AnnoReader::open(const string& filepath) {
  BgzReader bgz_reader;
  bgz_reader.open(filepath);
  int ncols = 0;
  if (!bgz_reader.eof()) {
    string line = bgz_reader.readline();
    vector<string> fields = util::str_split(line, "\t ");
    ncols = fields.size();
  }
  bgz_reader.close();

  this->tabixed_reader.open(filepath);
  if (ncols >= 5) {
    log_debug("reading tabixed annotation file");
    this->type = TABIXED;
  } else {
    log_debug("reading regenie annotation file");
    this->tabixed_reader.close();
    this->regenie_reader.open(filepath);
    this->type = REGENIE;
  }
}

void AnnoReader::close() {
  if (this->type == TABIXED) {
    this->tabixed_reader.close();
  } else if (this->type == REGENIE) {
    this->regenie_reader.close();
  }
}

annorec_t AnnoReader::readrec() {
  if (this->type == TABIXED) {
    return this->tabixed_reader.readrec();
  } else if (this->type == REGENIE) {
    return this->regenie_reader.readrec();
  } else {
    throw runtime_error("AnnoReader not opened");
  }
}

void AnnoReader::seek(const string& chrom, int pos) {
  if (this->type == TABIXED) {
    this->tabixed_reader.seek(chrom, pos);
  } else if (this->type == REGENIE) {
    this->regenie_reader.seek(chrom, pos);
  } else {
    throw runtime_error("AnnoReader not opened");
  }
}

bool AnnoReader::eof() {
  if (this->type == TABIXED) {
    return this->tabixed_reader.eof();
  } else if (this->type == REGENIE) {
    return this->regenie_reader.eof();
  } else {
    throw runtime_error("AnnoReader not opened");
  }
}

bool AnnoReader::closed() {
  if (this->type == TABIXED) {
    return this->tabixed_reader.closed();
  } else if (this->type == REGENIE) {
    return this->regenie_reader.closed();
  } else {
    throw runtime_error("AnnoReader not opened");
  }
}

bool AnnoReader::indexed() {
  if (this->type == TABIXED) {
    return this->tabixed_reader.indexed();
  } else if (this->type == REGENIE) {
    return this->regenie_reader.indexed();
  } else {
    throw runtime_error("AnnoReader not opened");
  }
}