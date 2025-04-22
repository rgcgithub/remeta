#ifndef ANNO_READER_H
#define ANNO_READER_H

#include "bgz_reader.hpp"
#include "regenie_anno_reader.hpp"
#include "tabixed_anno_reader.hpp"

class AnnoReader {
 public:
  AnnoReader() : type(NONE) {}

  void open(const string& filepath);

  void close();

  annorec_t readrec();

  void seek(const string& chrom, int pos);

  bool eof();

  bool closed();

  bool indexed();

 private:
  RegenieAnnoReader regenie_reader;
  TabixedAnnoReader tabixed_reader;

  enum { REGENIE, TABIXED, NONE } type;
};

#endif