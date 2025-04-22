#ifndef TABIXED_ANNO_READER_H
#define TABIXED_ANNO_READER_H

#include <string>

#include "regenie_anno_reader.hpp"
#include "bgz_reader.hpp"

class TabixedAnnoReader {
 public:
  TabixedAnnoReader() {};
  ~TabixedAnnoReader() { this->reader.close(); };
  TabixedAnnoReader(const TabixedAnnoReader& other) = delete;
  TabixedAnnoReader& operator=(TabixedAnnoReader other) = delete;

  void open(const std::string& filepath);

  void close();

  bool closed();

  bool eof();

  bool indexed();

  void seek(const std::string& chrom, int position);

  annorec_t readrec();

 private:
  std::string filepath;
  BgzReader reader;
};
#endif