#ifndef ANNOTATION_MAP_H
#define ANNOTATION_MAP_H

#include <string>
#include <unordered_map>
using namespace std;

#include "io/regenie_anno_reader.hpp"
#include "mask_set.hpp"

struct annomap_rec_t {
  string cpra;
  string chrom;
  string gene;
  int pos;
  int annotation;
};

class AnnotationMap {
 public:
  AnnotationMap(string filepath, MaskSet mask_set);
  int get_annotation(string cpra, string gene);

  static const int BUFFER_SIZE_MB = 1;

 private:
  string filepath;
  RegenieAnnoReader anno_reader;
  unordered_map<string, vector<annomap_rec_t> > anno_buffer;
  string buffer_chrom;
  int buffer_start;
  int buffer_end;
  MaskSet mask_set;

  annomap_rec_t next_rec;
  bool has_next_rec;

  void load_chunk(string chrom, int pos);
  // for unindexed files, load one chromosome at a time
  void load_chrom(string chrom);
  bool in_current_chunk(string chrom, int pos);
  bool in_next_chunk(string chrom, int pos);
};

#endif