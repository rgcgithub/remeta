#ifndef REGENIE_ANNO_READER_H
#define REGENIE_ANNO_READER_H

#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

/**
 * Data structure used to store index.
 *
 * TODO: Consider breaking out into own class.
 */
typedef int64_t rdr_size_t;
typedef struct anno_reader_idx {
  // Number of bases between each stored position in a file. A large stride
  // results in a smaller index but slower file access. A small stride results
  // in a larger index but faster file access.
  rdr_size_t stride{10000}; 

  // Store pointers to a position in each file in a vector, one vector per
  // chromosome. If i is an index and v a vector of pointers, then reading the
  // file starting from v[i] guarantees that all entries with position 
  // < min_pos + i*stride along the contig have been skipped.
  map<string,vector<rdr_size_t> > chr_ptr;

  // First position in the index per chromosome. If s is the starting 
  // position of the first entry of a chromosome chr in the annotation file,
  // then chr_to_min_pos[chr] = s - (s % stride)
  map<string,rdr_size_t> chr_to_min_pos;
} anno_reader_idx;


/**
 *  Data structure to store records.
 */
typedef struct annorec_t {
  string cpra, chrom, ref, alt, gene, annotation;
  int pos;
} annorec_t;


class RegenieAnnoReader {
 public:
  RegenieAnnoReader(string filepath);
  ~RegenieAnnoReader();
  RegenieAnnoReader(const RegenieAnnoReader& other);
  RegenieAnnoReader& operator=(RegenieAnnoReader other);

  friend void swap(RegenieAnnoReader& first, RegenieAnnoReader& second);

  string readline();
  annorec_t readrec();
  annorec_t parse_line(string line);
  void seek(string chrom, int pos);
  bool eof() { return this->at_eof; };
  bool indexed() { return this->has_index; };
  
  void build_index(int stride);
  void write_index(); // writes filepath.rgi index
  void load_index(); // loads filepath.rgi

 private:
  void check_eof();

  string filepath;
  bool has_index;
  string line;
  // buffer used for seek, if nonempty should be returned
  // before reading data from BGZF file
  string line_buffer; 
  bool at_eof;

  BGZF* bgzf;
  kstring_t buffer; // buffer to read data from bgzf

  anno_reader_idx index;
};

#endif