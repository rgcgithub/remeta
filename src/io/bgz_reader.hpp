/* bgz_reader.hpp
 * Author: Tyler Joseph
 *
 * HTSlib wrapper to read and query tabix indexed bgz files.
 *
 * Example:
 *   // load an Bgz file
 *   BgzReader reader = BgzReader("myfile.gz");
 *
 *   // read a line
 *   string line = reader.readline();
 *
 *   // if an index exists, skip to a region
 *   // otherwise throw a runtime_exception
 *   reader.seek("22", 15528158));
 *   line = reader.readline();
 *
 *   // reading a file
 *   while (!reader.eof()) {
 *      string line = reader.readline();
 *   }
 *
 */

#ifndef BGZ_READER_H
#define BGZ_READER_H

#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class BgzReader {
 public:
  BgzReader();
  BgzReader(string filepath);
  BgzReader(const BgzReader& other);
  BgzReader& operator=(BgzReader other);
  ~BgzReader();

  void open(string filepath);
  void close();

  friend void swap(BgzReader& first, BgzReader& second);

  string readline();
  virtual void seek(string chrom, hts_pos_t position);
  bool eof() { return this->at_eof; };
  bool closed() const { return this->is_closed; }
  bool indexed() const { return this->has_index; }
  string get_filepath() const { return this->filepath; }

  int64_t tell();
  void seek_addr(int64_t addr);

  template<typename T>
  T read() {
    ssize_t size = sizeof(T);
    if (this->closed()) {
      throw std::runtime_error("attempted to read closed file " + this->filepath);
    }
    if (this->at_eof) {
      throw std::runtime_error("read past end of file " + this->filepath);
    }

    T data;
    ssize_t bytes_read = bgzf_read(this->bgzf, &data, size);
    if (bytes_read != size) {
      throw std::runtime_error("failed to read " + this->filepath);
    }
    this->check_eof();
    return data;
  }

  string read_bytes(ssize_t size);

 private:
  void check_eof();
  void read_header();

  string filepath;
  BGZF* bgzf;
  tbx_t* bgzf_tbx;
  hts_itr_t* bgzf_itr;
  bool seek_set;
  string seek_buffer;
  kstring_t buffer;
  bool has_index;
  string line;
  bool at_eof;
  bool is_closed;
  int lines_read;
};

#endif
