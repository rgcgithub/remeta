#ifndef GENE_SET_READER_H
#define GENE_SET_READER_H

#include <memory>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;

#include "bgz_reader.hpp"

class Gene {
 public:
  Gene(
      string name,
      string chrom,
      const int start_pos,
      const int end_pos,
      shared_ptr<const unordered_set<string> > var_ids,
      shared_ptr<const vector<string> > var_ids_order
  );

  bool contains(const string& var_id) const;

  shared_ptr<const vector<string> > get_variants() const;

  size_t size() const;
  string get_name() const;
  string get_chrom() const;
  int get_start() const;
  int get_end() const;

 private:
  string name;
  string chrom;
  int start_pos;
  int end_pos;
  string raw_var_ids;
  bool var_ids_parsed;
  shared_ptr<const unordered_set<string> > var_ids;
  shared_ptr<const vector<string> > var_ids_order;
};

class GeneSetReader {
 public:
  GeneSetReader(const string& set_file);

  Gene read_gene();

  bool eof();

 private:
  string set_file;
  BgzReader reader;
};

#endif