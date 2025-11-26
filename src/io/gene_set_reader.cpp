#include "gene_set_reader.hpp"

Gene::Gene(string name,
           string chrom,
           const int start_pos,
           const int end_pos,
           shared_ptr<const unordered_set<string> > var_ids,
           shared_ptr<const vector<string> > var_ids_order)
 : name(name)
 , chrom(chrom)
 , start_pos(start_pos)
 , end_pos(end_pos)
 , var_ids(var_ids)
 , var_ids_order(var_ids_order) { };

bool Gene::contains(const string& var_id) const {
  return var_ids->find(var_id) != var_ids->end();
}

shared_ptr<const vector<string> > Gene::get_variants() const {
  return this->var_ids_order;
}

size_t Gene::size() const { return this->var_ids->size(); }
string Gene::get_name() const { return this->name; }
string Gene::get_chrom() const { return this->chrom; }
int Gene::get_start() const { return this->start_pos; }
int Gene::get_end() const { return this->end_pos; }

GeneSetReader::GeneSetReader(const string& set_file)
 : set_file(set_file)
 , reader(set_file) { }

bool GeneSetReader::eof() {
  return reader.eof();
}

Gene GeneSetReader::read_gene() {
  if (this->eof()) {
    throw runtime_error("read past end of file " + this->set_file);
  }

  string line = this->reader.readline();
  stringstream ss(line);
  
  string gene_name, chrom, var_ids;
  int start_pos;
  ss >> gene_name;
  ss >> chrom;
  ss >> start_pos;
  ss >> var_ids;

  string var_id;
  size_t i, j;
  int pos, end_pos = 0;
  shared_ptr<unordered_set<string> > var_ids_ptr = make_shared<unordered_set<string> >();
  shared_ptr<vector<string> > var_ids_order_ptr = make_shared<vector<string> >();
  ss = stringstream(var_ids);
  while (getline(ss, var_id, ',')) {
    if (var_ids_ptr->count(var_id) > 0) {
      continue;
    }
    var_ids_ptr->insert(var_id);
    var_ids_order_ptr->push_back(var_id);
    i = var_id.find(":");
    j = var_id.find(":", i+1);
    if (i != string::npos && j != string::npos) {
      pos = stoi(var_id.substr(i+1, j-i-1));
      start_pos = min(
        start_pos,
        pos
      );
      end_pos = max(
        end_pos,
        pos
      );
    }
  }

  return Gene(gene_name, chrom, start_pos, end_pos, var_ids_ptr, var_ids_order_ptr);
}