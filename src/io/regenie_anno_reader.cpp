#include "regenie_anno_reader.hpp"

#include <cstdio>

RegenieAnnoReader::RegenieAnnoReader()
 : has_index(false), at_eof(false), is_closed(true), bgzf(NULL) { }

RegenieAnnoReader::RegenieAnnoReader(string filepath) {
  this->filepath = filepath;
  this->bgzf = bgzf_open(filepath.c_str(), "r");
  this->buffer = KS_INITIALIZE;
  this->at_eof = false;
  this->is_closed = false;
  this->line_buffer = "";

  if (this->bgzf == NULL) {
    throw runtime_error("failed to open " + filepath);
  }

  ifstream index(filepath + ".rgi");
  if (index.good()) {
    this->has_index = true;
    this->load_index();
  } else {
    this->has_index = false;
  }

  this->check_eof();
  if (this->eof()) {
    throw runtime_error("annotation file " + filepath + " is empty");
  }
}

RegenieAnnoReader::~RegenieAnnoReader() {
  if (!is_closed) {
    bgzf_close(this->bgzf);
    ks_free(&this->buffer);
  }
}

RegenieAnnoReader::RegenieAnnoReader(const RegenieAnnoReader &other) {
  this->filepath = other.filepath;
  this->bgzf = bgzf_open(filepath.c_str(), "r");
  this->buffer = KS_INITIALIZE;
  this->at_eof = other.at_eof;
  this->line_buffer = other.line_buffer;

  if (this->bgzf == NULL) {
    throw runtime_error("failed to open " + filepath);
  }

  int ret = bgzf_seek(this->bgzf, bgzf_tell(other.bgzf), SEEK_SET);
  if (ret < 0) throw runtime_error("bgzf_seek failed in copy constructor");

  if (other.has_index) {
    this->load_index();
    this->has_index = true;
  } else {
    this->has_index = false;
  }
}

RegenieAnnoReader &RegenieAnnoReader::operator=(RegenieAnnoReader other) {
  swap(*this, other);
  return *this;
}

void swap(RegenieAnnoReader &first, RegenieAnnoReader &second) {
  std::swap(first.filepath, second.filepath);
  std::swap(first.has_index, second.has_index);
  std::swap(first.line, second.line);
  std::swap(first.line_buffer, second.line_buffer);
  std::swap(first.at_eof, second.at_eof);
  std::swap(first.is_closed, second.is_closed);
  std::swap(first.bgzf, second.bgzf);
  std::swap(first.buffer, second.buffer);
  std::swap(first.index, second.index);
}

void RegenieAnnoReader::open(const string& filepath) {
  RegenieAnnoReader other(filepath);
  swap(*this, other);
}

void RegenieAnnoReader::close() {
  if (!is_closed) {
    bgzf_close(this->bgzf);
    ks_free(&this->buffer);
    is_closed = true;
  }
}

string RegenieAnnoReader::readline() {
  if (this->eof()) {
    throw out_of_range("read past end of file " + this->filepath);
  }
  if (this->line_buffer != "") {
    string tmp(this->line_buffer);
    this->line_buffer = "";
    return tmp;
  }
  int ret = bgzf_getline(this->bgzf, '\n', &this->buffer);
  if (ret <= -2) {
    throw runtime_error("bgzf_getline failed when reading " + this->filepath);
  }
  stringstream ss(buffer.s);
  getline(ss, this->line);
  this->check_eof();
  return line;
}

annorec_t RegenieAnnoReader::readrec() {
  string line = this->readline();
  return parse_line(line);
}

annorec_t RegenieAnnoReader::parse_line(string line) {
  stringstream ss(line);
  annorec_t rec;
  try {
    ss >> rec.name;
    ss >> rec.gene;
    ss >> rec.annotation;

    size_t chr_end = rec.name.find(":");
    if (chr_end == string::npos) {
      throw runtime_error("could not parse cpra");
    }
    rec.chrom = rec.name.substr(0, chr_end);

    size_t pos_end = rec.name.find(":", chr_end + 1);
    rec.pos = stoi(rec.name.substr(chr_end + 1, pos_end - chr_end));  
  } catch(const std::exception& ex) {
    string msg = "could not parse line " + line + " in annotation file " + this->filepath;
    throw runtime_error(msg);
  }
  return rec;
}

void RegenieAnnoReader::check_eof() {
  if (bgzf_peek(this->bgzf) == -1) {
    this->at_eof = true;
  } else {
    this->at_eof = false;
  }
}

void RegenieAnnoReader::build_index(int stride) {
  this->index = anno_reader_idx();

  // open a new file for indexing
  BGZF *bgzf = bgzf_open(this->filepath.c_str(), "r");
  if (bgzf_compression(bgzf) != 2) {
    throw runtime_error("file is not BGZF-compressed");
  }
  this->index.stride = stride;

  kstring_t buffer = KS_INITIALIZE;
  int res = 0;
  string line;
  int last_chrom = -1;
  int last_pos = -1;
  while (bgzf_peek(bgzf) != -1) {
    rdr_size_t addr = bgzf_tell(bgzf);
    res = bgzf_getline(bgzf, '\n', &buffer);
    if (res < 0) throw runtime_error("index error: error reading file");

    stringstream ss(buffer.s);
    getline(ss, line);
    int chr_end = line.find(":");
    int pos_end = line.find(":", chr_end + 1);
    string chrom = line.substr(0, chr_end);
    int pos = stoi(line.substr(chr_end + 1, pos_end - chr_end));

    int chrom_int = (chrom != "X") ? stoi(chrom) : 23;
    if (chrom_int < last_chrom) {
      throw runtime_error("annotation file is not sorted by chromosome");
    } else if (chrom_int == last_chrom && pos < last_pos) {
      string last_chr_pos = to_string(last_chrom) + ":" + to_string(last_pos);
      string chr_pos = string(chrom) + ":" + to_string(pos);
      string msg = "annotation file is not sorted by position ( " +
                   last_chr_pos + " > " + chr_pos + " )";
      throw runtime_error(msg);
    } else {
      last_chrom = chrom_int;
      last_pos = pos;
    }

    if (index.chr_ptr.find(chrom) == index.chr_ptr.end()) {
      index.chr_to_min_pos[chrom] = pos - (pos % stride);
    }

    /*
     *  If the distance between the previous index and next
     *  index is > stride, this fills in the entries between.
     */
    size_t addr_idx = (pos - index.chr_to_min_pos[chrom]) / stride;
    for (auto i = index.chr_ptr[chrom].size(); i <= addr_idx; i++) {
      index.chr_ptr[chrom].push_back(addr);
    }
  }
  bgzf_close(bgzf);
  this->index = index;
  this->has_index = true;
}

/*
 *  TODO: Because duplicate addresses can exist in adjacent entries in an
 *  index, this creates a larger file than required. A better way to store
 *  the index would be to skip duplicates in the file, then fill them in
 *  when loading the index back into memory.
 */
void RegenieAnnoReader::write_index() {
  string outpath = this->filepath + ".rgi";

  ofstream out;
  out.open(outpath, ofstream::binary);

  // first write the size of the stride
  out.write((char *)&this->index.stride, sizeof(rdr_size_t));

  // write the number of chromosomes
  rdr_size_t nchrom = this->index.chr_ptr.size();
  out.write((char *)&nchrom, sizeof(nchrom));

  // for each chromosome:
  // write size of chrom string, chrom string, min pos, size of address vector
  // then write the addresses for the index of that chromosome
  for (auto it = this->index.chr_ptr.begin(); it != this->index.chr_ptr.end();
       it++) {
    rdr_size_t str_size = it->first.size() + 1;  // add one for null character
    rdr_size_t vec_size = it->second.size();
    rdr_size_t pos = this->index.chr_to_min_pos[it->first];
    out.write((char *)&str_size, sizeof(str_size));
    const char *chrom = it->first.c_str();
    out.write((char *)&chrom[0], str_size);
    out.write((char *)&pos, sizeof(pos));
    out.write((char *)&vec_size, sizeof(vec_size));
    for (size_t i = 0; i < it->second.size(); i++) {
      out.write((char *)&it->second[i], sizeof(rdr_size_t));
    }
  }
  out.close();
}

void RegenieAnnoReader::load_index() {
  anno_reader_idx index;
  this->index = index;

  string inpath = this->filepath + ".rgi";
  ifstream in(inpath, ios::binary);
  if (in.fail()) {
    throw runtime_error("failed to open index " + this->filepath + ".rgi");
  }

  in.read((char *)&this->index.stride, sizeof(this->index.stride));

  rdr_size_t nchrom;
  in.read((char *)&nchrom, sizeof(nchrom));

  for (rdr_size_t i = 0; i < nchrom; i++) {
    rdr_size_t str_size;
    in.read((char *)&str_size, sizeof(str_size));

    char chrom[str_size];
    in.read((char *)&chrom, str_size);

    rdr_size_t min_pos;
    in.read((char *)&min_pos, sizeof(min_pos));
    this->index.chr_to_min_pos[chrom] = min_pos;

    rdr_size_t naddr;
    in.read((char *)&naddr, sizeof(naddr));
    this->index.chr_ptr[chrom] = vector<rdr_size_t>(naddr);

    rdr_size_t addr;
    for (rdr_size_t j = 0; j < naddr; j++) {
      in.read((char *)&addr, sizeof(addr));
      this->index.chr_ptr[chrom][j] = addr;
    }
  }
  this->has_index = true;

  if (in.fail()) {
    throw runtime_error("failed to read index");
  }
}

void RegenieAnnoReader::seek(const string& chrom, int pos) {
  if (!this->has_index) {
    throw runtime_error("seek requires and index");
  }
  this->at_eof = false;

  if (this->index.chr_ptr.find(chrom) == this->index.chr_ptr.end()) {
    throw runtime_error("called seek on bad contig");
  }

  size_t addr_idx =
      (pos - this->index.chr_to_min_pos[chrom]) / this->index.stride;
  if (addr_idx >= this->index.chr_ptr[chrom].size()) {
    this->at_eof = true;
    return;
  }

  this->line_buffer = "";
  int ret =
      bgzf_seek(this->bgzf, this->index.chr_ptr[chrom][addr_idx], SEEK_SET);
  if (ret < 0) throw runtime_error("bgzf_seek error");
  string line_buffer = this->readline();
  annorec_t rec = this->parse_line(line_buffer);
  while (rec.pos < pos && rec.chrom == chrom) {
    line_buffer = this->readline();
    rec = this->parse_line(line_buffer);
  }
  this->line_buffer = line_buffer;
}
