#include "allele_frequency_map.hpp"

#include <string>
using std::string;

#include "logging.hpp"
#include "util.hpp"

AlleleFrequencyMap::AlleleFrequencyMap()
 : filepath("")
 , af_file()
 , af_buffer(unordered_map<string, af_info_t>())
 , buffer_chrom("")
 , buffer_start(0)
 , buffer_end(0) { }

void AlleleFrequencyMap::load(const string& filepath) {
  this->filepath = filepath;
  this->af_file.open(filepath);
}

void AlleleFrequencyMap::clear() {
  this->filepath = "";
  this->af_file.close();
  this->af_buffer.clear();
  this->buffer_chrom = "";
  this->buffer_start = 0;
  this->buffer_end = 0;
}

bool AlleleFrequencyMap::loaded() {
  return !this->af_file.closed();
}

af_info_t AlleleFrequencyMap::get_af_info(const string& cpra) {
  if (!this->loaded()) {
    log_error("AlleleFrequencyMap bad call to get_af_info (no allele frequency file loaded)", 1);
  }

  int i = cpra.find(":");
  int j = cpra.find(":", i+1);
  string chrom = cpra.substr(0, i);
  int pos = stoi(cpra.substr(i+1, j-i-1));

  if (this->af_file.indexed()) {
    this->load_chunk(chrom, pos);
  } else {
    this->load_chrom(chrom);
  }

  if (this->af_buffer.count(cpra) > 0) {
    return this->af_buffer.at(cpra);
  } else {
    return AF_INFO_NOT_FOUND;
  }
}

af_info_t AlleleFrequencyMap::parse_af_file_line(const string& line) {
  vector<string> line_split = util::str_split(line, "\t ");
  // If two columns are found, assumes the columns are CPRA and frequency.
  if (line_split.size() == 2) {
    string cpra = line_split[0];
    int i = cpra.find(":");
    int j = cpra.find(":", i+1);
    string chrom = cpra.substr(0, i);
    int pos = stoi(cpra.substr(i+1, j-i-1));

    af_info_t af_info {
      cpra,
      chrom,
      pos,
      stod(line_split[1]),
      UNSPECIFIED
    };
    return af_info;

  // If five columns found: ID FREQ IS_SINGLETON CHROM POS
  } else if (line_split.size() == 5) {
    return af_info_t {
      line_split[0],
      line_split[3],
      stoi(line_split[4]),
      stod(line_split[1]),
      line_split[2] == "1" ? YES : NO,
    };
  // If six columns are found, assumes a plink afreq file: CHROM ID REF ALT ALT_FREQS OBS_CT
  } else if (line_split.size() == 6) {
    string cpra = line_split[1];
    int i = cpra.find(":");
    int j = cpra.find(":", i+1);
    int pos = stoi(cpra.substr(i+1, j-i-1));
    return af_info_t {
      cpra,
      line_split[0],
      pos,
      stod(line_split[4]),
      UNSPECIFIED
    };
  } else {
    log_error("unrecognized allele frequency file format", 1);
    return AF_INFO_NOT_FOUND;
  }
}

// when a file is not indexed, load one chromosome at a time
void AlleleFrequencyMap::load_chrom(const string& chrom) {
  if (chrom == this->buffer_chrom) return;
  this->af_file.open(this->filepath);
  this->af_buffer.clear();
  this->buffer_chrom = chrom;

  string line;
  while (!this->af_file.eof()) {
    line = af_file.readline();
    if (line[0] == '#') continue;
    af_info_t af_info = this->parse_af_file_line(line);
    if (af_info.chrom == chrom) {
      this->af_buffer[af_info.cpra] = af_info;
    }
  }
}

void AlleleFrequencyMap::load_chunk(const string& chrom, const int& pos) {
  if (this->in_current_chunk(chrom, pos)) {
    return ;
  }

  if (!this->in_next_chunk(chrom, pos)) {
    this->af_file.seek(chrom, pos);
    this->buffer_chrom = chrom;
    this->buffer_start = pos;
    this->buffer_end = pos + 1e6*BUFFER_SIZE_MB;
  } else {
    this->buffer_chrom = chrom;
    this->buffer_start = this->buffer_end + 1;
    this->buffer_end = this->buffer_start + 1e6*BUFFER_SIZE_MB;
  }
  // clear the old buffer
  this->af_buffer.clear();

  // add any upcoming variants we have already read
  if (this->has_next_af_info) {
    this->af_buffer[this->next_af_info.cpra] = this->next_af_info;
    this->has_next_af_info = false;
  }

  // fill the buffer with the next chunk
  string line;
  while (!this->af_file.eof()) {
    line = this->af_file.readline();
    // skip header lines
    if (line[0] == '#') continue;
    af_info_t af_info = this->parse_af_file_line(line);
    this->af_buffer[af_info.cpra] = af_info;

    // if we read past the buffer, store the variant to add
    // when loading the next chunk
    if (af_info.pos > this->buffer_end) {
      this->next_af_info = af_info;
      this->has_next_af_info = true;
      break ;
    }
    if (af_info.chrom != this->buffer_chrom) {
      this->next_af_info = af_info;
      this->has_next_af_info = true;
      break ;
    }
  }
}

bool AlleleFrequencyMap::in_current_chunk(const string& chrom, const int& pos) {
  return chrom == this->buffer_chrom && pos >= this->buffer_start && pos <= this->buffer_end;
}

bool AlleleFrequencyMap::in_next_chunk(const string& chrom, const int& pos) {
  return chrom == this->buffer_chrom
                  && (pos >= this->buffer_start)
                  && (pos < this->buffer_end + 1e6*BUFFER_SIZE_MB);
}