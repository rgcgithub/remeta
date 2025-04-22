#include "allele_freq_writer.hpp"

#include <algorithm>
#include <string>
using std::string;

#include "../logging.hpp"
#include "../util.hpp"

AlleleFreqWriter::AlleleFreqWriter() {}

AlleleFreqWriter::~AlleleFreqWriter() {
  if (!this->writer.is_closed()) {
    this->close();
  }
}

AlleleFreqWriter::AlleleFreqWriter(AlleleFreqWriter&& other)
  : filepath(other.filepath) 
  , writer(std::move(other.writer))
  , variants_inserted(std::move(other.variants_inserted)) 
  , variants(std::move(other.variants)) {}

void AlleleFreqWriter::open(const string& filepath) {
  this->writer.open(filepath, "w");
}

void AlleleFreqWriter::close() {
  this->writer.close();
}

bool AlleleFreqWriter::closed() {
  return this->writer.is_closed();
}

void AlleleFreqWriter::add_freq(const string& vid, const string& chr, int pos, double aaf, bool is_singleton) {
  if (this->variants_inserted.count(vid) == 0) {
    this->variants_inserted.insert(vid);
    this->variants.push_back({vid, chr, pos, aaf, is_singleton, HTPv4Pos(chr, pos)});
  }
}

void AlleleFreqWriter::write_freqs_before(const std::string& chrom, int before_pos) {
  if (this->variants.size() == 0) return;
  std::sort(this->variants.begin(), this->variants.end(),
        [](const AlleleFreq& a, const AlleleFreq& b) {
          return a.htp_pos < b.htp_pos;
        });

  for (const AlleleFreq& freq : this->variants) {
    if (freq.chrom == chrom && freq.pos <= before_pos) {
      this->writer.write(freq.id + "\t" + htpv4_to_string(freq.aaf) + "\t" + (freq.is_singleton ? "1\t" : "0\t") + freq.chrom + "\t" + to_string(freq.pos) + "\n");
    }
  }
}

void AlleleFreqWriter::clear_freqs_before(const std::string& chrom, int before_pos) {
  auto it = this->variants.begin();
  while (it != this->variants.end()) {
    if (it->chrom == chrom && it->pos <= before_pos) {
      it = this->variants.erase(it);
      this->variants_inserted.erase(it->id);
    } else {
      ++it;
    }
  }
}