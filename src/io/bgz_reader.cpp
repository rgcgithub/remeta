#include "bgz_reader.hpp"

/*
 * Additional documentation on relevent HTSlib data structures and functions.
 *
 * struct tbx_t
 *      Data structure to store tabix index.
 *
 * struct hts_itr_t
 *      Data structure for iterator of HTS file types.
 *
 * struct kstring_t
 *      Data structure to store an entry from an HTS file. When an HTS
 *      function takes a void *data argument this is typically a kstring_t.
 *
 * tbx_t *tbx_index_load(const chart *fn)
 *      Takes path to a bgzipped file, returns tabix index of the file
 *
 * int tbx_name2id(tbx_t *tbx, const char *ss)
 *      Takes a tabix index and a sequence name, returns the integer
 *      id of that name.
 *
 * hts_itr_t *tbx_itr_queryi(tbx_t tbx, tid target_id, hts_pos_t beg, hts_pos_t end)
 *      Takes a tabix index, sequence id that is the result of tbx_name2id, and
 *      starting and ending regions (as integers). Returns a tabix iterator over
 *      those regions. htslib has several special values to make it easier to
 *      iterate over files
 *
 *      HTS_IDX_START : pass in place of a target id to interate from the start of the file
 *      HTS_POS_MIN   : value representing the smallest start position in a HTS file
 *      HTS_POS_MAX   : value representing the largest start position in a HTS file
 *
 *  hts_itr_t *tbx_itr_querys(tbx_t tbx, const char *reg)
 *      Takes a tabix index and a query string. Returns an iterator.
 *
 *  int *tbx_bgzf_itr_next(BGZF *fp, tbx_t tbx, hts_itr_t itr, void *data)
 *      Takes a tabix iterator, stores the current entry in data (typically a kstring_t) 
 *      and advances the iterator to the next entry.
 *
 */
BgzReader::BgzReader() {
  hts_set_log_level(HTS_LOG_OFF);
  this->filepath = filepath;
  this->bgzf = nullptr;
  this->buffer = KS_INITIALIZE;
  this->seek_set = false;
  this->seek_buffer = "";
  this->at_eof = true;
  this->is_closed = true;
  this->lines_read = 0;
  this->bgzf_tbx = nullptr;
  this->has_index = false;
}

BgzReader::BgzReader(string filepath) {
  hts_set_log_level(HTS_LOG_OFF);
  this->filepath = filepath;
  this->bgzf = bgzf_open(filepath.c_str(), "r");
  this->buffer = KS_INITIALIZE;
  this->seek_set = false;
  this->seek_buffer = "";
  this->at_eof = false;
  this->is_closed = false;
  this->lines_read = 0;

  if (this->bgzf == NULL) {
    throw runtime_error("failed to open " + filepath);
  }

  // Set up HTS lib iterator if the file is tabix indexed
  this->bgzf_tbx = tbx_index_load(filepath.c_str());
  if (this->bgzf_tbx == NULL) {
    this->has_index = false;
  } else {
    this->bgzf_itr = tbx_itr_queryi(this->bgzf_tbx, HTS_IDX_START, 0, 0);
    this->has_index = true;
  }
  this->check_eof();
}

BgzReader::~BgzReader() {
  ks_free(&this->buffer);

  if (!this->closed()) {
    bgzf_close(this->bgzf);
  }

  if (this->has_index) {
    tbx_destroy(this->bgzf_tbx);
    hts_itr_destroy(this->bgzf_itr);
  }
}

BgzReader::BgzReader(const BgzReader& other)
 : BgzReader() {
  if (other.closed()) {
    return ;
  }
  this->filepath = other.filepath;
  this->bgzf = bgzf_open(filepath.c_str(), "r");
  this->seek_set = other.seek_set;
  this->seek_buffer = other.seek_buffer;
  this->buffer = KS_INITIALIZE;
  this->at_eof = other.at_eof;
  this->is_closed = other.is_closed;
  this->lines_read = 0;

  if (this->bgzf == NULL) {
    throw runtime_error("failed to open " + filepath);
  }

  if (this->eof()) {
    return;
  }

  int compression = bgzf_compression(this->bgzf);
  // is a bgzip file
  if (compression == 2) {
    auto pos = bgzf_tell(other.bgzf);
    int res = bgzf_seek(this->bgzf, pos, SEEK_SET);
    if (res == -1 && !this->at_eof) throw runtime_error("BgzReader copy error in: " + this->filepath); 
  // If the file is uncompressed, scan to the current line of the file.
  // This can be slow for large files, so it's probably best to avoid
  // copying if possible.
  } else {
    for (int i = 0; i < other.lines_read; ++i) {
      this->readline();
    }
  }

  this->bgzf_tbx = tbx_index_load(filepath.c_str());
  if (this->bgzf_tbx == NULL) {
    this->has_index = false;
  } else {
    this->bgzf_itr = tbx_itr_queryi(
        this->bgzf_tbx, max(other.bgzf_itr->tid, other.bgzf_itr->curr_tid),
        max(other.bgzf_itr->beg, other.bgzf_itr->curr_end), HTS_POS_MAX);
    this->has_index = true;
  }
}

BgzReader& BgzReader::operator=(BgzReader other) {
  swap(*this, other);
  return *this;
}

void swap(BgzReader& first, BgzReader& second) {
  std::swap(first.filepath, second.filepath);
  std::swap(first.bgzf, second.bgzf);
  std::swap(first.bgzf_tbx, second.bgzf_tbx);
  std::swap(first.bgzf_itr, second.bgzf_itr);
  std::swap(first.seek_set, second.seek_set);
  std::swap(first.seek_buffer, second.seek_buffer);
  std::swap(first.buffer, second.buffer);
  std::swap(first.has_index, second.has_index);
  std::swap(first.line, second.line);
  std::swap(first.at_eof, second.at_eof);
  std::swap(first.is_closed, second.is_closed);
  std::swap(first.lines_read, second.lines_read);
}

void BgzReader::open(string filepath) {
  BgzReader other(filepath);
  swap(*this, other);
}

void BgzReader::close() {
  BgzReader other;
  swap(*this, other);
}

string BgzReader::readline() {
  if (this->eof()) {
    throw out_of_range("read past end of file " + this->filepath);
  } else if (this->closed()) {
    throw runtime_error("reading a closed file " + this->filepath);
  }

  if (this->seek_set) {
    line = this->seek_buffer;
    int ret = tbx_bgzf_itr_next(this->bgzf, this->bgzf_tbx, this->bgzf_itr,
                                &this->buffer);
    if (ret == -1) {
      this->seek_buffer = "";
    } else if (ret <= -2) {
      throw runtime_error("tbx_bgzf_itr_next error");
    } else {
      stringstream ss(buffer.s);
      getline(ss, this->seek_buffer);
    }
  } else {
    int ret = bgzf_getline(this->bgzf, '\n', &this->buffer);
    if (ret < 0) {
      throw runtime_error("bgzf_getline error " + to_string(ret) + " in file " + this->filepath);
    }
    line = string(this->buffer.s);
  }

  this->check_eof();
  ++this->lines_read;

  return line;
}

string BgzReader::read_bytes(ssize_t size) {
  if (this->closed()) {
    throw std::runtime_error("attempted to read closed file " + this->filepath);
  }
  if (this->at_eof) {
    throw std::runtime_error("read past end of file " + this->filepath);
  }

  char *data = new char[size];
  ssize_t bytes_read = bgzf_read(this->bgzf, data, size);
  if (bytes_read != size) {
    throw std::runtime_error("failed to read " + this->filepath);
  }
  string result(data, size);
  delete[] data;
  this->check_eof();
  return result;
}

void BgzReader::seek(string chrom, hts_pos_t position) {
  if (!this->has_index) {
    throw runtime_error("seek requires an index");
  } else if (this->closed()) {
    throw runtime_error("calling seek on a closed file");
  }

  int tid = tbx_name2id(this->bgzf_tbx, chrom.c_str());
  if (tid < 0) {
    throw runtime_error("called seek with invalid contig");
  }

  hts_itr_destroy(this->bgzf_itr);
  // hts iterators use 0-based indexing, but positions in Bgz use a 1-based index
  this->bgzf_itr = tbx_itr_queryi(this->bgzf_tbx, tid, position - 1, HTS_POS_MAX);
  int ret = tbx_bgzf_itr_next(this->bgzf, this->bgzf_tbx, this->bgzf_itr,
                              &this->buffer);
  if (ret == -1) {
    this->at_eof = true;
    this->seek_buffer = "";
  } else if (ret < -1) {
    throw runtime_error("tbx_bgzf_itr_next error");
  } else {
    stringstream ss(this->buffer.s);
    getline(ss, this->seek_buffer);
  }

  this->check_eof();
  this->seek_set = true;
  this->lines_read = 0;
}

void BgzReader::set_threads(int threads) {
  if (this->closed()) {
    throw runtime_error("calling set_threads on a closed file");
  }
  int ret = bgzf_mt(this->bgzf, threads, 0);
  if (ret < 0) {
    throw runtime_error("bgzf_mt error");
  }
}

int64_t BgzReader::tell() {
  if (this->closed()) {
    throw runtime_error("calling seek on a closed file");
  }
  return bgzf_tell(this->bgzf);
}

void BgzReader::seek_addr(int64_t addr) {
  if (this->closed()) {
    throw runtime_error("calling seek on a closed file");
  }

  if (bgzf_seek(this->bgzf, addr, SEEK_SET) == -1) {
    throw runtime_error("seek_addr failed on " + this->filepath);
  }
  this->check_eof();
}

void BgzReader::check_eof() {
  if (this->seek_set && this->seek_buffer == "") {
    this->at_eof = true;
  } else if (!this->seek_set && bgzf_peek(this->bgzf) == -1) {
    this->at_eof = true;
  } else {
    this->at_eof = false;
  }
}