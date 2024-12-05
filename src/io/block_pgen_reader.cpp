#include "block_pgen_reader.hpp"

#include <cmath>
#include <vector>
using namespace std;

#include <Eigen/Dense>
using Eigen::MatrixXd;
typedef Eigen::Array<bool,Eigen::Dynamic,1> ArrayXb;

#include <boost/algorithm/string.hpp>

#include "../util.hpp"

void BlockPgenReader::load(const string& pfile, const int& nthreads, const int& pvar_index_stride) {
  this->pfile = pfile;
  string pgen = pfile + ".pgen";
  string pvar = pfile + ".pvar";
  string psam = pfile + ".psam";

#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#endif 

  this->idx_stride = pvar_index_stride;
  this->load_psam(psam);
  this->load_pvar(pvar);

  vector<int> sample_subset_1based;
  string sample_id;
  size_t keep_list_size = this->keep_list.size();
  size_t remove_list_size = this->remove_list.size();
  for (size_t i = 1; i < this->samples.size() + 1; ++i) {
    sample_id = this->format_sample_id(this->samples[i-1]);
    if (keep_list_size > 0 && this->keep_list.count(sample_id) > 0) {
      sample_subset_1based.push_back(i);
    } else if (remove_list_size > 0 && this->remove_list.count(sample_id) == 0) {
      this->keep_list.insert(sample_id);
      sample_subset_1based.push_back(i);
    } else if (keep_list_size == 0 && remove_list_size == 0) {
      this->keep_list.insert(sample_id);
      sample_subset_1based.push_back(i);
    }
  }
  if (sample_subset_1based.size() != this->keep_list.size()) {
    throw runtime_error(
      "keep list specified " + to_string(this->keep_list.size())
      + " samples, but found "
      + to_string(sample_subset_1based.size())
      + " (if FID and IID are specified in the psam, did you use both in your keep list?)"
    );
  }

  if (sample_subset_1based.size() == 0) {
    throw runtime_error(
      "no samples remaining after filters"
    );
  }

  pgen_reader.Load(pgen, this->samples.size(), sample_subset_1based, nthreads);

  if (pgen_reader.GetRawSampleCt() != this->samples.size()) {
    throw runtime_error("sample number mismatch between psam and pgen files");
  }
  if (pgen_reader.GetVariantCt() != this->variants.size()) {
    throw runtime_error("variant number mismatch between pvar and pgen files");
  }
  if (pgen_reader.GetMaxAlleleCt() > 2) {
    throw runtime_error("pgen variants must be bi-allelic");
  }
}

void BlockPgenReader::load(const string& pfile,
                           const int& nthreads,
                           const int& pvar_index_stride,
                           const string& keep_file,
                           const string& remove_file) {
  if (keep_file != "" && remove_file != "") {
    throw runtime_error("cannot set both keep file and remove file");
  }
  if (keep_file != "") {
    this->load_keep_list(keep_file);
  } else if (remove_file != "") {
    this->load_remove_list(remove_file);
  }
  this->load(pfile, nthreads, pvar_index_stride);
}

void BlockPgenReader::load(const string& pfile,
                           const int& nthreads,
                           const int& pvar_index_stride,
                           const unordered_set<string>& keep_list,
                           const unordered_set<string>& remove_list) {
  this->keep_list = keep_list;
  this->remove_list = remove_list;
  this->load(pfile, nthreads, pvar_index_stride);
}

int BlockPgenReader::get_nblocks(const int& block_size) {
  if (this->read_mode == REGION) {
    return this->get_nblocks(block_size, this->region);
  } else {
    return ceil((double)this->variants.size() / (double)block_size);
  }
}

int BlockPgenReader::get_nblocks(const int& block_size, const block_pgen_region_t& region) {
  return ceil((double)(region.end_idx - region.start_idx + 1) / (double)block_size);
};

int BlockPgenReader::get_nsamples() {
  return (int)this->samples.size();
}

int BlockPgenReader::get_nkeep_samples() {
  return (int)this->keep_list.size();
}

int BlockPgenReader::get_nvariants() {
  if (this->read_mode == REGION) {
    return this->region.end_idx - this->region.start_idx + 1;
  } else if (this->extract_idx.size() > 0) {
    return (int)this->extract_idx.size();
  } else {
    return (int)this->variants.size();
  }
}

block_pgen_region_t BlockPgenReader::get_region(const string& chr, const int& start, const int& end) {
  if (start < 0 || end < 0 || end < start || this->chr_pos_idx.count(chr) == 0) {
    throw runtime_error(
      "bad call to BlockPgenReader::set_region (invalid region "
      + chr + ":" + to_string(start) + "-" + to_string(end) +  ")"
    );
  }
  block_pgen_region_t region {chr, start, end, 0, -1};

  if (start > this->chr_start_end[chr].second) {
    return region;
  }

  // first find the position in the index
  int pos = start - (start % this->idx_stride);
  int start_idx = 0;
  while (pos >= this->chr_start_end[chr].first && this->chr_pos_idx[chr].count(pos) == 0) {
    pos -= this->idx_stride;
  }
  if (this->chr_pos_idx[chr].count(pos) > 0) {
    start_idx = this->chr_pos_idx[chr][pos];
  } else {
    pos = this->chr_start_end[chr].first 
          - (this->chr_start_end[chr].first % this->idx_stride);
    start_idx = this->chr_pos_idx[chr][pos];
  }
  // then scan to the exact variant in variants
  while (start_idx < (int)this->variants.size()
         && this->variants[start_idx].chr == chr
         && this->variants[start_idx].pos < start) {
    ++start_idx;
  }

  // find the position in the index
  int end_idx = 0;
  pos = end - (end % this->idx_stride) + this->idx_stride;
  while (pos < this->chr_start_end[chr].second && this->chr_pos_idx[chr].count(pos) == 0) {
    pos += this->idx_stride;
  }

  if (this->chr_pos_idx[chr].count(pos) > 0) {
    end_idx = this->chr_pos_idx[chr][pos];
  } else {
    pos = this->chr_start_end[chr].second 
          - (this->chr_start_end[chr].second % this->idx_stride) + this->idx_stride;
    end_idx = this->chr_pos_idx[chr][pos];
  }
  // then scan to the exact variant in variants
  while (end_idx >= start_idx
         && this->variants[end_idx].chr == chr
         && this->variants[end_idx].pos > end) {
    --end_idx;
  }

  region.chr = chr;
  region.start_bp = start;
  region.end_bp = end;
  region.start_idx = start_idx;
  region.end_idx = end_idx;
  return region;
}

void BlockPgenReader::set_region(const string& chr, const int& start, const int& end) {
  block_pgen_region_t region = this->get_region(chr, start, end);
  this->read_mode = REGION;
  this->region = region;
}

void BlockPgenReader::set_extract_list(const unordered_set<string>& variant_ids) {
  this->extract_idx.clear();
  for ( const string& vid: variant_ids ) {
    this->extract_idx.insert(this->variant_idx[vid]);
  }
  this->variants.clear();
  this->load_pvar(this->pfile + ".pvar");
}

void BlockPgenReader::set_exclude_list(const unordered_set<string>& variant_ids) {
  // prioritize extract list over exclude list
  if (this->extract_idx.size() > 0) {
    return;
  }
  this->extract_idx.clear();
  for (size_t i = 0; i < this->variants.size(); ++i) {
    if (variant_ids.count(this->variants[i].id) == 0) {
      this->extract_idx.insert(this->variant_idx[this->variants[i].id]);
    }
  }
  this->variants.clear();
  this->load_pvar(this->pfile + ".pvar");
}

void BlockPgenReader::set_read_dosages() {
  if (!pgen_reader.DosagePresent()) {
    throw runtime_error("requested dosages but dosage data not present in pgen file");
  }
  this->read_dosages = true;
}

void BlockPgenReader::load_block(Eigen::MatrixXf& genotype_block,
                                 vector<pvar_variant_t>& variants,
                                 int& nvar_read,
                                 const int& block_index,
                                 const int& block_size) {
  if (genotype_block.rows() != (ssize_t)this->get_nkeep_samples()) {
    throw runtime_error("genotype block dimension mismatch");
  }

  int start_idx;
  int end_idx;
  if (this->read_mode == REGION) {
    start_idx = block_size * block_index + this->region.start_idx;
    end_idx = min(block_size * (block_index + 1) + this->region.start_idx, this->region.end_idx + 1);
  } else {
    start_idx = block_size * block_index;
    end_idx = min(block_size * (block_index + 1), (int)this->variants.size());
  }

  nvar_read = 0;
  variants.clear();
  variants.resize(end_idx - start_idx);

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = start_idx; i < end_idx; ++i) {
    int thread_num = 0;
    int variant_idx = this->variant_idx[this->variants[i].id];
#if defined(_OPENMP)
    thread_num = omp_get_thread_num();
#endif
    if (i - start_idx >= genotype_block.cols() || i - start_idx < 0) {
      throw runtime_error(
        "BlockPgenReader bad call to load_block (" 
         + to_string(i - start_idx) + "): this is a bug"
      );
    }

    if (this->read_dosages) {
      pgen_reader.ReadFloat(genotype_block.col(i - start_idx).data(),
                            this->get_nkeep_samples(),
                            thread_num,
                            variant_idx,
                            1);
    } else {
      pgen_reader.ReadFloatHardcalls(genotype_block.col(i - start_idx).data(),
                                    this->get_nkeep_samples(),
                                    thread_num,
                                    variant_idx,
                                    1);
    }

    ArrayXb missing = (genotype_block.col(i - start_idx).array() < 0) || (genotype_block.col(i - start_idx).array() > 2);
    double mu = missing.select(
      0,
      genotype_block.col(i - start_idx)
    ).sum() / this->get_nkeep_samples();
    for (ssize_t j = 0; j < genotype_block.rows(); ++j) {
      if (missing(j)) {
        genotype_block(j, i - start_idx) = mu;
      }
    }

    variants[i - start_idx] = this->variants[i];
  }
  //genotype_block = genotype_block.cwiseMin(2).cwiseMax(0);
  nvar_read = end_idx - start_idx;
}

void BlockPgenReader::load_block(Eigen::MatrixXf& genotype_block,
                                 vector<pvar_variant_t>& variants,
                                 int& nvar_read,
                                 const int& block_index,
                                 const int& block_size,
                                 const block_pgen_region_t& region) {
  read_mode_e current_read_mode = this->read_mode;
  block_pgen_region_t current_region = this->region;

  this->read_mode = REGION;
  this->region = region;
  this->load_block(genotype_block, variants, nvar_read, block_index, block_size);
  this->read_mode = current_read_mode;
  this->region = current_region;
}

bool BlockPgenReader::region_contains_variant(const string& vid) {
  if (this->variant_idx.count(vid) == 0 || (this->extract_idx.size() > 0 && this->extract_idx.count(this->variant_idx[vid]) == 0)) {
    return false;
  } else {
    pvar_variant_t var = this->variants[this->variant_idx[vid]];
    return var.chr == this->region.chr
           && var.pos >= this->region.start_bp
           && var.pos <= this->region.end_bp;
  }
}

void BlockPgenReader::load_pvar(const string& pvar) {
  this->chr_pos_idx.clear();
  this->chr_start_end.clear();

  BgzReader pvar_reader(pvar);
  string line;
  vector<string> row;
  string last_chr = "";
  int last_pos = 0;
  while (!pvar_reader.eof()) {
    line = pvar_reader.readline();
    if (line.substr(0, 2) == "##") {
      continue;
    }
    if (line == "") {
      break;
    }

    row = util::str_split(line, "\t ");
    if (row[0] == "#CHROM") {
      if (row.size() < 5
          || row[1] != "POS"
          || row[2] != "ID"
          || row[3] != "REF"
          || row[4] != "ALT") {
        throw runtime_error("improperly formated pvar file: expected columns #CHROM POS ID REF ALT");
      }
      continue;
    }

    pvar_variant_t variant {row[0], stoi(row[1]), row[2], (int)this->variants.size()};
    if (this->extract_idx.size() > 0 && this->extract_idx.count(this->variant_idx[variant.id]) != 0) {
      variant.index = this->variant_idx[variant.id];
      this->variants.push_back(variant);
    } else if (this->extract_idx.size() > 0 && this->extract_idx.count(this->variant_idx[variant.id]) == 0) {
      continue;
    } else if (this->extract_idx.size() == 0) {
      this->variant_idx[variant.id] = (int)this->variants.size();
      this->variants.push_back(variant);
    }

    if (last_chr == variant.chr && last_pos > variant.pos) {
      throw runtime_error("pfile variants must be sorted by position");
    }
    if (last_chr != variant.chr) {
      int last_index_bp = last_pos - (last_pos % this->idx_stride) + this->idx_stride;
      this->chr_pos_idx[last_chr][last_index_bp] = variants.size() - 1;
      this->chr_start_end[variant.chr] = make_pair(variant.pos, 0);
    }
    this->chr_start_end[variant.chr].second = max(variant.pos, this->chr_start_end[variant.chr].second);
    last_chr = variant.chr;
    last_pos = variant.pos;

    int index_bp = variant.pos - (variant.pos % this->idx_stride);
    if (this->chr_pos_idx.count(row[0]) == 0) {
      this->chr_pos_idx[variant.chr][index_bp] = variants.size() - 1;
    } else if (this->chr_pos_idx[variant.chr].count(index_bp) == 0) {
      this->chr_pos_idx[variant.chr][index_bp] = variants.size() - 1;
    }
  }

  int last_index_bp = last_pos - (last_pos % this->idx_stride) + this->idx_stride;
  this->chr_pos_idx[last_chr][last_index_bp] = variants.size() - 1;
}

void BlockPgenReader::load_psam(const string& psam) {
  BgzReader psam_reader(psam);
  string line;
  vector<string> row;

  // parse header
  int fid_col = -1;
  int iid_col = -1;
  int sex_col = -1;
  line = psam_reader.readline();
  row = util::str_split(line, "\t ");
  for (size_t i = 0; i < row.size(); ++i) {
    if (row[i] == "#FID" || row[i] == "FID") {
      fid_col = i;
    } else if (row[i] == "#IID" || row[i] == "IID") {
      iid_col = i;
    } else if (row[i] == "SEX") {
      sex_col = i;
    }
  }

  string fid;
  string iid;
  int8_t sex = -1;
  while (!psam_reader.eof()) {
    line = psam_reader.readline();
    if (line == "") {
      break;
    }

    row = util::str_split(line, "\t ");
    if (fid_col != -1) {
      fid = row[fid_col];
    } else {
      fid = "";
    }

    if (iid_col != -1) {
      iid = row[iid_col];
    } else {
      iid = "";
    }

    if (sex_col != -1 && row[sex_col] != "NA") {
      sex = (int8_t)stoi(row[sex_col]);
    } else {
      sex = -1;
    }

    this->samples.push_back({fid, iid, sex});
  }
}

void BlockPgenReader::load_keep_list(const string& keep_file) {
  BgzReader reader(keep_file);
  string line;
  vector<string> row;
  while (!reader.eof()) {
    line = reader.readline();
    row = util::str_split(line, "\t ");
    if (row.size() == 1) {
      this->keep_list.insert(row[0]);
    } else if (row.size() == 2) {
      this->keep_list.insert(row[0] + "_" + row[1]);
    } else {
      throw runtime_error(
        "found too many columns in keep file (expected <= 2, found "
        + to_string(row.size()) + ")");
    }
  }
}

void BlockPgenReader::load_remove_list(const string& remove_file) {
  BgzReader reader(remove_file);
  string line;
  vector<string> row;
  while (!reader.eof()) {
    line = reader.readline();
    row = util::str_split(line, "\t ");
    if (row.size() == 1) {
      this->remove_list.insert(row[0]);
    } else if (row.size() == 2) {
      this->remove_list.insert(row[0] + "_" + row[1]);
    } else {
      throw runtime_error(
        "found too many columns in remove file (expected <= 2, found "
        + to_string(row.size()) + ")");
    }
  }
}

string BlockPgenReader::format_sample_id(const psam_sample_t& sample) {
  if (sample.fid != "") {
    return sample.fid + "_" + sample.iid;
  } else {
    return sample.iid;
  }
}
