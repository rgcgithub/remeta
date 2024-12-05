#include "ref_ld_matrix_writer.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
using namespace std;

#include "../logging.hpp"

RefLDMatrixWriter::RefLDMatrixWriter(const string& file_prefix, const int& float_size) {
  this->is_closed = false;
  this->gene_ld_file.open(file_prefix + ".remeta.gene.ld", "w");
  this->buffer_ld_file.open(file_prefix + ".remeta.buffer.ld", "w");
  this->idx.open(file_prefix + ".remeta.ld.idx.gz", "w");
  this->float_size = float_size;

  this->gene_ld_file.write(string("remetaLD.v1.1"));
  this->buffer_ld_file.write(string("remetaLD.v1.1"));
  this->buffer_ld_file.write<int32_t>(this->float_size);

  if (float_size == 4 && sizeof(float) != 4) {
    throw runtime_error("bad float size: sizeof(float) != 4");
  }
}

RefLDMatrixWriter::~RefLDMatrixWriter() {
  if (!this->closed()) {
    this->gene_ld_file.close();
    this->buffer_ld_file.close();
    this->idx.close();
  }
}

void RefLDMatrixWriter::close() {
  if (this->is_closed) {
    return ;
  }

  this->gene_ld_file.close();
  this->buffer_ld_file.close();
  this->idx.close();
  this->is_closed = true;
}

void RefLDMatrixWriter::check_closed() {
  if (this->closed()) {
    throw runtime_error("RefLDMatrixWriter is attempting to write to closed files");
  }
}

void RefLDMatrixWriter::init_next_gene(const string& gene_name,
                                       const float& gene_sparsity_threshold,
                                       const float& buffer_sparsity_threshold,
                                       const vector<variant_id>& gene_variant_ids,
                                       const vector<float>& gene_variant_variances,
                                       const vector<variant_id>& buffer_variant_ids,
                                       const vector<float>& buffer_variant_variances,
                                       const int32_t& nbuffer_blocks) {
  this->check_closed();

  this->gene_gene_ld.nvariants = (int32_t)gene_variant_ids.size();
  this->gene_gene_ld.nentries = 0;
  this->gene_gene_ld.sparsity_threshold = gene_sparsity_threshold;
  this->gene_gene_ld.variances = gene_variant_variances;
  this->gene_gene_ld.data.clear();

  this->gene_buffer_ld.ngene_variants = (int32_t)gene_variant_ids.size();
  this->gene_buffer_ld.nbuffer_variants = (int32_t)buffer_variant_ids.size();
  this->gene_buffer_ld.nentries = 0;
  this->gene_buffer_ld.sparisty_threshold = buffer_sparsity_threshold;
  this->gene_buffer_ld.variances.clear();
  this->gene_buffer_ld.buffer_idx_map.clear();
  this->gene_buffer_ld.nbuffer_blocks = nbuffer_blocks;
  this->gene_buffer_ld.data.clear();
  for (const float& v : buffer_variant_variances) {
    this->gene_buffer_ld.variances.push_back(v);
  }

  this->idx_entries[gene_name].gene = gene_name;
  this->idx_entries[gene_name].gene_ld_addr = gene_ld_file.tell();
  this->idx_entries[gene_name].buffer_ld_addr = 0;
  this->idx_entries[gene_name].gene_variants = gene_variant_ids;
  this->idx_entries[gene_name].buffer_variants = buffer_variant_ids;

  for (const variant_id& id : buffer_variant_ids) {
    if (this->buffer_ld_variant_idx.count(id) == 0) {
      this->buffer_ld_variant_idx[id] = (int32_t)this->buffer_ld_variant_idx.size();
    }
    this->gene_buffer_ld.buffer_idx_map.push_back(this->buffer_ld_variant_idx[id]);
  }

  this->gene_variant_idx.clear();
  for (size_t i = 0; i < gene_variant_ids.size(); ++i) {
    this->gene_variant_idx[gene_variant_ids[i]] = (int32_t)i;
  }

  this->buffer_variant_idx.clear();
  for (size_t i = 0; i < buffer_variant_ids.size(); ++i) {
    if (this->gene_variant_idx.count(buffer_variant_ids[i]) > 0) {
      throw runtime_error(
        "duplicate variant " + buffer_variant_ids[i] +
        " found in gene variants and buffer variants: this is a bug"
      );
    }
    this->buffer_variant_idx[buffer_variant_ids[i]] = i;
  }
}

void RefLDMatrixWriter::add_gene_gene_ld_entry(const ref_ld_matrix_entry_t& data) {
  this->check_closed();

  ++this->gene_gene_ld.nentries;
  this->gene_gene_ld.data.push_back(
    matrix_entry_data_t {
      this->gene_variant_idx.at(data.row_id),
      this->gene_variant_idx.at(data.col_id),
      data.corr
    }
  );
}

void RefLDMatrixWriter::add_gene_buffer_ld_entry(const ref_ld_matrix_entry_t& data) {
  this->check_closed();

  int32_t row_id = -1;
  int32_t col_id = -1;

  if (this->gene_variant_idx.count(data.row_id) > 0) {
    row_id = this->gene_variant_idx[data.row_id];
    if (this->buffer_ld_variant_idx.count(data.col_id) == 0) {
      throw runtime_error("missing variant " + data.col_id + " from buffer variants");
    } else {
      col_id = this->buffer_ld_variant_idx[data.col_id];
    }
  } else if (this->buffer_variant_idx.count(data.row_id) > 0) {
    col_id = this->buffer_variant_idx[data.row_id];
    if (this->gene_variant_idx.count(data.col_id) == 0) {
      throw runtime_error("missing variant " + data.row_id + " from gene variants");
    } else {
      row_id = this->gene_variant_idx[data.col_id];
    }
  } else {
    throw runtime_error(
      "missing variants " + data.row_id + "," + data.col_id
      + " from current gene in RefLDMatrixWriter"  
    );
  }

  ++this->gene_buffer_ld.nentries;
  this->gene_buffer_ld.data.push_back(
    matrix_entry_data_t {
      row_id,
      col_id,
      data.corr
    }
  );
}

void RefLDMatrixWriter::write_gene_to_gene_ld_file() {
  this->check_closed();
  // this order matters - gene_gene_ld_t must be written to the LD file first
  this->write_gene_gene_ld();
  this->write_gene_buffer_ld();
}

void RefLDMatrixWriter::write_gene_to_index(const string& gene_name,
                                            const int64_t& first_block_block_ld_addr) {
  this->check_closed();
  if (this->idx_entries.count(gene_name) == 0) {
    throw runtime_error("RefLDMatrixWriter bad call to write_gene_to_index");
  }
  idx_entry_t idx_entry = this->idx_entries[gene_name];
  idx_entry.buffer_ld_addr = first_block_block_ld_addr;

  this->idx.write(
    idx_entry.gene + "\t" 
    + to_string(idx_entry.gene_ld_addr) + "\t"
    + to_string(idx_entry.buffer_ld_addr) + "\t"
  );

  for (size_t i = 0; i < idx_entry.gene_variants.size(); ++i) {
    this->idx.write(idx_entry.gene_variants[i]);
    if (i+1 < idx_entry.gene_variants.size()) {
      this->idx.write<char>(',');
    }
  }
  this->idx.write<char>('\t');

  for (size_t i = 0; i < idx_entry.buffer_variants.size(); ++i) {
    this->idx.write(idx_entry.buffer_variants[i]);
    if (i+1 < idx_entry.buffer_variants.size()) {
      this->idx.write<char>(',');
    }
  }
  this->idx.write<char>('\n');

  this->idx_entries.erase(gene_name);
}

void RefLDMatrixWriter::write_gene_gene_ld() {
  this->gene_ld_file.write<int32_t>(this->gene_gene_ld.nvariants);
  this->gene_ld_file.write<int32_t>(this->gene_gene_ld.nentries);
  this->gene_ld_file.write<int32_t>(this->gene_gene_ld.sparsity_threshold);
  for (const float& v : this->gene_gene_ld.variances) {
    this->gene_ld_file.write<float>(v);
  }
  for (const matrix_entry_data_t entry : this->gene_gene_ld.data) {
    this->gene_ld_file.write<matrix_entry_data_t>(entry);
  }
}

void RefLDMatrixWriter::write_gene_buffer_ld() {
  this->gene_ld_file.write<int32_t>(this->gene_buffer_ld.ngene_variants);
  this->gene_ld_file.write<int32_t>(this->gene_buffer_ld.nbuffer_variants);
  this->gene_ld_file.write<int32_t>(this->gene_buffer_ld.nentries);
  this->gene_ld_file.write<float>(this->gene_buffer_ld.sparisty_threshold);
  for (const float& v : this->gene_buffer_ld.variances) {
    this->gene_ld_file.write<float>(v);
  }
  for (const int32_t idx : this->gene_buffer_ld.buffer_idx_map) {
    this->gene_ld_file.write<int32_t>(idx);
  }
  this->gene_ld_file.write<int32_t>(this->gene_buffer_ld.nbuffer_blocks);
  for (const matrix_entry_data_t entry : this->gene_buffer_ld.data) {
    this->gene_ld_file.write<int32_t>(entry.row_idx);
    this->gene_ld_file.write<int32_t>(entry.col_idx);
    if (this->float_size == 1) {
      int8_t corr = entry.corr == 1 ? 0 : static_cast<int8_t>(round(entry.corr / pow_2_m7));
      this->gene_ld_file.write<int8_t>(corr);
    } else if (this->float_size == 2) {
      int16_t corr = entry.corr == 1 ? 0 : static_cast<int16_t>(round(entry.corr / pow_2_m15));
      this->gene_ld_file.write<int16_t>(corr);
    } else if (this->float_size == 4) {
      this->gene_ld_file.write<float>(entry.corr);
    } else {
      throw runtime_error("bad float size");
    }
  }
}

int64_t RefLDMatrixWriter::write_block_block_ld_to_buffer_ld_file(const block_block_ld_t& block_block_ld,
                                                                  const bool& within_block_ld) {
  int64_t addr = -1;
  if (within_block_ld) {
    this->buffer_ld_file.write<int32_t>(BLOCK_END_MARKER);
    addr = this->buffer_ld_file.tell();
    this->buffer_ld_file.write<int32_t>(BLOCK_START_MARKER);
  }

  std::unordered_map<int32_t, std::vector<matrix_row_entry_data_t> > block_block_ld_buffer;
  for (const ref_ld_matrix_entry_t& ref_entry : block_block_ld.data) {
    if (this->buffer_ld_variant_idx.count(ref_entry.row_id) == 0) {
      this->buffer_ld_variant_idx[ref_entry.row_id] = 
        (int32_t)this->buffer_ld_variant_idx.size();
    }
    if (this->buffer_ld_variant_idx.count(ref_entry.col_id) == 0) {
      this->buffer_ld_variant_idx[ref_entry.col_id] =
        (int32_t)this->buffer_ld_variant_idx.size();
    }

    int32_t row_idx = min(this->buffer_ld_variant_idx[ref_entry.row_id], this->buffer_ld_variant_idx[ref_entry.col_id]);
    int32_t col_idx = max(this->buffer_ld_variant_idx[ref_entry.row_id], this->buffer_ld_variant_idx[ref_entry.col_id]);
    block_block_ld_buffer[row_idx].push_back(
      matrix_row_entry_data_t {
        col_idx,
        ref_entry.corr
      }
    );
  }

  int32_t row;
  for ( const auto& it : block_block_ld_buffer ) {
    row = it.first;
    this->buffer_ld_file.write<int32_t>(ROW_START_MARKER);
    this->buffer_ld_file.write<int32_t>((int32_t)it.second.size());
    this->buffer_ld_file.write<int32_t>(row);
    for ( const matrix_row_entry_data_t& e : it.second ) {
      if (e.corr == 0) {
        throw runtime_error("found 0 correlation in block_block_ld: this is a bug");
      }
      this->buffer_ld_file.write<int32_t>(e.col_idx);
      if (this->float_size == 1) {
        int8_t corr = e.corr == 1 ? 0 : static_cast<int8_t>(round(e.corr / pow_2_m7));
        this->buffer_ld_file.write<int8_t>(corr);
      } else if (this->float_size == 2) {
        int16_t corr = e.corr == 1 ? 0 : static_cast<int16_t>(round(e.corr / pow_2_m15));
        this->buffer_ld_file.write<int16_t>(corr);
      } else if (this->float_size == 4) {
        this->buffer_ld_file.write<float>(e.corr);
      } else {
        throw runtime_error("bad float size");
      }
    }
  }

  return addr;
}