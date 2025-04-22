#include "annotation_map.hpp"

AnnotationMap::AnnotationMap(const string& filepath, MaskSet mask_set) 
 : filepath(filepath)
 , anno_reader()
 , buffer_chrom("")
 , buffer_start(0)
 , buffer_end(0)
 , mask_set(mask_set)
 , has_next_rec(false) {
  this->anno_reader.open(filepath);
  // if (!this->anno_reader.indexed()) {
  //   throw runtime_error("AnnotationMap requires " + filepath + " to be bgzipped and indexed with remeta index-anno.");
  // }
}


int AnnotationMap::get_annotation(const string& cpra, const string& gene) {
  size_t i = cpra.find(":");
  size_t j = cpra.find(":", i+1);
  if (i == string::npos || j == string::npos) {
    return -1;
  }

  string chrom = cpra.substr(0, i);
  int pos = stoi(cpra.substr(i+1, j-i-1));

  if (this->anno_reader.indexed()) {
    this->load_chunk(chrom, pos);
  } else {
    this->load_chrom(chrom);
  }

  for (size_t i = 0; i < anno_buffer[cpra].size(); i++) {
    if (anno_buffer[cpra][i].gene == gene) {
      return anno_buffer[cpra][i].annotation;
    }
  }
  return -1;
}


int AnnotationMap::get_annotation(const string& name, const string& gene, const string& chrom, int pos) {
  if (this->anno_reader.indexed()) {
    this->load_chunk(chrom, pos);
  } else {
    this->load_chrom(chrom);
  }

  for (size_t i = 0; i < anno_buffer[name].size(); i++) {
    if (anno_buffer[name][i].gene == gene) {
      return anno_buffer[name][i].annotation;
    }
  }
  return -1;
}


// when a file is not indexed, load one chromosome at a time
void AnnotationMap::load_chrom(const string& chrom) {
  if (chrom == this->buffer_chrom) return;
  this->anno_reader.close();
  this->anno_reader.open(this->filepath);
  this->anno_buffer = unordered_map<string, vector<annomap_rec_t> >();
  this->buffer_chrom = chrom;

  // fill the buffer with the next chunk
  while (!this->anno_reader.eof()) {
    annorec_t annorec = this->anno_reader.readrec();
    annomap_rec_t annomap_rec;
    annomap_rec.name = annorec.name;
    annomap_rec.chrom = annorec.chrom;
    annomap_rec.pos = annorec.pos;
    annomap_rec.gene = annorec.gene;
    annomap_rec.annotation = this->mask_set.anno_to_int(annorec.annotation);

    if (annomap_rec.chrom == chrom) {
      this->anno_buffer[annomap_rec.name].push_back(annomap_rec);
    }
  }
}


void AnnotationMap::load_chunk(const string& chrom, int pos) {
  if (this->in_current_chunk(chrom, pos)) {
    return ;
  }

  if (!this->in_next_chunk(chrom, pos)) {
    this->anno_reader.seek(chrom, pos);
    this->buffer_chrom = chrom;
    this->buffer_start = pos;
    this->buffer_end = pos + 1e6*BUFFER_SIZE_MB;
  } else {
    this->buffer_chrom = chrom;
    this->buffer_start = this->buffer_end + 1;
    this->buffer_end = this->buffer_start + 1e6*BUFFER_SIZE_MB;
  }
  // clear the old buffer
  this->anno_buffer = unordered_map<string, vector<annomap_rec_t> >();

  // add any upcoming variants we have already read
  if (this->has_next_rec) {
    this->anno_buffer[this->next_rec.name].push_back(this->next_rec);
    this->has_next_rec = false;
  }

  // fill the buffer with the next chunk
  while (!this->anno_reader.eof()) {
    annorec_t annorec = this->anno_reader.readrec();
    annomap_rec_t annomap_rec;
    annomap_rec.name = annorec.name;
    annomap_rec.chrom = annorec.chrom;
    annomap_rec.pos = annorec.pos;
    annomap_rec.gene = annorec.gene;
    annomap_rec.annotation = this->mask_set.anno_to_int(annorec.annotation);
    this->anno_buffer[annomap_rec.name].push_back(annomap_rec);

    // if we read past the buffer, store the variant to add
    // when loading the next chunk
    if (annomap_rec.pos > this->buffer_end) {
      this->next_rec = annomap_rec;
      this->has_next_rec = true;
      break ;
    }
    if (annomap_rec.chrom != this->buffer_chrom) {
      this->next_rec = annomap_rec;
      this->has_next_rec = true;
      break ;
    }
  }
}


bool AnnotationMap::in_current_chunk(const string& chrom, int pos) {
  return chrom == buffer_chrom && pos >= buffer_start && pos <= buffer_end;
}


bool AnnotationMap::in_next_chunk(const string& chrom, int pos) {
  return chrom == buffer_chrom
                  && (pos >= buffer_start)
                  && (pos < buffer_end + 1e6*BUFFER_SIZE_MB);
}