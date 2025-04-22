#include "source_map.hpp"

#include <string>
#include <unordered_map>
using std::string;
using std::unordered_map;

#include "io/bgz_reader.hpp"
#include "util.hpp"

SourceMap::SourceMap()
 : filepath("")
 , is_loaded(false) { }

void SourceMap::load(const string& filepath) {
  this->source_map.clear();

  BgzReader reader;
  reader.open(filepath);
  string line;
  vector<string> split;
  while (!reader.eof()) {
    line = reader.readline();
    split = util::str_split(line, "\t ");
    this->source_map[split[0]] = split[1];

    if (split.size() >= 3) {
      this->category_map[split[0]] = split[2];
    }
  }
  this->is_loaded = true;
}

void SourceMap::clear() {
  this->source_map.clear();
  this->is_loaded = false;
}

bool SourceMap::loaded() {
  return this->is_loaded;
}

string SourceMap::get_source_short(const string& source) {
  if (this->source_map.count(source) > 0) {
    return this->source_map.at(source);
  } else {
    return source;
  }
}

string SourceMap::get_category(const string& source) {
  if (this->category_map.count(source) > 0) {
    return this->category_map.at(source);
  } else {
    return "";
  }
}