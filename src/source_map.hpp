#ifndef SOURCE_MAP_H
#define SOURCE_MAP_H

#include <string>
#include <unordered_map>

class SourceMap {
 public:
  SourceMap();

  void load(const std::string& filepath);

  void clear();

  bool loaded();

  std::string get_source_short(const std::string& source);

  std::string get_category(const std::string& source);

 private:
  std::string filepath;
  bool is_loaded;
  std::unordered_map<std::string, std::string> source_map;
  std::unordered_map<std::string, std::string> category_map;
};

#endif