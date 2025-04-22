#ifndef REMETA_UTIL_H
#define REMETA_UTIL_H

#include <string>
#include <vector>

#include "io/htpv4_reader.hpp"

namespace util {
  std::vector<std::string> str_split(const std::string& str, const std::string& delimiters);
  bool is_cpra(const std::string& variant_name);
  bool file_exists(std::string filename);
  std::string format_cohort_meta(const std::vector<std::string>& cohorts, std::vector<int> cohorts_to_include);
  std::vector<std::vector<htpv4_record_t> > get_htp_variants(const std::string& variant_file, const std::vector<std::string>& htp_files, const string& chr = "");
}

#endif