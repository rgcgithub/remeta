#include "util.hpp"

#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

#include "logging.hpp"

namespace util {
  vector<string> str_split(const string& str, const string& delimiters) {
    vector<string> split;
    size_t beg = 0;
    size_t end = str.find_first_of(delimiters, beg);
    while (end != string::npos) {
      if (delimiters.find(str.substr(beg, 1)) == string::npos) {
        split.push_back(str.substr(beg, end - beg));
      }
      beg = end+1;
      end = str.find_first_of(delimiters, beg);
    }
    if (delimiters.find(str, beg) == string::npos) {
      split.push_back(str.substr(beg, end - beg));
    }
    return split;
  }

  bool is_cpra(const string& name) {
    return std::count(name.begin(), name.end(), ':') == 3;
  }

  bool file_exists(string filename) {
    if (FILE *file = fopen(filename.c_str(), "r")) {
      fclose(file);
      return true;
    } else {
      return false;
    }
  }

  string format_cohort_meta(const vector<string>& cohorts, vector<int> cohorts_to_include) {
    string cohort_meta = "";
    sort(cohorts_to_include.begin(), cohorts_to_include.end());
    for (size_t i = 0; i < cohorts_to_include.size(); ++i) {
      cohort_meta += cohorts.at(cohorts_to_include[i]);
      if (i < cohorts_to_include.size() - 1) {
        cohort_meta += "__";
      }
    }
    return cohort_meta;
  }

  vector<vector<htpv4_record_t> > get_htp_variants(const string& variant_file, const vector<string>& htp_files, const string& chr) {
    BgzReader var_reader(variant_file);
    set<string> variants_set;
    vector<string> variants_vec;
    string line;
    while (!var_reader.eof()) {
      line = var_reader.readline();
      variants_set.insert(line);
      variants_vec.push_back(line);
    }

    vector<vector<htpv4_record_t> > htp_entries(htp_files.size());
    htpv4_record_t rec;
    for (size_t i = 0; i < htp_files.size(); ++i) {;
      HTPv4Reader htp_reader(htp_files[i]);
      if (htp_reader.indexed()) {
        for ( const string& vid: variants_vec ) {
          vector<string> cpra = util::str_split(vid, ":");
          if (cpra.size() != 4) {
            log_warning("invalid variant ID " + vid);
            continue;
          }

          hts_pos_t pos = (hts_pos_t)stoi(cpra[1]);
          htp_reader.seek(cpra[0], pos);
          rec.pos = 0;
          while (rec.pos <= pos && !htp_reader.eof()) {
            rec = htp_reader.readrec();
            if (rec.name == vid && (chr == "" || rec.chr == chr)) {
              htp_entries[i].push_back(rec);
            }
          }
        }
      } else {
        while (!htp_reader.eof()) {
          rec = htp_reader.readrec();
          if (variants_set.count(rec.name) > 0 && (chr == "" || rec.chr == chr)) {
            htp_entries[i].push_back(rec);
          }
        }
      }
    }
    return htp_entries;
  }
}