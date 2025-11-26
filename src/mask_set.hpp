#ifndef MASK_SET_H
#define MASK_SET_H

#include <htslib/khash.h>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include "io/bgz_reader.hpp"
#include "logging.hpp"

using namespace std;

constexpr double SINGLETON = -1;

class Mask {
 public:
  Mask(string name, double freq, vector<int> annotations)
   : name(name)
   , freq(freq)
   , annotations(annotations) {
    if (freq == SINGLETON) {
      alt = name + ".singleton";
    } else {
      alt = name + "." + format_freq(freq);
    }
  };

  bool contains(int annotation, double aaf, bool is_singleton) const {
    return (is_singleton && freq == SINGLETON && annotations.at(annotation) == 1)
      || (aaf <= freq && annotations.at(annotation) == 1);
  }

  double freq_bin() const {
    return freq;
  }

  vector<int> get_annotations() const { 
    return annotations;
  }

  string name;
  string alt;
 private:
  double freq;
  vector<int> annotations;

  string format_freq(double freq) const {
    string f = to_string(freq);
    while (f.back() == '0') {
      f.pop_back();
    }
    return f;
  }
};


class MaskSet {
 public:
  MaskSet(string mask_def_file, vector<double> freq_bins) {
    BgzReader reader(mask_def_file);
    string line;
    stringstream ss;
    vector<vector<string> > mask_annotations;
    vector<string> mask_names_vec;
    while (!reader.eof()) {
      line = reader.readline();
      ss = stringstream(line);
      string mask_name, anno_str;

      ss >> mask_name;
      if (mask_name.find(".") != string::npos) {
        log_error("in " + mask_def_file + ": mask names cannot contain '.'", 1);
      }

      mask_names_vec.push_back(mask_name);
      mask_names.insert(mask_name);

      ss >> anno_str;
      ss = stringstream(anno_str);
      string anno;
      vector<string> annotations;
      while(getline(ss, anno, ',')) {
        if (anno_int.find(anno) == anno_int.end()) {
          anno_int[anno] = anno_int.size();
        }
        annotations.push_back(anno);
      }
      mask_annotations.push_back(annotations);
    }

    for (size_t i = 0; i < mask_annotations.size(); ++i) {
      vector<int> annotations_int(anno_int.size(), 0);
      for (string annotation : mask_annotations[i]) {
        annotations_int[anno_int[annotation]] = 1;
      }
      for (double freq : freq_bins) {
        masks.push_back(Mask(mask_names_vec[i], freq, annotations_int));
      }
    }
  }

  MaskSet(vector<Mask> masks, unordered_map<string, int> anno_int)
   : masks(masks)
   , anno_int(anno_int) {
    for (const Mask& mask : masks) {
      mask_names.insert(mask.name);
    }
  }

  int n_masks() const { return masks.size(); }

  // TODO: returning a mask creates a copy, which is very slow
  // this function prevents accidentally copying a mask, but a
  // better solution would be to improve the mask data structure
  // to avoid an expensive copy
  bool mask_contains(int mask, string annotation, double aaf, bool is_singleton) const {
    return masks[mask].contains(anno_int.at(annotation), aaf, is_singleton);
  }

  bool mask_contains(int mask, int annotation, double aaf, bool is_singleton) const {
    return annotation != -1 && masks.at(mask).contains(annotation, aaf, is_singleton);
  }

  bool mask_name_contains(string mask_name) const {
    return mask_names.find(mask_name) != mask_names.end(); 
  }

  string get_mask_name(int mask) const {
    return masks.at(mask).name;
  }

  string get_mask_alt(int mask) const {
    return masks.at(mask).alt;
  }

  double get_mask_bin(int mask) const {
    return masks.at(mask).freq_bin();
  }

  vector<int> get_mask_annotations(int mask) const {
    return masks.at(mask).get_annotations();
  }

  int anno_to_int(string anno) const {
    if (anno_int.find(anno) == anno_int.end()) {
      return -1;
    } else {
      return anno_int.at(anno);
    }
  }

  unordered_map<string, int> get_anno_map() const {
    return anno_int;
  }

  vector<Mask> get_masks() const {
    return masks;
  }

 private:
  vector<Mask> masks;
  unordered_map<string, int> anno_int;
  std::unordered_set<std::string> mask_names;
};

#endif
