#include "htpv4_pos.hpp"

int HTPv4Pos::chr_to_int(const string& chr) {
  string chr_only = chr;
  if (chr.find("chr") != string::npos) {
    chr_only = chr.substr(3, chr.length());
  }

  if (chr_only == "X") {
    return 23;
  } else if (chr_only == "y") {
    return 24;
  } else {
    int t = stoi(chr_only);
    if (t < 1 || t > 23) {
      throw runtime_error("unable to convert chromosome " + 
        chr + " (" + chr_only + ") "
        " to integer index");
    }
    return t;
  }
}

pair<int, HTPv4Pos> get_min_pos(const vector<HTPv4Pos>& positions) {
  HTPv4Pos min = HTPV4_EOF;
  int index = -1;
  for (size_t i = 0; i < positions.size(); ++i) {
    if (min == HTPV4_EOF && positions[i] != HTPV4_EOF) {
      index = i;
      min = positions[i];
    } else if (positions[i] != HTPV4_EOF && positions[i] < min) {
      index = i;
      min = positions[i];
    }
  }
  return make_pair(index, min);
}