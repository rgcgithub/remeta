#ifndef GCOORDINATE_H
#define GCOORDINATE_H

#include <stdexcept>
#include <string>
#include <vector>
using namespace std;

class HTPv4Pos {
 public:
  HTPv4Pos() : int_chr(-2), pos(-2) {}
  HTPv4Pos(int chr, int pos) : int_chr(chr), pos(pos) {}
  HTPv4Pos(string chr, int pos) : int_chr(chr_to_int(chr)), pos(pos) {}

  int chr_to_int(const string& chr);

  bool operator==(const HTPv4Pos& other) const {
    return int_chr == other.int_chr && pos == other.pos;
  }

  bool operator!=(const HTPv4Pos&  other) const {
    return !(*this == other);
  }

  bool operator<(const HTPv4Pos&  other) const {
    return int_chr < other.int_chr || (other.int_chr == int_chr && pos < other.pos);
  }

  bool operator<=(const HTPv4Pos&  other) const {
    return (*this < other) || (other == *this);
  }

  bool operator>(const HTPv4Pos&  other) const {
    return !(*this <= other);
  }

  bool operator>=(const HTPv4Pos& other) const {
    return !(*this < other);
  }

  int int_chr;
  int pos;
};

const HTPv4Pos HTPV4_EOF(-1, -1);
const HTPv4Pos HTPV4_NULL_POS(-2, -2);

pair<int, HTPv4Pos> get_min_pos(const vector<HTPv4Pos>& positions);

#endif