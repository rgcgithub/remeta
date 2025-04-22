/* htpv4_reader.hpp
 * Author: Tyler Joseph
 *
 * Read a raw, gzipped, or bgzipped HTPv4 file using HTSlib. Can rapidly query
 * tabixed HTPv4 files. Skips header lines starting with # and the line with
 * column names.
 *
 * Inherits from the BgzReader class.
 *
 * To index an HTPv4 file:
 *   $ bgzip myfile
 *   $ tabix -s 2 -b 3 -e 3 -S 1 myfile.gz
 *
 * Example:
 *   // load an HTPv4 file
 *   HTPv4Reader reader = HTPvReader("myfile.gz");
 *
 *   // read lines or records
 *   string line = reader.readline();
 *   htpv4_record_t rec = reader.readrec();
 *
 *   // if an index exists, skip to a region
 *   // otherwise throws a runtime_error
 *   reader.seek("22", 15528158));
 *   htv4_record_t rec = reader.readrec();
 *
 *   // reading a file
 *   while (!reader.eof()) {
 *      string line = reader.readline();
 *   }
 */

#ifndef HTPv4Reader_H
#define HTPv4Reader_H

#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <string>
using namespace std;

#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

#include "bgz_reader.hpp"

/**
 * Data structure to store HTPv4 records.
 *
 * Numeric types with NA values are set to HTPv4_NA.
 *
 * Info fields are stored in a map<string,string>. Fields without
 * values (e.g. NO_BETA) map to the empty string "".
 *
 */
const int HTPv4_NA = -2147483648; // for numeric types
const string HTPv4_ESTIMATED_GENOTYPE_COUNTS_FLAG = "EGC_FLAG";
const string HTPv4_PVAL_STRING = "HTPv4_PVAL_STRING";
typedef struct {
  string name;
  string chr;
  int pos;
  string ref;
  string alt;
  string trait;
  string cohort;
  string model;
  double effect;
  double lci_effect;
  double uci_effect;
  double pval;
  double aaf;
  int num_cases;
  double cases_ref;
  double cases_het;
  double cases_alt;
  int num_controls;
  double controls_ref;
  double controls_het;
  double controls_alt;
  map<string, string> info;
} htpv4_record_t;

string htpv4_to_string(int i);
string htpv4_to_string(double d);
string htpv4_rounded_double_to_string(double d, bool keep_as_float);
string htpv4_make_info_string(const map<string, string>& info);
string htpv4_pval_to_string(double pval, const string& mlog10p, const string& pval_string);
double htpv4_pval_string_to_log10p(string pval);

class HTPv4Reader : public BgzReader {
 public:
  HTPv4Reader() : BgzReader(){};

  HTPv4Reader(string filepath);

  ~HTPv4Reader(){};

  htpv4_record_t readrec();

  void open(string filepath);

  void seek(string chrom, hts_pos_t position);

  bool eof();

  static double get_beta(const htpv4_record_t& rec);

  static double get_se(const htpv4_record_t& rec);

  static bool has_beta(const htpv4_record_t& rec);

  static bool has_se(const htpv4_record_t& rec);

 private:
  string first_line;
};



#endif
