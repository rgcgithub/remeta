#include "ref_ld_genetic_map.hpp"

#include <cmath>
#include <string>
#include <map>
using namespace std;

#include <boost/algorithm/string.hpp>
#include <Eigen/Dense>
using Eigen::VectorXd;
using Eigen::MatrixXd;

#include "bgz_reader.hpp"
#include "../util.hpp"

RefLDGeneticMap::RefLDGeneticMap(const string& map_file) {
  this->map_file = map_file;
  VectorXd cm_beg(this->lookback_size);
  VectorXd cm_end(this->lookback_size);
  VectorXd bp_beg(this->lookback_size);
  VectorXd bp_end(this->lookback_size);
  int entries_read = 0;

  BgzReader map(map_file);
  while (!map.eof()) {
    genetic_map_entry_t map_entry = this->parse_map_entry(map.readline());
    if (map_entry == NULL_GENETIC_MAP_ENTRY) {
      continue;
    }

    if (entries_read == 0) {
      this->start_bp = map_entry.bp;
      this->start_cm = map_entry.cm;
    }
    this->end_bp = map_entry.bp;
    this->end_cm = map_entry.cm;

    if (entries_read < this->lookback_size) {
      cm_beg[entries_read] = map_entry.cm;
      bp_beg[entries_read] = map_entry.bp;
    }
    cm_end[entries_read % this->lookback_size] = map_entry.cm;
    bp_end[entries_read % this->lookback_size] = map_entry.bp;
    ++entries_read;

    this->bp_cm[map_entry.bp] = map_entry.cm;
    this->cm_bp[map_entry.cm] = map_entry.bp;
  }

  if (entries_read < this->lookback_size) {
    cm_beg = VectorXd(cm_beg(Eigen::seqN(0, entries_read)));
    cm_end = VectorXd(cm_end(Eigen::seqN(0, entries_read)));
    bp_beg = VectorXd(bp_beg(Eigen::seqN(0, entries_read)));
    bp_end = VectorXd(bp_end(Eigen::seqN(0, entries_read)));
  }

  this->slope_cm_bp_beg = this->compute_slope_leastsq(cm_beg, bp_beg);
  this->slope_bp_cm_beg = this->compute_slope_leastsq(bp_beg, cm_beg);
  this->slope_cm_bp_end = this->compute_slope_leastsq(cm_end, bp_end);
  this->slope_bp_cm_end = this->compute_slope_leastsq(bp_end, cm_end);
}

double RefLDGeneticMap::bp_to_cm(double bp) {
  map<double, double>::iterator it_lb = this->bp_cm.lower_bound(bp);
  if (it_lb->first == bp) {
    return it_lb->second;
  } else if (it_lb == this->bp_cm.begin()) {
    return max(
      0.0, 
      this->start_cm + this->slope_cm_bp_beg * (bp - this->start_bp)
    );
  } else if (it_lb == this->bp_cm.end()) {
    return this->end_cm + this->slope_cm_bp_end * (bp - this->end_bp);
  } else {
    double bp2 = it_lb->first;
    double cm2 = it_lb->second;
    --it_lb;
    double bp1 = it_lb->first;
    double cm1 = it_lb->second;
    double slope = (cm2 - cm1) / (bp2 - bp1);
    return cm1 + slope * (bp - bp1);
  }
}

double RefLDGeneticMap::cm_to_bp(double cm) {
  map<double, double>::iterator it_lb = this->cm_bp.lower_bound(cm);
  if (it_lb->first == cm) {
    return it_lb->second;
  } else if (it_lb == this->cm_bp.begin()) {
    return max(
      0.0,
      this->start_bp + this->slope_bp_cm_beg * (cm - this->start_cm)
    );
  } else if (it_lb == this->cm_bp.end()) {
    return this->end_bp + this->slope_bp_cm_end * (cm - this->end_cm);
  } else {
    double cm2 = it_lb->first;
    double bp2 = it_lb->second;
    --it_lb;
    double cm1 = it_lb->first;
    double bp1 = it_lb->second;

    double slope = (bp2 - bp1) / (cm2 - cm1);
    return bp1 + slope * (cm - cm1);
  }
}

genetic_map_entry_t RefLDGeneticMap::parse_map_entry(const string& line) {
  vector<string> cols = util::str_split(line, "\t ");
  if (cols.size() != 3) {
    throw runtime_error(
      "unable to parse map file " + this->map_file 
      + ": incorrect number of columns for line "
      + line + " (expected 3, found "
      + to_string(cols.size()) + ")"
    );
  }

  genetic_map_entry_t map_entry;
  if (boost::algorithm::to_lower_copy(cols[0]) == "pos") {
    return NULL_GENETIC_MAP_ENTRY;
  } else {
    return genetic_map_entry_t {
      stod(cols[0]),
      cols[1],
      stod(cols[2])
    };
  }
}

double RefLDGeneticMap::compute_slope_leastsq(VectorXd y, VectorXd x) {
  // add an intercept
  MatrixXd X(x.rows(), 2);
  X.col(0) = x;
  X.col(1) = VectorXd::Zero(x.rows()).array() + 1;

  // solve and return the slope
  MatrixXd solution = X.colPivHouseholderQr().solve(y);
  return solution(0);
}