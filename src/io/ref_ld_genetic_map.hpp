#ifndef REF_LD_GENETIC_MAP_H
#define REF_LD_GENETIC_MAP_H

#include <string>
#include <map>

#include "../lapack_complex.hpp"
#include <Eigen/Dense>

struct genetic_map_entry_t {
	double bp;
	std::string chr;
	double cm;

	bool operator==(const genetic_map_entry_t& other) {
		return this->bp == other.bp
			&& this->chr == other.chr
			&& this->cm == other.cm;
	}
};
const genetic_map_entry_t NULL_GENETIC_MAP_ENTRY = {-1, "", -1};

class RefLDGeneticMap {
 public:
  RefLDGeneticMap(const std::string& map_file);

  // position in bases to position in cm
  double bp_to_cm(double bp);

  // position in cm to position in bases
  double cm_to_bp(double cm);

 private:
	const int lookback_size = 100;

	std::string map_file;
	std::map<double, double> cm_bp;
  std::map<double, double> bp_cm;
  double start_bp;
  double end_bp;
  double start_cm;
  double end_cm;
  double slope_cm_bp_beg;
  double slope_bp_cm_beg;
  double slope_cm_bp_end;
  double slope_bp_cm_end;

  genetic_map_entry_t parse_map_entry(const std::string& line);

  double compute_slope_leastsq(Eigen::VectorXd y, Eigen::VectorXd x);
};

#endif