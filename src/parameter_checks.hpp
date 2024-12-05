#ifndef REMETA_PARAMETER_CHECKS_H
#define REMETA_PARAMETER_CHECKS_H

#include <string>
#include <vector>
using namespace std;

#include <stdio.h>

#include "enums.hpp"
#include "logging.hpp"

namespace parameter_checks {

  trait_type_e check_trait_type(string trait_type);

  void check_file_exists(string filename);

  void check_htp_files(vector<string> htp_files);

  pvma_method_e check_pvma_method(string method);

  void check_cohorts(vector<string> cohorts, vector<string> htp_files);

  void check_ld_files(vector<string> ld_prefixes, const size_t& ncohorts);

  void check_rho_values(vector<double> rho_values);

  void check_aac(double aac);

  void check_af(double af);

  weight_strategy_e check_weight_strategy(string weight_strategy);

  af_strategy_e check_af_strategy(string af_strategy);

  void check_cf_files(vector<string> cf_files, vector<string> htp_files);

  void check_anno_file_indexed(string file);

  void check_block_mapper_params(const double& buffer_mb, 
                                 const double& buffer_cm, 
                                 const string& genetic_map_file);

  void check_conditional_analysis_files(const string& condition_file,
                                        const vector<string>& conditional_htp_files,
                                        const int& n_cohorts);

  void check_aaf_bins(const vector<double>& aaf_bins);

  singleton_def_e check_singleton_def(const string& singleton_def);

  void check_pval(double pval);

  void check_case_control_ratio(double case_control_ratio);

  void check_float_size(int float_size, double buffer_r2);
}


#endif