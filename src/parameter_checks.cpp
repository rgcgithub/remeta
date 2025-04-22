#include "parameter_checks.hpp"

#include <boost/algorithm/string.hpp>
#include <iostream>

#include "util.hpp"

namespace parameter_checks {
  trait_type_e check_trait_type(string trait_type) {
    boost::algorithm::to_lower(trait_type);
    if (trait_type == "qt") {
      return QT;
    } else if (trait_type == "bt") {
      return BT;
    } else {
      log_error("trait must be one of qt or bt", 1);
    }
    return QT;
  }

  void check_file_exists(string filename) {
    if (!util::file_exists(filename)) {
      log_error(filename + " does not exist", 1);
    }
  }

  void check_htp_files(vector<string> htp_files) {
    for (string file : htp_files) {
      check_file_exists(file);
    }
  }

  pvma_method_e check_pvma_method(string method, bool unweighted, bool two_sided) {
    pvma_method_e pvma_method;
    if (method == "stouffers") {
      pvma_method = STOUFFERS;
    } else if (method == "fishers") {
      pvma_method = FISHERS;
    } else {
      log_error("method must be one of: stouffers, fishers", 1);
    }
    
    if (pvma_method != STOUFFERS && unweighted) {
      log_error("--unweighted is only valid with --method stouffers", 1);
    }
    if (pvma_method != STOUFFERS && two_sided) {
      log_error("--two-sided is only valid with --method stouffers", 1);
    }

    return pvma_method;
  }

  void check_cohorts(vector<string> cohorts, vector<string> htp_files) {
    if (cohorts.size() != htp_files.size()) {
      log_error("--cohorts and --htp have different number of parameters", 1);
    }
  }

  void check_ld_files(vector<string> ld_prefixes, const size_t& ncohorts) {
    if (ncohorts != ld_prefixes.size()) {
      log_error("number of LD files does not match number of cohorts", 1);
    }

    for (const string& ld_prefix : ld_prefixes) {
      bool is_rg_format = util::file_exists(ld_prefix + ".rg.ld")
                          && util::file_exists(ld_prefix + ".rg.ld.idx.gz");
      bool is_rm_format = util::file_exists(ld_prefix + ".remeta.buffer.ld")
                          && util::file_exists(ld_prefix + ".remeta.gene.ld")
                          && util::file_exists(ld_prefix + ".remeta.ld.idx.gz");

      if (!is_rg_format && !is_rm_format) {
        log_error("could not find LD matrices with prefix " + ld_prefix, 1);
      }
    }
  }

  void check_rho_values(vector<double> rho_values) {
    for (const double& rho : rho_values) {
      if (rho < 0 || rho > 1) {
        log_error("rho-values must be between 0 and 1.", 1);
      }
    }
  }

  void check_aac(double aac) {
    if (aac < 0) {
      log_error("aac must be >= 0", 1);
    }
  }

  void check_af(double af) {
    if (af < 0 || af > 1) {
      log_error("af must be between 0 and 1", 1);
    }
  }

  weight_strategy_e check_weight_strategy(string weight_strategy) {
    boost::algorithm::to_lower(weight_strategy);
    if (weight_strategy == "beta") {
      return USE_BETA_WEIGHTS;
    } else if (weight_strategy == "uniform") {
      return USE_UNIFORM_WEIGHTS;
    } else {
      log_error("weight-strategy must be one of 'beta' or 'uniform', found: " + weight_strategy, 1);
    }
    return USE_BETA_WEIGHTS;
  }

  af_strategy_e check_af_strategy(string af_strategy) {
    boost::algorithm::to_lower(af_strategy);
    if (af_strategy == "overall") {
      return USE_OVERALL_AF;
    } else if (af_strategy == "max") {
      return USE_MAX_AF;
    } else {
      log_error("af-strategy must be one of 'overall' or 'max', found: " + af_strategy, 1);
    }
    return USE_OVERALL_AF;
  }

  void check_cf_files(vector<string> cf_files, vector<string> htp_files) {
    if (cf_files.size() != htp_files.size()) {
      log_error("number of corr-factors files does not match number of htp files");
    }
    for (const string& cf_file : cf_files) {
      check_file_exists(cf_file);
    }
  }

  void check_anno_file_indexed(string file) {
    if (!util::file_exists(file + ".rgi") && !util::file_exists(file + ".tbi")) {
      log_warning(file + " is not indexed: this may impact speed and memory usage");
    }
  }

  void check_block_mapper_params(const double& buffer_mb, 
                                 const double& buffer_cm, 
                                 const string& genetic_map_file) {
    if (buffer_mb != 0 && buffer_cm != 0) {
      log_error("only one of --buffer-mb and --buffer-cm can be set", 1);
    } else if (buffer_cm != 0 && genetic_map_file == "") {
      log_error("--buffer-cm requires --genetic-map to be set", 1);
    } else if (genetic_map_file != "" && buffer_cm == 0) {
      log_error("--genetic-map requires --buffer-cm to be set", 1);
    } else if (genetic_map_file != "") {
      check_file_exists(genetic_map_file);
    }
  }

  void check_conditional_analysis_files(const string& condition_file,
                                        const vector<string>& conditional_htp_files,
                                        const int& n_cohorts) {
    if(!util::file_exists(condition_file)) {
      log_error(condition_file + " does not exist", 1);
    }
    for ( const string& htp_file : conditional_htp_files ) {
      if (!util::file_exists(htp_file)) {
        log_error(htp_file + " does not exist", 1);
      }
    }

    if ((int)conditional_htp_files.size() != n_cohorts) {
      log_error(
        "number of --conditional-files does not match number of --htp files ("
        + to_string(conditional_htp_files.size()) + " != " + to_string(n_cohorts) + ")"
      , 1);
    }
  }

  void check_aaf_bins(const vector<double>& aaf_bins) {
    for (const double& aaf : aaf_bins) {
      if (aaf <= 0 || aaf > 1) {
        log_error("allele frequencies passed to --aaf-bins must be >= 0 and < 1", 1);
      }
    }
  }

  singleton_def_e check_singleton_def(const string& singleton_def) {
    if (singleton_def == "within") {
      return WITHIN;
    } else if (singleton_def == "across") {
      return ACROSS;
    } else if (singleton_def == "omit") {
      return OMIT;
    } else {
      log_error("unrecognized --singleton-def " + singleton_def, 1);
      return WITHIN;
    }
  }

  void check_pval(double pval) {
    if (pval < 0 || pval > 1) {
      log_error("pvalues must be in the interval [0, 1]", 1);
    }
  }

  void check_case_control_ratio(double case_control_ratio) {
    if (case_control_ratio < 0) {
      log_error("case control ratio must be non-negative", 1);
    }
  }

  void check_float_size(int float_size, double buffer_r2) {
    if (float_size != 1 && float_size != 2 && float_size != 4) {
      log_error("float size must be either 1, 2 or 4", 1);
    }
    if (buffer_r2 < 0.0078125*0.0078125 && float_size == 1) {
      log_error("1 byte floats do not have enough precision for buffer_r2 " + to_string(buffer_r2), 1);
    }
  }
}