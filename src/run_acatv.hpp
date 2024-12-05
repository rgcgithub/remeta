#ifndef RUN_ACATV_H
#define RUN_ACATV_H

#include "lapack_complex.hpp"

#include <string>
#include <vector>

void run_acatv(const std::vector<std::string>& htp_files,
               const std::vector<std::string>& cohort_meta,
               const std::string& anno_file,
               const std::string& set_list_file,
               const std::string& mask_def_file,
               const std::string& trait_name,
               const std::string& trait_type,
               const std::string& out_prefix,
               const std::string& af_strategy,
               const std::string& af_file,
               const std::string& weight_strategy,
               const int& min_aac,
               const double& max_aaf,
               const std::string& chr,
               const std::string& gene,
               const std::string& extract_file,
               const std::string& exclude_file,
               const std::vector<std::string>& sources,
               const std::vector<std::string>& ld_prefixes,
               const std::string& condition_list_file,
               const std::vector<std::string>& condition_htp_files,
               const int& threads,
	       const bool& write_mask_snplist);

#endif
