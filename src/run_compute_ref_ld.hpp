#ifndef RUN_COMPUTE_REF_LD_H
#define RUN_COMPUTE_REF_LD_H

#include <string>
using namespace std;

void run_compute_ref_ld(const std::string& target_pfile,
                        const std::string& gene_file,
                        const std::string& chrom,
                        const std::string& out_prefix,
                        const std::string& buffer_pfile,
                        const double& buffer_mb,
                        const double& buffer_cm,
                        const std::string& genetic_map_file,
                        const double& target_r2,
                        const double& buffer_r2,
                        const int& float_size,
                        const int& block_size,
                        const int& nthreads,
                        const std::string& target_extract,
                        const std::string& target_exclude,
                        const std::string& target_keep,
                        const std::string& target_remove,
                        const std::string& buffer_extract,
                        const std::string& buffer_exclude,
                        const bool& skip_buffer,
                        const bool& use_dosages);

#endif