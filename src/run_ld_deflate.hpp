#ifndef RUN_LD_DEFLATE_H
#define RUN_LD_DEFLATE_H

#include <string>

void run_ld_deflate(const std::string& ld_file,
                    const int& sample_size,
                    const double& sparsity_threshold,
                    const std::string& out);

#endif