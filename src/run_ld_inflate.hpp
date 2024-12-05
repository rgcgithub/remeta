#ifndef RUN_LD_INFLATE_H
#define RUN_LD_INFLATE_H

#include "lapack_complex.hpp"

#include <string>
using namespace std;

#include "logging.hpp"

void run_ld_inflate(const string& mat_prefix, 
                    const string& out_prefix,
                    const string& gene_name,
                    const string& extract);

#endif