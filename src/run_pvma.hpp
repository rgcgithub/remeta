#ifndef REMETA_RUN_PVMA_H
#define REMETA_RUN_PVMA_H

#include<string>
#include <vector>
using namespace std;

void run_pvma(const vector<string>& htp_files,
              const vector<string>& cohorts,
              const string& trait_name,
              const string& trait_type,
              const string& out_prefix,
              const string& method,
              const bool& unweighted,
              const bool& two_sided,
              const string& chr,
              const string& extract_file,
              const string& exclude_file,
              const bool& skip_beta,
              const bool& skip_genep,
              const vector<string>& sources,
              const string& source_def);

#endif