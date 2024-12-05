#ifndef RUN_ESMA_H
#define RUN_ESMA_H

#include <string>
#include <vector>
using namespace std;

void run_esma(const vector<string>& htp_files,
              const vector<string>& cohorts,
              const string& trait_name,
              const string& trait_type,
              const string& out_prefix,
              const string& chr,
              const string& extract_file,
              const string& exclude_file,
              const vector<string>& sources,
              const string& source_def);

#endif