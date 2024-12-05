#ifndef RUN_GENEP_H
#define RUN_GENEP_H

#include <string>
#include <vector>

void run_genep(const std::vector<std::string>& htp_files,
               const std::string& genep_def_file,
               const std::string& chr,
               const std::string& burden_model,
               const std::string& acatv_model,
               const std::string& skato_model,
               const bool& include_sbat,
               const std::string& extract,
               const std::string& exclude,
               const bool& drop_regenie_gene,
               const std::string& out_prefix);

#endif