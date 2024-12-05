#include "run_index_anno.hpp"

#include "io/regenie_anno_reader.hpp"
#include "logging.hpp"
#include "parameter_checks.hpp"
namespace pc = parameter_checks;

void run_index_anno(const string& anno_file, const int& stride) {
  pc::check_file_exists(anno_file);

  RegenieAnnoReader reader(anno_file);
  try {
    log_info("indexing " + anno_file + " (stride : " + to_string(stride) + ")");
    reader.build_index(stride);
    log_info("writing index...");
    reader.write_index();
  } catch (const std::exception& ex) {
    log_error(ex.what(), 1);
  }
}