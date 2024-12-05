#include "run_esma.hpp"

#include <chrono>

#include "io/htpv4_reader.hpp"
#include "io/htpv4_writer.hpp"
#include "meta/es_meta_analyzer.hpp"
#include "meta/variant_meta_analyzer.hpp"
#include "enums.hpp"
#include "htpv4_pos.hpp"
#include "logging.hpp"
#include "parameter_checks.hpp"
namespace pc = parameter_checks;
#include "variant_filter.hpp"

void run_esma(const vector<string>& htp_files,
              const vector<string>& cohorts,
              const string& trait_name,
              const string& trait_type,
              const string& out_prefix,
              const string& chr,
              const string& extract_file,
              const string& exclude_file,
              const vector <string>& sources,
              const string& source_def) {
  pc::check_htp_files(htp_files);
  pc::check_cohorts(cohorts, htp_files);
  trait_type_e trait = pc::check_trait_type(trait_type);

  VariantFilter vf;
  if (extract_file != "") {
    pc::check_file_exists(extract_file);
    vf.set_extract_file(extract_file);
  }
  if (exclude_file != "") {
    pc::check_file_exists(exclude_file);
    vf.set_exclude_file(exclude_file);
  }
  vf.set_effect_not_na();
  vf.set_info_has_se();
  if (sources.size() > 0) {
    vf.set_info_source_is_one_of(sources);
  }

  ESMetaAnalyzer meta(
    trait,
    trait_name,
    cohorts,
    (int)htp_files.size()
  );

  if (source_def != "") {
    pc::check_file_exists(source_def);
    meta.set_source_definitions(source_def);
  }

  HTPv4Writer out(out_prefix + ".esma.remeta.gz", "w");
  out.writeheader();
  vector<HTPv4Reader> readers;
  for (const string& htp_file : htp_files) {
    readers.push_back(HTPv4Reader(htp_file));
  }

  run_variant_meta_analysis(meta, vf, readers, out, chr);
}