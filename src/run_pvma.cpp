#include "run_pvma.hpp"

#include <chrono>

#include "io/htpv4_reader.hpp"
#include "io/htpv4_writer.hpp"
#include "meta/pv_meta_analyzer.hpp"
#include "meta/variant_meta_analyzer.hpp"
#include "enums.hpp"
#include "htpv4_pos.hpp"
#include "logging.hpp"
#include "parameter_checks.hpp"
namespace pc = parameter_checks;
#include "variant_filter.hpp"

bool htp_record_has_beta(const htpv4_record_t& rec) {
  return rec.effect != HTPv4_NA;
}

bool htp_record_is_genep(const htpv4_record_t& rec) {
  return rec.model.find("GENEP") != string::npos;
}

bool skip_htp_record(const htpv4_record_t& rec,
                     const bool& skip_beta,
                     const bool& skip_genep) {
  return (skip_genep && htp_record_is_genep(rec))
         || (skip_beta && htp_record_has_beta(rec));
}

void run_pvma(const vector<string>& htp_files,
              const vector<string>& cohorts,
              const string& trait_name,
              const string& trait_type,
              const string& out_prefix,
              const string& method,
              const bool& unweighted,
              const string& chr,
              const string& extract_file,
              const string& exclude_file,
              const bool& skip_beta,
              const bool& skip_genep,
              const vector<string>& sources,
              const string& source_def) {
  pc::check_htp_files(htp_files);
  pc::check_cohorts(cohorts, htp_files);
  trait_type_e trait = pc::check_trait_type(trait_type);
  pvma_method_e pvma_method = pc::check_pvma_method(method);
  
  VariantFilter vf;
  if (extract_file != "") {
    pc::check_file_exists(extract_file);
    vf.set_extract_file(extract_file);
  }
  if (exclude_file != "") {
    pc::check_file_exists(exclude_file);
    vf.set_exclude_file(exclude_file);
  }
  if (skip_beta) {
    vf.set_effect_is_na();
  }
  if (skip_genep) {
    vf.set_skip_genep();
  }
  if (sources.size() > 0) {
    vf.set_info_source_is_one_of(sources);
  }
  
  PVMetaAnalyzer meta(
    pvma_method,
    trait,
    trait_name,
    cohorts,
    unweighted
  );

  if (source_def != "") {
    pc::check_file_exists(source_def);
    meta.set_source_definitions(source_def);
  }

  HTPv4Writer out(out_prefix + ".pvma.remeta.gz", "w");
  out.writeheader();
  vector<HTPv4Reader> readers;
  for (const string& htp_file : htp_files) {
    readers.push_back(HTPv4Reader(htp_file));
  }

  run_variant_meta_analysis(meta, vf, readers, out, chr);
}