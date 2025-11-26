#include "run_genep.hpp"

#include "meta/genep_meta_analyzer.hpp"
#include "parameter_checks.hpp"
namespace pc = parameter_checks;

using namespace std;

void run_genep(const vector<string>& htp_files,
               const string& genep_def_file,
               const string& chr,
               const std::string& burden_model,
               const std::string& acatv_model,
               const std::string& skato_model,
               const std::string& sbat_model,
               const std::string& extract,
               const std::string& exclude,
               const bool& drop_regenie_gene,
               const string& out_prefix) {
  pc::check_htp_files(htp_files);
  if (genep_def_file != "") {
    pc::check_file_exists(genep_def_file);
  }

  GenePMetaAnalyzer meta(genep_def_file,
                         burden_model,
                         acatv_model,
                         skato_model,
                         sbat_model);

  VariantFilter vf;
  vf.set_skip_genep();

  if (extract != "") {
    pc::check_file_exists(extract);
    vf.set_extract_file(extract);
  }
  if (exclude != "") {
    pc::check_file_exists(exclude);
    vf.set_exclude_file(exclude);
  }

  if (drop_regenie_gene && (sbat_model == "ADD-WGR-BURDEN-SBAT" || sbat_model == "ADD-WGR-BURDEN-SBAT-META")) {
    vf.set_keep_remeta_gene_and_sbat();
  } else if (drop_regenie_gene) {
    vf.set_keep_remeta_gene_only();
  }

  HTPv4Writer out(out_prefix + ".merge.remeta.gz", "w");
  out.writeheader();
  vector<HTPv4Reader> readers;
  for (const string& htp_file : htp_files) {
    readers.push_back(HTPv4Reader(htp_file));
  }

  run_variant_meta_analysis(meta, vf, readers, out, chr);
}
