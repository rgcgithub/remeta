#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <htslib/hts_log.h>

#include "logging.hpp"
#include "run_acatv.hpp"
#include "run_compute_ref_ld.hpp"
#include "run_esma.hpp"
#include "run_genep.hpp"
#include "run_htp.hpp"
#include "run_index_anno.hpp"
#include "run_ld_deflate.hpp"
#include "run_ld_inflate.hpp"
#include "run_pvma.hpp"

void log_options_in_effect(const vector<string>& options) {
  string opt_fmt = " ";
  for (const string& option : options) {
    if (opt_fmt != " " && option[0] == '-') {
      opt_fmt += "\n  " + option;
    } else {
      opt_fmt += " " + option;
    }
  }
  log_info("Options in effect:\n" + opt_fmt + "\n\n");
}

void log_program_info(const string& command, const vector<string>& options) {
  string version = "remeta " VERSION_NUMBER;
  #if defined(WITH_MKL)
    version += " + Intel MKL\n";
  #elif defined(WITH_OPENBLAS)
    version += " + OpenBLAS\n";
  #endif
  log_info(version);
  log_info("Command: remeta " + command + "\n");
  log_options_in_effect(options);
}

void pvma(const vector<string>& options) {
  po::options_description opts(R""""(remeta: p-value meta-analysis

Usage: remeta pvma [OPTIONS]

Options (* Mandatory))"""");

  opts.add_options()
    ("htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->required()->multitoken(),
      "List of HTP input files. (*)")
    ("cohorts", po::value<vector<string> >()->value_name("NAME1 NAME2 ...")->required()->multitoken(),
      "List of cohort names per file. (*)")
    ("trait-name", po::value<string>()->value_name("NAME")->required(),
      "Name of trait. (*)")
    ("trait-type", po::value<string>()->value_name("TYPE")->required(),
      "One of BT or QT. (*)")
    ("out", po::value<string>()->value_name("PREFIX")->required(),
     "Prefix for output file. (*)")
    ("method", po::value<string>()->value_name("METHOD")->default_value("stouffers"),
      "One of stouffers or fishers")
    ("unweighted", "Omit sample size weighting (affects stouffers only).")
    ("chr", po::value<string>()->value_name("CHR"),
     "Run only on chromosome CHR.")
    ("extract", po::value<string>()->value_name("FILE"),
      "Include only the variants with IDs listed in this file (one per line).")
    ("exclude", po::value<string>()->value_name("FILE"),
      "Exclude variants with IDs listed in this file (one per line).")
    ("skip-beta", "Skip entries with an effect size estimate.")
    ("sources", po::value<vector<string> >()->value_name("SOURCE1 SOURCE2 ...")->multitoken(),
      "Only meta-analyze variants where the info field SOURCE is one of SOURCE1 SOURCE2 ...")
    ("source-def", po::value<string>()->value_name("FILE"),
      "Two column file mapping long SOURCE info fields to shorter SOURCE info fields")
    ("log-level", po::value<string>()->value_name("LEVEL"), "Set logging level.")
    ("help", "Print this message and exit.");

  po::options_description rgc_opts("RGC options");
  rgc_opts.add_options()
    ("skip-genep", "Skip GENE_P entries (i.e. where Name ends in .GENE and the test is not SBAT)");

  po::options_description all_opts("All options");
  all_opts.add(opts).add(rgc_opts);

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try {
    po::parsed_options parsed = po::command_line_parser(options)
      .options(all_opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  if (vm.count("log-level") && vm["log-level"].as<string>() == "debug") {
    init_logging("debug", vm["out"].as<string>() + ".pvma.log");
  } else if (vm.count("log-level")) {
    cerr << "error: unrecognized --log-level: " + vm["log-level"].as<string>()
         << endl;
    exit(1);
  } else {
    init_logging("info", vm["out"].as<string>() + ".pvma.log");
  }

  bool unweighted = false;
  if (vm.count("unweighted") > 0) {
    unweighted =true;
  }
  string chr = "";
  if (vm.count("chr") > 0) {
    chr = vm["chr"].as<string>();
  }
  string extract = "";
  if (vm.count("extract") > 0) {
    extract = vm["extract"].as<string>();
  }
  string exclude = "";
  if (vm.count("exclude") > 0) {
    exclude = vm["exclude"].as<string>();
  }
  vector<string> sources;
  if (vm.count("sources") > 0) {
    sources = vm["sources"].as<vector<string> >();
  }
  bool skip_beta = false;
  if (vm.count("skip-beta") > 0) {
    skip_beta = true;
  }
  bool skip_genep = false;
  if (vm.count("skip-genep") > 0) {
    skip_genep = true;
  }
  string source_def = "";
  if (vm.count("source-def") > 0) {
    source_def = vm["source-def"].as<string>(); 
  }

  log_program_info("pvma", options);

  run_pvma(
    vm["htp"].as<vector<string> >(),
    vm["cohorts"].as<vector<string> >(),
    vm["trait-name"].as<string>(),
    vm["trait-type"].as<string>(),
    vm["out"].as<string>(),
    vm["method"].as<string>(),
    unweighted,
    chr,
    extract,
    exclude,
    skip_beta,
    skip_genep,
    sources,
    source_def
  );
}

void esma(const vector<string>& options) {
  po::options_description opts(R""""(remeta: effect size meta-analysis

Usage: remeta esma [OPTIONS]

Options (* Mandatory))"""");

  opts.add_options()
    ("htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->required()->multitoken(),
      "List of HTP input files. (*)")
    ("cohorts", po::value<vector<string> >()->value_name("NAME1 NAME2 ...")->required()->multitoken(),
      "List of cohort names per file. (*)")
    ("trait-name", po::value<string>()->value_name("TYPE")->required(),
      "Name of trait. (*)")
    ("trait-type", po::value<string>()->value_name("TYPE")->required(),
      "One of BT or QT. (*)")
    ("out", po::value<string>()->value_name("PREFIX")->required(),
     "Prefix for output file. (*)")
    ("chr", po::value<string>()->value_name("CHR"),
     "Run only on chromosome CHR.")
    ("extract", po::value<string>()->value_name("FILE"),
      "Include only the variants with IDs listed in this file (one per line).")
    ("exclude", po::value<string>()->value_name("FILE"),
      "Exclude variants with IDs listed in this file (one per line).")
    ("sources", po::value<vector<string> >()->value_name("SOURCE1 SOURCE2 ...")->multitoken(),
      "Only meta-analyze variants where the info field SOURCE is one of SOURCE1 SOURCE2 ...")
    ("source-def", po::value<string>()->value_name("FILE"),
      "Two column file mapping long SOURCE info fields to shorter SOURCE info fields")
    ("log-level", po::value<string>()->value_name("LEVEL"), "Set logging level.")
    ("help", "Print this message and exit.");

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try {
    po::parsed_options parsed = po::command_line_parser(options)
      .options(opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  if (vm.count("log-level") && vm["log-level"].as<string>() == "debug") {
    init_logging("debug", vm["out"].as<string>() + ".esma.log");
  } else if (vm.count("log-level")) {
    cerr << "error: unrecognized --log-level: " + vm["log-level"].as<string>()
         << endl;
    exit(1);
  } else {
    init_logging("info", vm["out"].as<string>() + ".esma.log");
  }

  string chr = "";
  if (vm.count("chr") > 0) {
    chr = vm["chr"].as<string>();
  }
  string extract = "";
  if (vm.count("extract") > 0) {
    extract = vm["extract"].as<string>();
  }
  string exclude = "";
  if (vm.count("exclude") > 0) {
    exclude = vm["exclude"].as<string>();
  }
  vector<string> sources;
  if (vm.count("sources") > 0) {
    sources = vm["sources"].as<vector<string> >();
  }
  string source_def = "";
  if (vm.count("source-def") > 0) {
    source_def = vm["source-def"].as<string>(); 
  }

  log_program_info("esma", options);

  run_esma(vm["htp"].as<vector<string> >(),
           vm["cohorts"].as<vector<string> >(),
           vm["trait-name"].as<string>(),
           vm["trait-type"].as<string>(),
           vm["out"].as<string>(),
           chr,
           extract,
           exclude,
           sources,
           source_def         
  );
}

void acatv(const vector<string>& options) {
  po::options_description opts(R""""(remeta: aggregated cauchy association test meta-analysis

Usage: remeta acatv [OPTIONS]

Options (* Mondatory))"""");

  opts.add_options()
    ("htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->required()->multitoken(),
      "List of HTP input files. (*)")
    ("cohorts", po::value<vector<string> >()->value_name("NAME1 NAME2 ...")->required()->multitoken(),
      "List of cohort names per file. (*)")
    ("anno-file", po::value<string>()->value_name("FILE")->required(),
      "File with variant annotations. Bgzipped and indexed with 'index-anno'. (*)")
    ("set-list", po::value<string>()->value_name("FILE")->required(),
      "Regenie set-list file. (*)")
    ("mask-def", po::value<string>()->value_name("FILE")->required(),
      "Regenie mask-def file. (*)")
    ("trait-name", po::value<string>()->value_name("NAME")->required(),
      "Name of trait. (*)")
    ("trait-type", po::value<string>()->value_name("TYPE")->required(),
      "One of bt or qt. (*)")
    ("out", po::value<string>()->value_name("PREFIX")->required(),
      "Prefix for output file. (*)")
    ("weight-strategy", po::value<string>()->value_name("STR")->default_value("beta"),
      "Strategy to compute variant weights. One of 'beta' or 'uniform'.")
    ("af-strategy", po::value<string>()->value_name("STR")->default_value("overall"),
      "Strategy to compute variant allele frequences. One of 'overall' or 'max'.")
    ("aaf-file", po::value<string>()->value_name("FILE"),
      "Use precomputed alternate allele frequencies from an external file.")
    ("max-aaf", po::value<double>()->value_name("FLOAT")->default_value(0.01),
      "Maximum allele frequency for a variant to be included in mask.")
    ("min-aac", po::value<int>()->value_name("INT")->default_value(5),
      "Minimum AAC across cohorts for a variant to be included in a mask.")
    ("chr", po::value<string>()->value_name("CHR"),
      "Run only on chromosome CHR.")
    ("gene", po::value<string>()->value_name("GENE"),
      "If set, run only on this gene.")
    ("extract", po::value<string>()->value_name("FILE"),
      "Include only the variants with IDs listed in this file (one per line).")
    ("exclude", po::value<string>()->value_name("FILE"),
      "Exclude variants with IDs listed in this file (one per line).")
    ("sources", po::value<vector<string> >()->value_name("SOURCE1 SOURCE2 ...")->multitoken(),
      "Only include variants where the info field SOURCE is one of SOURCE1 SOURCE2 ...")
    ("ld-prefixes", po::value<vector<string> >()->value_name("FILE1 FILE2 ...")->required()->multitoken(),
      "Prefixes to LD files per cohort.")
    ("condition-list", po::value<string>()->value_name("FILE"),
      "File with variants to condition on (one per line). Requires --ld-prefixes option.")
    ("condition-htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->multitoken(),
      "List of HTP files with summary statistics of conditional variants per study. Requires --ld-prefixes option.") 
    ("write-mask-snplist", "Write file with list of variants included in each mask.")
    ("threads", po::value<int>()->value_name("INT")->default_value(1),
      "Number of threads to use.")
    ("log-level", po::value<string>()->value_name("LEVEL"), "Set logging level.")
    ("help", "Print this message and exit.");

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try { 
    po::parsed_options parsed = po::command_line_parser(options)
      .options(opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  if (vm.count("log-level") && vm["log-level"].as<string>() == "debug") {
    init_logging("debug", vm["out"].as<string>() + ".acatv.log");
  } else if (vm.count("log-level")) {
    cerr << "error: unrecognized --log-level: " + vm["log-level"].as<string>()
         << endl;
    exit(1);
  } else {
    init_logging("info", vm["out"].as<string>() + ".acatv.log");
  }

  string af_file;
  if (vm.count("aaf-file") > 0) {
    af_file = vm["aaf-file"].as<string>();
  }
  string chr;
  if (vm.count("chr") > 0) {
    chr = vm["chr"].as<string>();
  }
  string gene;
  if (vm.count("gene") > 0) {
    gene = vm["gene"].as<string>();
  }
  string extract;
  if (vm.count("extract") > 0) {
    extract = vm["extract"].as<string>();
  }
  string exclude;
  if (vm.count("exclude") > 0) {
    exclude = vm["exclude"].as<string>();
  }
  vector<string> sources;
  if (vm.count("sources") > 0) {
    sources = vm["sources"].as<vector<string> >();
  }
  vector<string> ld_prefixes;
  if (vm.count("ld-prefixes") > 0) {
    ld_prefixes = vm["ld-prefixes"].as<vector<string> >();
  }
  string condition_file;
  if (vm.count("condition-list") > 0) {
    condition_file = vm["condition-list"].as<string>();
  }
  vector<string> condition_htp;
  if (vm.count("condition-htp") > 0) {
    condition_htp = vm["condition-htp"].as<vector<string> >();
  }

  log_program_info("acatv", options);

  run_acatv(
    vm["htp"].as<vector<string> >(),
    vm["cohorts"].as<vector<string> >(),
    vm["anno-file"].as<string>(),
    vm["set-list"].as<string>(),
    vm["mask-def"].as<string>(),
    vm["trait-name"].as<string>(),
    vm["trait-type"].as<string>(),
    vm["out"].as<string>(),
    vm["af-strategy"].as<string>(),
    af_file,
    vm["weight-strategy"].as<string>(),
    vm["min-aac"].as<int>(),
    vm["max-aaf"].as<double>(),
    chr,
    gene,
    extract,
    exclude,
    sources,
    ld_prefixes,
    condition_file,
    condition_htp,
    vm["threads"].as<int>(),
    vm.count("write-mask-snplist") > 0
  );
}

void skato(const vector<string>& options) {
  po::options_description opts(R""""(remeta: omnibus sequence kernal association test meta-analysis

Usage: remeta skato [OPTIONS]

Options (* Mandatory))"""");

  opts.add_options()
    ("htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->required()->multitoken(),
      "List of HTP input files. (*)")
    ("ld-prefixes", po::value<vector<string> >()->value_name("FILE1 FILE2 ...")->required()->multitoken(),
      "Prefixes to LD files per cohort. (*)")
    ("cohorts", po::value<vector<string> >()->value_name("NAME1 NAME2 ...")->required()->multitoken(),
      "List of cohort names per file. (*)")
    ("anno-file", po::value<string>()->value_name("FILE")->required(),
      "File with variant annotations. Bgzipped and indexed with 'index-anno'. (*)")
    ("set-list", po::value<string>()->value_name("FILE")->required(),
      "Regenie set-list file. (*)")
    ("mask-def", po::value<string>()->value_name("FILE")->required(),
      "Regenie mask-def file. (*)")
    ("trait-name", po::value<string>()->value_name("NAME")->required(),
      "Name of trait. (*)")
    ("trait-type", po::value<string>()->value_name("TYPE")->required(),
      "One of bt or qt. (*)")
    ("out", po::value<string>()->value_name("PREFIX")->required(),
      "Prefix for output file. (*)")
    ("max-aaf", po::value<double>()->value_name("FLOAT")->default_value(0.01),
      "Maximum allele frequency for a variant to be included in mask.")
    ("min-aac", po::value<int>()->value_name("INT")->default_value(1),
      "Minimum AAC across cohorts for a variant to be included in a mask.")
    ("rho-values",
      po::value<vector<double> >()->value_name("r1 r2 ...")->default_value(
      vector<double>{0, 0.1*0.1, 0.2*0.2, 0.3*0.3, 0.4*0.4, 0.5*0.5, 0.5, 1},
      "0 0.01 0.04 0.09 0.16 0.25 0.5 1")->multitoken(),
      "Rho values for SKATO.")
    ("weight-strategy", po::value<string>()->value_name("STR")->default_value("beta"),
      "Strategy to compute variant weights. One of 'beta' or 'uniform'.")
    ("condition-list", po::value<string>()->value_name("FILE"),
      "File with variants to condition on (one per line).")
    ("condition-htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->multitoken(),
      "List of HTP files with summary statistics of conditional variants per study.")
    ("af-strategy", po::value<string>()->value_name("STR")->default_value("overall"),
      "Strategy to compute variant allele frequences. One of 'overall' or 'max'.")
    ("aaf-file", po::value<string>()->value_name("FILE"),
      "Use precomputed alternate allele frequencies from an external file.")
    ("spa-pval", po::value<double>()->value_name("FLOAT")->default_value(0.05, "0.05"),
      "Apply SPA when the burden p-value is below spa-pval (BTs only, not applied to ACATV).")
    ("spa-ccr", po::value<double>()->value_name("FLOAT")->default_value(0.01, "0.01"),
      "Apply SPA when # cases / # controls < spa-ccr (BTs only, not applied to ACATV).")
    ("chr", po::value<string>()->value_name("CHR"),
      "Run only on chromosome CHR.")
    ("gene", po::value<string>()->value_name("GENE"),
      "If set, run only on this gene.")
    ("extract", po::value<string>()->value_name("FILE"),
      "Include only the variants with IDs listed in this file (one per line).")
    ("exclude", po::value<string>()->value_name("FILE"),
      "Exclude variants with IDs listed in this file (one per line).")
    ("sources", po::value<vector<string> >()->value_name("SOURCE1 SOURCE2 ...")->multitoken(),
      "Only include variants where the info field SOURCE is one of SOURCE1 SOURCE2 ...")
    ("write-mask-snplist", "Write file with list of variants included in each mask.")
    ("recompute-score", "Recompute score statistics from betas and standard errors when missing in input.")
    ("keep-variants-not-in-ld-mat", "Keep variants absent from the LD matrix instead of dropping them.")
    ("threads", po::value<int>()->value_name("INT")->default_value(1),
      "Number of threads to use.")
    ("log-level", po::value<string>()->value_name("LEVEL"), "Set logging level.")
    ("help", "Print this message and exit.");

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try { 
    po::parsed_options parsed = po::command_line_parser(options)
      .options(opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  if (vm.count("log-level") && vm["log-level"].as<string>() == "debug") {
    init_logging("debug", vm["out"].as<string>() + ".skato.log");
  } else if (vm.count("log-level")) {
    cerr << "error: unrecognized --log-level: " + vm["log-level"].as<string>()
         << endl;
    exit(1);
  } else {
    init_logging("info", vm["out"].as<string>() + ".skato.log");
  }

  string condition_file;
  if (vm.count("condition-list") > 0) {
    condition_file = vm["condition-list"].as<string>();
  }
  vector<string> condition_htp;
  if (vm.count("condition-htp") > 0) {
    condition_htp = vm["condition-htp"].as<vector<string> > ();
  }
  string af_file;
  if (vm.count("aaf-file") > 0) {
    af_file = vm["aaf-file"].as<string>();
  }
  string chr;
  if (vm.count("chr") > 0) {
    chr = vm["chr"].as<string>();
  }
  string gene;
  if (vm.count("gene") > 0) {
    gene = vm["gene"].as<string>();
  }
  string extract = "";
  if (vm.count("extract") > 0) {
    extract = vm["extract"].as<string>();
  }
  string exclude = "";
  if (vm.count("exclude") > 0) {
    exclude = vm["exclude"].as<string>();
  }
  vector<string> sources;
  if (vm.count("sources") > 0) {
    sources = vm["sources"].as<vector<string> >();
  }

  log_program_info("skato", options);
 
  run_htp(
    vm["htp"].as<vector<string> >(),
    vm["ld-prefixes"].as<vector<string> >(),
    vm["cohorts"].as<vector<string> >(),
    vm["anno-file"].as<string>(),
    vm["set-list"].as<string>(),
    vm["mask-def"].as<string>(),
    vm["trait-name"].as<string>(),
    vm["trait-type"].as<string>(),
    vm["out"].as<string>() + ".skato",
    {},
    "",
    "",
    true,
    vm["max-aaf"].as<double>(),
    vm["rho-values"].as<vector<double> >(),
    vm["min-aac"].as<int>(),
    vm["weight-strategy"].as<string>(),
    false,
    0,
    0,
    "",
    true,
    condition_file,
    condition_htp,
    vm["af-strategy"].as<string>(),
    af_file,
    vm["spa-pval"].as<double>(),
    vm["spa-ccr"].as<double>(),
    chr,
    gene,
    extract,
    exclude,
    sources,
    vm["threads"].as<int>(),
    false,
    vm.count("write-mask-snplist") > 0,
    false,
    vm.count("recompute-score") > 0,
    vm.count("keep-variants-not-in-ld-mat") > 0
  );
}

void burden(const vector<string>& options) {
  po::options_description opts(R""""(remeta: weighted sum test meta-analysis

Usage: remeta burden [OPTIONS]

Options (* Mandatory))"""");

  opts.add_options()
    ("htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->required()->multitoken(),
      "List of HTP input files. (*)")
    ("ld-prefixes", po::value<vector<string> >()->value_name("FILE1 FILE2 ...")->required()->multitoken(),
      "Prefixes to LD files per cohort. (*)")
    ("cohorts", po::value<vector<string> >()->value_name("NAME1 NAME2 ...")->required()->multitoken(),
      "List of cohort names per file. (*)")
    ("anno-file", po::value<string>()->value_name("FILE")->required(),
      "File with variant annotations. Bgzipped and indexed with 'index-anno'. (*)")
    ("set-list", po::value<string>()->value_name("FILE")->required(),
      "Regenie set-list file. (*)")
    ("mask-def", po::value<string>()->value_name("FILE")->required(),
      "Regenie mask-def file. (*)")
    ("trait-name", po::value<string>()->value_name("NAME")->required(),
      "Name of trait. (*)")
    ("trait-type", po::value<string>()->value_name("TYPE")->required(),
      "One of bt or qt. (*)")
    ("out", po::value<string>()->value_name("PREFIX")->required(),
      "Prefix for output file. (*)")
    ("aaf-bins", po::value<vector<double> >()->value_name("FLOAT FLOAT ...")->multitoken()->default_value(
      vector<double>{1e-4, 1e-3, 5e-3, 1e-2},
      "0.0001 0.001 0.005 0.01"),
      "Allele frequency cutoffs for building masks for burden testing.")
    ("singleton-def", po::value<string>()->value_name("STRING")->default_value("within"),
      "Define singletons for the singleton mask within cohorts or across cohorts. One of 'within', 'across' or 'omit'.")
    ("weight-strategy", po::value<string>()->value_name("STR")->default_value("uniform"),
      "Strategy to compute variant weights for burden testing. One of 'beta' or 'uniform'.")
    ("condition-list", po::value<string>()->value_name("FILE"),
      "File with variants to condition on (one per line).")
    ("condition-htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->multitoken(),
      "List of HTP files with summary statistics of conditional variants per study.")
    ("af-strategy", po::value<string>()->value_name("STR")->default_value("overall"),
      "Strategy to compute variant allele frequences. One of 'overall' or 'max'.")
    ("aaf-file", po::value<string>()->value_name("FILE"),
      "Use precomputed alternate allele frequencies from an external file.")
    ("spa-pval", po::value<double>()->value_name("FLOAT")->default_value(0.05, "0.05"),
      "Apply SPA when the burden p-value is below spa-pval (BTs only, not applied to ACATV).")
    ("spa-ccr", po::value<double>()->value_name("FLOAT")->default_value(0.01, "0.01"),
      "Apply SPA when # cases / # controls < spa-ccr (BTs only, not applied to ACATV).")
    ("chr", po::value<string>()->value_name("CHR"),
      "Run only on chromosome CHR.")
    ("gene", po::value<string>()->value_name("GENE"),
      "If set, run only on this gene.")
    ("extract", po::value<string>()->value_name("FILE"),
      "Include only the variants with IDs listed in this file (one per line).")
    ("exclude", po::value<string>()->value_name("FILE"),
      "Exclude variants with IDs listed in this file (one per line).")
    ("sources", po::value<vector<string> >()->value_name("SOURCE1 SOURCE2 ...")->multitoken(),
      "Only include variants where the info field SOURCE is one of SOURCE1 SOURCE2 ...")
    ("write-cohort-burden-tests", "Compute and store per cohort burden tests.")
    ("write-mask-snplist", "Write file with list of variants included in each mask.")
    ("recompute-score", "Recompute score statistics from betas and standard errors when missing in input.")
    ("keep-variants-not-in-ld-mat", "Keep variants absent from the LD matrix instead of dropping them.")
    ("threads", po::value<int>()->value_name("INT")->default_value(1),
      "Number of threads to use.")
    ("log-level", po::value<string>()->value_name("LEVEL"), "Set logging level.")
    ("help", "Print this message and exit.");

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try { 
    po::parsed_options parsed = po::command_line_parser(options)
      .options(opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  if (vm.count("log-level") && vm["log-level"].as<string>() == "debug") {
    init_logging("debug", vm["out"].as<string>() + ".burden.log");
  } else if (vm.count("log-level")) {
    cerr << "error: unrecognized --log-level: " + vm["log-level"].as<string>()
         << endl;
    exit(1);
  } else {
    init_logging("info", vm["out"].as<string>() + ".burden.log");
  }

  string condition_file;
  if (vm.count("condition-list") > 0) {
    condition_file = vm["condition-list"].as<string>();
  }
  vector<string> condition_htp;
  if (vm.count("condition-htp") > 0) {
    condition_htp = vm["condition-htp"].as<vector<string> > ();
  }
  string af_file;
  if (vm.count("aaf-file") > 0) {
    af_file = vm["aaf-file"].as<string>();
  }
  string chr;
  if (vm.count("chr") > 0) {
    chr = vm["chr"].as<string>();
  }
  string gene;
  if (vm.count("gene") > 0) {
    gene = vm["gene"].as<string>();
  }
  string extract = "";
  if (vm.count("extract") > 0) {
    extract = vm["extract"].as<string>();
  }
  string exclude = "";
  if (vm.count("exclude") > 0) {
    exclude = vm["exclude"].as<string>();
  }
  vector<string> sources;
  if (vm.count("sources") > 0) {
    sources = vm["sources"].as<vector<string> >();
  }

  log_program_info("burden", options);

  run_htp(
    vm["htp"].as<vector<string> >(),
    vm["ld-prefixes"].as<vector<string> >(),
    vm["cohorts"].as<vector<string> >(),
    vm["anno-file"].as<string>(),
    vm["set-list"].as<string>(),
    vm["mask-def"].as<string>(),
    vm["trait-name"].as<string>(),
    vm["trait-type"].as<string>(),
    vm["out"].as<string>() + ".burden",
    vm["aaf-bins"].as<vector<double> >(),
    vm["singleton-def"].as<string>(),
    vm["weight-strategy"].as<string>(),
    false,
    0,
    {},
    0,
    "",
    true,
    0,
    0,
    "",
    true,
    condition_file,
    condition_htp,
    vm["af-strategy"].as<string>(),
    af_file,
    vm["spa-pval"].as<double>(),
    vm["spa-ccr"].as<double>(),
    chr,
    gene,
    extract,
    exclude,
    sources,
    vm["threads"].as<int>(),
    vm.count("write-cohort-burden-tests") > 0,
    vm.count("write-mask-snplist") > 0,
    false,
    vm.count("recompute-score") > 0,
    vm.count("keep-variants-not-in-ld-mat") > 0
  );
}

void htp(const vector<string>& options) {
  po::options_description opts(R""""(remeta: perform burden, SKATO, and ACATV meta-analysis.

Usage: remeta gene [OPTIONS]

Options (* Mandatory))"""");

  opts.add_options()
    ("htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->required()->multitoken(),
      "List of HTP input files. (*)")
    ("ld-prefixes", po::value<vector<string> >()->value_name("FILE1 FILE2 ...")->multitoken(),
      "Prefixes to LD files per cohort. (*)")
    ("cohorts", po::value<vector<string> >()->value_name("NAME1 NAME2 ...")->required()->multitoken(),
      "List of cohort names per file. (*)")
    ("anno-file", po::value<string>()->value_name("FILE")->required(),
      "File with variant annotations. Bgzipped and indexed with 'index-anno'. (*)")
    ("set-list", po::value<string>()->value_name("FILE")->required(),
      "Regenie set-list file. (*)")
    ("mask-def", po::value<string>()->value_name("FILE")->required(),
      "Regenie mask-def file. (*)")
    ("trait-name", po::value<string>()->value_name("NAME")->required(),
      "Name of trait. (*)")
    ("trait-type", po::value<string>()->value_name("TYPE")->required(),
      "One of bt or qt. (*)")
    ("out", po::value<string>()->value_name("PREFIX")->required(),
      "Prefix for output file. (*)")
    ("burden-aaf-bins", po::value<vector<double> >()->value_name("FLOAT FLOAT ...")->multitoken()->default_value(
      vector<double>{1e-4, 1e-3, 5e-3, 1e-2},
      "0.0001 0.001 0.005 0.01"),
      "Allele frequency cutoffs for building masks for burden testing.")
    ("burden-singleton-def", po::value<string>()->value_name("STRING")->default_value("within"),
      "Define singletons for the singleton mask within cohorts or across cohorts. One of 'within', 'across' or 'omit'.")
    ("burden-weight-strategy", po::value<string>()->value_name("STR")->default_value("uniform"),
      "Strategy to compute variant weights for burden testing. One of 'beta' or 'uniform'.")
    ("skip-burden", "Do not run burden testing.")
    ("skato-max-aaf", po::value<double>()->value_name("FLOAT")->default_value(0.01),
      "Maximum allele frequency for a variant to be included in mask for SKATO.")
    ("skato-rho-values",
      po::value<vector<double> >()->value_name("r1 r2 ...")->default_value(
      vector<double>{0, 0.1*0.1, 0.2*0.2, 0.3*0.3, 0.4*0.4, 0.5*0.5, 0.5, 1},
      "0 0.01 0.04 0.09 0.16 0.25 0.5 1")->multitoken(),
      "Rho values for SKATO.")
    ("skato-min-aac", po::value<int>()->value_name("INT")->default_value(1),
      "Minimum AAC across cohorts for a variant to be included in a mask for SKATO.")
    ("skato-weight-strategy", po::value<string>()->value_name("STR")->default_value("beta"),
      "Strategy to compute variant weights for SKATO. One of 'beta' or 'uniform'.")
    ("skip-skato", "Do not run SKATO.")
    ("acatv-max-aaf", po::value<double>()->value_name("FLOAT")->required()->default_value(0.01),
      "Maximum allele frequency for a variant to be included in mask for ACATV.")
    ("acatv-min-aac", po::value<int>()->value_name("INT")->default_value(5),
      "Minimum AAC across cohorts for a variant to be included in a mask for ACATV.")
    ("acatv-weight-strategy", po::value<string>()->value_name("STR")->default_value("beta"),
      "Strategy to compute variant weights for ACATV. One of 'beta' or 'uniform'.")
    ("skip-acatv", "Do not run ACATV.")
    ("condition-list", po::value<string>()->value_name("FILE"),
      "File with variants to condition on (one per line).")
    ("condition-htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->multitoken(),
      "List of HTP files with summary statistics of conditional variants per cohort .")
    ("af-strategy", po::value<string>()->value_name("STR")->default_value("overall"),
      "Strategy to compute variant allele frequences. One of 'overall' or 'max'.")
    ("aaf-file", po::value<string>()->value_name("FILE"),
      "Use precomputed alternate allele frequencies from an external file.")
    ("spa-pval", po::value<double>()->value_name("FLOAT")->default_value(0.05, "0.05"),
      "Apply SPA when the burden p-value is below spa-pval (BTs only, not applied to ACATV).")
    ("spa-ccr", po::value<double>()->value_name("FLOAT")->default_value(0.01, "0.01"),
      "Apply SPA when # cases / # controls < spa-ccr (BTs only, not applied to ACATV).")
    ("chr", po::value<string>()->value_name("CHR"),
      "Run only on chromosome CHR.")
    ("gene", po::value<string>()->value_name("GENE"),
      "If set, run only on this gene.")
    ("extract", po::value<string>()->value_name("FILE"),
      "Include only the variants with IDs listed in this file (one per line).")
    ("exclude", po::value<string>()->value_name("FILE"),
      "Exclude variants with IDs listed in this file (one per line).")
    ("sources", po::value<vector<string> >()->value_name("SOURCE1 SOURCE2 ...")->multitoken(),
      "Only include variants where the info field SOURCE is one of SOURCE1 SOURCE2 ...")
    ("write-cohort-burden-tests", "Compute and store per cohort burden tests (ignores changes to --burden-weight-strategy).")
    ("write-mask-snplist", "Write file with list of variants included in each mask.")
    ("recompute-score", "Recompute score statistics from betas and standard errors when missing in input.")
    ("keep-variants-not-in-ld-mat", "Keep variants absent from the LD matrix instead of dropping them.")
    ("ignore-mask-ld", "Ignore LD between variants in a mask.")
    ("threads", po::value<int>()->value_name("INT")->default_value(1),
      "Number of threads to use.")
    ("log-level", po::value<string>()->value_name("LEVEL"), "Set logging level.")
    ("help", "Print this message and exit.");

  // po::options_description dev_opts("Development options");
  // dev_opts.add_options();

  po::options_description all_opts("All options");
  all_opts.add(opts);//.add(dev_opts);

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try { 
    po::parsed_options parsed = po::command_line_parser(options)
      .options(all_opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  if (vm.count("log-level") && vm["log-level"].as<string>() == "debug") {
    init_logging("debug", vm["out"].as<string>() + ".gene.log");
  } else if (vm.count("log-level")) {
    cerr << "error: unrecognized --log-level: " + vm["log-level"].as<string>()
         << endl;
    exit(1);
  } else {
    init_logging("info", vm["out"].as<string>() + ".gene.log");
  }

  vector<string> ld_prefixes;
  if (vm.count("ld-prefixes") > 0) {
    ld_prefixes = vm["ld-prefixes"].as<vector<string> >();
  }

  string condition_file;
  if (vm.count("condition-list") > 0) {
    condition_file = vm["condition-list"].as<string>();
  }
  vector<string> condition_htp;
  if (vm.count("condition-htp") > 0) {
    condition_htp = vm["condition-htp"].as<vector<string> > ();
  }
  string af_file;
  if (vm.count("aaf-file") > 0) {
    af_file = vm["aaf-file"].as<string>();
  }
  string chr;
  if (vm.count("chr") > 0) {
    chr = vm["chr"].as<string>();
  }
  string gene;
  if (vm.count("gene") > 0) {
    gene = vm["gene"].as<string>();
  }
  string extract = "";
  if (vm.count("extract") > 0) {
    extract = vm["extract"].as<string>();
  }
  string exclude = "";
  if (vm.count("exclude") > 0) {
    exclude = vm["exclude"].as<string>();
  }
  vector<string> sources;
  if (vm.count("sources") > 0) {
    sources = vm["sources"].as<vector<string> >();
  }

  log_program_info("gene", options);

  run_htp(
    vm["htp"].as<vector<string> >(),
    ld_prefixes,
    vm["cohorts"].as<vector<string> >(),
    vm["anno-file"].as<string>(),
    vm["set-list"].as<string>(),
    vm["mask-def"].as<string>(),
    vm["trait-name"].as<string>(),
    vm["trait-type"].as<string>(),
    vm["out"].as<string>() + ".gene",
    vm["burden-aaf-bins"].as<vector<double> >(),
    vm["burden-singleton-def"].as<string>(),
    vm["burden-weight-strategy"].as<string>(),
    vm.count("skip-burden") > 0,
    vm["skato-max-aaf"].as<double>(),
    vm["skato-rho-values"].as<vector<double> >(),
    vm["skato-min-aac"].as<int>(),
    vm["skato-weight-strategy"].as<string>(),
    vm.count("skip-skato") > 0,
    vm["acatv-max-aaf"].as<double>(),
    vm["acatv-min-aac"].as<int>(),
    vm["acatv-weight-strategy"].as<string>(),
    vm.count("skip-acatv") > 0,
    condition_file,
    condition_htp,
    vm["af-strategy"].as<string>(),
    af_file,
    vm["spa-pval"].as<double>(),
    vm["spa-ccr"].as<double>(),
    chr,
    gene,
    extract,
    exclude,
    sources,
    vm["threads"].as<int>(),
    vm.count("write-cohort-burden-tests") > 0,
    vm.count("write-mask-snplist") > 0,
    vm.count("ignore-mask-ld") > 0,
    vm.count("recompute-score") > 0,
    vm.count("keep-variants-not-in-ld-mat") > 0
  );
}

void genep(const vector<string>& options) {
  po::options_description opts(R""""(remeta: merge results and compute meta-analysis GENE_P

Usage: remeta merge [OPTIONS]

Options (* Mandatory))"""");

  opts.add_options()
    ("htp", po::value<vector<string> >()->value_name("HTP1 HTP2 ...")->required()->multitoken(),
      "List of HTP input files. (*)")
    ("out", po::value<string>()->value_name("PREFIX")->required(),
      "Prefix for output file. (*)")
    ("genep-def", po::value<string>()->value_name("FILE"),
      "File with masks to group for GENE_P.")
    ("chr", po::value<string>()->value_name("CHR"),
      "Run only on chromosome CHR.")
    ("burden-model", po::value<string>()->value_name("MODEL")->default_value("REMETA-BURDEN-META"),
      "Model column to collapse burden p-values.")
    ("acatv-model", po::value<string>()->value_name("MODEL")->default_value("REMETA-ACATV-META"),
      "Model column to collapse ACATV p-values.")
    ("skato-model", po::value<string>()->value_name("MODEL")->default_value("REMETA-SKATO-ACAT-META"),
      "Model column to collapse SKATO p-values.")
    ("extract", po::value<string>()->value_name("FILE"),
      "Include only the variants with IDs listed in this file (one per line).")
    ("exclude", po::value<string>()->value_name("FILE"),
      "Exclude variants with IDs listed in this file (one per line).")
    ("include-sbat", "Include SBAT PVMA in GENE_P if available.")
    ("log-level", po::value<string>()->value_name("LEVEL"), "Set logging level.")
    ("help", "Print this message and exit.");

  po::options_description rgc_opts("RGC options");
  rgc_opts.add_options()
    ("drop-regenie-gene", "Filter out regenie's gene-based tests (except possibly SBAT).");

  po::options_description all_opts("All options");
  all_opts.add(opts).add(rgc_opts);

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try {
    po::parsed_options parsed = po::command_line_parser(options)
      .options(all_opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  if (vm.count("log-level") && vm["log-level"].as<string>() == "debug") {
    init_logging("debug", vm["out"].as<string>() + ".merge.log");
  } else if (vm.count("log-level")) {
    cerr << "error: unrecognized --log-level: " + vm["log-level"].as<string>()
         << endl;
    exit(1);
  } else {
    init_logging("info", vm["out"].as<string>() + ".merge.log");
  }

  string genep_file;
  if (vm.count("genep-def") > 0) {
    genep_file = vm["genep-def"].as<string>();
  }
  string chr;
  if (vm.count("chr") > 0) {
    chr = vm["chr"].as<string>();
  }
  string extract = "";
  if (vm.count("extract") > 0) {
    extract = vm["extract"].as<string>();
  }
  string exclude = "";
  if (vm.count("exclude") > 0) {
    exclude = vm["exclude"].as<string>();
  }
  bool include_sbat = false;
  if (vm.count("include-sbat")) {
    include_sbat = true;
  }
  bool drop_regenie_gene = false;
  if (vm.count("drop-regenie-gene")) {
    drop_regenie_gene = true;
  }

  log_program_info("merge", options);

  run_genep(
    vm["htp"].as<vector<string> >(),
    genep_file,
    chr,
    vm["burden-model"].as<string>(),
    vm["acatv-model"].as<string>(),
    vm["skato-model"].as<string>(),
    include_sbat,
    extract,
    exclude,
    drop_regenie_gene,
    vm["out"].as<string>()
  );
}

void ld_inflate(const vector<string>& options) {
  po::options_description opts(R""""(remeta: decompress LD matrices for specific genes

Usage: remeta ld-inflate [OPTIONS]

Options (* Mandatory))"""");

  opts.add_options()
    ("mat-prefix", po::value<string>()->value_name("PREFIX")->required(),
      "Prefix to LD files. (*)")
    ("out", po::value<string>()->value_name("PREFIX")->required(),
      "Prefix to output files. (*)")
    ("gene", po::value<string>()->value_name("GENE")->required(),
      "Only inflate specified gene. (*)")
    ("extract", po::value<string>()->value_name("FILE")->required(),
      "File with list of variants IDs to include in the LD matrix. (*)")
    ("log-level", po::value<string>()->value_name("LEVEL"), "Set logging level.")
    ("help", "Print this message and exit.");

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try {
    po::parsed_options parsed = po::command_line_parser(options)
      .options(opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  if (vm.count("log-level") && vm["log-level"].as<string>() == "debug") {
    init_logging("debug", vm["out"].as<string>() + ".ld_inflate.log");
  } else if (vm.count("log-level")) {
    cerr << "error: unrecognized --log-level: " + vm["log-level"].as<string>()
         << endl;
    exit(1);
  } else {
    init_logging("info", vm["out"].as<string>() + ".ld_inflate.log");
  }

  log_program_info("ld-inflate", options);

  run_ld_inflate(vm["mat-prefix"].as<string>(),
                 vm["out"].as<string>(),
                 vm["gene"].as<string>(),
                 vm["extract"].as<string>());
}

void ld_deflate(const vector<string>& options) {
  po::options_description opts(R""""(remeta: convert LD matrices to remeta format

Usage: remeta ld-deflate [OPTIONS]

Options (* Mandatory))"""");

  opts.add_options()
    ("ld-file", po::value<string>()->value_name("FILE")->required(),
      "Prefix to LD file. (*)")
    ("sample-size", po::value<int>()->value_name("INT")->required(),
      "Sample size used to compute LD. (*)")
    ("sparsity-threshold", po::value<double>()->value_name("FLOAT")->required(),
      "Drop corr where abs(corr^2) < sparsity threshold used for sparsity. (*)")
    ("out", po::value<string>()->value_name("PREFIX")->required(),
      "Prefix to output files. (*)")
    ("log-level", po::value<string>()->value_name("LEVEL"), "Set logging level.")
    ("help", "Print this message and exit.");

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try {
    po::parsed_options parsed = po::command_line_parser(options)
      .options(opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  if (vm.count("log-level") && vm["log-level"].as<string>() == "debug") {
    init_logging("debug", vm["out"].as<string>() + ".ld_deflate.log");
  } else if (vm.count("log-level")) {
    cerr << "error: unrecognized --log-level: " + vm["log-level"].as<string>()
         << endl;
    exit(1);
  } else {
    init_logging("info", vm["out"].as<string>() + ".ld_deflate.log");
  }

  log_program_info("ld-deflate", options);

  run_ld_deflate(vm["ld-file"].as<string>(),
                 vm["sample-size"].as<int>(),
                 vm["sparsity-threshold"].as<double>(),
                 vm["out"].as<string>());
}

void compute_ref_ld(const vector<string>& options) {
  po::options_description opts(R""""(remeta: compute reference LD matrices from plink2 pgen/pvar/psam files.

Usage: remeta compute-ref-ld [OPTIONS]

Options (* Mandatory))"""");

  opts.add_options()
    ("target-pfile", po::value<string>()->value_name("PREFIX")->required(),
      "Prefix to pgen/pvar/psam files in target regions (typically WES). (*)")
    ("gene-list", po::value<string>()->value_name("FILE")->required(),
      "List of genes to include in LD matrix in four column format: GENE_NAME CHR START END. (*)")
    ("chr", po::value<string>()->value_name("CHR")->required(),
      "Chromosome to run. (*)")
    ("out", po::value<string>()->value_name("PREFIX")->required(),
      "Prefix to output files. (*)")
    ("buffer-pfile", po::value<string>()->value_name("PREFIX"),
      "Prefix to pgen/pvar/psam files to use for buffer regions (typically array or imputed genotypes).")
    ("buffer-mb", po::value<double>()->value_name("FLOAT"),
      "Buffer in Mb around each gene to search for variants to include LD file.")
    ("buffer-cm", po::value<double>()->value_name("FLOAT"),
      "Buffer in cM around each gene to search for variants to include in LD file (requires a genetic map).")
    ("genetic-map", po::value<string>()->value_name("FILE"),
      "Path to genetic map in three column format: POS CHR CM. Required for --buffer-cm option.")
    ("target-r2", po::value<double>()->value_name("FLOAT")->default_value(1e-4),
      "Drop target (gene) LD matrix entries where r2 < target_r2.")
    ("buffer-r2", po::value<double>()->value_name("FLOAT")->default_value(1e-4),
      "Drop buffer (conditional) LD matrix entries where r2 < buffer_r2.")
    ("float-size", po::value<int>()->value_name("INT")->default_value(1),
      "Size of float in bytes to store LD of buffer variants. Possible values: 1, 2, or 4")
    ("block-size", po::value<int>()->value_name("INT")->default_value(2048),
      "Number of genotypes loaded into memory = 2*block_size.")
    ("threads", po::value<int>()->value_name("INT")->default_value(1),
      "Number of threads for computation.")
    ("target-extract", po::value<string>()->value_name("FILE"),
      "Extract list of variants to include in target regions (e.g. exonic variants)")
    ("target-exclude", po::value<string>()->value_name("FILE"),
      "Exclude list of variants from target file (e.g. non-coding variants)")
    ("target-keep", po::value<string>()->value_name("FILE"),
      "File with list of samples to keep (one sample per line, matching columns of psam file)")
    ("target-remove", po::value<string>()->value_name("FILE"),
      "File with list of samples to remove (one per line, matching columns of psam file)")
    ("buffer-extract", po::value<string>()->value_name("FILE"),
      "Extract list of variants to include in buffer regions (e.g. imputed variants)")
    ("buffer-exclude", po::value<string>()->value_name("FILE"),
      "Exclude list of variants from buffer file (e.g. low quality variants)")
    ("skip-buffer", "Exclude all buffer variants from LD calculation.")
    ("log-level", po::value<string>()->value_name("LEVEL"), "Set logging level.")
    ("help", "Print this message and exit.");

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try {
    po::parsed_options parsed = po::command_line_parser(options)
      .options(opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  if (vm.count("log-level") && vm["log-level"].as<string>() == "debug") {
    init_logging("debug", vm["out"].as<string>() + ".compute_ref_ld.log");
  } else if (vm.count("log-level")) {
    cerr << "error: unrecognized --log-level: " + vm["log-level"].as<string>()
         << endl;
    exit(1);
  } else {
    init_logging("info", vm["out"].as<string>() + ".compute_ref_ld.log");
  }

  log_program_info("compute-ref-ld", options);

  string buffer_pfile = "";
  if (vm.count("buffer-pfile") > 0) {
    buffer_pfile = vm["buffer-pfile"].as<string>();
  }
  double buffer_mb = 0;
  if (vm.count("buffer-mb") > 0) {
    buffer_mb = vm["buffer-mb"].as<double>();
  }
  double buffer_cm = 0;
  if (vm.count("buffer-cm") > 0) {
    buffer_cm = vm["buffer-cm"].as<double>();
  }
  string genetic_map;
  if (vm.count("genetic-map")) {
    genetic_map = vm["genetic-map"].as<string>();
  }
  string target_extract = "";
  if (vm.count("target-extract") > 0) {
    target_extract = vm["target-extract"].as<string>();
  }
  string target_exclude = "";
  if (vm.count("target-exclude")) {
    target_exclude = vm["target-exclude"].as<string>();
  }
  string target_keep = "";
  if (vm.count("target-keep") > 0) {
    target_keep = vm["target-keep"].as<string>();
  }
  string target_remove = "";
  if (vm.count("target-remove")) {
    target_remove = vm["target-remove"].as<string>();
  }
  string buffer_extract = "";
  if (vm.count("buffer-extract") > 0) {
    buffer_extract = vm["buffer-extract"].as<string>();
  }
  string buffer_exclude = "";
  if (vm.count("buffer-exclude")) {
    buffer_exclude = vm["buffer-exclude"].as<string>();
  }
  bool skip_buffer = false;
  if (vm.count("skip-buffer")) {
    skip_buffer = true;
  }

  try {
    run_compute_ref_ld(
      vm["target-pfile"].as<string>(),
      vm["gene-list"].as<string>(),
      vm["chr"].as<string>(),
      vm["out"].as<string>(),
      buffer_pfile,
      buffer_mb,
      buffer_cm,
      genetic_map,
      vm["target-r2"].as<double>(),
      vm["buffer-r2"].as<double>(),
      vm["float-size"].as<int>(),
      vm["block-size"].as<int>(),
      vm["threads"].as<int>(),
      target_extract,
      target_exclude,
      target_keep,
      target_remove,
      buffer_extract,
      buffer_exclude,
      skip_buffer
    );
  } catch (std::exception& e) {
    log_error(e.what(), 1);
  }
}

void index_anno(const vector<string>& options) {
  po::options_description opts(R""""(remeta: index annotation files.

Usage: remeta index-anno [OPTIONS]

Options (* Mandatory))"""");

  opts.add_options()
    ("file", po::value<string>()->value_name("FILE")->required(),
      "Path to annotation file. (*)")
    ("stride", po::value<int>()->value_name("INT")->default_value(100000),
      "Length in bases between each index pointer.")
    ("help", "Print this message and exit.");

  if (options.size() < 1) {
    cerr << opts;
    exit(1);
  }

  po::variables_map vm;
  try {
    po::parsed_options parsed = po::command_line_parser(options)
      .options(opts)
      .run();
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    exit(1);
  }

  init_logging("info", "");

  log_program_info("index-anno", options);

  run_index_anno(vm["file"].as<string>(), vm["stride"].as<int>());
}

int main(int argc, char *argv[]) {
  // supress htslib logging
  hts_set_log_level(HTS_LOG_OFF);

  string version = "remeta " VERSION_NUMBER;
  string description = version + "\n" +
      "Copyright (c) 2024 Tyler Joseph, Joelle Mbatchou, Jonathan Marchini.\nDistributed under the MIT License.\n";

  #if defined(WITH_MKL)
    description += "Using Intel MKL with Eigen.\n";
  #elif defined(WITH_OPENBLAS)
    description += "Using BLAS/LAPACK from OpenBLAS with Eigen.\n";
  #endif
  description += "\nUsage: remeta [OPTIONS] COMMAND [ARGS]...\n\nOptions";

  po::options_description main_opts(description);
  main_opts.add_options()
    ("help,h", "Print this message and exit.")
    ("version,v", "Print version.");

  string commands =
R""""(
Commands:

* Meta-Analyses
  pvma         Perform p-value meta-analysis.
  esma         Perform effect size meta-analysis.
  gene         Perform burden, SKATO, and ACATV meta-analysis.
  merge        Merge results and compute meta-analysis GENEP.

* Utilities
  compute-ref-ld   Compute a reference LD matrix from plink2 pgen/pvar/psam files.
  index-anno       Create an index from a bgzipped REGENIE annotation file.
)"""";

  // these need to mirror positional arguments for
  // positional arguments to be parsed correctly
  po::options_description hidden_opts("");
  hidden_opts.add_options()
    ("command", po::value<string>(), "command to execute")
    ("subargs", po::value<vector<string> >(), "arguments for command");

  po::options_description all_opts;
  all_opts.add(main_opts).add(hidden_opts);

  po::positional_options_description pos;
  pos.add("command", 1)
     .add("subargs", -1);

  po::variables_map vm;

  po::parsed_options parsed = po::command_line_parser(argc, argv)
    .options(all_opts)
    .positional(pos)
    .allow_unregistered()
    .run();

  try {
    po::store(parsed, vm);
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    cerr << "error: " << e.what() << endl;
    return 1;
  }

  if (vm.count("help") && argc <= 2) {
    cerr << main_opts
         << commands
         << endl;
    return 1;
  }
  if (vm.count("version")) {
    cout << VERSION_NUMBER << endl;
    return 0;
  }

  if (!vm.count("command")) {
    cerr << "error: missing command"
         << endl << endl
         << main_opts
         << commands;
    return 1;
  }

  string cmd = vm["command"].as<string>();
  vector<string> subargs = po::collect_unrecognized(parsed.options, po::include_positional);
  subargs.erase(subargs.begin());
  if (cmd == "pvma") {
    pvma(subargs);
  } else if (cmd == "esma") {
    esma(subargs);
  } else if (cmd == "acatv") {
    acatv(subargs);
  } else if (cmd == "skato") {
    skato(subargs);
  } else if (cmd == "burden") {
    burden(subargs);
  } else if (cmd == "gene" || cmd == "htp") {
    htp(subargs);
  } else if (cmd == "merge" || cmd == "genep") {
    genep(subargs);
  } else if (cmd == "ld-deflate") {
    ld_deflate(subargs);
  } else if (cmd == "ld-inflate") {
    ld_inflate(subargs);
  } else if (cmd == "compute-ref-ld") {
    compute_ref_ld(subargs);
  } else if (cmd == "index-anno") {
    index_anno(subargs);
  } else {
    cerr << "error: unrecognized command " << cmd << endl;
    return 1;
  }
  return 0;
}
