#include "set_meta_analyzer.hpp"

#include "../io/htpv4_reader.hpp"
#include "../io/htpv4_writer.hpp"
#include "../htpv4_pos.hpp"
#include "../logging.hpp"
#include "../variant_filter.hpp"

void run_set_meta_analysis(SetMetaAnalyzer& meta,
                           VariantFilter& vf,
                           vector<HTPv4Reader>& in_files,
                           const string& set_list_file,
                           HTPv4Writer& out,
                           const string& chr,
                           const string& gene) {
  log_info("running meta-analysis");
  if (chr != "") {
    for (size_t i = 0; i < in_files.size(); ++i) {
      if (in_files[i].eof()) {
        continue;
      }
      try {
        in_files[i].seek(chr, 0);
      } catch (const runtime_error& err) {
        continue;
      }
    }
  }

  vector<HTPv4Pos> htp_positions;
  for (size_t i = 0; i < in_files.size(); ++i) {
    if (!in_files[i].eof()) {
      htpv4_record_t rec = in_files[i].readrec();
      htp_positions.push_back(HTPv4Pos(rec.chr, rec.pos));
      if (vf.include_htp_record(rec)) {
        meta.add_line(rec, i);
      }
    } else {
      htp_positions.push_back(HTPV4_EOF);
    }
  }

  GeneSetReader gene_reader(set_list_file);
  HTPv4Pos gene_start = HTPV4_NULL_POS;
  HTPv4Pos gene_end = HTPV4_NULL_POS;
  HTPv4Pos prv_gene_start = HTPV4_NULL_POS;
  string prv_gene_chr = "";
  htpv4_record_t rec;
  while (!gene_reader.eof()) {
    Gene g = gene_reader.read_gene();

    // if gene option is set, skip all genes other than the one requested
    if (gene != "" && g.get_name() != gene) {
      continue;
    }

    // if chr option is set, skip all genes that are not on the requested chromosome
    if (chr != "" && g.get_chrom() != chr) {
      continue;
    }
    gene_start = HTPv4Pos(g.get_chrom(), g.get_start());
    gene_end = HTPv4Pos(g.get_chrom(), g.get_end());

    if (prv_gene_start != HTPV4_NULL_POS && gene_start < prv_gene_start) {
      log_error("set-list file " + set_list_file + " must be sorted by start position", 1);
    }

    // clear htpv4 records that are no longer required
    if (prv_gene_chr != "" && g.get_chrom() != prv_gene_chr) {
      meta.clear_before(prv_gene_chr, 1000000000);
    } else if (g.get_chrom() == prv_gene_chr && gene_start > prv_gene_start) {
      meta.clear_before(prv_gene_chr, prv_gene_start.pos);
    }

    // add htpv4 records required for the current gene
    pair<int, HTPv4Pos> arg_min = get_min_pos(htp_positions);
    if (gene_start.int_chr < arg_min.second.int_chr) {
      continue;
    }

    while (arg_min.second != HTPV4_EOF && arg_min.second <= gene_end) {
      if (in_files[arg_min.first].eof()) {
        htp_positions[arg_min.first] = HTPV4_EOF;
        arg_min = get_min_pos(htp_positions);
        continue;
      }

      rec = in_files[arg_min.first].readrec();
      HTPv4Pos rec_pos = HTPv4Pos(rec.chr, rec.pos);
      if (htp_positions[arg_min.first] != HTPV4_EOF && rec_pos < htp_positions[arg_min.first]) {
        log_error("HTP files must be sorted by chromosome and position: HTP file " + to_string(arg_min.first) + " is not sorted", 1);
      }

      if (rec_pos >= gene_start && vf.include_htp_record(rec)) {
        meta.add_line(rec, arg_min.first);
      }

      htp_positions[arg_min.first] = HTPv4Pos(rec.chr, rec.pos);
      arg_min = get_min_pos(htp_positions);
    }

    // write results to the output file
    log_info("meta-analyzing " + g.get_name());
    for (const htpv4_record_t& rec: meta.meta_analyze_gene(g)) {
      out.writerec(rec);
    }

    prv_gene_start = gene_start;
    prv_gene_chr = g.get_chrom();
  }
  out.close();
  log_info("finished!");
}