#include "variant_meta_analyzer.hpp"

#include <chrono>

#include "../io/htpv4_reader.hpp"
#include "../io/htpv4_writer.hpp"
#include "../htpv4_pos.hpp"
#include "../logging.hpp"
#include "../variant_filter.hpp"

void run_variant_meta_analysis(VariantMetaAnalyzer& meta,
                               VariantFilter& vf,
                               vector<HTPv4Reader>& in_files,
                               HTPv4Writer& out,
                               const string& chr) {
  if (chr != "") {
    for (size_t i = 0; i < in_files.size(); ++i) {
      try {
        in_files[i].seek(chr, 0);
      } catch (const runtime_error& err) {
        continue;
      }
    }
  }

  std::chrono::duration<double> processing_time(0);
  std::chrono::time_point<std::chrono::steady_clock> t;

  vector<HTPv4Pos> htp_positions;
  vector<htpv4_record_t> htp_records;
  for (size_t i = 0; i < in_files.size(); ++i) {
    htpv4_record_t rec;
    if (!in_files[i].eof()) {
      rec = in_files[i].readrec();
      htp_positions.push_back(HTPv4Pos(rec.chr, rec.pos));
    } else {
      htp_positions.push_back(HTPV4_EOF);
    }
    htp_records.push_back(rec);
  }
  pair<int, HTPv4Pos> min_htp_pos = get_min_pos(htp_positions);

  htpv4_record_t *rec;
  int next_record_id;
  string last_chr = "";
  int records_added = 0;
  t = std::chrono::steady_clock::now();
  while (min_htp_pos.second != HTPV4_EOF) {
    next_record_id = min_htp_pos.first;
    rec = &htp_records[next_record_id];

    if (last_chr != "" && rec->chr != last_chr) {
      for (const htpv4_record_t& rec: meta.meta_analyze_before(last_chr, 1000000000)) {
        out.writerec(rec);
      }
      meta.clear_before(last_chr, 1000000000);

      processing_time = std::chrono::steady_clock::now() - t;
      t = std::chrono::steady_clock::now();
      log_info("\t* completed in " + to_string(processing_time.count()) + " seconds");
      log_info("processing chromosome " + rec->chr);
    } else if (last_chr != rec->chr) {
      log_info("processing chromosome " + rec->chr);   
    }
    last_chr = rec->chr;

    if (chr == "" || rec->chr == chr) {
      if (vf.include_htp_record(*rec)) {
        meta.add_line(htp_records[next_record_id], next_record_id);;
        ++records_added;
      }
      for (const htpv4_record_t& rec: meta.meta_analyze_before(rec->chr, rec->pos)) {
        out.writerec(rec);
      }
      meta.clear_before(rec->chr, rec->pos);
    }

    if (in_files[next_record_id].eof()) {
      htp_positions[next_record_id] = HTPV4_EOF;
    } else {
      htp_records[next_record_id] = in_files[next_record_id].readrec();
      HTPv4Pos rec_pos = HTPv4Pos(htp_records[next_record_id].chr, htp_records[next_record_id].pos);
      if (htp_positions[next_record_id] != HTPV4_EOF && rec_pos < htp_positions[next_record_id]) {
        log_error("HTP files must be sorted by chromosome and position: HTP file " + to_string(next_record_id) + " is not sorted", 1);
      }
      htp_positions[next_record_id] = rec_pos;
    }
    min_htp_pos = get_min_pos(htp_positions);
  }

  for (const htpv4_record_t& rec: meta.meta_analyze_before(last_chr, 1000000000)) {
    out.writerec(rec);
  }

  processing_time = std::chrono::steady_clock::now() - t;
  log_info("\t* completed in " + to_string(processing_time.count()) + " seconds");
  t = std::chrono::steady_clock::now();

  out.close();
  log_info("finished!");
}