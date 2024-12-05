#include "ref_ld_block_mapper.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;

#include "bgz_reader.hpp"
#include "ref_ld_genetic_map.hpp"
#include "../logging.hpp"
#include "../util.hpp"

void parse_genelist_entry(string& gene,
                          string& chr,
                          int& start,
                          int& end,
                          const string& line,
                          const string& genelist_file) {
  vector<string> row = util::str_split(line, "\t ");
  if (row.size() < 4) {
    throw runtime_error("genelist_file " + genelist_file + " has fewer than 3 columns");
  }

  gene = row[0];
  chr = row[1];
  start = stoi(row[2]);
  end = stoi(row[3]);

  if (to_string(start) != row[2] || to_string(end) != row[3]) {
    log_error("bad format for gene-list file (does your file have columns GENE CHR START END ?)", 1);
  }
}

LDBlockPair::LDBlockPair(const ld_block_t& ld_block1,
                         const ld_block_t& ld_block2) {
  if (ld_block1 < ld_block2) {
    this->ld_block1 = ld_block1;
    this->ld_block2 = ld_block2;
  } else {
    this->ld_block1 = ld_block2;
    this->ld_block2 = ld_block1;
  }
}

bool LDBlockPair::operator<(const LDBlockPair& other) const {
  return (this->ld_block1 < other.ld_block1)
    || (this->ld_block1 == other.ld_block1 && this->ld_block2 < other.ld_block2);
}

bool LDBlockPair::operator==(const LDBlockPair& other) const {
  return this->ld_block1 == other.ld_block1 
    && this->ld_block2 == other.ld_block2;
}

ld_block_map_t RefLDBlockMapper::map_ld_blocks_mb(const string& genelist_file, 
                                                  const string& chrom,
                                                  const double& buffer_mb) {
  BgzReader genelist(genelist_file);
  string line;
  string gene;
  string chr;
  int start;
  int end;
  int window = (int)ceil(1E6*buffer_mb);
  vector<ld_block_breakpoint_t> breakpoints;
  while (!genelist.eof()) {
    line = genelist.readline();
    parse_genelist_entry(gene, chr, start, end, line, genelist_file);
    if (chr != chrom) {
      continue;
    }
    breakpoints.push_back({gene, max(0, start - window), true});
    breakpoints.push_back({gene, end + window, false});
  }
  if (breakpoints.size() == 0) {
    throw runtime_error(
      "failed to find LD blocks in " + genelist_file  
      + " for chr " + chrom
    );
  }
  return this->map_ld_blocks(breakpoints);
}

ld_block_map_t RefLDBlockMapper::map_ld_blocks_cm(const std::string& genelist_file,
                                                  const std::string& genetic_map_file,
                                                  const std::string& chrom,
                                                  const double& buffer_cm) {
  BgzReader genelist(genelist_file);
  RefLDGeneticMap genetic_map(genetic_map_file);

  string line;
  string gene;
  string chr;
  int start;
  int end;
  vector<ld_block_breakpoint_t> breakpoints;
  while (!genelist.eof()) {
    line = genelist.readline();
    parse_genelist_entry(gene, chr, start, end, line, genelist_file);
    if (chr != chrom) {
      continue;
    }
    int start_pos = round(
      max(0.,
        genetic_map.cm_to_bp(
          genetic_map.bp_to_cm(start) - buffer_cm
        )
      )
    );
    int end_pos = round(
      genetic_map.cm_to_bp(
        genetic_map.bp_to_cm(end) + buffer_cm
      )
    );
    if (end_pos < 0) {
      throw runtime_error(
        "overflow encountered mapping cm to bp"
      );
    }

    breakpoints.push_back({gene, start_pos, true});
    breakpoints.push_back({gene, end_pos, false});
  }
  if (breakpoints.size() == 0) {
    throw runtime_error(
      "failed to find LD blocks in " + genelist_file  
      + " for chr " + chrom
    );
  }
  return this->map_ld_blocks(breakpoints);
}

ld_block_map_t RefLDBlockMapper::map_ld_blocks(vector<ld_block_breakpoint_t>& breakpoints) {
  std::sort(breakpoints.begin(), breakpoints.end());
  set<LDBlockPair> ld_block_pairs;
  ld_block_map_t ld_block_map;
  vector<gene_id> current_genes;
  for (size_t i = 0; i < breakpoints.size()-1; ++i) {
    if (breakpoints[i].is_start) {
      current_genes.push_back(breakpoints[i].gene);
    } else if (!breakpoints[i].is_start) {
      current_genes.erase(
        std::find(
          current_genes.begin(),
          current_genes.end(),
          breakpoints[i].gene
        )
      );
    }

    if (current_genes.size() > 0 && breakpoints[i].pos != breakpoints[i+1].pos) {
      ld_block_t ld_block{ breakpoints[i].pos, breakpoints[i+1].pos };
      ld_block_pairs.insert(LDBlockPair(ld_block, ld_block));
      
      for (const gene_id& gene : current_genes) {
        for (const ld_block_t& cross_ld_block : ld_block_map.gene_ld_blocks[gene]) {
          ld_block_pairs.insert(LDBlockPair(cross_ld_block, ld_block));
        }

        ld_block_map.ld_block_genes[ld_block].push_back(gene);
        ld_block_map.gene_ld_blocks[gene].push_back(ld_block);
      }
    }
  }

  for (const LDBlockPair& ld_block_pair : ld_block_pairs) {
    ld_block_map.sorted_ld_block_pairs.push_back(ld_block_pair);
  }
  return ld_block_map;
}