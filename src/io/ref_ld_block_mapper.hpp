#ifndef REF_LD_BLOCK_MAPPER_H
#define REF_LD_BLOCK_MAPPER_H

#include <string>
#include <unordered_map>
#include <vector>

typedef std::string gene_id;

struct ld_block_t {
  int start_bp;
  int end_bp;

  inline bool operator<(const ld_block_t& other) const {
    return this->start_bp < other.start_bp;
  }

  inline bool operator==(const ld_block_t& other) const {
    return this->start_bp == other.start_bp
           && this->end_bp == other.end_bp;
  }
};

class LDBlockPair {
 public:
  LDBlockPair(const ld_block_t& ld_block1, const ld_block_t& ld_block2);

  bool operator<(const LDBlockPair& other) const;

  bool operator==(const LDBlockPair& other) const;

  int min_pos() { return this->ld_block1.start_bp; }

  int max_pos() { return this->ld_block2.end_bp; }

  ld_block_t ld_block1;
  ld_block_t ld_block2;
};

// https://en.cppreference.com/w/cpp/utility/hash
template<>
struct std::hash<ld_block_t> {
  size_t operator()(ld_block_t const& ld_block) const noexcept {
    size_t h1 = std::hash<int>{}(ld_block.start_bp);
    size_t h2 = std::hash<int>{}(ld_block.end_bp);
    return h1 ^ (h2 << 1); 
  }
};

template<>
struct std::hash<LDBlockPair> {
  size_t operator()(LDBlockPair const& block_pair) const noexcept {
    size_t h1 = std::hash<ld_block_t>{}(block_pair.ld_block1);
    size_t h2 = std::hash<ld_block_t>{}(block_pair.ld_block2);
    return h1 ^ (h2 << 1);
  }
};

struct ld_block_map_t {
  // sorted list of LD block pairs along the chromosome
  std::vector<LDBlockPair> sorted_ld_block_pairs;

  // maps an ld_block_t to the genes with overlapping LD windows
  std::unordered_map<ld_block_t, std::vector<gene_id> > ld_block_genes;

  // maps a gene to a list of LD blocks for the LD window of that gene
  std::unordered_map<gene_id, std::vector<ld_block_t> > gene_ld_blocks;

  int nblock_pairs(const gene_id& gene) {
    int n = (int)this->gene_ld_blocks[gene].size();
    if (n == 1) {
      return 1;
    } else {
      return n*(n-1)/2;
    }
  }
};

class RefLDBlockMapper {
 public:
  RefLDBlockMapper() {};

  // The genelist_file has four columns
  //     GENE_NAME    CHROM     START_BP     END_BP
  ld_block_map_t map_ld_blocks_mb(const std::string& genelist_file,
                                  const std::string& chrom,
                                  const double& buffer_mb);

  ld_block_map_t map_ld_blocks_cm(const std::string& genelist_file,
                                  const std::string& genetic_map_file,
                                  const std::string& chrom,
                                  const double& buffer_cm);

 private:
  struct ld_block_breakpoint_t {
    gene_id gene;
    int pos;
    bool is_start;

    inline bool operator<(const ld_block_breakpoint_t& other) {
      return this->pos < other.pos;
    }
  };

  ld_block_map_t map_ld_blocks(std::vector<ld_block_breakpoint_t>& breakpoints);

};

#endif