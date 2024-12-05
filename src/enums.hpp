#ifndef REMETA_ENUMS_H
#define REMETA_ENUMS_H

enum trait_type_e {
  QT,
  BT
};

// used for burden testing
enum singleton_def_e {
  WITHIN,
  ACROSS,
  OMIT
};

// used in skato_meta_analyzer.hpp
enum af_strategy_e {
  USE_OVERALL_AF,
  USE_MAX_AF,
  USE_EXTERNAL_AF
};

// used in skato_meta_analyzer.hpp
enum weight_strategy_e {
  USE_BETA_WEIGHTS,
  USE_UNIFORM_WEIGHTS
};

enum pvma_method_e {
  STOUFFERS,
  FISHERS
};

#endif