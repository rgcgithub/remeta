#ifndef SPA_INCLUDED_H
#define SPA_INCLUDED_H

namespace stat::spa {
  const double SPA_FAILED = -9;

  double compute_spa_pval_from_sum_stats(const double& score,
                                         const double& std,
                                         const double& burden_aaf,
                                         const double& num_cases,
                                         const double& num_controls,
                                         const bool& one_sided);

  double compute_spa_chival_from_sum_stats(const double& score,
                                           const double& std,
                                           const double& burden_aaf,
                                           const double& num_cases,
                                           const double& num_controls,
                                           const bool& one_sided);

  double compute_spa_pval_from_geno_counts(const double& score,
                                           const double& std,
                                           const double& nhom_ref,
                                           const double& nhet,
                                           const double& nhom_alt,
                                           const double& num_cases,
                                           const double& num_controls,
                                           const bool& one_sided);

  double compute_spa_chival_from_geno_counts(const double& score,
                                             const double& std,
                                             const double& nhom_ref,
                                             const double& nhet,
                                             const double& nhom_alt,
                                             const double& num_cases,
                                             const double& num_controls,
                                             const bool& one_sided);
}

#endif