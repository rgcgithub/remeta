#ifndef SPA_INCLUDED_H
#define SPA_INCLUDED_H

#include "../lapack_complex.hpp"

#include <vector>

#include <Eigen/Dense>

namespace stat::spa {
  const double SPA_FAILED = -2147483648;

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

  double compute_spa_pval_from_gc_per_study(const Eigen::VectorXd& scores,
                                            const Eigen::VectorXd& stds,
                                            const Eigen::VectorXd& nhom_ref,
                                            const Eigen::VectorXd& nhet,
                                            const Eigen::VectorXd& nhom_alt,
                                            const Eigen::VectorXd& ncases,
                                            const Eigen::VectorXd& ncontrols,
                                            bool one_sided);

  double compute_spa_chival_from_gc_per_study(const Eigen::VectorXd& scores,
                                              const Eigen::VectorXd& stds,
                                              const Eigen::VectorXd& nhom_ref,
                                              const Eigen::VectorXd& nhet,
                                              const Eigen::VectorXd& nhom_alt,
                                              const Eigen::VectorXd& ncases,
                                              const Eigen::VectorXd& ncontrols,
                                              bool one_sided);
}

#endif