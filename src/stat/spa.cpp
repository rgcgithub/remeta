#include "spa.hpp"

#include <cmath>
#include <utility>
using std::pair;

#include <iostream>
using namespace std;

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/roots.hpp>

namespace stat::spa {

struct sum_stat_spa_params {
  double N;
  double mu;

  double Ghom_ref;
  double Ghet;
  double Ghom_alt;

  double nhom_ref;
  double nhet;
  double nhom_alt;
};

inline double K(const double& t, const sum_stat_spa_params& params) {
  return params.nhom_ref*(log(1 - params.mu + params.mu*exp(t*params.Ghom_ref)) - t*params.mu*params.Ghom_ref)
    + params.nhet*(log(1 - params.mu + params.mu*exp(t*params.Ghet)) - t*params.mu*params.Ghet)
    + params.nhom_alt*(log(1 - params.mu + params.mu*exp(t*params.Ghom_alt)) - t*params.mu*params.Ghom_alt);
}

inline double K1(const double& t, const sum_stat_spa_params& params) {
  return params.nhom_ref*params.Ghom_ref*params.mu / (params.mu + (1-params.mu)*exp(-t*params.Ghom_ref))
    + params.nhet*params.Ghet*params.mu / (params.mu + (1-params.mu)*exp(-t*params.Ghet))
    + params.nhom_alt*params.Ghom_alt*params.mu / (params.mu + (1-params.mu)*exp(-t*params.Ghom_alt))
    - params.nhom_ref*params.Ghom_ref*params.mu
    - params.nhet*params.Ghet*params.mu
    - params.nhom_alt*params.Ghom_alt*params.mu;
}

inline double K2(const double& t, const sum_stat_spa_params& params) {
  return params.nhom_ref*pow(params.Ghom_ref, 2)*params.mu*(1-params.mu)*exp(-t*params.Ghom_ref) / pow(params.mu + (1-params.mu)*exp(-t*params.Ghom_ref), 2)
    + params.nhet*pow(params.Ghet, 2)*params.mu*(1-params.mu)*exp(-t*params.Ghet) /  pow(params.mu + (1-params.mu)*exp(-t*params.Ghet), 2)
    + params.nhom_alt*pow(params.Ghom_alt, 2)*params.mu*(1-params.mu)*exp(-t*params.Ghom_alt) /  pow(params.mu + (1-params.mu)*exp(-t*params.Ghom_alt), 2);
}

// inline double solve_for_delta(const double& t, const sum_stat_spa_params& params) {
//   std::uintmax_t max_iter = 100;
//   return boost::math::tools::newton_raphson_iterate(
//     [&t, &params](const double& d) {
//       return std::make_pair(K1(d, params) - t, K2(d, params));
//     },
//     t,
//     -100.0,
//     100.0,
//     10,
//     max_iter
//   );
// }

inline double solve_for_delta(const double& t, const sum_stat_spa_params& params) {
  int max_iter = 100;
  double tol = 1e-6;
  double t_min = 0;
  double t_max = 0;
  if (params.Ghom_ref > 0) {
    t_min -= params.nhom_ref * params.Ghom_ref;
    t_max += params.nhom_ref * params.Ghom_ref;
  }
  if (params.Ghet > 0) {
    t_min -= params.nhet * params.Ghet;
    t_max += params.nhet * params.Ghet;
  }
  if (params.Ghom_alt > 0) {
    t_min -= params.nhom_alt * params.Ghom_alt;
    t_max += params.nhom_alt * params.Ghom_alt;
  }

  if (t < t_min || t > t_max) {
    return SPA_FAILED;
  }

  int it = 0;
  double step = 0;
  double delta = t;
  double diff = K1(delta, params) - t;
  while (abs(diff) > tol && it < max_iter) {
    step = -(K1(delta, params) - t) / K2(delta, params);
    if (t+ step > t_max or t + step < t_min) {
      delta = 0.5*(t_max + t_min);
    } else {
      delta = t + step;
    }

    diff = K1(delta, params) - t;
    if (diff < 0) {
      t_min = delta;
    } else {
      t_max = delta;
    }
    ++it;
  }
  return delta;
}

inline int sgn(const double& d) {
  if (d < 0) {
    return -1;
  } else {
    return 1;
  }
}

double compute_spa_pval_from_sum_stats(const double& score,
                                       const double& std,
                                       const double& burden_aaf,
                                       const double& num_cases,
                                       const double& num_controls,
                                       const bool& one_sided) {
  double N = num_cases + num_controls;
  double nhom_ref = pow(1 - burden_aaf, 2)*N;
  double nhet = 2*burden_aaf*(1 - burden_aaf)*N;
  double nhom_alt = pow(burden_aaf, 2)*N;
  return compute_spa_pval_from_geno_counts(score, std, nhom_ref, nhet, nhom_alt, num_cases, num_controls, one_sided);
}

double compute_spa_pval_from_geno_counts(const double& score,
                                         const double& std,
                                         const double& nhom_ref,
                                         const double& nhet,
                                         const double& nhom_alt,
                                         const double& num_cases,
                                         const double& num_controls,
                                         const bool& one_sided) {
  double t = score / std;
  double N = nhom_ref + nhet + nhom_alt;
  double mu = num_cases / (num_cases + num_controls);
  double burden_aaf = (nhet + 2*nhom_alt) / (2*N);
  double Ghom_ref = -2*burden_aaf / std;
  double Ghet = (1 - 2*burden_aaf) / std;
  double Ghom_alt = (2 - 2*burden_aaf) / std;

  const sum_stat_spa_params params {
    N, mu, Ghom_ref, Ghet, Ghom_alt, nhom_ref, nhet, nhom_alt
  };

  double delta = solve_for_delta(t, params);
  if (abs(K1(delta, params) - t) > 1e-2) {
    return SPA_FAILED;
  }

  boost::math::normal dist(0, 1);
  double w = sgn(delta)*sqrt(2*(t*delta - K(delta, params)));
  double v = delta*sqrt(K2(delta, params));
  double z = w + (1/w)*log(v / w);
  double pval = 0;
  // compute the upper tail
  try {
    if (z > 0) {
    pval = cdf(complement(dist, z));
    } else {
      pval = cdf(dist, z);
    }
  } catch (const std::domain_error& e) {
    return SPA_FAILED;
  }

  if (!one_sided) {
    t = -t;
    delta = -delta;
    if (abs(-K1(-delta, params) - t) > 1e-2) {
      return SPA_FAILED;  
    }

    w = sgn(delta)*sqrt(2*(t*delta + K(-delta, params)));
    v = delta*sqrt(K2(delta, params));
    z = w + (1/w)*log(v / w);

    try {
      if (z > 0) {
        pval += cdf(complement(dist, z));
      } else {
        pval += cdf(dist, z);
      }
    } catch (const std::domain_error& e) {
      return SPA_FAILED;
    }
  }

  if (pval <= 0 || pval >= 1 || boost::math::isnan(pval) || boost::math::isinf(pval)) {
    return SPA_FAILED;
  } else {
    return pval;
  }
}

double compute_spa_chival_from_sum_stats(const double& score,
                                         const double& std,
                                         const double& burden_aaf,
                                         const double& num_cases,
                                         const double& num_controls,
                                         const bool& one_sided) {
  double N = num_cases + num_controls;
  double nhom_ref = pow(1 - burden_aaf, 2)*N;
  double nhet = 2*burden_aaf*(1 - burden_aaf)*N;
  double nhom_alt = pow(burden_aaf, 2)*N;
  return compute_spa_chival_from_geno_counts(score, std, nhom_ref, nhet, nhom_alt, num_cases, num_controls, one_sided);
}

double compute_spa_chival_from_geno_counts(const double& score,
                                           const double& std,
                                           const double& nhom_ref,
                                           const double& nhet,
                                           const double& nhom_alt,
                                           const double& num_cases,
                                           const double& num_controls,
                                           const bool& one_sided) {
  double pval = compute_spa_pval_from_geno_counts(score,
                                                  std,
                                                  nhom_ref,
                                                  nhet,
                                                  nhom_alt,
                                                  num_cases,
                                                  num_controls,
                                                  one_sided);
  if (pval == SPA_FAILED) {
    return pval;
  }

  try {
    return boost::math::quantile(boost::math::complement(boost::math::chi_squared(1), pval));
  } catch (const std::overflow_error& e) {
    return SPA_FAILED;
  }
}

inline double Kvec(double t, const vector<sum_stat_spa_params>& params) {
  double r = 0;
  for (size_t i = 0; i < params.size(); ++i) {
    if (params[i].N == 0) continue;
    r += K(t, params[i]);
  }
  return r;
}

inline double K1vec(double t, const vector<sum_stat_spa_params>& params) {
  double r = 0;
  for (size_t i = 0; i < params.size(); ++i) {
    if (params[i].N == 0) continue;
    r += K1(t, params[i]);
  }
  return r;
}

inline double K2vec(double t, const vector<sum_stat_spa_params>& params) {
  double r = 0;
  for (size_t i = 0; i < params.size(); ++i) {
    if (params[i].N == 0) continue;
    r += K2(t, params[i]);
  }
  return r;
}

inline double solve_for_delta_vec(double t, const vector<sum_stat_spa_params>& params) {
  int max_iter = 100;
  double tol = 1e-6;
  double t_min = 0;
  double t_max = 0;
  for (size_t i = 0; i < params.size(); ++i) {
    if (params[i].N == 0) continue;
    if (params[i].Ghom_ref > 0) {
      t_min -= params[i].nhom_ref * params[i].Ghom_ref;
      t_max += params[i].nhom_ref * params[i].Ghom_ref;
    }
    if (params[i].Ghet > 0) {
      t_min -= params[i].nhet * params[i].Ghet;
      t_max += params[i].nhet * params[i].Ghet;
    }
    if (params[i].Ghom_alt > 0) {
      t_min -= params[i].nhom_alt * params[i].Ghom_alt;
      t_max += params[i].nhom_alt * params[i].Ghom_alt;
    }
  }

  if (t < t_min || t > t_max) {
    return SPA_FAILED;
  }

  int it = 0;
  double step = 0;
  double delta = t;
  double diff = K1vec(delta, params) - t;
  while (abs(diff) > tol && it < max_iter) {
    step = -(K1vec(delta, params) - t) / K2vec(delta, params);
    if (t+ step > t_max or t + step < t_min) {
      delta = 0.5*(t_max + t_min);
    } else {
      delta = t + step;
    }

    diff = K1vec(delta, params) - t;
    if (diff < 0) {
      t_min = delta;
    } else {
      t_max = delta;
    }
    ++it;
  }
  return delta;
}

double compute_spa_pval_from_gc_per_study(const Eigen::VectorXd& scores,
                                          const Eigen::VectorXd& stds,
                                          const Eigen::VectorXd& nhom_ref,
                                          const Eigen::VectorXd& nhet,
                                          const Eigen::VectorXd& nhom_alt,
                                          const Eigen::VectorXd& ncases,
                                          const Eigen::VectorXd& ncontrols,
                                          bool one_sided) {
  double t = 0;
  vector<sum_stat_spa_params> params;
  for (ssize_t i = 0; i < scores.rows(); ++i) {
    double aaf = (nhet[i] + 2*nhom_alt[i]) / (2*(ncases[i] + ncontrols[i]));
    params.push_back({
      ncases[i] + ncontrols[i],
      ncases[i] / (ncases[i] + ncontrols[i]),
      -2*aaf / stds[i],
      (1 - 2*aaf) / stds[i],
      (2 - 2*aaf) / stds[i],
      nhom_ref[i],
      nhet[i],
      nhom_alt[i]
    });
    if (params[i].N != 0) {
      t += scores[i] / stds[i];
    }
  }

  double delta = solve_for_delta_vec(t, params);
  if (abs(K1vec(delta, params) - t) > 1e-2) {
    return SPA_FAILED;
  }

  boost::math::normal dist(0, 1);
  double w = sgn(delta)*sqrt(2*(t*delta - Kvec(delta, params)));
  double v = delta*sqrt(K2vec(delta, params));
  double z = w + (1/w)*log(v / w);
  double pval = 0;
  // compute the upper tail
  try {
    if (z > 0) {
    pval = cdf(complement(dist, z));
    } else {
      pval = cdf(dist, z);
    }
  } catch (const std::domain_error& e) {
    return SPA_FAILED;
  }

  if (!one_sided) {
    t = -t;
    delta = -delta;
    if (abs(-K1vec(-delta, params) - t) > 1e-2) {
      return SPA_FAILED;  
    }

    w = sgn(delta)*sqrt(2*(t*delta + Kvec(-delta, params)));
    v = delta*sqrt(K2vec(delta, params));
    z = w + (1/w)*log(v / w);

    try {
      if (z > 0) {
        pval += cdf(complement(dist, z));
      } else {
        pval += cdf(dist, z);
      }
    } catch (const std::domain_error& e) {
      return SPA_FAILED;
    }
  }

  if (pval <= 0 || pval >= 1 || boost::math::isnan(pval) || boost::math::isinf(pval)) {
    return SPA_FAILED;
  } else {
    return pval;
  }
}

double compute_spa_chival_from_gc_per_study(const Eigen::VectorXd& scores,
                                            const Eigen::VectorXd& stds,
                                            const Eigen::VectorXd& nhom_ref,
                                            const Eigen::VectorXd& nhet,
                                            const Eigen::VectorXd& nhom_alt,
                                            const Eigen::VectorXd& ncases,
                                            const Eigen::VectorXd& ncontrols,
                                            bool one_sided) {
  double pval = compute_spa_pval_from_gc_per_study(scores,
                                                   stds,
                                                   nhom_ref,
                                                   nhet,
                                                   nhom_alt,
                                                   ncases,
                                                   ncontrols,
                                                   one_sided);
  if (pval == SPA_FAILED) {
    return pval;
  }

  try {
    return boost::math::quantile(boost::math::complement(boost::math::chi_squared(1), pval));
  } catch (const std::overflow_error& e) {
    return SPA_FAILED;
  }
}

}