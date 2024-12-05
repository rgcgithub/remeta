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

inline double solve_for_delta(const double& t, const sum_stat_spa_params& params) {
  std::uintmax_t max_iter = 100;
  return boost::math::tools::newton_raphson_iterate(
    [&t, &params](const double& d) {
      return std::make_pair(K1(d, params) - t, K2(d, params));
    },
    t,
    -100.0,
    100.0,
    10,
    max_iter
  );
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
    mu, Ghom_ref, Ghet, Ghom_alt, nhom_ref, nhet, nhom_alt
  };

  double delta = solve_for_delta(t, params);
  if (abs(K1(delta, params) - t) > 1e-2 || 100 - abs(delta) < 1) {
    return SPA_FAILED;
  }

  boost::math::normal dist(0, 1);
  double w = sgn(delta)*sqrt(2*(t*delta - K(delta, params)));
  double v = delta*sqrt(K2(delta, params));
  double z = w + (1/w)*log(v / w);
  double pval = 0;
  // compute the upper tail
  try {
    if (t > 0) {
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
      if (t > 0) {
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

}