#ifndef STAT_TESTS_DEFINED_H
#define STAT_TESTS_DEFINED_H

#include "../lapack_complex.hpp"

#include <vector>

#include <boost/math/distributions.hpp>
#include <Eigen/Dense>

#include "qfc.h"


namespace stat::tests {
  struct test_result_t {
    double pval;
    double log10p;

    bool operator==(const test_result_t& other) {
      return this->pval == other.pval && this->log10p == other.log10p;
    }

    bool operator!=(const test_result_t& other) {
      return !this->operator==(other);
    }
  };

  const double MIN_Z = boost::math::quantile(boost::math::normal(), std::numeric_limits<double>::min());

  const test_result_t TEST_FAILED {-1, -1};
  
  // Computes log10 of the CDF (lower tail) of a standard normal distribution
  double log10p_normal(const double& stat, const bool& two_sided);

  test_result_t chisq_df1(const double& stat);
  test_result_t chisq_dfk(const double& stat, const double& df);
  test_result_t normal(const double& stat, const bool& two_sided);

  test_result_t skato(const double& Qskat, const double& Qburden, double rho, const Eigen::MatrixXd& LD);
  test_result_t wst_burden(const double& Qburden, const Eigen::MatrixXd& LD);
  test_result_t wst_burden(const double& Qburden, const double& ld_sum);
  test_result_t acat(const std::vector<double>& log10_pvals, double pmax_thr = 0.999);
  test_result_t weighted_acat(const std::vector<double>& log10_pvals, const std::vector<double>& weights, double pmax_thr = 0.999);

  test_result_t fishers(const std::vector<double>& log10_pvals);
  test_result_t stouffers(const std::vector<double>& log10_pvals, const std::vector<double>& weights, double adj_var = 0);
  test_result_t stouffers_two_sided(const std::vector<double>& log10_pvals, const std::vector<double>& weights, const std::vector<double>& signs, double adj_var = 0);

  test_result_t mixture_of_chisq(const double& stat, const Eigen::VectorXd& weights);
  test_result_t mixture_of_chisq_davies(const double& stat, Eigen::VectorXd& weights, const bool& force_stringent);
  test_result_t mixture_of_chisq_kuonen(const double& stat, const Eigen::VectorXd& weights);
  test_result_t mixture_of_chisq_liu(const double& stat, const Eigen::VectorXd& weights);
}

#endif