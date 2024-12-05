#ifndef MISC_INCLUDED_H
#define MISC_INCLUDED_H

#include "../lapack_complex.hpp"

#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace stat::misc {

  std::pair<double, double> wst_get_effect_size_unconditional(const Eigen::VectorXd& betas, const Eigen::VectorXd& ses, const Eigen::Ref<const Eigen::MatrixXd> cov);
  std::pair<double, double> wst_get_effect_size_unconditional(const Eigen::VectorXd& betas, const Eigen::VectorXd& ses, const Eigen::Ref<const Eigen::MatrixXd> full_cov, const std::vector<int>& mask_indices);
  std::pair<double, double> meta_analyze_effects(const Eigen::VectorXd& betas, const Eigen::VectorXd& ses);
  double estimate_burden_aaf(const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> ld);
  double estimate_burden_aaf(const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> full_ld, const std::vector<int>& mask_indices);
  double estimate_nhets(const Eigen::VectorXd& variant_nhets, const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> ld);
  double estimate_nhets(const Eigen::VectorXd& variant_nhets, const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> full_ld, const std::vector<int>& mask_indices);
  double estimate_nalts(const Eigen::VectorXd& variant_nalts, const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> ld);
  double estimate_nalts(const Eigen::VectorXd& variant_nalts, const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> full_ld, const std::vector<int>& mask_indices);
  double get_se_from_beta_pval(const double& beta, const double& pval);

  void cov_to_corr(Eigen::VectorXd& variances, Eigen::MatrixXd& corr, const Eigen::MatrixXd& cov);
  void rescale_cov(Eigen::MatrixXd& cov, const Eigen::VectorXd& by_diag);
  void rescale_cov_block(Eigen::MatrixXd& G,
                         Eigen::MatrixXd& C,
                         Eigen::MatrixXd& G_C,
                         const Eigen::VectorXd& by_G_diag,
                         const Eigen::VectorXd& by_C_diag);

};

#endif
