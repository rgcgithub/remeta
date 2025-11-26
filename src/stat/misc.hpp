#ifndef MISC_INCLUDED_H
#define MISC_INCLUDED_H

#include "../lapack_complex.hpp"

#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace stat::misc {

  std::pair<double, double> wst_get_effect_size_unconditional(const Eigen::VectorXd& betas, const Eigen::VectorXd& ses, const Eigen::Ref<const Eigen::MatrixXd> cov);
  std::pair<double, double> wst_get_effect_size_unconditional(const Eigen::VectorXd& betas, const Eigen::VectorXd& ses, const Eigen::Ref<const Eigen::MatrixXd> full_cov, const std::vector<int>& mask_indices);
  std::pair<double, double> meta_analyze_effects(const Eigen::VectorXd& betas, const Eigen::VectorXd& ses);
  double estimate_burden_aaf(const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> ld);
  double estimate_burden_aaf(const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> full_ld, const std::vector<int>& mask_indices);
  double estimate_nhets(const Eigen::VectorXd& variant_nhets, const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> ld);
  double estimate_nhets(const Eigen::VectorXd& variant_nhets, const Eigen::VectorXd& variant_aafs, const Eigen::SparseMatrix<double>& ld_sp);
  double estimate_nhets(const Eigen::VectorXd& variant_nhets, const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> full_ld, const std::vector<int>& mask_indices);
  double estimate_nalts(const Eigen::VectorXd& variant_nalts, const Eigen::VectorXd& variant_aafs, const Eigen::SparseMatrix<double>& ld_sp);
  double estimate_nalts(const Eigen::VectorXd& variant_nalts, const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> ld);
  double estimate_nalts(const Eigen::VectorXd& variant_nalts, const Eigen::VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> full_ld, const std::vector<int>& mask_indices);
  double get_se_from_beta_log10p(double beta, double log10pval);
  double get_chisq_stat_from_log10p(double log10p);

  void cov_to_corr(Eigen::VectorXd& variances, Eigen::MatrixXd& corr, const Eigen::MatrixXd& cov);
  void rescale_cov(Eigen::MatrixXd& cov, const Eigen::VectorXd& by_diag);
  void rescale_cov(Eigen::SparseMatrix<double>& cov, const Eigen::VectorXd& cov_diag, const Eigen::VectorXd& by_diag);
  void rescale_cov_block(Eigen::MatrixXd& G,
                         Eigen::MatrixXd& C,
                         Eigen::MatrixXd& G_C,
                         const Eigen::VectorXd& by_G_diag,
                         const Eigen::VectorXd& by_C_diag);
  void rescale_cov_block(Eigen::SparseMatrix<double>& Gsp,
                         Eigen::MatrixXd& C,
                         Eigen::MatrixXd& G_C,
                         const Eigen::VectorXd& G_diag,
                         const Eigen::VectorXd& by_G_diag,
                         const Eigen::VectorXd& by_C_diag);

  void collapse_anno(Eigen::SparseMatrix<double>& Gsp,
                     Eigen::Ref<Eigen::VectorXd> scores,
                     Eigen::Ref<Eigen::VectorXd> scorev,
                     Eigen::Ref<Eigen::VectorXd> ses,
                     const Eigen::VectorXd& anno_indicators);

};

#endif
