#include "misc.hpp"

#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <utility>
#include <vector>
using std::max;
using std::min;
using std::pair;
using std::unordered_set;
using std::vector;
using std::cout;
using std::endl;

#include <boost/math/distributions.hpp>

#include <Eigen/Dense>
using Eigen::DiagonalMatrix;
using Eigen::MatrixXf;
using Eigen::MatrixXd;
using Eigen::Ref;
using Eigen::VectorXf;
using Eigen::VectorXd;
#include <Eigen/Sparse>
using Eigen::SparseMatrix;
using Eigen::Triplet;


namespace stat::misc {

pair<double, double> wst_get_effect_size_unconditional(const Eigen::VectorXd& betas,
                                                       const Eigen::VectorXd& ses,
                                                       const Eigen::Ref<const Eigen::MatrixXd> cov) {
  vector<double> weights;
  vector<ssize_t> indices;
  double beta = 0;
  double weight_sum = 0;
  double w = 0;
  for (ssize_t i = 0; i < betas.size(); ++i) {
    if (betas(i) != 0 && ses(i) != 0 && cov(i, i) != 0) {
      indices.push_back(i);
      w = 1/(ses(i)*ses(i));
      weights.push_back(w);
      weight_sum += w;
      beta += betas(i) * w;
    }
  }
  beta /= weight_sum;

  for (size_t i = 0; i < weights.size(); ++i) {
    weights[i] /= weight_sum;
  }

  double var = 0;
  double s;
  for (size_t i = 0; i < indices.size(); ++i) {
    for (size_t j = 0; j < indices.size(); ++j) {
      if ( (cov(indices[i], indices[i]) == 0) || (cov(indices[j], indices[j]) == 0) ) {
        continue;
      }
      s = ses(indices[i]) * ses(indices[j]) / sqrt(cov(indices[i], indices[i])*cov(indices[j], indices[j]));
      var += weights[i] * weights[j] * cov(indices[i], indices[j]) * s;
    }
  }

  pair<double, double> result;
  result.first = beta;
  result.second = sqrt(var);
  return result;
}

pair<double, double> wst_get_effect_size_unconditional(const Eigen::VectorXd& betas,
                                                       const Eigen::VectorXd& ses,
                                                       const Eigen::Ref<const Eigen::MatrixXd> full_cov,
                                                       const vector<int>& mask_indices) {
  vector<double> weights;
  vector<ssize_t> indices;
  double beta = 0;
  double weight_sum = 0;
  double w = 0;
  int i0, j0;
  for (ssize_t i = 0; i < betas.size(); ++i) {
    i0 = mask_indices[i];
    if (betas(i) != 0 && ses(i) != 0 && full_cov(i0, i0) != 0) {
      indices.push_back(i);
      w = 1/(ses(i)*ses(i));
      weights.push_back(w);
      weight_sum += w;
      beta += betas(i) * w;
    }
  }
  beta /= weight_sum;

  for (size_t i = 0; i < weights.size(); ++i) {
    weights[i] /= weight_sum;
  }

  double var = 0;
  double s;
  for (size_t i = 0; i < indices.size(); ++i) {
    i0 = mask_indices[indices[i]];
    for (size_t j = 0; j < indices.size(); ++j) {
      j0 = mask_indices[indices[j]];
      if ( (full_cov(i0, i0) == 0) || (full_cov(j0, j0) == 0) ) {
        continue;
      }
      s = ses(indices[i]) * ses(indices[j]) / sqrt(full_cov(i0, i0)*full_cov(j0, j0));
      var += weights[i] * weights[j] * full_cov(i0, j0) * s;
    }
  }

  pair<double, double> result;
  result.first = beta;
  result.second = sqrt(var);
  return result;
}


double get_se_from_beta_log10p(double beta, double log10p) {
  double chival = get_chisq_stat_from_log10p(log10p);
  double z = sqrt(chival);
  return abs(beta / z);
}


double get_chisq_stat_from_log10p(double log10p) {
  boost::math::chi_squared chisq1(1);
  if (log10p < log10(std::numeric_limits<double>::min())) {
    double val = -log10p * log(100) + log(2/M_PI);
    return val - log(val);
  } else {
    return boost::math::quantile(
      boost::math::complement(
        boost::math::chi_squared(1),
        pow(10, log10p)    
      )
    );
  }
}

// Let: 
//   * G=X_1 + X_2 be a genotype, where X_i \in {0, 1}
//   * cov be the covariance between the genotypes of two variants G_1 and G_2
//   * f1 and f2 be the allele frequency of variants 1 and 2
//
// Computes E[X_1 X_2] from cov, f1, and g2
double cov_to_prob(const double& cov, const double& f1, const double& f2) {
  return max(0.5*cov + f1*f2, 0.0);
}

// Let X=(X_1, ..., X_p) be the (haploid) genotypes of a variants in a max.
// The burden genotype is Y=max(X_1, ..., X_p) \in {0, 1}
// We want to compute E[Y] = P(Y = 1) = 1 - P(Y = 0)
// If we model the X_{1..p} as a Markov chain, we can estimate:
// P(Y = 0) = P(X_1 = 0) \prod_{t=2}^T P(X_t = 0 | X_{t-1}=0)
// P(X_t = 0) = 1 -f_t
// P(X_t = 0 | X_{t-1}=0) can be computed from the genotype covariance matrix
double estimate_burden_aaf(const VectorXd& variant_aafs, const Eigen::Ref<const MatrixXd> cov) {
  unordered_set<int> variants_remaining;
  int cur_idx = 0;
  double cur_aaf = variant_aafs(0);
  for (int i = 0; i < variant_aafs.size(); ++i) {
    variants_remaining.insert(i);
    if (variant_aafs(i) > cur_aaf) {
      cur_idx = i;
      cur_aaf = variant_aafs(i);
    }
  }
  double max_aaf = cur_aaf;

  double logp = log(1 - variant_aafs(cur_idx));
  variants_remaining.erase(cur_idx);
  while (variants_remaining.size() > 0) {
    int next_idx = 0;
    double next_aaf = 0;
    double next_prob = 0;
    double prob_x1x2;
    for ( const int& i : variants_remaining ) {
      prob_x1x2 = cov_to_prob(cov(cur_idx, i), variant_aafs(cur_idx), variant_aafs(i));
      if (prob_x1x2 >= next_prob) {
        next_idx = i;
        next_aaf = variant_aafs(i);
        next_prob = prob_x1x2;
      }
    }
    logp += log(1 - cur_aaf - next_aaf + next_prob) - log(1 - cur_aaf);

    variants_remaining.erase(next_idx);
    cur_idx = next_idx;
    cur_aaf = next_aaf;
    next_prob = 0;
  }

  double aaf = 1 - exp(logp);
  return max(aaf, max_aaf);
}

double estimate_burden_aaf(const VectorXd& variant_aafs, const Eigen::Ref<const MatrixXd> full_cov, const vector<int>& mask_indices) {
  unordered_set<int> variants_remaining;
  int cur_idx = 0;
  double cur_aaf = variant_aafs(0);
  for (int i = 0; i < variant_aafs.size(); ++i) {
    variants_remaining.insert(i);
    if (variant_aafs(i) > cur_aaf) {
      cur_idx = i;
      cur_aaf = variant_aafs(i);
    }
  }
  double max_aaf = cur_aaf;

  double logp = log(1 - variant_aafs(cur_idx));
  variants_remaining.erase(cur_idx);
  int cur_idx0, i0;
  while (variants_remaining.size() > 0) {
    cur_idx0 =  mask_indices[cur_idx];
    int next_idx = 0;
    double next_aaf = 0;
    double next_prob = 0;
    double prob_x1x2;
    for ( const int& i : variants_remaining ) {
      i0 = mask_indices[i];
      prob_x1x2 = cov_to_prob(full_cov(cur_idx0, i0), variant_aafs(cur_idx), variant_aafs(i));
      if (prob_x1x2 >= next_prob) {
        next_idx = i;
        next_aaf = variant_aafs(i);
        next_prob = prob_x1x2;
      }
    }
    logp += log(1 - cur_aaf - next_aaf + next_prob) - log(1 - cur_aaf);

    variants_remaining.erase(next_idx);
    cur_idx = next_idx;
    cur_aaf = next_aaf;
    next_prob = 0; 
  }

  double aaf = 1 - exp(logp);
  return max(aaf, max_aaf);
}

double estimate_nhets(const VectorXd& variant_nhets, const VectorXd& variant_aafs, const Eigen::Ref<const MatrixXd> ld) {
  double nhets = variant_nhets(0);
  double max_hets = nhets;
  double f = 0;
  double c = 0;
  bool skip = false;
  for (ssize_t i = 1; i < variant_nhets.size(); ++i) {
    if (ld(i, i) == 0 || variant_nhets(i) == 0) continue;
    c = 0;
    skip = false;
    for (ssize_t j = 0; j < i; ++j) {
      // if there is a low mac in perfect LD, then skip variant i to not overcount
      if (variant_nhets(j) > 0
          && variant_nhets(j) < 20
          && variant_aafs(i) == variant_aafs(j)
          && ld(i, j) == ld(i, i) && ld(i, i) == ld(j, j)) {
        skip = true;
        break;
      }
      if (ld(i, j) <= 0 || variant_nhets(j) == 0) {
        f = 0;
      } else {
        f = (cov_to_prob(ld(i, j), variant_aafs(i), variant_aafs(j)) / variant_aafs(i)) * (1 - variant_aafs(i))
          + (1 - cov_to_prob(ld(i, j), variant_aafs(i), variant_aafs(j)) / variant_aafs(i)) * variant_aafs(i);
      }
      f = max(0., min(1., f));
      c += log(1 - f);
    }

    if (!skip) {
      nhets += variant_nhets(i) * exp(c);
      max_hets = max(max_hets, variant_nhets(i));
    }
  }
  return max(nhets, max_hets);
}

double estimate_nhets(const VectorXd& variant_nhets, const VectorXd& variant_aafs, const SparseMatrix<double>& ld_sp) {
  double nhets = variant_nhets(0);
  double max_hets = nhets;
  double f = 0;
  double c = 0;
  bool skip = false;
  int j = 0;
  VectorXd ld_diag = ld_sp.diagonal();
  for (int i = 1; i < ld_sp.outerSize(); ++i) {
    if (ld_diag(i) == 0 || variant_nhets(i) == 0) continue;
    c = 0;
    skip = false;
    for (SparseMatrix<double>::InnerIterator it(ld_sp, i); it.index() < i; ++it) {
      j = it.index();
      // if there is a low mac in perfect LD, then skip variant i to not overcount
      if (variant_nhets(j) > 0
          && variant_nhets(j) < 20
          && variant_aafs(i) == variant_aafs(j)
          && it.value() == ld_diag(i) && ld_diag(i) == ld_diag(j)) {
        skip = true;
        break;
      }
      if (it.value() <= 0 || variant_nhets(j) == 0) {
        f = 0;
      } else {
        f = (cov_to_prob(it.value(), variant_aafs(i), variant_aafs(j)) / variant_aafs(i)) * (1 - variant_aafs(i))
          + (1 - cov_to_prob(it.value(), variant_aafs(i), variant_aafs(j)) / variant_aafs(i)) * variant_aafs(i);
      }
      f = max(0., min(1., f));
      c += log(1 - f);
    }

    if (!skip) {
      nhets += variant_nhets(i) * exp(c);
      max_hets = max(max_hets, variant_nhets(i));
    }
  }
  return max(nhets, max_hets);
}

double estimate_nhets(const VectorXd& variant_nhets, const VectorXd& variant_aafs, const Eigen::Ref<const MatrixXd> full_ld, const vector<int>& mask_indices) {
  double nhets = variant_nhets(0);
  double max_hets = nhets;
  double f = 0;
  double c = 0;
  bool skip = false;
  int i0, j0;
  for (ssize_t i = 1; i < variant_nhets.size(); ++i) {
    i0 = mask_indices[i];
    if (full_ld(i0, i0) == 0 || variant_nhets(i) == 0) continue;
    c = 0;
    skip = false;
    for (ssize_t j = 0; j < i; ++j) {
      j0 = mask_indices[j];
      // if there is a low mac in perfect LD, then skip variant i to not overcount
      if (variant_nhets(j) > 0
          && variant_nhets(j) < 20
          && variant_aafs(i) == variant_aafs(j)
          && full_ld(i0, j0) == full_ld(i0, i0) && full_ld(i0, i0) == full_ld(j0, j0)) {
        skip = true;
        break;
      }
      if (full_ld(i0, j0) <= 0 || variant_nhets(j) == 0) {
        f = 0;
      } else {
        f = (cov_to_prob(full_ld(i0, j0), variant_aafs(i), variant_aafs(j)) / variant_aafs(i)) * (1 - variant_aafs(i))
          + (1 - cov_to_prob(full_ld(i0, j0), variant_aafs(i), variant_aafs(j)) / variant_aafs(i)) * variant_aafs(i);
      }
      f = max(0., min(1., f));
      c += log(1 - f);
    }

    if (!skip) {
      nhets += variant_nhets(i) * exp(c);
      max_hets = max(max_hets, variant_nhets(i));
    }
  }
  return max(nhets, max_hets);
}

double estimate_nalts(const VectorXd& variant_nalts, const VectorXd& variant_aafs, const Eigen::Ref<const Eigen::MatrixXd> ld) {
  double nalts = variant_nalts(0);
  double max_alts = nalts;
  double f = 0;
  double c = 0;
  for (ssize_t i = 1; i < variant_nalts.size(); ++i) {
    if (ld(i, i) == 0 || variant_nalts(i) == 0) continue;
    c = 0;
    for (ssize_t j = 0; j < i; ++j) {
      if (ld(i, j) <= 0 || variant_nalts(j) == 0) {
        f = 0;
      } else {
        f = cov_to_prob(ld(i, j), variant_aafs(i), variant_aafs(j)) / variant_aafs(i);
      }
      f = max(0., min(1., f));
      c += log(1 - f*f);
    }
    nalts += variant_nalts(i) * exp(c);
    max_alts = max(max_alts, variant_nalts(i));
  }
  return max(nalts, max_alts);
}

double estimate_nalts(const VectorXd& variant_nalts, const VectorXd& variant_aafs, const SparseMatrix<double>& ld_sp) {
  double nalts = variant_nalts(0);
  double max_alts = nalts;
  double f = 0;
  double c = 0;
  int j = -1;
  VectorXd ld_diag = ld_sp.diagonal();
  for (int i = 1; i < ld_sp.outerSize(); ++i) {
    if (ld_diag(i) == 0 || variant_nalts(i) == 0) continue;
    c = 0;
    for (SparseMatrix<double>::InnerIterator it(ld_sp, i); it.index() < i; ++it) {
      j = it.index();
      if (it.value() <= 0 || variant_nalts(j) == 0) {
        f = 0;
      } else {
        f = cov_to_prob(it.value(), variant_aafs(i), variant_aafs(j)) / variant_aafs(i);
      }
      f = max(0., min(1., f));
      c += log(1 - f*f);
    }
    nalts += variant_nalts(i) * exp(c);
    max_alts = max(max_alts, variant_nalts(i));
  }
  return max(nalts, max_alts);
}

double estimate_nalts(const VectorXd& variant_nalts, const VectorXd& variant_aafs, const Eigen::Ref<const MatrixXd> full_ld, const vector<int>& mask_indices) {
  double nalts = variant_nalts(0);
  double max_alts = nalts;
  double f = 0;
  double c = 0;
  int i0, j0;
  for (ssize_t i = 1; i < variant_nalts.size(); ++i) {
    i0 = mask_indices[i];
    if (full_ld(i0, i0) == 0 || variant_nalts(i) == 0) continue;
    c = 0;
    for (ssize_t j = 0; j < i; ++j) {
      j0 = mask_indices[j];
      if (full_ld(i0, j0) <= 0 || variant_nalts(j) == 0) {
        f = 0;
      } else {
        f = cov_to_prob(full_ld(i0, j0), variant_aafs(i), variant_aafs(j)) / variant_aafs(i);
      }
      f = max(0., min(1., f));
      c += log(1 - f*f);
    }
    nalts += variant_nalts(i) * exp(c);
    max_alts = max(max_alts, variant_nalts(i));
  }
  return max(nalts, max_alts);
}

pair<double, double> meta_analyze_effects(const VectorXd& betas, const VectorXd& ses) {
  vector<double> weights;
  double weight_sum = 0;
  double beta = 0;
  double weight = 0;
  double sse = 0;
  for (int i = 0; i < betas.size(); ++i) {
    if (ses[i] != 0) {
      sse = ses[i];
      weight = 1 / pow(sse, 2);
      weights.push_back(weight);
      weight_sum += weight;
      beta += weight * betas[i]; 
    }
  }

  beta /= weight_sum;
  double se = sqrt(1 / weight_sum);
  //double z = beta / se;

  //double pval = 2 * cdf(this->s, -abs(z));
  pair<double, double> result;
  result.first = beta;
  result.second = se;
  
  return result;
}

void cov_to_corr(VectorXd& variances, MatrixXd& corr, const MatrixXd& cov) {
  variances = cov.diagonal();
  MatrixXd tmp = (variances.array() > 0).select(
                    (variances.array() > 0).select(variances, 1)
                                          .array()
                                          .sqrt()
                                          .inverse()
                                          .matrix(),
                    0).asDiagonal();
  corr = tmp * cov * tmp;
}

void rescale_cov(MatrixXd& cov, const VectorXd& by_diag) {
  DiagonalMatrix<double, Eigen::Dynamic> D(by_diag.size());
  for (ssize_t i = 0; i < by_diag.size(); ++i) {
    if (cov(i, i) == 0) {
      D.diagonal()[i] = 0;
    } else {
      D.diagonal()[i] = sqrt(abs(by_diag(i) / cov(i, i)));
    }
  }
  cov = D * cov * D;
}

void rescale_cov(SparseMatrix<double>& cov, const Eigen::VectorXd& cov_diag, const Eigen::VectorXd& by_diag) {
  DiagonalMatrix<double, Eigen::Dynamic> D(by_diag.size());
  for (ssize_t i = 0; i < by_diag.size(); ++i) {
    if (cov_diag(i) == 0) {
      D.diagonal()[i] = 0;
    } else {
      D.diagonal()[i] = sqrt(abs(by_diag(i) / cov_diag(i)));
    }
  }
  cov = D * cov * D;
}

void rescale_cov_block(MatrixXd& G,
                       MatrixXd& C,
                       MatrixXd& G_C,
                       const VectorXd& by_G_diag,
                       const VectorXd& by_C_diag) {
  VectorXd G_diag = G.diagonal();
  for (ssize_t i = 0; i < G_diag.size(); ++i) {
    if (G(i, i) == 0) {
      G_diag[i] = 1;
    } else {
      G_diag[i] = G(i, i);
    }
  }
  VectorXd C_diag = C.diagonal();
  for (ssize_t i = 0; i < C_diag.size(); ++i) {
    if (C(i, i) == 0) {
      C_diag[i] = 1;
    } else {
      C_diag[i] = C(i, i);
    }
  }

  MatrixXd D1 = (by_G_diag.array() / G_diag.array())
                .sqrt()
                .matrix()
                .asDiagonal();
  MatrixXd D2 = (by_C_diag.array() / C_diag.array())
                .sqrt()
                .matrix()
                .asDiagonal(); 

  G = D1 * G * D1;
  C = D2 * C * D2;
  G_C = D1 * G_C * D2;
}

void rescale_cov_block(SparseMatrix<double>& Gsp,
                       MatrixXd& C,
                       MatrixXd& G_C,
                       const VectorXd& G_diag,
                       const VectorXd& by_G_diag,
                       const VectorXd& by_C_diag) {
  VectorXd G_diag_nonzero = G_diag;
  for (ssize_t i = 0; i < G_diag.size(); ++i) {
    if (G_diag[i] == 0) {
      G_diag_nonzero[i] = 1;
    } else {
      G_diag_nonzero[i] = G_diag[i];
    }
  }
  VectorXd C_diag = C.diagonal();
  for (ssize_t i = 0; i < C_diag.size(); ++i) {
    if (C(i, i) == 0) {
      C_diag[i] = 1;
    } else {
      C_diag[i] = C(i, i);
    }
  }

  SparseMatrix<double> D1 = SparseMatrix<double>(
                              (by_G_diag.array() / G_diag_nonzero.array())
                              .sqrt()
                              .matrix()
                              .asDiagonal()
                            );
  MatrixXd D2 = (by_C_diag.array() / G_diag_nonzero.array())
                .sqrt()
                .matrix()
                .asDiagonal(); 

  SparseMatrix<double> Gtmp = D1 * Gsp * D1;
  std::swap(Gsp, Gtmp);
  C = D2 * C * D2;
  G_C = D1 * G_C * D2;
}

void collapse_anno(SparseMatrix<double>& Gsp,
                   Eigen::Ref<VectorXd> scores,
                   Eigen::Ref<VectorXd> scorev,
                   Eigen::Ref<VectorXd> ses,
                   const VectorXd& anno_indicators) {
  ssize_t anno_idx = 0;
  for (ssize_t v = 0; v < scores.size(); ++v) {
    if (anno_indicators(v) == 1) {
      anno_idx = v;
      break;
    }
  }

  vector<Triplet<double> > anno_collapsor_triplets;
  for (ssize_t v = 0; v < scores.size(); ++v) {
    if (anno_indicators(v) == 1) {
      anno_collapsor_triplets.push_back(Triplet<double>(anno_idx, v, 1));
    } else {
      anno_collapsor_triplets.push_back(Triplet<double>(v, v, 1));
    }
  }
  SparseMatrix<double> anno_collapsor(Gsp.rows(), Gsp.cols());
  anno_collapsor.setFromTriplets(anno_collapsor_triplets.begin(), anno_collapsor_triplets.end());

  double burden_score = (scores.transpose() * anno_indicators).value();
  double burden_var = (anno_indicators.transpose() * Gsp * anno_indicators).value();
  double z = burden_score / sqrt(burden_var);
  if (burden_score == 0 || burden_var == 0) {
    return ;
  }

  VectorXd G_diag = Gsp.diagonal();
  SparseMatrix<double> beta_cov;
  VectorXd D = (G_diag.array() > 0).select(
    ses.array() / G_diag.array().sqrt(),
    VectorXd::Zero(ses.rows())
  ).matrix();
  beta_cov = D.asDiagonal() * Gsp * D.asDiagonal();
  VectorXd beta_weights = (beta_cov.diagonal().array() > 0).select(
    beta_cov.diagonal().array().inverse(),
    0
  );

  double weight_sum = (anno_indicators.transpose() * beta_weights).value();
  double beta_var = (
    (anno_indicators.array()*beta_weights.array()).matrix().transpose()/weight_sum 
    * beta_cov
    * (anno_indicators.array()*beta_weights.array()).matrix()/weight_sum
  ).value();
  double anno_se = sqrt( beta_var );

  Gsp = anno_collapsor * Gsp * anno_collapsor.transpose();
  G_diag = Gsp.diagonal();
  for (ssize_t i = 0; i < scores.size(); ++i) {
    if (i == anno_idx) {
      ses(i) = anno_se;
      scores(i) = z / anno_se;
      scorev(i) = 1 / (anno_se * anno_se);
      D(i) = sqrt(scorev(i) / G_diag(i));
    } else if (anno_indicators(i) == 1) {
      ses(i) = 0;
      scores(i) = 0;
      scorev(i) = 0;
      D(i) = 0;
    } else {
      D(i) = 1;
    }
  }

  Gsp = D.asDiagonal() * Gsp * D.asDiagonal();
}

};
