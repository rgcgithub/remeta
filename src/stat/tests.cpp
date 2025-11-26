#include "tests.hpp"

#include <algorithm>
#include <cmath>
#include <vector>
using std::max;
using std::min;
using std::vector;

#include <boost/math/distributions.hpp>
using boost::math::cauchy;
using boost::math::cdf;
using boost::math::complement;

#include <Eigen/Dense>
using Eigen::ArrayXd;
using Eigen::ArrayXi;
using Eigen::Map;
using Eigen::Ref;
using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::Array<bool,Eigen::Dynamic,1> ArrayXb;
typedef Eigen::Map<const Eigen::ArrayXd > MapcArXd;

#include <iostream>
using std::cout;
using std::endl;

#include "Faddeeva.hh"
#include "qfc.h"

// the implementations of many of these methods are originally from regenie
namespace stat::tests {

// based on https://github.com/scipy/scipy/blob/v1.13.0/scipy/special/_faddeeva.cxx#L85-L105
double log10p_normal(const double& stat, const bool& two_sided) {
  double logp = 0;
  if (stat < -1) {
    logp = log(Faddeeva::erfcx(-stat*M_SQRT1_2)/2) - stat*stat/2;
  } else {
    logp = log1p(-Faddeeva::erfc(stat*M_SQRT1_2)/2);
  }

  if (two_sided) {
    return logp / log(10) + log10(2);
  } else {
    return logp / log(10);
  }
}

// https://ritscm.regeneron.com/projects/RGC-STATML/repos/regenie/browse/src/Regenie.cpp#1698
Eigen::ArrayXi get_true_indices(const Ref<const ArrayXb>&  bool_arr){

  ArrayXi v_indices ( bool_arr.count() );
  for(int i = 0, j = 0; i < bool_arr.size(); i++)
    if(bool_arr(i)) v_indices(j++) = i;

  return v_indices;
}

// https://ritscm.regeneron.com/projects/RGC-STATML/repos/regenie/browse/src/Joint_Tests.cpp#283
// logpvals are assumed to be -log10 pvalues
// robust to low pvalues
double get_acat_robust(const Eigen::Ref<const ArrayXd>& logpvals, const Eigen::Ref<const ArrayXd>& weights, double pmax_thr) {
  // if single pval, return log10p
  int n_pv = ((weights!=0) && (logpvals >= 0)).count();
  if(n_pv == 0) return 1;
  else if(n_pv == 1) return -1*logpvals.maxCoeff();

  cauchy dc(0,1);
  double lpv_thr = 15, lpval_out;

  // split pvals by thr
  int n_A = ((weights!=0) && (logpvals >= lpv_thr)).count(); // very small pvals
  int n_B = ((weights!=0) && (logpvals >= 0) && (logpvals < lpv_thr)).count();
  double wsum = (logpvals >= 0).select(weights, 0).sum();
  double l_TA = 0, TB = 0;

  // T_A
  if(n_A > 0){ // compute on log scale to handle the very small pvalues
    ArrayXi vind = get_true_indices((weights!=0) && (logpvals >= lpv_thr));
    ArrayXd lp = logpvals( vind ), ws = weights( vind ) / wsum;
    ArrayXd zvec = lp * log(10) + ws.log() - log(M_PI);
    double zmax = zvec.maxCoeff();
    l_TA = zmax + log( (zvec - zmax).exp().sum() );
  }
  // T_B (can be negative)
  if(n_B > 0){
    ArrayXi vind = get_true_indices((weights!=0) && (logpvals >= 0) && (logpvals < lpv_thr));
    ArrayXd pv = pow(10, -logpvals(vind)).min(pmax_thr); // avoid pvalues of 1
    ArrayXd ws = weights( vind ) / wsum;
    TB = ( ws * tan( M_PI * (0.5 - pv)) ).sum();
  }

  // T_ACAT = TA + TB
  if(n_A == 0){ // avoid computing log(TB) as TB can be negative
    lpval_out = ( TB >= 8886111 ? -log(TB) - log(M_PI) : log(cdf(complement(dc, TB))) );
  } else if ((n_B == 0) || (TB == 0)){
    lpval_out = ( l_TA >= 16 ? -l_TA - log(M_PI) : log(cdf(complement(dc, exp(l_TA)))) );
  } else {
    double lsum; // get sum on log scale
    if(TB < 0){
      double l_abs_TB = log(fabs(TB));
      if(l_abs_TB < l_TA)
        lsum = l_TA + log1p(-exp(l_abs_TB - l_TA));
      else { // compute log(-Tacat)
        lsum = l_abs_TB + log1p(-exp(l_TA - l_abs_TB)); 
        lpval_out = ( lsum >= 16 ? log1p(-exp(-lsum-log(M_PI))) : log(cdf(complement(dc, -exp(lsum)))) );
        return lpval_out/log(10);
      }
    } else {
      double l_TB = log(TB);
      lsum = fmax(l_TA, l_TB) + log1p(exp(-fabs(l_TB - l_TA)));
    } 
    lpval_out = ( lsum >= 16 ? -lsum - log(M_PI) : log(cdf(complement(dc, exp(lsum) ))) );
  }

  // return log10P
  return lpval_out/log(10);
}

// https://ritscm.regeneron.com/projects/RGC-STATML/repos/regenie/browse/src/Regenie.cpp#1772
double get_chisq_stat_from_logp(const double& log10p) {
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

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

test_result_t chisq_df1(const double& stat) {
  if (boost::math::isnan(stat) || boost::math::isinf(stat)) {
    return stat::tests::TEST_FAILED;
  }

  boost::math::chi_squared chisq1(1);
  double pval = boost::math::cdf(boost::math::complement(chisq1, stat));
  if (pval == 0) {
    return test_result_t {
      std::numeric_limits<double>::min(),
      log10(2) - 0.5 * log10( 2 * M_PI * stat ) - 0.5 * stat * M_LOG10E
    };
  } else {
    return test_result_t {
      pval,
      log10(pval)
    };
  }
}

test_result_t normal(const double& stat, const bool& two_sided) {
  if (stat <= MIN_Z) {
    return test_result_t {
      std::numeric_limits<double>::min(),
      log10p_normal(stat, two_sided)
    };
  } else if (-stat <= MIN_Z) {
      return test_result_t {
        1-std::numeric_limits<double>::min(),
        log10p_normal(stat, two_sided)
      };
  } else {
    double pval = boost::math::cdf(boost::math::normal(), stat);
    test_result_t result {pval, log10(pval)};
    if (two_sided) {
      result.pval *= 2;
      result.log10p += log10(2);
    }
    return result;
  }
}

test_result_t chisq_dfk(const double& stat, const double& df) {
  boost::math::chi_squared chisqK(df);
  double pval = boost::math::cdf(boost::math::complement(chisqK, stat));
  if (pval == 0) {
    return test_result_t {
      std::numeric_limits<double>::min(),
      log10(2) - 0.5 * df * log10(2)
               - boost::math::lgamma(df * 0.5) / log(10) 
               + 0.5 * (df-2) * log10(stat)
               - 0.5 * stat * M_LOG10E
    };
  } else {
    return test_result_t {
      pval,
      log10(pval)
    };
  }
}

test_result_t skato(const double& Qskat, const double& Qburden, double rho, const Eigen::MatrixXd& LD) {
  if (rho == 1) return wst_burden(Qburden, LD.sum());
  double tol = 1e-5;
  double Qrho = (1 - rho)*Qskat + rho*Qburden;

  int m = LD.rows();
  double c1 = sqrt(1 - rho);
  double c2 = sqrt(1 - rho + m * rho), gamma1;
  Eigen::VectorXd b = LD.rowwise().sum();
  gamma1 = b.sum();
  Eigen::MatrixXd Psi = ((1-rho) * LD.array() + 
      c1 * (c2-c1)/m * (b.rowwise().replicate(m) + b.transpose().colwise().replicate(m)).array() + // last term is outer sum of b
      (c2-c1)/m * (c2-c1)/m * gamma1
      ).matrix();

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Psi, Eigen::EigenvaluesOnly);
  Psi.resize(0, 0);

  double eigvals_max = es.eigenvalues().tail(1)(0);
  int nng = (es.eigenvalues().array() >= 0).count();
  int nonzero = (es.eigenvalues().array() > 
                  ( (es.eigenvalues().array() >= 0).select(es.eigenvalues().array(), 0).sum() / nng * tol) ).count();
  Eigen::VectorXd eigvals = es.eigenvalues().tail(nonzero);

  if (eigvals.size() == 0) {
    return stat::tests::TEST_FAILED;
  }

  eigvals /= eigvals_max;
  Qrho /= eigvals_max;

  if (rho == 1 || eigvals.size() == 1) {
    return chisq_df1(Qrho);
  } else {
    return mixture_of_chisq(Qrho, eigvals);
  }
}

test_result_t wst_burden(const double& Qburden, const Eigen::MatrixXd& LD) {
  double c = LD.sum();
  if (c > 0) {
    return chisq_df1(Qburden / c);
  } else {
    return stat::tests::TEST_FAILED;
  }
}

test_result_t wst_burden(const double& Qburden, const double& ld_sum) {
  if (ld_sum > 0) {
    return chisq_df1(Qburden / ld_sum);
  } else {
    return stat::tests::TEST_FAILED;
  }
}

test_result_t acat(const vector<double>& log10_pvals, double pmax_thr) {
  ArrayXd wts = ArrayXd::Constant(log10_pvals.size(), 1); // uniform weights
  double log10p = get_acat_robust(
    -1*Map<const ArrayXd>(log10_pvals.data(), (int)log10_pvals.size()), 
    wts,
    pmax_thr
  );
  if (log10p == 1) {
    return TEST_FAILED;
  } else {
    return test_result_t {
      max(pow(10, log10p), std::numeric_limits<double>::min()),
      log10p
    };
  }
}

test_result_t weighted_acat(const vector<double>& log10_pvals, const vector<double>& weights, double pmax_thr) {
  double log10p = get_acat_robust(
    -1*Map<const ArrayXd>(log10_pvals.data(), log10_pvals.size()),
    Map<const ArrayXd>(weights.data(), weights.size()),
    pmax_thr
  );
  if (log10p == 1) {
    return TEST_FAILED;
  } else {
    return test_result_t {
      max(pow(10, log10p), std::numeric_limits<double>::min()),
      log10p
    };
  }
}

test_result_t fishers(const vector<double>& log10_pvals) {
  double chi2 = 0;
  for (double log10p : log10_pvals) {
    chi2 += log10p / log10(exp(1));
  }
  chi2 *= -2;
  return chisq_dfk(chi2, 2*log10_pvals.size());
}

test_result_t stouffers(const vector<double>& log10_pvals, const vector<double>& weights, double adj_var) {
  boost::math::normal s;
  double numer = 0;
  double denom = 0;
  double stat = 0;
  double pv = 0;
  for (size_t i = 0; i < log10_pvals.size(); ++i) {
    if (log10_pvals[i] < log10(std::numeric_limits<double>::min())) {
      stat = sqrt(get_chisq_stat_from_logp(log10_pvals[i]));
    } else {
      pv = min(pow(10, log10_pvals[i]), 1-1e-7);
      stat = boost::math::quantile(boost::math::complement(s, pv));
    }
    numer += weights[i] * stat;
    denom += weights[i] * weights[i];
  }
  denom += adj_var;
  denom = sqrt(denom);
  double z = -numer / denom;
  return normal(z, false);
}

test_result_t stouffers_two_sided(const vector<double>& log10_pvals, const vector<double>& weights, const vector<double>& signs, double adj_var) {  
  boost::math::normal s;
  double numer = 0;
  double denom = 0;
  double stat = 0;
  double pv = 0;
  for (size_t i = 0; i < log10_pvals.size(); ++i) {
    if (log10_pvals[i] - log10(2) < log10(std::numeric_limits<double>::min())) {
      stat = sqrt(get_chisq_stat_from_logp(log10_pvals[i]) - log10(2));
    } else {
      pv = min(pow(10, log10_pvals[i]), 1-1e-7);
      stat = boost::math::quantile(boost::math::complement(s, pv/2));
    }
    numer += signs[i] * weights[i] * stat;
    denom += weights[i] * weights[i];
  }
  denom += adj_var;
  denom = sqrt(denom);
  double z = -abs(numer) / denom;
  return normal(z, true);
}

test_result_t mixture_of_chisq(const double& stat, const Eigen::VectorXd& weights) {
  double pv_davies_thr = 1e-5; // davies can be unreliable if pv is too small

  // re-scale so that max lambda is 1 (lambda is sorted)
  double newQ = stat / weights.tail(1)(0);
  VectorXd newL = weights / weights.tail(1)(0);
  test_result_t result = mixture_of_chisq_davies(newQ, newL, false);

  // if failed or is very low, use SPA
  if (result.pval <= pv_davies_thr || result == TEST_FAILED){ 
    result = mixture_of_chisq_kuonen(newQ, newL); // SPA
  }
  // if SPA failed
  if (result == TEST_FAILED) {
    result = mixture_of_chisq_davies(newQ, newL, true);
  }
  // only use mod Liu if Davies/SPA failed
  if (result == TEST_FAILED) {
    result = mixture_of_chisq_liu(newQ, newL);
  }

  return result;
}

test_result_t mixture_of_chisq_davies(const double& stat, Eigen::VectorXd& weights, const bool& force_stringent) {
  // use default lim/acc values from CompQuadForm R package and SKAT resp.
  int k = weights.size(), ifault = 0, lim = 1e4; // p & error
  double cdf, pv, acc1 = 1e-6;
  if (force_stringent) { lim=1e6; acc1 = 1e-9; }
  ArrayXd nc = ArrayXd::Constant(k, 0); // ncp
  ArrayXi df = ArrayXi::Constant(k, 1); // df
  ArrayXd tr = ArrayXd::Constant(7, 0); // params for qf

  try {
    cdf = qf(weights.data(), nc.data(), df.data(), k, 0, stat, lim, acc1, tr.data(), &ifault); 
    pv = 1 - cdf;
  } catch (...){
    return TEST_FAILED;
  }

  if((ifault != 0) || (pv <= 0) || (pv > 1))
    return TEST_FAILED;
  else {
    return test_result_t {pv, log10(pv)};
  }
}

// used by mixture_of_chisq_kounen
double get_tmin_lambda(const double& q, const ArrayXd& lambdas) {
  if (lambdas(0) < 0)  // not applicable here since matrix is psd
    return 1 / (2 * lambdas(0));
  else if (q > lambdas.sum())
    return 0;
  else
    return -0.5 * lambdas.size() / q;
}

double get_tmax_lambda(const ArrayXd& lambdas) {
  // return 1 / (2 * lambdas.tail(1)(0));
  return 0.5 - 1e-8;  // lambdas are re-scaled so max=1
}

double K_lambda(const double& t, const ArrayXd& lambdas) {
  return -0.5 * (1 - 2 * t * lambdas).log().sum();
}

double Kp_lambda(const double& t, const ArrayXd& lambdas) {
  return (lambdas / (1 - 2 * t * lambdas)).sum();
}

double Kpp_lambda(const double& t, const ArrayXd& lambdas) {
  return ((2 * lambdas.square()) / (1 - 2 * t * lambdas).square()).sum();
}

bool valid_bounds(double& fmin, double const& tmin, double const& tmax,
                  const double& q, const ArrayXd& lambdas) {
  fmin = Kp_lambda(tmin, lambdas) - q;
  double fmax = Kp_lambda(tmax, lambdas) - q;

  return ((fmin <= 0) && (fmax >= 0));
}

void solve_kp(bool& success, double& t_new, const double& q, const double& tmin,
              const double& tmax, const ArrayXd& lambdas) {
  int niter_cur = 0, niter_max = 1e3;
  double min_x, max_x, t_old, f_old, f_new, hess, tol = 1e-8;

  min_x = tmin, max_x = tmax;
  t_old = min_x;
  // check sign switches
  if (!valid_bounds(f_old, min_x, tmax, q, lambdas)) {
    success = false;
    return;
  }

  while (niter_cur++ < niter_max) {
    hess = Kpp_lambda(t_old, lambdas);
    t_new = t_old - f_old / hess;
    f_new = Kp_lambda(t_new, lambdas) - q;

    // cerr << "#" << niter_cur << ": t=" << t_old << "->" << t_new << " f(t)="
    // << f_new << "; bounds = (" << min_x << "," << max_x << ")\n";
    if (fabs(f_new) < tol) break;

    // update bounds on root
    if ((t_new > min_x) && (t_new < max_x)) {
      if (f_new > 0)
        max_x = t_new;
      else
        min_x = t_new;
    } else {  // bisection method if t_new went out of bounds and re-compute
              // f_new
      t_new = (min_x + max_x) * 0.5;
      f_new = Kp_lambda(t_new, lambdas) - q;
      if (f_new <= 0)
        min_x = t_new;  // reduce interval
      else
        max_x = t_new;
    }

    t_old = t_new;
    f_old = f_new;
  }

  // If didn't converge
  success = niter_cur <= niter_max;
  // cerr << "#iterations = " << niter_cur << "; f= " << f_new << endl;
}

double get_spa_pv(const double& root, const double& q, const ArrayXd& lambdas) {
  double u, w, r;
  boost::math::normal nd(0, 1);

  w = sgn(root) * sqrt(2 * (q * root - K_lambda(root, lambdas)));
  u = root * sqrt(Kpp_lambda(root, lambdas));
  if (fabs(u) < 1e-4) return -1;
  // cout << root << " " << q << " " << w << " " << u << " " <<  K_lambda(root,
  // lambdas) << endl;
  r = w + log(u / w) / w;
  return boost::math::cdf(boost::math::complement(nd, r));
}

test_result_t mixture_of_chisq_kuonen(const double& stat, const Eigen::VectorXd& weights) {
  bool success = false;
  double pv, t_root = -1;
  MapcArXd lambdas(weights.data(), weights.size(), 1);

  // lambdas are sorted in increasing order (from eigen)
  double tmin = get_tmin_lambda(stat, lambdas);
  double tmax = get_tmax_lambda(lambdas);
  if(tmax < tmin) return TEST_FAILED;

  solve_kp(success, t_root, stat, tmin, tmax, lambdas);
  if(!success) return TEST_FAILED;

  try {
    pv = get_spa_pv(t_root, stat, lambdas);
  } catch (...) {
    return TEST_FAILED;
  }

  if((pv <= 0) || (pv > 1)) return TEST_FAILED;

  return test_result_t {pv, log10(pv)};
}

// used by mixture_of_chisq_liu
void get_cvals(ArrayXd& cvals, const VectorXd& lambdas) {
  // cvals = [muQ, invsQ, muX, sX, df, ncp]
  double c1, c2, c3, c4, s1, s1_sq, s2, df, ncp, a;

  c1 = lambdas.sum();
  c2 = lambdas.squaredNorm();
  c3 = lambdas.array().pow(3).sum();
  c4 = lambdas.array().pow(4).sum();
  s1 = c3 / c2 / sqrt(c2);
  s1_sq = s1 * s1;
  s2 = c4 / (c2 * c2);
  if (s1_sq <= s2) {
    df = 1 / s2;
    a = sqrt(df);
    ncp = 0;
  } else {
    a = 1 / (s1 - sqrt(s1_sq - s2));
    ncp = (s1 * a - 1) * a * a;
    df = a * a - 2 * ncp;
  }

  cvals(0) = c1;                // muQ
  cvals(1) = 1 / sqrt(2 * c2);  // invsQ
  cvals(2) = df + ncp;          // muX
  cvals(3) = sqrt(2) * a;       // sX
  cvals(4) = df;
  cvals(5) = ncp;
}

test_result_t mixture_of_chisq_liu(const double& stat, const Eigen::VectorXd& weights) {
  ArrayXd cvals(6);
  get_cvals(cvals, weights);
  
  double pv;
  double tstar = (stat - cvals(0)) * cvals(1);
  double val = tstar * cvals(3) + cvals(2);

  if (val < 0) {
    return TEST_FAILED;
  }

  // 0 ncp gives strange behavior with non_central_chi_squared (returns -cdf instead of 1-cdf)
  if (boost::math::isnan(cvals(4)) || boost::math::isinf(cvals(4))) {
    return TEST_FAILED;
  } else if (cvals(5) == 0) {
    return chisq_dfk(val, cvals(4));
  } else if (boost::math::isnan(cvals(5)) || boost::math::isinf(cvals(5))) {
    return TEST_FAILED;
  } else {
    pv = boost::math::cdf(
      boost::math::complement(
        boost::math::non_central_chi_squared(cvals(4), cvals(5)), val
      )
    );
    if (pv <= 0 || pv > 1) {
      return TEST_FAILED;
    } else {
      return test_result_t {pv, log10(pv)};
    }
  }
}

}