#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppNumerical.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <iostream>
#include <cmath>
#include <Rdefines.h>

using namespace arma;
using namespace Numer;
using namespace std;
using namespace Rcpp;


// =================================================================================================
// The Slash Distribution
// =================================================================================================

// Probability density function --------------------------------------------------------------------

// [[Rcpp::export]]
arma::vec slash_PDF(const arma::vec& x, const arma::vec& mu, const arma::vec& sigma,
                    arma::vec& nu, const bool log_p = false) {

  arma::uword n = x.n_rows;
  arma::vec result(n);
  arma::vec m = (x - mu) / sigma;

  for (arma::uword i = 0; i < n; ++i) {

    if (x(i) == mu(i)){

      result[i] = -log(sigma(i)) + log(2 * nu(i)) - log(2 * nu(i) + 1) - log(sqrt(2 * M_PI));

    }else{

      result[i] = -log(sigma(i)) + log(nu(i)) + nu(i) * log(2) +
        R::lgammafn(nu(i) + 0.5) +
        R::pgamma(arma::as_scalar(0.5 * pow(m(i), 2)), nu(i) + 0.5, 1.0, true, true) -
        log(sqrt(M_PI)) - (nu(i) + 0.5) * log(pow(m(i), 2));
    }

  }

  if (log_p) {
    return result;
  }else{
    return exp(result);
  }

}

// Cumulative distribution function ----------------------------------------------------------------
class slashIntegrand: public Func
{
private:
  double qtl;
  double nu;
public:
  slashIntegrand(double qtl_, double nu_) : qtl(qtl_), nu(nu_) {}

  double operator()(const double& x) const
  {
    return 2.0 * nu * std::pow(x, 2.0 * nu - 1.0) * R::pnorm(x * qtl, 0, 1, true, false);
  }
};


// [[Rcpp::export]]
arma::vec slash_CDF(const arma::vec& qtl, const arma::vec& mu, const arma::vec& sigma, arma::vec& nu)
{
  int n = qtl.n_rows;
  const double lower = 0, upper = 1;
  double err_est;
  int err_code;

  arma::vec m = (qtl - mu) / sigma;

  arma::vec result(n);

  for (arma::uword i = 0; i < n; ++i) {

    slashIntegrand f(m(i), nu(i));

    result[i] = arma::as_scalar(integrate(f, lower, upper, err_est, err_code));
  }

  return result;

}

// Quantile distribution function ------------------------------------------------------------------

class slash_AUX: public MFuncGrad
{
private:
  double p;
  double nu;
public:

  slash_AUX(double p_, double nu_) : p(p_), nu(nu_) {}

  double f_grad(Constvec& x, Refvec grad)
  {

    const double lower = 0, upper = 1;
    double err_est;
    int err_code;

    slashIntegrand f(x[0], nu);

    double aux1;

    if (x[0] == 0.0){

      aux1 = log(2 * nu) - log(2 * nu + 1) - log(sqrt(2 * M_PI));

    }else{

      aux1 = log(nu) + nu * log(2) + R::lgammafn(nu + 0.5) +
        R::pgamma(arma::as_scalar(0.5 * pow(x[0], 2)), nu + 0.5, 1.0, true, true) -
        log(sqrt(M_PI)) - (nu + 0.5) * log(pow(x[0], 2));
    }


    double aux2 = arma::as_scalar(integrate(f, lower, upper, err_est, err_code)) - p;

    grad[0] = 2 * aux2 * aux1;


    return aux2 * aux2;
  }
};


// [[Rcpp::export]]
arma::vec slash_QTF(const arma::vec& p, const arma::vec& nu)
{

  int n = p.n_rows;
  arma::vec result(n);

  Eigen::VectorXd x(1);

  for (arma::uword i = 0; i < n; ++i) {

    x[0] = R::qnorm5(p(i), 0, 1, true, false);

    slash_AUX obj(p(i), nu(i));

    double fopt;
    int res = optim_lbfgs(obj, x, fopt);

    result[i] = x[0];
  }

  return result;

}



// double slashCDF_dbl(double& qtl, double& nu)
// {
//   const double lower = 0, upper = 1;
//   double err_est;
//   int err_code;
//
//   slashIntegrand f(qtl, nu);
//
//   return arma::as_scalar(integrate(f, lower, upper, err_est, err_code));
//
// }
//
// // [[Rcpp::export]]
// NumericVector slash_QTF(NumericVector p, NumericVector nu){
//
//   int n = p.length();
//   double a; double b; double t = 1e-10;
//   double c; double d; double e; double fa; double fb; double fc;
//   double m; double tol; double macheps = 2.220446049250313e-016;
//   double pp; double q; double r; double s; double sa; double sb;
//
//   NumericVector result(n);
//
//   for ( arma::uword i = 0; i < n; ++i ) {
//
//     if (p(i) > 0.5){
//
//       a = 0;
//       b = 50;
//       fa = 0.5 - p(i);
//       fb = slashCDF_dbl(b, nu(i));
//       fb = fb - p(i);
//
//     } else{
//
//       a = -50;
//       b = 0;
//       fa = slashCDF_dbl(a, nu(i));
//       fa = fa - p(i);
//       fb = 0.5 - p(i);
//     }
//
//     sa = a;
//     sb = b;
//     c = sa;
//     fc = fa;
//     e = sb - sa;
//     d = e;
//
//     if (fa * fb > 0.0) stop("f() values at end points not of opposite sign.");
//
//     for ( ; ; )
//     {
//       if (abs(fc) < abs(fb)){
//         sa = sb;
//         sb = c;
//         c = sa;
//         fa = fb;
//         fb = fc;
//         fc = fa;
//       }
//       tol = 2.0*macheps*abs(sb) + t;
//       m = 0.5*(c - sb);
//
//       if (abs(m) <= tol || fb == 0.0 ){
//         break;
//       }
//
//       if (abs(e) < tol || abs(fa) <= abs(fb)){
//         e = m; d = e;
//       } else {
//         s = fb/fa;
//         if (sa == c){
//           pp = 2.0*m*s;
//           q = 1.0 - s;
//         } else {
//           q = fa/fc;
//           r = fb/fc;
//           pp = s*(2.0*m*q*(q - r) - (sb - sa)*(r - 1.0));
//           q = (q - 1.0)*(r - 1.0)*(s - 1.0);
//         }
//         if ( 0.0 < pp ){
//           q = - q;
//         } else {
//           pp = - pp;
//         }
//         s = e;
//         e = d;
//         if (2.0*pp < 3.0*m*q - abs(tol*q) && pp < abs(0.5*s*q)){
//           d = pp/q;
//         } else {
//           e = m;
//           d = e;
//         }
//       }
//       sa = sb;
//       fa = fb;
//       if (tol < abs(d)){
//         sb = sb + d;
//       } else if ( 0.0 < m ){
//         sb = sb + tol;
//       } else {
//         sb = sb - tol;
//       }
//       fb = slashCDF_dbl(sb,nu(i));
//       fb = fb - p(i);
//       if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0)){
//         c = sa;
//         fc = fa;
//         e = sb - sa;
//         d = e;
//       }
//     }
//
//     result[i] = sb;
//
//   }
//
//   return result;
//
// }

// double aux_slash(double& qtl, double& p, double& nu)
// {
//   const double lower = 0.0, upper = 1.0;
//   double err_est;
//   int err_code;
//
//   slashIntegrand f(qtl, nu);
//
//   return arma::as_scalar(integrate(f, lower, upper, err_est, err_code)) - p;
//
// }
//
// // [[Rcpp::export]]
// arma::vec slash_QTF2(const arma::vec& p, arma::vec& nu)
// {
//   int n = p.n_rows;
//   int j;
//   arma::vec rg(2);
//
//   // Extract R's uniroot function
//   Rcpp::Environment stats("package:stats");
//   Rcpp::Function uniroot = stats["uniroot"];
//
//   arma::vec result(n);
//
//   for(arma::uword i = 0; i < n; ++i) {
//
//     if(0.5 < p(i)){
//       rg = {0, 1};
//       j = 2;
//       while(slashCDF_dbl(rg(1), nu(i)) < p(i)){
//         rg(1) = j;
//         j++;
//       }
//     }else{
//       rg = {-1, 0};
//       j = 2;
//       while(slashCDF_dbl(rg(0), nu(i)) > p(i)){
//         rg(0) = -j;
//         j++;
//       }
//     }
//
//     Rcpp::List opt_results = uniroot(Rcpp::_["f"] = Rcpp::InternalFunction(&aux_slash),
//                                      Rcpp::_["lower"] = rg(0),
//                                      Rcpp::_["upper"] = rg(1),
//                                      // Pass in the other parameters as everything
//                                      // is scoped environmentally
//                                      Rcpp::_["p"] = p(i),
//                                      Rcpp::_["nu"] = nu(i));
//
//     result[i] = opt_results[0];
//   }
//
//   return result;
//
// }







