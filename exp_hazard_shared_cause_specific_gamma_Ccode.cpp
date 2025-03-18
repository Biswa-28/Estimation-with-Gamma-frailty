// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#include <vector>
#include <iostream>
using namespace Numer;

class MargSubDist: public Func{
public:
  MargSubDist(std::vector<double> mu,std::vector<double> tau,int j){
    mu_vec = mu;
    tau_vec = tau;
    j_1 = j;
  }
public:
  std::vector<double> mu_vec;
  std::vector<double> tau_vec;
  int j_1;
  double operator()(const double& u) const
  {
    double m = 1;
    for(int i=0;i<mu_vec.size();i++){
      m = m*pow((1 + tau_vec[i]*u*(1/mu_vec[i])),-1/tau_vec[i]);
    }
    m = m*pow((1 + tau_vec[j_1]*u*(1/mu_vec[j_1])),-1);
    m = m*(1/mu_vec[j_1]);
    return m;
  }
};

// [[Rcpp::export]]
double marg_subdis_shared_cause_cpp(const double t,const int l,
    const std::vector<double> mu_vec_1,const std::vector<double> tau_vec_1)
{
  MargSubDist f(mu_vec_1,tau_vec_1,l - 1);
  double err_est;
  int err_code;
  const double value = integrate(f,0,t,err_est,err_code);
  return value;
};