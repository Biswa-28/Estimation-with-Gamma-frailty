// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <math.h>
#include <vector>
#include <iostream>
#include <RcppNumerical.h>
using namespace Numer;

class JointSubDist: public MFunc
{
  public:
    JointSubDist(double k_1, double k_2, std::vector<double> mu_1,
    std::vector<double> mu_2,std::vector<double> tau){
      j_1 = k_1;
      j_2 = k_2;
      mu_vec_1 = mu_1;
      mu_vec_2 = mu_2;
      tau_vec = tau;
    }
public:
   double j_1;
   double j_2;
   std::vector<double> mu_vec_1;
   std::vector<double> mu_vec_2;
   std::vector<double> tau_vec;
   // doubleegrand of jodouble sub-distribution function
   double operator()(Constvec& x)
   {
     double c_1 = x[0]*(1/mu_vec_1[0]);
     double c_2 = x[1]*(1/mu_vec_2[0]);
     double c_3 = x[0]*(1/mu_vec_1[1]);
     double c_4 = x[1]*(1/mu_vec_2[1]);
     double c_5 = c_1 + c_2;
     double c_6 = c_3 + c_4;
     double b_1 = 1.00000 + (tau_vec[0]*c_5);
     double b_2 = 1.00000 + (tau_vec[1]*c_6);
     if(j_1 == j_2)
     {
       double j_3 = j_1;
       if(j_3 == 0)
       {
         double d_1 = pow(b_1,-2.00000 -double(1/tau_vec[0]));
         double d_2 = pow(b_2,-double(1/tau_vec[1]));
         double d_3 = (1/mu_vec_1[0])*(1/mu_vec_2[0])*(1.00000 + tau_vec[0])*d_1*d_2;
         return d_3;
       }
       else
       {
         double e_1 = pow(b_1,-double(1/tau_vec[0]));
         double e_2 = pow(b_2,-2.00000 -double(1/tau_vec[1]));
         double e_3 = (1/mu_vec_1[1])*(1/mu_vec_2[1])*(1.00000 + tau_vec[1])*e_1*e_2;
         return e_3;
       }
     }
     else
     {
       double f_1 = pow(b_1,-1.00000 -double(1/tau_vec[0]));
       double f_2 = pow(b_2,-1.00000 -double(1/tau_vec[1]));
       double f_3 = (1/mu_vec_1[j_1])*(1/mu_vec_2[j_2])*f_1*f_2;
       return f_3;
     }
   }
};
class jointsubdist_in: public Func
{
public:
  jointsubdist_in(double k_1, double k_2, double y_1,std::vector<double> mu_1,
  std::vector<double> mu_2,std::vector<double> tau){
    j_1 = k_1;
    j_2 = k_2;
    y = y_1;
    mu_vec_1 = mu_1;
    mu_vec_2 = mu_2;
    tau_vec = tau;
  }
public:
  double j_1;
  double j_2;
  double y;
  std::vector<double> mu_vec_1;
  std::vector<double> mu_vec_2;
  std::vector<double> tau_vec;
  // doubleegrand of jodouble sub-distribution function
  double operator()(const double& x) const
  {
    double c_1 = x*(1/mu_vec_1[0]);
    double c_2 = y*(1/mu_vec_2[0]);
    double c_3 = x*(1/mu_vec_1[1]);
    double c_4 = y*(1/mu_vec_2[1]);
    double c_5 = c_1 + c_2;
    double c_6 = c_3 + c_4;
    double b_1 = 1 + (tau_vec[0]*c_5);
    double b_2 = 1 + (tau_vec[1]*c_6);
    if(j_1 == j_2)
    {
      double j_3 = j_1;
      if(j_3 == 0)
      {
        double d_1 = pow(b_1,-2.00000 -double(1/tau_vec[0]));
        double d_2 = pow(b_2,-double(1/tau_vec[1]));
        double d_3 = (1/mu_vec_1[0])*(1/mu_vec_2[0])*(1.00000 + tau_vec[0])*d_1*d_2;
        return d_3;
      }
      else
      {
        double e_1 = pow(b_1,-double(1/tau_vec[0]));
        double e_2 = pow(b_2,-2.00000 -double(1/tau_vec[1]));
        double e_3 = (1/mu_vec_1[1])*(1/mu_vec_2[1])*(1.00000 + tau_vec[1])*e_1*e_2;
        return e_3;
      }
    }
    else
    {
      double f_1 = pow(b_1,-1.00000 -double(1/tau_vec[0]));
      double f_2 = pow(b_2,-1.00000 -double(1/tau_vec[1]));
      double f_3 = (1/mu_vec_1[j_1])*(1/mu_vec_2[j_2])*f_1*f_2;
      return f_3;
    }
  }
};
class jointsubdist: public Func
{
  public:
  jointsubdist(double k_3,double k_4,double t_1,std::vector<double> mu_3,
  std::vector<double> mu_4,std::vector<double> tau_1)
  {
    k_5 = k_3;
    k_6 = k_4;
    t_3 = t_1;
    mu_vec_3 = mu_3;
    mu_vec_4 = mu_4;
    tau_vec_2 = tau_1;
  }
  double k_5;
  double k_6;  
  double t_3;
  std::vector<double> mu_vec_3;
  std::vector<double> mu_vec_4;
  std::vector<double> tau_vec_2;
  double operator()(const double& v) const
  {
    jointsubdist_in f(k_5,k_6,v,mu_vec_3,mu_vec_4,tau_vec_2);
    double err_est;
    int err_code;
    const double subdiv = 100; 
    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod201;
    const double value = integrate(f,0,t_3,err_est,err_code,subdiv,rule);
    return value;
  }
};
// [[Rcpp::export]]
Rcpp::List joint_sub_dis_func_cpp_1(const double t_5, const double t_6,
const double l_3, const double l_4,const std::vector<double> mu_vec_5,
const std::vector<double> mu_vec_6,const std::vector<double> tau_vec_3)
{
  jointsubdist f(l_3 - 1,l_4 - 1,t_5,mu_vec_5,mu_vec_6,tau_vec_3);
  double err_est;
  int err_code;
  const double upper = t_6;
  const double subdiv = 100;
  const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod201;
  const double res = integrate(f, 0, upper, err_est, err_code, subdiv, rule);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}
// [[Rcpp::export]]
Rcpp::List joint_sub_dis_func_cpp(const double t_1, const double t_2,
const double l_1, const double l_2,const std::vector<double> mu_vec_3,
const std::vector<double> mu_vec_4,const std::vector<double> tau_vec_1){
  JointSubDist f(l_1 - 1,l_2 - 1,mu_vec_3,mu_vec_4,tau_vec_1);
  Eigen::VectorXd lower(2);
  lower << 0, 0;
  Eigen::VectorXd upper(2);
  upper << t_1, t_2;
  double err_est;
  int err_code;
  const double maxeval = 100000;
  const double value = integrate(f, lower, upper, err_est, err_code, maxeval);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = value,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code);
}
