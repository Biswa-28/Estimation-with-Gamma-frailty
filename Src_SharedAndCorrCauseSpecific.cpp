// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#include <vector>
using namespace Numer;

class margsubdis_exp_hazard_shared_cause_specific: public Func{
public:
  margsubdis_exp_hazard_shared_cause_specific(int j,
                std::vector<double> mu_vec, std::vector<double> tau_vec)
  {
    mu_vec_1 = mu_vec;
    j_1 = j;
    tau_vec_1 = tau_vec;
  }
  std::vector<double> mu_vec_1;
  int j_1;
  std::vector<double> tau_vec_1;
  double operator()(const double& u) const
  {
    double prod = 1.00000;
    for(int i = 0;i < mu_vec_1.size();i++)
    {
      double d_1 = log(tau_vec_1[i]) + log(u) + log(mu_vec_1[i]);
      double d_2 = exp(d_1);
      double c = 1.00000 + d_2;
      double d_3 = -1.000000*log(tau_vec_1[i]);
      double d_4 = exp(d_3);
      double d_5 = -d_4*log(c);
      prod = prod*exp(d_5);
    }
    double c_1 = 1.00000 + (tau_vec_1[j_1]*u*mu_vec_1[j_1]);
    double c_2 = pow(c_1,-1.000000);
    double integrand = mu_vec_1[j_1]*c_2*prod;
    return integrand;
  }
};

class margsubdis_weibull_hazard_shared_cause_specific: public Func{
public:
  margsubdis_weibull_hazard_shared_cause_specific(int j, double beta,
                                                  std::vector<double> mu_vec, std::vector<double> tau_vec)
  {
    mu_vec_1 = mu_vec;
    beta_1 = beta;
    j_1 = j;
    tau_vec_1 = tau_vec;
  }
  std::vector<double> mu_vec_1;
  int j_1;
  std::vector<double> tau_vec_1;
  double beta_1;
  double operator()(const double& u) const
  {
    double prod = 1.00000;
    for(int i = 0;i < mu_vec_1.size();i++)
    {
      double c = pow(u*mu_vec_1[i],beta_1);
      double c_1 = 1.00000 + (tau_vec_1[i]*c);
      prod = prod*pow(c_1,-pow(tau_vec_1[i],-1.000000));
    }
    double c_2 = pow(u*mu_vec_1[j_1],beta_1);
    double c_3 = pow((1.00000 + tau_vec_1[j_1]*c_2),-1.00000);
    double c_4 = beta_1*pow(mu_vec_1[j_1],beta_1)*pow(u,beta_1 - 1.00000);
    double integrand = c_3*c_4*prod;
    return integrand;
  }
};

class marg_subdis_exp_corr_cause: public Func{
public:
  marg_subdis_exp_corr_cause(double t_3,std::vector<double> mu_3,
                             int j_3,std::vector<double> tau_3)
  {
    t_5 = t_3;
    mu_vec_3 = mu_3;
    j_5 = j_3;
    tau_vec_3 = tau_3;
  }
  double t_5;
  std::vector<double> mu_vec_3;
  int j_5;
  std::vector<double> tau_vec_3;
  double operator()(const double& u) const
  {
    double prod = 1.00000;
    for(int i = 0;i < mu_vec_3.size();i++)
    {
      prod = prod*pow((1.00000 + (tau_vec_3[i]*tau_vec_3[i]*u*(1.00000/mu_vec_3[i]))),-pow(tau_vec_3[i]*tau_vec_3[i],-1.000000));
    }
    double integrand = (1.00000/mu_vec_3[j_5])*prod*pow((1.00000 + (tau_vec_3[j_5]*tau_vec_3[j_5]*u*(1.00000/mu_vec_3[j_5]))),-1.00000);
    return integrand;
  }
};

class marg_subdis_weibull_corr_cause: public Func{
public:
  marg_subdis_weibull_corr_cause(int j,double alpha,
                                 std::vector<double> mu,std::vector<double> tau)
  {
    j_1 = j;
    beta = alpha;
    mu_vec = mu;
    tau_vec = tau;
  }
  int j_1;
  double beta;
  std::vector<double> mu_vec;
  std::vector<double> tau_vec;
  double operator()(const double& u) const
  {
    double log_marg_surv = 0.00000;
    int L = mu_vec.size();
    std::vector<double> tau_vec_sq(L);
    std::vector<double> tau_vec_sq_inv(L);
    std::vector<double> a_1(L);
    for(int i = 0;i < L;i++)
    {
      tau_vec_sq[i] = tau_vec[i]*tau_vec[i];
      tau_vec_sq_inv[i] = -1.000000*std::log(tau_vec_sq[i]);
      tau_vec_sq_inv[i] = std::exp(tau_vec_sq_inv[i]);
      a_1[i] = beta*(std::log(mu_vec[i]) + std::log(u));
      a_1[i] = std::exp(a_1[i]);
      a_1[i] = 1.000000 + (tau_vec_sq[i]*a_1[i]);
      log_marg_surv = log_marg_surv - (tau_vec_sq_inv[i]*std::log(a_1[i]));
    }
    double log_con = std::log(beta) + (beta*std::log(mu_vec[j_1]));
    double log_var = ((beta - 1.00000)*std::log(u)) - (tau_vec_sq_inv[j_1]*std::log(a_1[j_1]));
    double log_integrand = log_con + log_var + log_marg_surv;
    double integrand = std::exp(log_integrand);
    return integrand;
  }
};

class joint_sub_dis_exp_shared_cause_specific: public MFunc{
public:
  joint_sub_dis_exp_shared_cause_specific(int j_1,int j_2,
                                          std::vector<double> mu_vec_1,std::vector<double> mu_vec_2,std::vector<double> tau,int l)
  {
    mu_vec_3 = mu_vec_1;
    mu_vec_4 = mu_vec_2;
    j_3 = j_1;
    j_4 = j_2;
    tau_vec = tau;
    L = l;
  }
  std::vector<double> mu_vec_3;
  std::vector<double> mu_vec_4;
  int j_3;
  int j_4;
  std::vector<double> tau_vec;
  int L;
  double operator()(Constvec& u)
  {
    std::vector<double> a(L);
    std::vector<double> b(L);
    std::vector<double> c(L);
    // first defining the joint survival function
    double log_joint_surv = 0.0000000;
    for(int i = 0;i < L;i++)
    {
      a[i] = 1.00000 + (tau_vec[i]*((u[0]*mu_vec_3[i]) + (u[1]*mu_vec_4[i])));
      b[i] = -1.000000*std::log(tau_vec[i]);
      b[i] = std::exp(b[i]);
      c[i] = -1.00000*b[i]*std::log(a[i]);
      log_joint_surv = log_joint_surv + c[i];
    }
    if(j_3 == j_4)
    {
      int j = j_3;
      double log_con = std::log(1 + tau_vec[j]) + std::log(mu_vec_3[j]) + std::log(mu_vec_4[j]);
      double log_f = -2.00000*std::log(a[j]);
      double log_integrand = log_con + log_joint_surv + log_f;
      double integrand = std::exp(log_integrand);
      return integrand;
    }
    else
    {
      double log_con_1 = std::log(mu_vec_3[j_3]) + std::log(mu_vec_4[j_4]);
      double log_g = -1.00000*std::log(a[j_3]) - 1.00000*std::log(a[j_4]);
      double log_integrand = log_con_1 + log_joint_surv + log_g;
      double integrand = std::exp(log_integrand);
      return integrand;
    }
  }
};

class joint_sub_dis_weibull_shared_cause_specific: public MFunc{
public:
  joint_sub_dis_weibull_shared_cause_specific(int j_1,int j_2,double alpha_1, double alpha_2,
                                              std::vector<double> mu_vec_1,std::vector<double> mu_vec_2,std::vector<double> tau,int l)
  {
    j_3 = j_1;
    j_4 = j_2;
    beta_1 = alpha_1;
    beta_2 = alpha_2;
    mu_vec_3 = mu_vec_1;
    mu_vec_4 = mu_vec_2;
    tau_vec = tau;
    p = l;
  }
  int j_3;
  int j_4;
  double beta_1;
  double beta_2;
  std::vector<double> mu_vec_3;
  std::vector<double> mu_vec_4;
  std::vector<double> tau_vec;
  int p;
  double operator()(Constvec& u)
  {
    std::vector<double> a_1(p);
    std::vector<double> a_2(p);
    std::vector<double> a_3(p);
    std::vector<double> tau_vec_inv(p);
    double log_joint_surv = 0.000000;
    for(int i = 0;i < p;i++)
    {
      tau_vec_inv[i] = -1.000000*std::log(tau_vec[i]);
      tau_vec_inv[i] = std::exp(tau_vec_inv[i]);
      a_1[i] = (beta_1*std::log(mu_vec_3[i])) + (beta_1*std::log(u[0]));
      a_1[i] = std::exp(a_1[i]);
      a_2[i] = (beta_2*std::log(mu_vec_4[i])) + (beta_2*std::log(u[1]));
      a_2[i] = std::exp(a_2[i]);
      a_3[i] = 1.000000 + tau_vec[i]*(a_1[i] + a_2[i]);
      log_joint_surv = log_joint_surv - (tau_vec_inv[i]*std::log(a_3[i]));
    }
    if(j_3 == j_4)
    {
      int j = j_3;
      double log_con = std::log(1.00000 + tau_vec[j]) + std::log(beta_1) + std::log(beta_2)
        + (beta_1*std::log(mu_vec_3[j])) + (beta_2*std::log(mu_vec_4[j]));
      double log_var = ((beta_1 - 1.00000)*std::log(u[0])) + ((beta_2 - 1.00000)*std::log(u[1]));
      double log_var_1 = -2.00000*std::log(a_3[j]);
      double log_integrand = log_con + log_var + log_var_1 + log_joint_surv;
      double integrand = std::exp(log_integrand);
      return integrand;
    }
    else
    {
      double log_con = std::log(beta_1) + std::log(beta_2) + (beta_1*std::log(mu_vec_3[j_3])) +
        (beta_2*std::log(mu_vec_4[j_4]));
      double log_var = ((beta_1 - 1.00000)*std::log(u[0])) + ((beta_2 - 1.00000)*std::log(u[1]));
      double log_var_1 = -std::log(a_3[j_3]) - std::log(a_3[j_4]);
      double log_integrand = log_con + log_var + log_var_1 + log_joint_surv;
      double integrand = std::exp(log_integrand);
      return integrand;
    }
  }
};

class joint_sub_dis_exp_corr_cause: public MFunc{
public:
  joint_sub_dis_exp_corr_cause(int j_3,int j_4,
                               std::vector<double> mu_3,std::vector<double> mu_4,std::vector<double> tau,
                               std::vector<double> rho,int l)
  {
    mu_vec_3 = mu_3;
    mu_vec_4 = mu_4;
    j_5 = j_3;
    j_6 = j_4;
    tau_vec = tau;
    rho_vec = rho;
    p = l;
  }
  std::vector<double> mu_vec_3;
  std::vector<double> mu_vec_4;
  int j_5;
  int j_6;
  std::vector<double> tau_vec;
  std::vector<double> rho_vec;
  int p;
  double operator()(Constvec& u)
  {
    std::vector<double> tau_vec_sq(p);
    std::vector<double> tau_vec_sq_inv(p);
    std::vector<double> z(p);
    std::vector<double> w(p);
    std::vector<double> d_1(p);
    std::vector<double> d_2(p);
    std::vector<double> d_3(p);
    // first defining the joint survival function
    double log_joint_surv = 0.0000000;
    for(int i = 0;i < p;i++)
    {
      tau_vec_sq[i] = tau_vec[i]*tau_vec[i];
      tau_vec_sq_inv[i] = -1.000000*std::log(tau_vec_sq[i]);
      tau_vec_sq_inv[i] = std::exp(tau_vec_sq_inv[i]);
      z[i] = rho_vec[i]*tau_vec_sq_inv[i];
      w[i] = (1.00000 - rho_vec[i])*tau_vec_sq_inv[i];
      d_1[i] = 1.00000 + (tau_vec_sq[i]*u[0]*mu_vec_3[i]);
      d_2[i] = 1.00000 + (tau_vec_sq[i]*u[1]*mu_vec_4[i]);
      d_3[i] = d_1[i] + d_2[i] - 1.00000;
      log_joint_surv = log_joint_surv - (z[i]*std::log(d_3[i]));
      log_joint_surv = log_joint_surv - (w[i]*std::log(d_1[i]));
      log_joint_surv = log_joint_surv - (w[i]*std::log(d_2[i]));
    }
    if(j_5 == j_6)
    {
      int j = j_5;
      double log_con = (2.00000*std::log(tau_vec_sq[j])) + std::log(mu_vec_3[j]) + std::log(mu_vec_4[j]);
      double e_1 = z[j];
      double log_f_1 = std::log(e_1) + std::log(1 + e_1) - (2.00000*std::log(d_3[j]));
      double f_1 = std::exp(log_f_1);
      double log_f_2 = std::log(e_1) + std::log(w[j]) - (1.00000*std::log(d_2[j])) - (1.00000*std::log(d_3[j]));
      double f_2 = std::exp(log_f_2);
      double log_f_3 = std::log(e_1) + std::log(w[j]) - (1.00000*std::log(d_1[j])) - (1.00000*std::log(d_3[j]));
      double f_3 = std::exp(log_f_3);
      double log_f_4 = (2.00000*std::log(w[j])) - (1.00000*std::log(d_1[j])) - (1.00000*std::log(d_2[j]));
      double f_4 = std::exp(log_f_4);
      double log_integrand = log_con + log_joint_surv + std::log(f_1 + f_2 + f_3 + f_4);
      double integrand = std::exp(log_integrand);
      return integrand;
    }
    else
    {
      double log_con_1 = std::log(tau_vec_sq[j_5]) + std::log(tau_vec_sq[j_6]) + std::log(mu_vec_3[j_5]) + std::log(mu_vec_4[j_6]);
      double e_3 = z[j_5];
      double e_4 = z[j_6];
      double log_g_1 = std::log(e_3) + std::log(e_4) - (1.00000*std::log(d_3[j_5])) - (1.00000*std::log(d_3[j_6]));
      double g_1 = std::exp(log_g_1);
      double log_g_2 = std::log(e_3) + std::log(w[j_6]) - (1.00000*std::log(d_2[j_6])) - (1.00000*std::log(d_3[j_5]));
      double g_2 = std::exp(log_g_2);
      double log_g_3 = std::log(e_4) + std::log(w[j_5]) - (1.00000*std::log(d_1[j_5])) - (1.00000*std::log(d_3[j_6]));
      double g_3 = std::exp(log_g_3);
      double log_g_4 = std::log(w[j_5]) + std::log(w[j_6]) - (1.00000*std::log(d_1[j_5])) - (1.00000*std::log(d_2[j_6]));
      double g_4 = std::exp(log_g_4);
      double log_integrand = log_con_1 + log_joint_surv + std::log(g_1 + g_2 + g_3 + g_4);
      double integrand = std::exp(log_integrand);
      return integrand;
    }
  }
};

class joint_sub_dis_exp_corr_cause_full: public MFunc{
public:
  joint_sub_dis_exp_corr_cause_full(int j_1,int j_2,
                                    std::vector<double> mu_1,std::vector<double> mu_2,std::vector<double> tau_1,
                                    std::vector<double> tau_2,std::vector<double> rho,int l)
  {
    mu_vec_1 = mu_1;
    mu_vec_2 = mu_2;
    j_3 = j_1;
    j_4 = j_2;
    tau_vec_1 = tau_1;
    tau_vec_2 = tau_2;
    rho_vec = rho;
    p = l;
  }
  std::vector<double> mu_vec_1;
  std::vector<double> mu_vec_2;
  int j_3;
  int j_4;
  std::vector<double> tau_vec_1;
  std::vector<double> tau_vec_2;
  std::vector<double> rho_vec;
  int p;
  double operator()(Constvec& u)
  {
    std::vector<double> tau_vec_1_inv(p);
    std::vector<double> tau_vec_2_inv(p);
    std::vector<double> z(p);
    std::vector<double> w_1(p);
    std::vector<double> w_2(p);
    std::vector<double> d_1(p);
    std::vector<double> d_2(p);
    std::vector<double> d_3(p);
    // first defining the joint survival function
    double log_joint_surv = 0.000000;
    for(int i = 0;i < p;i++)
    {
      tau_vec_1_inv[i] = -1.000000*std::log(tau_vec_1[i]);
      tau_vec_1_inv[i] = std::exp(tau_vec_1_inv[i]);
      tau_vec_2_inv[i] = -1.000000*std::log(tau_vec_2[i]);
      tau_vec_2_inv[i] = std::exp(tau_vec_2_inv[i]);
      z[i] = rho_vec[i]*tau_vec_1_inv[i]*tau_vec_2_inv[i];
      w_1[i] = (tau_vec_1_inv[i]*tau_vec_1_inv[i]) - z[i];
      w_2[i] = (tau_vec_2_inv[i]*tau_vec_2_inv[i]) - z[i];
      d_1[i] = 1.000000 + (tau_vec_1[i]*tau_vec_1[i]*u[0]*mu_vec_1[i]);
      d_2[i] = 1.000000 + (tau_vec_2[i]*tau_vec_2[i]*u[1]*mu_vec_2[i]);
      d_3[i] = d_1[i] + d_2[i] - 1.000000;
      log_joint_surv = log_joint_surv - (z[i]*std::log(d_3[i]));
      log_joint_surv = log_joint_surv - (w_1[i]*std::log(d_1[i]));
      log_joint_surv = log_joint_surv - (w_2[i]*std::log(d_2[i]));
    }
    if(j_3 == j_4)
    {
      int j = j_3;
      double log_con = (2.00000*(std::log(tau_vec_1[j]) + std::log(tau_vec_2[j]))) + std::log(mu_vec_1[j]) + std::log(mu_vec_2[j]);
      double e_1 = z[j];
      double log_f_1 = std::log(e_1) + std::log(1 + e_1) - (2.00000*std::log(d_3[j]));
      double f_1 = std::exp(log_f_1);
      double log_f_2 = std::log(e_1) + std::log(w_2[j]) - (1.00000*std::log(d_2[j])) - (1.00000*std::log(d_3[j]));
      double f_2 = exp(log_f_2);
      double log_f_3 = std::log(e_1) + std::log(w_1[j]) - (1.00000*std::log(d_1[j])) - (1.00000*std::log(d_3[j]));
      double f_3 = exp(log_f_3);
      double log_f_4 = (std::log(w_1[j]) + std::log(w_2[j])) - (1.00000*std::log(d_1[j])) - (1.00000*std::log(d_2[j]));
      double f_4 = std::exp(log_f_4);
      double log_integrand = log_con + log_joint_surv + std::log(f_1 + f_2 + f_3 + f_4);
      double integrand = std::exp(log_integrand);
      return integrand;
    }
    else
    {
      double log_con_1 = (2.000000*std::log(tau_vec_1[j_3])) + (2.000000*std::log(tau_vec_2[j_4])) + std::log(mu_vec_1[j_3]) + std::log(mu_vec_2[j_4]);
      double e_3 = z[j_3];
      double e_4 = z[j_4];
      double log_g_1 = std::log(e_3) + std::log(e_4) - (1.00000*std::log(d_3[j_3])) - (1.00000*std::log(d_3[j_4]));
      double g_1 = std::exp(log_g_1);
      double log_g_2 = std::log(e_3) + std::log(w_2[j_4]) - (1.00000*std::log(d_2[j_4])) - (1.00000*std::log(d_3[j_3]));
      double g_2 = std::exp(log_g_2);
      double log_g_3 = std::log(e_4) + std::log(w_1[j_3]) - (1.00000*std::log(d_1[j_3])) - (1.00000*std::log(d_3[j_4]));
      double g_3 = std::exp(log_g_3);
      double log_g_4 = std::log(w_1[j_3]) + std::log(w_2[j_4]) - (1.00000*std::log(d_1[j_3])) - (1.00000*std::log(d_2[j_4]));
      double g_4 = std::exp(log_g_4);
      double log_integrand = log_con_1 + log_joint_surv + std::log(g_1 + g_2 + g_3 + g_4);
      double integrand = std::exp(log_integrand);
      return integrand;
    }
  }
};


class joint_sub_dis_weibull_corr_cause: public MFunc{
public:
  joint_sub_dis_weibull_corr_cause(int j_1, int j_2, std::vector<double> mu_vec_1,
                                   std::vector<double> mu_vec_2, double alpha_1, double alpha_2, std::vector<double> tau,
                                   std::vector<double> rho,int L)
  {
    j_3 = j_1;
    j_4 = j_2;
    mu_vec_3= mu_vec_1;
    mu_vec_4= mu_vec_2;
    beta_1 = alpha_1;
    beta_2 = alpha_2;
    tau_vec = tau;
    rho_vec = rho;
    p = L;
  }
  int j_3;
  int j_4;
  std::vector<double> mu_vec_3;
  std::vector<double> mu_vec_4;
  double beta_1;
  double beta_2;
  std::vector<double> tau_vec;
  std::vector<double> rho_vec;
  int p;
  double operator()(Constvec& u)
  {
    std::vector<double> tau_vec_sq(p);
    std::vector<double> tau_vec_sq_inv(p);
    std::vector<double> z(p);
    std::vector<double> w(p);
    std::vector<double> d_1(p);
    std::vector<double> d_2(p);
    std::vector<double> d_3(p);
    double log_joint_surv = 0.000000;
    for(int i = 0;i < p;i++)
    {
      tau_vec_sq[i] = tau_vec[i]*tau_vec[i];
      tau_vec_sq_inv[i] = -1.000000*std::log(tau_vec_sq[i]);
      tau_vec_sq_inv[i] = std::exp(tau_vec_sq_inv[i]);
      z[i] = rho_vec[i]*tau_vec_sq_inv[i];
      w[i] = (1.00000 - rho_vec[i])*tau_vec_sq_inv[i];
      d_1[i] = beta_1*(std::log(mu_vec_3[i]) + std::log(u[0]));
      d_1[i] = std::exp(d_1[i]);
      d_1[i] = 1.000000 + (tau_vec_sq[i]*d_1[i]);
      d_2[i] = beta_2*(std::log(mu_vec_4[i]) + std::log(u[1]));
      d_2[i] = std::exp(d_2[i]);
      d_2[i] = 1.000000 + (tau_vec_sq[i]*d_2[i]);
      d_3[i] = d_1[i] + d_2[i] - 1.000000;
      log_joint_surv = log_joint_surv - (z[i]*std::log(d_3[i]));
      log_joint_surv = log_joint_surv - (w[i]*std::log(d_1[i]));
      log_joint_surv = log_joint_surv - (w[i]*std::log(d_2[i]));
    }
    if(j_3 == j_4)
    {
      int j = j_3;
      double log_con = (2.00000*std::log(tau_vec_sq[j])) + std::log(beta_1) +
        std::log(beta_2) + (beta_1*std::log(mu_vec_3[j])) + (beta_2*std::log(mu_vec_4[j]));
      double log_var = ((beta_1 - 1.00000)*std::log(u[0])) + ((beta_2 - 1.00000)*std::log(u[1]));
      double e_1 = z[j];
      double e_2 = w[j];
      double log_g_1 = std::log(e_1) + std::log(1.00000 + e_1) - (2.00000*std::log(d_3[j]));
      double g_1 = std::exp(log_g_1);
      double log_g_2 = std::log(e_1) + std::log(e_2) - std::log(d_3[j]) - std::log(d_1[j]);
      double g_2 = std::exp(log_g_2);
      double log_g_3 = std::log(e_1) + std::log(e_2) - std::log(d_3[j]) - std::log(d_2[j]);
      double g_3 = std::exp(log_g_3);
      double log_g_4 = (2.00000*std::log(e_2)) - std::log(d_1[j]) - std::log(d_2[j]);
      double g_4 = std::exp(log_g_4);
      double log_integrand = log_con + log_var + log_joint_surv + std::log(g_1 + g_2 + g_3 + g_4);
      double integrand = std::exp(log_integrand);
      return integrand;
    }
    else
    {
      double log_con_1 = std::log(tau_vec_sq[j_3]) + std::log(tau_vec_sq[j_4]) +
        std::log(beta_1) + std::log(beta_2) + (beta_1*std::log(mu_vec_3[j_3])) +
        (beta_2*std::log(mu_vec_4[j_4]));
      double log_var = ((beta_1 - 1.00000)*std::log(u[0])) + ((beta_2 - 1.00000)*std::log(u[1]));
      double e_3 = z[j_3];
      double e_4 = z[j_4];
      double log_h_1 = std::log(e_3) + std::log(e_4) - std::log(d_3[j_3]) - std::log(d_3[j_4]);
      double h_1 = std::exp(log_h_1);
      double log_h_2 = std::log(e_3) + std::log(w[j_4]) - std::log(d_3[j_3]) - std::log(d_2[j_4]);
      double h_2 = std::exp(log_h_2);
      double log_h_3 = std::log(e_4) + std::log(w[j_3]) - std::log(d_3[j_4]) - std::log(d_1[j_3]);
      double h_3 = std::exp(log_h_3);
      double log_h_4 = std::log(w[j_3]) + std::log(w[j_4]) - std::log(d_1[j_3]) - std::log(d_2[j_4]);
      double h_4 = std::exp(log_h_4);
      double log_integrand = log_con_1 + log_var + log_joint_surv + std::log(h_1 + h_2 + h_3 + h_4);
      double integrand = std::exp(log_integrand);
      return integrand;
    }
  }
};

// [[Rcpp::export]]
double marg_subdis_exp_shared_cause_specific_cpp(const double t,const int j,
                                                 const std::vector<double> mu_vec,const std::vector<double> tau_vec){
  margsubdis_exp_hazard_shared_cause_specific f(j - 1,mu_vec,tau_vec);
  double err_est;
  int err_code;
  const double value = integrate(f,0.000000,t,err_est,err_code);
  return value;
}

// [[Rcpp::export]]
double marg_subdis_weibull_shared_cause_specific_cpp(const double t,const int j,
                                                     const std::vector<double> mu_vec,const double beta,const std::vector<double> tau_vec){
  margsubdis_weibull_hazard_shared_cause_specific f(j - 1,beta,mu_vec,tau_vec);
  double err_est;
  int err_code;
  const double value = integrate(f,0.000000,t,err_est,err_code);
  return value;
}

// [[Rcpp::export]]
double marg_subdis_exp_corr_cause_specific_cpp(const double t,const int j,const std::vector<double> mu_vec,
                                    const std::vector<double> tau_vec) {
  marg_subdis_exp_corr_cause f(t,mu_vec,j - 1,tau_vec);
  double err_est;
  int err_code;
  const double value = integrate(f,0.00000,t,err_est,err_code);
  return value;
}

// [[Rcpp::export]]
double marg_subdis_weibull_corr_cause_specific_cpp(const double t,const int j,const std::vector<double> mu_vec,
                                          const double beta, const std::vector<double> tau_vec) {
  marg_subdis_weibull_corr_cause f(j - 1,beta,mu_vec,tau_vec);
  double err_est;
  int err_code;
  const double value = integrate(f,0.00000,t,err_est,err_code);
  return value;
}

// [[Rcpp::export]]
Rcpp::List joint_sub_dis_exp_shared_cause_specific_cpp(const double t_1, const double t_2,
                                                       const int j_1,const int j_2,const std::vector<double> mu_vec_1,const std::vector<double> mu_vec_2,
                                                       const std::vector<double> tau_vec,const int L)
{
  joint_sub_dis_exp_shared_cause_specific f(j_1 - 1,j_2 - 1,mu_vec_1,mu_vec_2,tau_vec,L);
  Eigen::VectorXd lower(2);
  lower << 0.000000, 0.000000;
  Eigen::VectorXd upper(2);
  upper << t_1, t_2;
  double err_est;
  int err_code;
  const int maxeval = 200000;
  const double eps_rel = 1e-13;
  const double res = integrate(f, lower, upper, err_est, err_code, maxeval, eps_rel);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}

// [[Rcpp::export]]
Rcpp::List joint_sub_dis_weibull_shared_cause_specific_cpp(const double t_1, const double t_2,
                                                           const int j_1,const int j_2,const std::vector<double> mu_vec_1,const std::vector<double> mu_vec_2,
                                                           const double beta_1,const double beta_2,const std::vector<double> tau_vec,const int L)
{
  joint_sub_dis_weibull_shared_cause_specific f(j_1 - 1,j_2 - 1,beta_1,beta_2,mu_vec_1,mu_vec_2,tau_vec,L);
  Eigen::VectorXd lower(2);
  lower << 0.000000, 0.000000;
  Eigen::VectorXd upper(2);
  upper << t_1, t_2;
  double err_est;
  int err_code;
  const int maxeval = 200000;
  const double eps_rel = 1e-13;
  const double res = integrate(f, lower, upper, err_est, err_code, maxeval, eps_rel);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}

// [[Rcpp::export]]
Rcpp::List joint_sub_dis_exp_corr_cause_specific_cpp(const double t_1, const double t_2,
                                          const int j_1,const int j_2,const std::vector<double> mu_vec_1,const std::vector<double> mu_vec_2,
                                          const std::vector<double> tau_vec_new,const std::vector<double> phi_vec,const int L) {
  joint_sub_dis_exp_corr_cause f(j_1 - 1,j_2 - 1,mu_vec_1,mu_vec_2,tau_vec_new,phi_vec,L);
  Eigen::VectorXd lower(2);
  lower << 0.000000, 0.000000;
  Eigen::VectorXd upper(2);
  upper << t_1, t_2;
  double err_est;
  int err_code;
  const int maxeval = 200000;
  const double eps_rel = 1e-13;
  const double res = integrate(f, lower, upper, err_est, err_code, maxeval, eps_rel);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}

// [[Rcpp::export]]
Rcpp::List joint_sub_dis_exp_corr_cause_specific_full_cpp(const double t_1, const double t_2,
                                                          const int j_1,const int j_2,const std::vector<double> mu_vec_1,const std::vector<double> mu_vec_2,
                                                          const std::vector<double> tau_vec_new_1,const std::vector<double> tau_vec_new_2,const std::vector<double> phi_vec,const int L) {
  joint_sub_dis_exp_corr_cause_full f(j_1 - 1,j_2 - 1,mu_vec_1,mu_vec_2,tau_vec_new_1,tau_vec_new_2,phi_vec,L);
  Eigen::VectorXd lower(2);
  lower << 0.000000, 0.000000;
  Eigen::VectorXd upper(2);
  upper << t_1, t_2;
  double err_est;
  int err_code;
  const int maxeval = 200000;
  const double eps_rel = 1e-13;
  const double res = integrate(f, lower, upper, err_est, err_code, maxeval, eps_rel);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}


// [[Rcpp::export]]
Rcpp::List joint_sub_dis_weibull_corr_cause_specific_cpp(const double t_1, const double t_2,
                                                const int j_1,const int j_2,const double beta_1,
                                                const double beta_2,const std::vector<double> mu_vec_1,const std::vector<double> mu_vec_2,
                                                const std::vector<double> tau_vec_new,const std::vector<double> phi_vec,const int L) {
  joint_sub_dis_weibull_corr_cause f(j_1 - 1,j_2 - 1,mu_vec_1,mu_vec_2,beta_1,beta_2,
                                     tau_vec_new,phi_vec,L);
  Eigen::VectorXd lower(2);
  lower << 0.000000, 0.000000;
  Eigen::VectorXd upper(2);
  upper << t_1, t_2;
  double err_est;
  int err_code;
  const int maxeval = 200000;
  const double eps_rel = 1e-13;
  const double res = integrate(f, lower, upper, err_est, err_code, maxeval, eps_rel);
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}
