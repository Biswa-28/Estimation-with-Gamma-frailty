lambda_vec_1 <- c(0.2,0.25)
lambda_vec_2 <- c(0.15,0.1)
sigma_1 <- 0.95
sigma_2 <- 0.85
rho <- 0.8
###------------------------------###
### Gamma correlated frailty model with
### exponential cause specific hazards
### joint sub-density function of correlated Gamma frailty
joint_sub_density_correlated_gamma <- function(
    t_1,t_2,j_1,j_2,mu_vec_1,mu_vec_2,tau_1,tau_2,phi)
{
  mu_1 <- sum(mu_vec_1)
  mu_2 <- sum(mu_vec_2)
  tau_sq_1 <- (tau_1)^2
  tau_sq_2 <- (tau_2)^2
  c_1 <- phi/(tau_1*tau_2)
  c_2 <- -c_1 + 1/tau_sq_1
  c_3 <- -c_1 + 1/tau_sq_2
  a_1 <- 1 + (tau_sq_1*mu_1*t_1) + (tau_sq_2*mu_2*t_2)
  a_2 <- 1 + (tau_sq_1*mu_1*t_1)
  a_3 <- 1 + (tau_sq_2*mu_2*t_2)
  b_1 <- (-c_1)*log(a_1)
  b_2 <- (-c_2)*log(a_2)
  b_3 <- (-c_3)*log(a_3)
  b <- exp(b_1 + b_2 + b_3)
  A_1 <- (c_1*(c_1 + 1))/(a_1*a_1)
  A_2 <- (c_1*c_3)/(a_1*a_3)
  A_3 <- (c_1*c_2)/(a_1*a_2)
  A_4 <- (c_2*c_3)/(a_2*a_3)
  f <- tau_sq_1*tau_sq_2*mu_vec_1[j_1]*mu_vec_2[j_2]*
    b*(A_1 + A_2 + A_3 + A_4)
  return(f)
}
### joint survival function of correlated Gamma frailty
joint_surv_func_correlated_gamma <- function(t_1,t_2,
                                             mu_vec_1,mu_vec_2,tau_1,tau_2,phi)
{
  mu_1 <- sum(mu_vec_1)
  mu_2 <- sum(mu_vec_2)
  c_1 <- phi/(tau_1*tau_2)
  c_2 <- -c_1 + 1/(tau_1^2) 
  c_3 <- -c_1 + 1/(tau_2^2)  
  b_1 <- 1 + ((tau_1^2)*mu_1*t_1) + ((tau_2^2)*mu_2*t_2)
  b_2 <- 1 + ((tau_1^2)*mu_1*t_1)
  b_3 <- 1 + ((tau_2^2)*mu_2*t_2)
  f_1 <- (-c_1)*log(b_1) 
  f_2 <- (-c_2)*log(b_2)
  f_3 <- (-c_3)*log(b_3)
  f <- exp(f_1 + f_2 + f_3)
  return(f)
}
### marginal sub-density function of correlated Gamma frailty
marg_subdensity_correlated_gamma <- function(t_1,j_1,
                                             mu_vec_1,tau_1)
{
  mu_1 <- sum(mu_vec_1)
  return(mu_vec_1[j_1]*((1 + (tau_1^2)*mu_1*t_1)^(-1 - 1/(tau_1^2))))
}
## one integrand appearing in denominator
integrand_1 <- function(t_1,t_2,j_2,mu_vec_1,mu_vec_2,tau_1,tau_2,phi)
{
  mu_1 <- sum(mu_vec_1)
  mu_2 <- sum(mu_vec_2)
  tau_sq_1 <- (tau_1)^2
  tau_sq_2 <- (tau_2)^2
  c_1 <- phi/(tau_1*tau_2)
  c_2 <- -c_1 + 1/tau_sq_1
  c_3 <- -c_1 + 1/tau_sq_2
  a_1 <- ((1 + tau_sq_1*mu_1*t_1 + tau_sq_2*mu_2*t_2))
  a_2 <- ((1 + tau_sq_1*mu_1*t_1))
  a_3 <- ((1 + tau_sq_2*mu_2*t_2))
  b_1 <- (-c_1)*log(a_1)
  b_2 <- (-c_2)*log(a_2)
  b_3 <- (-c_3)*log(a_3)
  b <- exp(b_1 + b_2 + b_3)
  A_1 <- (c_1*(c_1 + 1))/(a_1*a_1)
  A_2 <- (c_1*c_3)/(a_1*a_3)
  A_3 <- (c_1*c_2)/(a_1*a_2)
  A_4 <- (c_2*c_3)/(a_2*a_3)
  f <- tau_sq_1*tau_sq_2*mu_1*mu_vec_2[j_2]*
    b*(A_1 + A_2 + A_3 + A_4)
  return(f)
}
## another one integrand appearing in denominator
integrand_2 <- function(t_1,t_2,j_1,mu_vec_1,mu_vec_2,tau_1,tau_2,phi)
{
  mu_1 <- sum(mu_vec_1)
  mu_2 <- sum(mu_vec_2)
  tau_sq_1 <- (tau_1)^2
  tau_sq_2 <- (tau_2)^2
  c_1 <- phi/(tau_1*tau_2)
  c_2 <- -c_1 + 1/tau_sq_1
  c_3 <- -c_1 + 1/tau_sq_2
  a_1 <- ((1 + tau_sq_1*mu_1*t_1 + tau_sq_2*mu_2*t_2))
  a_2 <- ((1 + tau_sq_1*mu_1*t_1))
  a_3 <- ((1 + tau_sq_2*mu_2*t_2))
  b_1 <- (-c_1)*log(a_1)
  b_2 <- (-c_2)*log(a_2)
  b_3 <- (-c_3)*log(a_3)
  b <- exp(b_1 + b_2 + b_3)
  A_1 <- (c_1*(c_1 + 1))/(a_1*a_1)
  A_2 <- (c_1*c_3)/(a_1*a_3)
  A_3 <- (c_1*c_2)/(a_1*a_2)
  A_4 <- (c_2*c_3)/(a_2*a_3)
  f <- tau_sq_1*tau_sq_2*mu_2*mu_vec_1[j_1]*
    b*(A_1 + A_2 + A_3 + A_4)
  return(f)
}
### Cross hazard ratio function of correlated Gamma frailty
cr_function_exp_hazard_correlated_gamma <- function(
    x_1,x_2,i_1,i_2,mu_vec_1,mu_vec_2,tau_1,tau_2,phi)
{
  joint_sub_den <- joint_sub_density_correlated_gamma(x_1,x_2,i_1,i_2,
                                                      mu_vec_1,mu_vec_2,tau_1,tau_2,phi)
  joint_surv <- joint_surv_func_correlated_gamma(x_1,x_2,
                                                 mu_vec_1,mu_vec_2,tau_1,tau_2,phi)
  I_1 <- integrate(integrand_1,
                   lower=x_1,upper=Inf,t_2=x_2,j_2=i_2,mu_vec_1=mu_vec_1,
                   mu_vec_2=mu_vec_2,tau_1=tau_1,tau_2=tau_2,phi=phi)$value
  I_2 <- integrate(integrand_2,
                   lower=x_2,upper=Inf,t_1=x_1,j_1=i_1,mu_vec_1=mu_vec_1,
                   mu_vec_2=mu_vec_2,tau_1=tau_1,tau_2=tau_2,phi=phi)$value
  cross_ratio <- (joint_sub_den*joint_surv)/(I_1*I_2)
  return(cross_ratio)
}
mon_1 <- seq(0.01,100,by=0.1)
mon_2 <- c(0.05,2,10,50,100)
## for j_{1} = 1 and j_{2} = 1
cr_correlated_1_11 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=0.05,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_2_11 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=2,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_3_11 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=10,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_4_11 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=50,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_5_11 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=100,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
plot(mon_1,cr_correlated_1_11,type = "l",ylab = latex2exp::TeX(r'($CR_{1,1}(t_1,t_2;\xi,\theta)$)'),
     main = "correlated Gamma frailty",lty="solid",ylim = c(1,1.65),xlab = latex2exp::TeX(r'($t_1$)'))
lines(mon_1,cr_correlated_2_11,lty="dashed")
lines(mon_1,cr_correlated_3_11,lty="dotted")
lines(mon_1,cr_correlated_4_11,lty="dotdash")
lines(mon_1,cr_correlated_5_11,lty="longdash")
legend("topright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),lty=c("solid","dashed",
"dotted","dotdash","longdash"),cex=0.5)
## for j_{1} = 1 and j_{2} = 2
cr_correlated_1_12 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                                x_2=0.05,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_2_12 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                                x_2=2,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                          tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_3_12 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                                x_2=10,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                          tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_4_12 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                                x_2=50,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                          tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_5_12 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                                x_2=100,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                          tau_1=sigma_1,tau_2=sigma_2,phi=rho)
plot(mon_1,cr_correlated_1_12,type = "l",ylab = latex2exp::TeX(r'($CR_{1,2}(t_1,t_2;\xi,\theta)$)'),
     main = "correlated Gamma frailty",lty="solid",ylim = c(1,2),xlab = latex2exp::TeX(r'($t_1$)'))
lines(mon_1,cr_correlated_2_12,lty="dashed")
lines(mon_1,cr_correlated_3_12,lty="dotted")
lines(mon_1,cr_correlated_4_12,lty="dotdash")
lines(mon_1,cr_correlated_5_12,lty="longdash")
legend("topright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),lty=c("solid","dashed",
"dotted","dotdash","longdash"),cex=0.5)
## for j_{1} = 2 and j_{2} = 1
cr_correlated_1_21 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=0.05,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_2_21 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=2,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_3_21 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=10,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_4_21 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=50,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_5_21 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=100,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
plot(mon_1,cr_correlated_1_21,type = "l",ylab = latex2exp::TeX(r'($CR_{2,1}(t_1,t_2;\xi,\theta)$)'),
     main = "correlated Gamma frailty",lty="solid",ylim = c(1,2),xlab = latex2exp::TeX(r'($t_1$)'))
lines(mon_1,cr_correlated_2_21,lty="dashed")
lines(mon_1,cr_correlated_3_21,lty="dotted")
lines(mon_1,cr_correlated_4_21,lty="dotdash")
lines(mon_1,cr_correlated_5_21,lty="longdash")
legend("topright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
                            latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),lty=c("solid","dashed",
                                                                                                                                          "dotted","dotdash","longdash"),cex=0.5)
## for j_{1} = 2 and j_{2} = 2
cr_correlated_1_22 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=0.05,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_2_22 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=2,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_3_22 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=10,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_4_22 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=50,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
cr_correlated_5_22 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
                             x_2=100,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                             tau_1=sigma_1,tau_2=sigma_2,phi=rho)
plot(mon_1,cr_correlated_1,type = "l",ylab = latex2exp::TeX(r'($CR_{j_{2},j_{2}}(t_1,t_2;\xi,\theta)$)'),
     main = "correlated Gamma frailty",lty="solid",ylim = c(1,7),xlab = latex2exp::TeX(r'($t_1$)'), cex = 0.8)
lines(mon_1,cr_correlated_2,lty="dashed")
lines(mon_1,cr_correlated_3,lty="dotted")
lines(mon_1,cr_correlated_4,lty="dotdash")
lines(mon_1,cr_correlated_5,lty="longdash")
legend(80,4.6, legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 0.2$)'),
latex2exp::TeX(r'($t_{2} = 0.5$)'),latex2exp::TeX(r'($t_{2} = 0.9$)'),latex2exp::TeX(r'($t_{2} = 2$)')),lty=c("solid","dashed",
                "dotted","dotdash","longdash"),cex=0.5)
###------------------------------###
### Shared gamma cause specific frailty model with
### exponential cause specific hazards
### unconditional joint sub-density function for shared cause specific Gamma frailty
joint_sub_density_shared_cause_specific_gamma <- function(
    t_1,t_2,j_1,j_2,mu_vec_1,mu_vec_2,tau_vec,M)
{
  a <- sum((-1/(tau_vec^2))*log(replicate(M,1) + ((tau_vec^2)*((mu_vec_1*t_1) + (mu_vec_2*t_2)))))
  if(j_1 == j_2)
  {
    j <- j_1
    b <- -2*log(1 + (tau_vec[j]^2)*((mu_vec_1[j]*t_1) + (mu_vec_2[j]*t_2)))
    f <- (1 + (tau_vec[j]^2))*mu_vec_1[j]*mu_vec_2[j]*exp(a + b)
    return(f)
  }else
  {
    b <- -log(1 + (tau_vec[j_1]^2)*((mu_vec_1[j_1]*t_1) + (mu_vec_2[j_1]*t_2)))
    c <- -log(1 + (tau_vec[j_2]^2)*((mu_vec_1[j_2]*t_1) + (mu_vec_2[j_2]*t_2)))
    f <- mu_vec_1[j_1]*mu_vec_2[j_2]*exp(a + b + c)
    return(f)
  }
}
### unconditional joint survival function for shared cause specific Gamma frailty
joint_surv_shared_cause_specific_gamma <- function(t_1,t_2,mu_vec_1,
                                                   mu_vec_2,tau_vec,M)
{
  a <- sum((-1/(tau_vec)^2)*log(replicate(M,1) +
       (tau_vec^2)*((mu_vec_1*t_1) + (mu_vec_2*t_2))))
  return(exp(a))
}
## integrand_1_shared_cause for different cause of failure for two individual 
integrand_1_shared_cause_diff_cause <- function(t_1,t_2,j_2,mu_vec_1,
                                     mu_vec_2,tau_vec,M)
{
  a <- sum((-1/(tau_vec^2))*log(replicate(M,1) + (tau_vec^2)*((mu_vec_1*t_1) + (mu_vec_2*t_2))))
  b <- mu_vec_1*((replicate(M,1) + (tau_vec^2)*((mu_vec_1*t_1) + (mu_vec_2*t_2)))^(-1))
  c <- mu_vec_2[j_2]*(1 + (tau_vec[j_2]^2)*((mu_vec_1[j_2]*t_1) + (mu_vec_2[j_2]*t_2)))^(-1)
  f <- exp(a)*sum(b)*c
  return(f)
}
## integrand_2_shared_cause for different cause of failure for two individual 
integrand_2_shared_cause_diff_cause <- function(t_1,t_2,j_1,mu_vec_1,
                                     mu_vec_2,tau_vec,M)
{
  a <- sum((-1/(tau_vec^2))*log(replicate(M,1) + (tau_vec^2)*(mu_vec_1*t_1 + mu_vec_2*t_2)))
  b <- mu_vec_2*((replicate(M,1) + (tau_vec^2)*((mu_vec_1*t_1) + (mu_vec_2*t_2)))^(-1))
  c <- mu_vec_1[j_1]*(1 + (tau_vec[j_1]^2)*((mu_vec_1[j_1]*t_1) + (mu_vec_2[j_1]*t_2)))^(-1)
  f <- exp(a)*sum(b)*c
  return(f)
}
cr_function_exp_hazard_shared_cause_specific_gamma <- function(
    x_1,x_2,i_1,i_2,mu_vec_1,mu_vec_2,tau_vec,M)
{
  joint_surv <- joint_surv_shared_cause_specific_gamma(x_1,x_2,
                                mu_vec_1,mu_vec_2,tau_vec,M)
  if(i_1 == i_2)
  {
    i <- i_1
    joint_sub_den <- joint_sub_density_shared_cause_specific_gamma(x_1,x_2,i,i,
                                                                   mu_vec_1,mu_vec_2,tau_vec,M)
    integrand_1 <- function(u)
    {
      return(integrand_1_shared_cause_diff_cause(u,x_2,i,
                                      mu_vec_1,mu_vec_2,tau_vec,M))
    }
    integrand_2 <- function(v)
    {
      return(integrand_2_shared_cause_diff_cause(x_1,v,i,
                                      mu_vec_1,mu_vec_2,tau_vec,M))
    }
    I_1 <- integrate(Vectorize(integrand_1),
                     lower=x_1,upper=Inf)$value
    I_2 <- integrate(Vectorize(integrand_2),
                     lower=x_2,upper=Inf)$value
    cross_ratio <- (joint_sub_den*joint_surv)/(I_1*I_2)
    return(cross_ratio)
  }else
  {
    joint_sub_den  <- joint_sub_density_shared_cause_specific_gamma(x_1,x_2,i_1,i_2,
                            mu_vec_1,mu_vec_2,tau_vec,M)
    integrand_1 <- function(u)
    {
      return(integrand_1_shared_cause_diff_cause(u,x_2,i_2,
                                      mu_vec_1,mu_vec_2,tau_vec,M))
    }
    integrand_2 <- function(v)
    {
      return(integrand_2_shared_cause_diff_cause(x_1,v,i_1,
                                      mu_vec_1,mu_vec_2,tau_vec,M))
    }
    I_1 <- integrate(Vectorize(integrand_1),
                     lower=x_1,upper=Inf)$value
    I_2 <- integrate(Vectorize(integrand_2),
                     lower=x_2,upper=Inf)$value
    cross_ratio <- (joint_sub_den*joint_surv)/(I_1*I_2)
    return(cross_ratio)
  }
}
mon_1 <- seq(0.1,100,by=0.1)
mon_2 <- c(0.05,0.2,0.5,0.9,2)
lambda_vec_1 <- c(0.2,0.25)
lambda_vec_2 <- c(0.15,0.1)
sigma_vec <- c(1.18,0.671)             #1.28
## for j_{1} = 1 and j_{2} = 1
cr_shared_cause_1_11 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
x_2=0.05,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                          tau_vec=sigma_vec,M=2)
cr_shared_cause_2_11 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                            x_2=2,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                            tau_vec=sigma_vec,M=2)
cr_shared_cause_3_11 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                            x_2=10,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                            tau_vec=sigma_vec,M=2)
cr_shared_cause_4_11 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                            x_2=50,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                            tau_vec=sigma_vec,M=2)
cr_shared_cause_5_11 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                            x_2=100,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                            tau_vec=sigma_vec,M=2)
plot(mon_1,cr_shared_cause_1_11,type = "l", xlab = latex2exp::TeX(r'($t_1$)'),ylim = c(4,6.5),
main = latex2exp::TeX(r'($j_{1} = 1$ and $j_{2} = 1$)'),lty="solid",
ylab = latex2exp::TeX(r'($CR_{1,1}(t_1,t_2;\xi,\theta)$)'))
lines(mon_1,cr_shared_cause_2_11,lty="dashed")
lines(mon_1,cr_shared_cause_3_11,lty="dotted")
lines(mon_1,cr_shared_cause_4_11,lty="dotdash")
lines(mon_1,cr_shared_cause_5_11,lty="longdash")
legend("topright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),
lty=c("solid","dashed","dotted","dotdash","longdash"),cex=0.8)
## for j_{1} = 1 and j_{2} = 2
cr_shared_cause_1_12 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=0.05,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_2_12 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=2,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_3_12 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=10,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_4_12 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=50,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_5_12 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=100,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
plot(mon_1,cr_shared_cause_1_12,type = "l", xlab = latex2exp::TeX(r'($t_1$)'),ylim = c(1.85,2.2),
     main = latex2exp::TeX(r'($j_{1} = 1$ and $j_{2} = 2$)'),lty="solid",
     ylab = latex2exp::TeX(r'($CR_{1,2}(t_1,t_2;\xi,\theta)$)'))
lines(mon_1,cr_shared_cause_2_12,lty="dashed")
lines(mon_1,cr_shared_cause_3_12,lty="dotted")
lines(mon_1,cr_shared_cause_4_12,lty="dotdash")
lines(mon_1,cr_shared_cause_5_12,lty="longdash")
legend("topright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),
lty=c("solid","dashed","dotted","dotdash","longdash"),cex=0.8)
## for j_{1} = 2 and j_{2} = 1
cr_shared_cause_1_21 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=0.05,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_2_21 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=2,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_3_21 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=10,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_4_21 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=50,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_5_21 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=100,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
plot(mon_1,cr_shared_cause_1_21,type = "l", xlab = latex2exp::TeX(r'($t_1$)'),ylim = c(1.65,1.9),
     main = latex2exp::TeX(r'($j_{1} = 2$ and $j_{2} = 1$)'),lty="solid",
     ylab = latex2exp::TeX(r'($CR_{2,1}(t_1,t_2;\xi,\theta)$)'))
lines(mon_1,cr_shared_cause_2_21,lty="dashed")
lines(mon_1,cr_shared_cause_3_21,lty="dotted")
lines(mon_1,cr_shared_cause_4_21,lty="dotdash")
lines(mon_1,cr_shared_cause_5_21,lty="longdash")
legend("topright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),
       lty=c("solid","dashed","dotted","dotdash","longdash"),cex=0.8)
## for j_{1} = 2 and j_{2} = 2
cr_shared_cause_1_22 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=0.05,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_2_22 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=2,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_3_22 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=10,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_4_22 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=50,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
cr_shared_cause_5_22 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
                               x_2=100,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                               tau_vec=sigma_vec,M=2)
plot(mon_1,cr_shared_cause_1_22,type = "l", xlab = latex2exp::TeX(r'($t_1$)'),ylim = c(2.2,2.7),
     main = latex2exp::TeX(r'($j_{1} = 2$ and $j_{2} = 2$)'),lty="solid",
     ylab = latex2exp::TeX(r'($CR_{2,2}(t_1,t_2;\xi,\theta)$)'))
lines(mon_1,cr_shared_cause_2_22,lty="dashed")
lines(mon_1,cr_shared_cause_3_22,lty="dotted")
lines(mon_1,cr_shared_cause_4_22,lty="dotdash")
lines(mon_1,cr_shared_cause_5_22,lty="longdash")
legend("bottomright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),
       lty=c("solid","dashed","dotted","dotdash","longdash"),cex=0.8)
###------------------------------###
### Correlated gamma cause specific frailty model with
### exponential cause specific hazards
sigma_vec_1 <- c(1.2,1.8)
sigma_vec_2 <- c(0.8,0.4)
rho_vec <- c(0.7,0.25)
joint_sub_density_correlated_cause_specific_gamma <- function(
    t_1,t_2,j_1,j_2,mu_vec_1,mu_vec_2,tau_vec_1,tau_vec_2,phi_vec)
{
  c_vec_1 <- phi_vec/(tau_vec_1*tau_vec_2)
  c_vec_2 <- -c_vec_1 + (tau_vec_1)^(-2)
  c_vec_3 <- -c_vec_1 + (tau_vec_2)^(-2)
  a_1 <- ((1 + (tau_vec_1^2)*mu_vec_1*t_1 + (tau_vec_2^2)*mu_vec_2*t_2))
  a_2 <- ((1 + (tau_vec_1^2)*mu_vec_1*t_1))
  a_3 <- ((1 + (tau_vec_2^2)*mu_vec_2*t_2))
  b_1 <- -c_vec_1*log(a_1)
  b_2 <- -c_vec_2*log(a_2)
  b_3 <- -c_vec_3*log(a_3)
  joint_surv <- exp(sum(b_1 + b_2 + b_3))
  if(j_1 == j_2)
  {
    j <- j_1
    A_1 <- (c_vec_1[j]*(1 + c_vec_1[j]))/(a_1[j]^2)
    A_2 <- (c_vec_1[j]*c_vec_2[j])/(a_1[j]*a_2[j])
    A_3 <- (c_vec_1[j]*c_vec_3[j])/(a_1[j]*a_3[j])
    A_4 <- (c_vec_2[j]*c_vec_3[j])/(a_2[j]*a_3[j])
    f <- (tau_vec_1[j]^2)*(tau_vec_2[j]^2)*mu_vec_1[j]*mu_vec_2[j]*
      joint_surv*(A_1 + A_2 + A_3 + A_4)
    return(f)
  }
  else
  {
    A_1 <- (c_vec_1[j_1]*c_vec_1[j_2])/(a_1[j_1]*a_1[j_2])
    A_2 <- (c_vec_1[j_1]*c_vec_3[j_2])/(a_1[j_1]*a_3[j_2])
    A_3 <- (c_vec_1[j_2]*c_vec_2[j_1])/(a_1[j_2]*a_2[j_1])
    A_4 <- (c_vec_2[j_1]*c_vec_3[j_2])/(a_2[j_1]*a_3[j_2])
    f <- (tau_vec_1[j_1]^2)*(tau_vec_2[j_2]^2)*mu_vec_1[j_1]*
      mu_vec_2[j_2]*joint_surv*(A_1 + A_2 + A_3 + A_4)
    return(f)
  }
  
}
joint_surv_corr_cause_specific_gamma <- function(t_1,t_2,mu_vec_1,
                                                 mu_vec_2,tau_vec_1,tau_vec_2,phi_vec)
{
  a_1 <- phi_vec/(tau_vec_1*tau_vec_2)
  a_2 <- (1/tau_vec_1)^2 - a_1
  a_3 <- (1/tau_vec_2)^2 - a_1
  b_1 <- 1 + ((tau_vec_1^2)*mu_vec_1*t_1) + ((tau_vec_2^2)*mu_vec_2*t_2)
  b_2 <- 1 + ((tau_vec_1^2)*mu_vec_1*t_1)
  b_3 <- 1 + ((tau_vec_2^2)*mu_vec_2*t_2)
  c_1 <- (-a_1)*log(b_1)
  c_2 <- (-a_2)*log(b_2)
  c_3 <- (-a_3)*log(b_3)
  f <- exp(sum(c_1 + c_2 + c_3))
  return(f)
}
integrand_1_correlated_cause <- function(
    t_1,t_2,j_2,mu_vec_1,mu_vec_2,tau_vec_1,tau_vec_2,phi_vec)
{
  c_vec_1 <- phi_vec/(tau_vec_1*tau_vec_2)
  c_vec_2 <- -c_vec_1 + (tau_vec_1)^(-2)
  c_vec_3 <- -c_vec_1 + (tau_vec_2)^(-2)
  a_1 <- ((1 + (tau_vec_1^2)*mu_vec_1*t_1 + (tau_vec_2^2)*mu_vec_2*t_2))
  a_2 <- ((1 + (tau_vec_1^2)*mu_vec_1*t_1))
  a_3 <- ((1 + (tau_vec_2^2)*mu_vec_2*t_2))
  c_1 <- (-c_vec_1)*log(a_1)
  c_2 <- (-c_vec_2)*log(a_2)
  c_3 <- (-c_vec_3)*log(a_3)
  joint_surv <- exp(sum(c_1 + c_2 + c_3))
  d <- (tau_vec_1^2)*mu_vec_1
  A_1 <- ((c_vec_1[j_2])/(a_1[j_2]))*sum(d*(c_vec_1/a_1))
  A_2 <- ((c_vec_3[j_2])/(a_3[j_2]))*sum(d*(c_vec_1/a_1))
  A_3 <- ((c_vec_1[j_2])/(a_1[j_2]))*sum(d*(c_vec_2/a_2))
  A_4 <- ((c_vec_3[j_2])/(a_3[j_2]))*sum(d*(c_vec_2/a_2))
  f <- (tau_vec_2[j_2]^2)*mu_vec_2[j_2]*joint_surv*(A_1 + A_2 + A_3 + A_4)
  return(f)
}
integrand_2_correlated_cause <- function(
    t_1,t_2,j_1,mu_vec_1,mu_vec_2,tau_vec_1,tau_vec_2,phi_vec)
{
  c_vec_1 <- phi_vec/(tau_vec_1*tau_vec_2)
  c_vec_2 <- -c_vec_1 + (tau_vec_1)^(-2)
  c_vec_3 <- -c_vec_1 + (tau_vec_2)^(-2)
  a_1 <- ((1 + (tau_vec_1^2)*mu_vec_1*t_1 + (tau_vec_2^2)*mu_vec_2*t_2))
  a_2 <- ((1 + (tau_vec_1^2)*mu_vec_1*t_1))
  a_3 <- ((1 + (tau_vec_2^2)*mu_vec_2*t_2))
  c_1 <- (-c_vec_1)*log(a_1)
  c_2 <- (-c_vec_2)*log(a_2)
  c_3 <- (-c_vec_3)*log(a_3)
  joint_surv <- exp(sum(c_1 + c_2 + c_3))
  d <- (tau_vec_2^2)*mu_vec_2
  A_1 <- ((c_vec_1[j_1])/(a_1[j_1]))*sum(d*(c_vec_1/a_1))
  A_2 <- ((c_vec_2[j_1])/(a_2[j_1]))*sum(d*(c_vec_1/a_1))
  A_3 <- ((c_vec_1[j_1])/(a_1[j_1]))*sum(d*(c_vec_3/a_3))
  A_4 <- ((c_vec_2[j_1])/(a_2[j_1]))*sum(d*(c_vec_3/a_3))
  f <- (tau_vec_1[j_1]^2)*mu_vec_1[j_1]*joint_surv*(A_1 + A_2 + A_3 + A_4)
  return(f)
}
cr_function_exp_hazard_correlated_cause_specific_gamma <- function(
    x_1,x_2,i_1,i_2,mu_vec_1,mu_vec_2,tau_vec_1,tau_vec_2,phi_vec)
{
  joint_surv <- joint_surv_corr_cause_specific_gamma(x_1,x_2,mu_vec_1,
                                  mu_vec_2,tau_vec_1,tau_vec_2,phi_vec)
  if(i_1 == i_2)
  {
    i <- i_1
    joint_sub_den <- joint_sub_density_correlated_cause_specific_gamma(x_1,x_2,i,i,mu_vec_1,mu_vec_2,tau_vec_1,tau_vec_2,phi_vec)
    integrand_1 <- function(u)
    {
      return(integrand_1_correlated_cause(u,x_2,i,
                                          mu_vec_1,mu_vec_2,tau_vec_1,tau_vec_2,phi_vec))
    }
    integrand_2 <- function(v)
    {
      return(integrand_2_correlated_cause(x_1,v,i,
                                          mu_vec_1,mu_vec_2,tau_vec_1,tau_vec_2,phi_vec))
    }
    I_1 <- integrate(Vectorize(integrand_1),
                     lower=x_1,upper=Inf)$value
    I_2 <- integrate(Vectorize(integrand_2),
                     lower=x_2,upper=Inf)$value
    cross_ratio <- (joint_sub_den*joint_surv)/(I_1*I_2)
    return(cross_ratio)
  }
  else
  {
    joint_sub_den <- joint_sub_density_correlated_cause_specific_gamma(x_1,x_2,i_1,i_2,mu_vec_1,mu_vec_2,tau_vec_1,tau_vec_2,phi_vec)
    
    integrand_1 <- function(u)
    {
      return(integrand_1_correlated_cause(u,x_2,i_2,
                                          mu_vec_1,mu_vec_2,tau_vec_1,tau_vec_2,phi_vec))
    }
    integrand_2 <- function(v)
    {
      return(integrand_2_correlated_cause(x_1,v,i_1,
                                          mu_vec_1,mu_vec_2,tau_vec_1,tau_vec_2,phi_vec))
    }
    I_1 <- integrate(Vectorize(integrand_1),
                     lower=x_1,upper=Inf)$value
    I_2 <- integrate(Vectorize(integrand_2),
                     lower=x_2,upper=Inf)$value
    cross_ratio <- (joint_sub_den*joint_surv)/(I_1*I_2)
    return(cross_ratio)
  }
}
mon_1 <- seq(0.1,100,by=0.1)
mon_2 <- c(0.05,0.2,0.5,0.9,2)
### parameter vectors
lambda_vec_1 <- c(0.2,0.25)
lambda_vec_2 <- c(0.15,0.1)
sigma_vec_1 <- c(1.2,1.8)
sigma_vec_2 <- c(0.8,0.4)
rho_vec <- c(0.7,0.25)
## for j_{1} = 1 and j_{2} = 1
cr_correlated_cause_1_11 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=0.05,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_2_11 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=2,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_3_11 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=10,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_4_11 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=50,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_5_11 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=100,i_1=1,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
plot(mon_1,cr_correlated_cause_1_11,type = "l",xlab=latex2exp::TeX(r'($t_1$)'),
     ylab = latex2exp::TeX(r'($CR_{1,1}(t_1,t_2;\xi,\theta)$)'),
     main = latex2exp::TeX(r'($j_{1} = 1$ and $j_{2} = 1$)'),lty="solid",ylim = c(1,4))
lines(mon_1,cr_correlated_cause_2_11,lty="dashed")
lines(mon_1,cr_correlated_cause_3_11,lty="dotted")
lines(mon_1,cr_correlated_cause_4_11,lty="dotdash")
lines(mon_1,cr_correlated_cause_5_11,lty="longdash")
legend("topright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),
lty=c("solid","dashed","dotted","dotdash","longdash"),cex=0.8)
## for j_{1} = 1 and j_{2} = 2
cr_correlated_cause_1_12 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
x_2=0.05,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_2_12 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
x_2=2,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_3_12 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
x_2=10,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_4_12 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
x_2=50,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_5_12 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
x_2=100,i_1=1,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
plot(mon_1,cr_correlated_cause_1_12,type = "l",xlab=latex2exp::TeX(r'($t_1$)'),
     ylab = latex2exp::TeX(r'($CR_{1,2}(t_1,t_2;\xi,\theta)$)'),
     main = latex2exp::TeX(r'($j_{1} = 1$ and $j_{2} = 2$)'),lty="solid",ylim = c(1,1.5))
lines(mon_1,cr_correlated_cause_2_12,lty="dashed")
lines(mon_1,cr_correlated_cause_3_12,lty="dotted")
lines(mon_1,cr_correlated_cause_4_12,lty="dotdash")
lines(mon_1,cr_correlated_cause_5_12,lty="longdash")
legend("topright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
                            latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),
       lty=c("solid","dashed","dotted","dotdash","longdash"),cex=0.8)
## for j_{1} = 2 and j_{2} = 1
cr_correlated_cause_1_21 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=0.05,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_2_21 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=2,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_3_21 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=10,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_4_21 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=50,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_5_21 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=100,i_1=2,i_2=1,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
plot(mon_1,cr_correlated_cause_1_21,type = "l",xlab=latex2exp::TeX(r'($t_1$)'),
     ylab = latex2exp::TeX(r'($CR_{2,1}(t_1,t_2;\xi,\theta)$)'),
     main = latex2exp::TeX(r'($j_{1} = 2$ and $j_{2} = 1$)'),lty="solid",ylim = c(1,1.4))
lines(mon_1,cr_correlated_cause_2_21,lty="dashed")
lines(mon_1,cr_correlated_cause_3_21,lty="dotted")
lines(mon_1,cr_correlated_cause_4_21,lty="dotdash")
lines(mon_1,cr_correlated_cause_5_21,lty="longdash")
legend("topright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
                            latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),
       lty=c("solid","dashed","dotted","dotdash","longdash"),cex=0.8)
## for j_{1} = 2 and j_{2} = 2
cr_correlated_cause_1_22 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=0.05,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_2_22 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=2,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_3_22 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=10,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_4_22 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=50,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
cr_correlated_cause_5_22 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
                                   x_2=100,i_1=2,i_2=2,mu_vec_1=lambda_vec_1,mu_vec_2=lambda_vec_2,
                                   tau_vec_1=sigma_vec_1,tau_vec_2=sigma_vec_2,phi_vec=rho_vec)
plot(mon_1,cr_correlated_cause_1_22,type = "l",xlab=latex2exp::TeX(r'($t_1$)'),
     ylab = latex2exp::TeX(r'($CR_{2,2}(t_1,t_2;\xi,\theta)$)'),
     main = latex2exp::TeX(r'($j_{1} = 2$ and $j_{2} = 2$)'),lty="solid",ylim = c(1,1.5))
lines(mon_1,cr_correlated_cause_2_22,lty="dashed")
lines(mon_1,cr_correlated_cause_3_22,lty="dotted")
lines(mon_1,cr_correlated_cause_4_22,lty="dotdash")
lines(mon_1,cr_correlated_cause_5_22,lty="longdash")
legend("topright", legend=c(latex2exp::TeX(r'($t_{2} = 0.05$)'),latex2exp::TeX(r'($t_{2} = 2$)'),
                            latex2exp::TeX(r'($t_{2} = 10$)'),latex2exp::TeX(r'($t_{2} = 50$)'),latex2exp::TeX(r'($t_{2} = 100$)')),
       lty=c("solid","dashed","dotted","dotdash","longdash"),cex=0.8)
################## Association in different fitted models in dataset #######################
## call the data in r environment
data <- read.csv("data.csv")
## checking the structure and missingness of the data
str(data)
d <- attach(data)
Age <- Age[-c(380)]
## transforming the monitoring time in months
chng_year_to_mths <- function(x)
{
  chac_vec <- unlist(strsplit(x,split = ""))
  chac_vec <- chac_vec[which(chac_vec != " ")]
  pos_y <- which((chac_vec == "y")|(chac_vec == "Y"))
  pos_r <- which(chac_vec == "r")
  pos_s <- which(chac_vec == "s")
  pos_m <- which(chac_vec == "m")
  l_y <- length(pos_y)
  l_r <- length(pos_r)
  if(length(pos_s) == 0)
  {
    l_s <- 0
  }else
  {
    l_s <- length(pos_s[1])
  }
  l_m <- length(pos_m)
  if((l_y == 0)&&(l_r == 0)&&(l_m == 1))
  {
    return(as.numeric(paste(chac_vec[1:(pos_m - 1)],collapse = "")))
  }
  if((l_y == 1)&&(l_r == 1)&&(l_m == 0))
  {
    return(12*as.numeric(paste(chac_vec[1:(pos_y - 1)],collapse = "")))
  }
  if((l_y == 1)&&(l_r == 1)&&(l_s == 0)&&(l_m == 1))
  {
    t_1 <- (12*as.numeric(paste(chac_vec[1:(pos_y - 1)],collapse = "")))
    t_2 <-  as.numeric(paste(chac_vec[(pos_r + 1):(pos_m - 1)],collapse = ""))
    return(t_1 + t_2)
  }
  if((l_y == 1)&&(l_r == 1)&&(pos_s[1] > pos_m)&&(l_m == 1))
  {
    t_1 <- (12*as.numeric(paste(chac_vec[1:(pos_y - 1)],collapse = "")))
    t_2 <-  as.numeric(paste(chac_vec[(pos_r + 1):(pos_m - 1)],collapse = ""))
    return(t_1 + t_2)
  }
  if((l_y == 1)&&(l_r == 1)&&(pos_s[1] < pos_m)&&(l_m == 1))
  {
    t_1 <- (12*as.numeric(paste(chac_vec[1:(pos_y - 1)],collapse = "")))
    t_2 <-  as.numeric(paste(chac_vec[(pos_s[1] + 1):(pos_m - 1)],collapse = ""))
    return(t_1 + t_2)
  }
}
x <- unlist(sapply(Age,chng_year_to_mths))
Age <- c(x[1:379],12.5,x[380:879])
data_1 <- cbind(Age,data[,c(3,4)])
delete_rows <- as.numeric(rownames(subset(data_1,Loss == "Y" & Type == "")))
data_1 <- data_1[-c(delete_rows),]
failure_cause_gen <- function(x)
{
  loss <- x[1]
  type <- x[2]
  if((loss == "N") & (type == ""))     ## 1
  {
    return(c(0,0))
  }
  if((loss == "Y") & (type == "SNHL"))  ## 2
  {
    return(c(1,1))
  }
  if((loss == "Y") & (type == "Conductive"))  ## 3
  {
    return(c(2,2))
  }
  if((loss == "Y") & (type == "Mixed"))  ## 4
  {
    return(c(3,3))
  }
  if((loss == "Y") & (type == "R-Mixed,L-SNHL"))  ##5
  {
    return(c(1,3))
  }
  if((loss == "Y") & (type == "R-SNHL,L-Mixed"))  ##6
  {
    return(c(3,1))
  }
  if((loss == "Y") & (type == "conductive"))  ## 7
  {
    return(c(2,2))
  }
  if((loss == "Y") & (type == "R-Mixed,L-Conductive"))  ##8
  {
    return(c(2,3))
  }
  if((loss == "Y") & (type == "R-Conductive, L-Mixed"))  ##9
  {
    return(c(3,2))
  }
  if((loss == "Y") & (type == "R-Mixed, L-SNHL")) ## 10 
  {
    return(c(1,3))
  }
  if((loss == "Y") & (type == "R-Conductive "))  ## 11
  {
    return(c(0,2))
  }
  if((loss == "Y") & (type == "L-Conductive, R-Mixed"))  ## 12
  {
    return(c(2,3))
  }  
  if((loss == "Y") & (type == "L-Conductive "))  ## 13
  {
    return(c(2,0))
  }  
}
data_failure_cause <- t(apply(data_1[,c(2,3)],1,failure_cause_gen))
hearing_data <- cbind(data_1[,1],data_failure_cause)
############################_______________________________##################################
########### Exponential cause specific hazard with correlated Gamma frailty model ##########
lambda_vec_1_mle <- c(0.6652,0.0431,0.0921)
lambda_vec_2_mle <- c(0.6392,0.0442,0.0931)
sigma_1_mle <- 2.3903
sigma_2_mle <- 2.4003
rho_mle <- 0.9902
mon_1 <- unique(sort(hearing_data[,1]))[1:100]
mon_2 <- c(4,53,106,252,444)
cr_correlated_fitted_1 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
x_2=mon_2[1],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_1=sigma_1_mle,tau_2=sigma_2_mle,phi=rho_mle)
cr_correlated_fitted_2 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
x_2=mon_2[2],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_1=sigma_1_mle,tau_2=sigma_2_mle,phi=rho_mle)
cr_correlated_fitted_3 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
x_2=mon_2[3],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_1=sigma_1_mle,tau_2=sigma_2_mle,phi=rho_mle)
cr_correlated_fitted_4 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
x_2=mon_2[4],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_1=sigma_1_mle,tau_2=sigma_2_mle,phi=rho_mle)
cr_correlated_fitted_5 <- sapply(mon_1,cr_function_exp_hazard_correlated_gamma,
x_2=mon_2[5],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_1=sigma_1_mle,tau_2=sigma_2_mle,phi=rho_mle)
plot(mon_1,cr_correlated_fitted_1,type = "l",ylab="Cross ratio function",
     main = "Correlated Gamma frailty",lty="solid",ylim = c(3,8),xlab="Monitring time")
lines(mon_1,cr_correlated_fitted_2,lty="dashed")
lines(mon_1,cr_correlated_fitted_3,lty="dotted")
lines(mon_1,cr_correlated_fitted_4,lty="dotdash")
lines(mon_1,cr_correlated_fitted_5,lty="longdash")
##############################_____________________________#####################################
###### Exponential cause specific hazard with shared cause specific Gamma frailty ###### 
lambda_vec_1_mle <- c(35.95820,0.00605,0.00338)
lambda_vec_2_mle <- c(39.21070,0.00550,0.00345)
sigma_vec_mle <- c(12.92450,45.58919,7.65704)
mon_1 <- unique(sort(hearing_data[,1]))[1:100]
mon_2 <- c(4,53,106,252,444)
cr_shared_cause_fitted_1 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
x_2=mon_2[1],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_vec=sigma_vec_mle,M=3)
cr_shared_cause_fitted_2 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
x_2=mon_2[2],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_vec=sigma_vec_mle,M=3)
cr_shared_cause_fitted_3 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
x_2=mon_2[3],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_vec=sigma_vec_mle,M=3)
cr_shared_cause_fitted_4 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
x_2=mon_2[4],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_vec=sigma_vec_mle,M=3)
cr_shared_cause_fitted_5 <- sapply(mon_1,cr_function_exp_hazard_shared_cause_specific_gamma,
x_2=mon_2[5],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_vec=sigma_vec_mle,M=3)
plot(mon_1,cr_shared_cause_fitted_1,type = "l",xlab = "Monitoring time",ylim = c(30,60),
     main = "Shared cause specific Gamma",lty="solid",ylab="Cross ratio function")
lines(mon_1,cr_shared_cause_fitted_2,lty="dashed")
lines(mon_1,cr_shared_cause_fitted_3,lty="dotted")
lines(mon_1,cr_shared_cause_fitted_4,lty="dotdash")
lines(mon_1,cr_shared_cause_fitted_5,lty="longdash")
##############################_____________________________#####################################
###### Exponential cause specific hazard with correlated cause specific Gamma frailty ###### 
lambda_vec_1_mle <- c(2485.23120,0.02021,0.01821)
lambda_vec_2_mle <- c(2645.33701,0.01867,0.01881)
sigma_vec_mle <- c(4.65674,5.07337,1.99483)
rho_vec_mle <- c(0.999,0.995,0.991)
mon_1 <- unique(sort(hearing_data[,1]))[1:100]
mon_2 <- c(4,53,106,252,444)
cr_correlated_cause_fitted_1 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
x_2=mon_2[1],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_vec_1=sigma_vec_mle,tau_vec_2=sigma_vec_mle,phi_vec=rho_vec_mle)
cr_correlated_cause_fitted_2 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
x_2=mon_2[2],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_vec_1=sigma_vec_mle,tau_vec_2=sigma_vec_mle,phi_vec=rho_vec_mle)
cr_correlated_cause_fitted_3 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
x_2=mon_2[3],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_vec_1=sigma_vec_mle,tau_vec_2=sigma_vec_mle,phi_vec=rho_vec_mle)
cr_correlated_cause_fitted_4 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
x_2=mon_2[4],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_vec_1=sigma_vec_mle,tau_vec_2=sigma_vec_mle,phi_vec=rho_vec_mle)
cr_correlated_cause_fitted_5 <- sapply(mon_1,cr_function_exp_hazard_correlated_cause_specific_gamma,
x_2=mon_2[5],i_1=1,i_2=2,mu_vec_1=lambda_vec_1_mle,mu_vec_2=lambda_vec_2_mle,
tau_vec_1=sigma_vec_mle,tau_vec_2=sigma_vec_mle,phi_vec=rho_vec_mle)
plot(mon_1,cr_correlated_cause_fitted_1,type = "l", xlab = "x1",ylab = "Cross Ratio (x1,x2)",
main = "Correlated cause specific Gamma",lty="solid")
lines(mon_1,cr_correlated_cause_fitted_2,lty="dashed")
lines(mon_1,cr_correlated_cause_fitted_3,lty="dotted")
lines(mon_1,cr_correlated_cause_fitted_4,lty="dotdash")
lines(mon_1,cr_correlated_cause_fitted_5,lty="longdash")

