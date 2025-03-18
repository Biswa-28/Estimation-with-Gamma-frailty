lambda_vec_1 <- c(3,2)      # c(3,2)
lambda_vec_2 <- c(2.5,0.8)  # c(2.5,0.8)
sigma_1 <- 0.9           # 0.9
sigma_2 <- 0.7           # 0.7
rho <- 0.65 #0 <= rho <= min(sigma_1/sigma_2,sigma_2/sigma_1) = 0.75  # 0.65
lambda_1 <- sum(lambda_vec_1)
lambda_2 <- sum(lambda_vec_2)
p_cen <- 0.8
integrand_1 <- function(x,mu,mu_1,mu_2,tau_1,tau_2,phi)
{
  d_1 <- mu_1*(tau_1^2)*x
  d_2 <- mu_2*(tau_2^2)*x
  b_1 <- 1 + d_1 + d_2
  b_2 <- 1 + d_1
  b_3 <- 1 + d_2
  c_1 <- phi/(tau_1*tau_2)
  c_2 <- c_1 - (1/(tau_1^2))
  c_3 <- c_1 - (1/(tau_2^2))
  return(dexp(x,rate = mu)*(b_1^(-c_1))*(b_2^c_2)*(b_3^c_3))
}
integrand_2 <- function(x,mu,mu_1,tau_1)
{
  a_1 <- 1 + (mu_1*(tau_1^2)*x)
  return(dexp(x,rate = mu)*(a_1^(-1/(tau_1^2))))
}
## equation to solve for mu for a given censoring prob
f_1 <- function(mu)
{
  return((integrate(integrand_1,0,Inf,mu=mu,mu_1=lambda_1,
                    mu_2=lambda_2,tau_1=sigma_1,tau_2=sigma_2,phi=rho)$value - p_cen))
}
## solution of mu_hat
mu_hat <- nleqslv::nleqslv(0.004,f_1)$x
L_1 <- length(lambda_vec_1)
L_2 <- length(lambda_vec_2)
p_1 <- (lambda_vec_1)/sum(lambda_vec_1)
p_2 <- (lambda_vec_2)/sum(lambda_vec_2)
p_3 <- outer(p_1,p_2)
## function to generate failure causes for two individuals
failure_cause_gen <- function(x,mu_1,mu_2,tau_1,tau_2,phi,M_1,M_2,q_1,q_2,q_3)
{
  d_1 <- mu_1*(tau_1^2)*x
  d_2 <- mu_2*(tau_2^2)*x
  b_1 <- 1 + d_1 + d_2
  b_2 <- 1 + d_1
  b_3 <- 1 + d_2
  c_1 <- phi/(tau_1*tau_2)
  c_2 <- c_1 - (1/tau_1^2)
  c_3 <- c_1 - (1/tau_2^2)
  log_marg_surv_1 <- (-(1/tau_1^2))*log(b_2)
  log_marg_surv_2 <- (-(1/tau_2^2))*log(b_3)
  log_joint_surv <-  ((-c_1))*log(b_1) + (c_2*log(b_2)) + (c_3*log(b_3))
  joint_surv <- exp(log_joint_surv)
  marg_surv_1 <- exp(log_marg_surv_1)
  marg_surv_2 <- exp(log_marg_surv_2)
  joint_failure_prob <- 1 - marg_surv_1 - marg_surv_2 + joint_surv
  prob_1 <- q_3*joint_failure_prob
  prob_1 <- rbind(q_2*(marg_surv_1 - joint_surv),prob_1)
  prob_1 <- cbind(c(joint_surv,q_1*(marg_surv_2 - joint_surv)),
                  prob_1)
  s <- which(rmultinom(1,1,prob_1)!=0)
  s_1 <- floor(s/(M_1 + 1))
  s_2 <- s - ((M_1 + 1)*s_1)
  if(s_2 != 0)
  {
    return(c(s_2 - 1,s_1))
  }else
  {
    return(c(M_1,s_1 - 1))
  }
}
## function to generate the data for a given sample size n
corr_data_gen <- function(n,mu_1,mu_2,tau_1,tau_2,phi,M_1,M_2,q_4,q_5,q_6)
{
  mon_time <- rexp(n,rate=mu_hat)
  F <- sapply(mon_time,failure_cause_gen,tau_1=tau_1,tau_2=tau_2,mu_1=mu_1,
              mu_2=mu_2,phi=phi,M_1=M_1,M_2=M_2,q_1=q_4,q_2=q_5,q_3=q_6)
  return(cbind(mon_time,t(F)))
}
##Calculation of true cell probabilities
surv_1 <- integrate(integrand_2,0,Inf,mu=mu_hat,
                    mu_1=lambda_1,tau_1=sigma_1)$value
surv_2 <- integrate(integrand_2,0,Inf,mu=mu_hat,
                    mu_1=lambda_2,tau_1=sigma_2)$value
surv_joint <- integrate(integrand_1,0,Inf,mu=mu_hat,
          mu_1=lambda_1,mu_2=lambda_2,tau_1=sigma_1,tau_2=sigma_2,phi=rho)$value
subdist_func_1 <- (p_1)*(1 - surv_1)
subdist_func_2 <- (p_2)*(1 - surv_2)
c <- (1 - (surv_1) - (surv_2) + surv_joint)
joint_subdist_func <- (p_3)*c
P_cell_1 <- rbind(p_2*(surv_1 - p_cen),joint_subdist_func)
P_cell_1 <- cbind(c(p_cen,p_1*(surv_2 - p_cen)),P_cell_1)
freq_cause <- function(D,M_1,M_2)
{
  n <- nrow(D)
  M <- replicate(n,0)
  for(i in 1:n)
  {
    M[i] = ((M_1 + 1)*D[i,][2]) + D[i,][3] + 1
  }
  freq <- table(M)
  observed_cat <- as.numeric(names(freq))
  all_possible_cat <- 1:((M_1 + 1)*(M_2 + 1))
  x <- replicate((M_1 + 1)*(M_2 + 1),0)
  if(length(observed_cat) == (M_1 + 1)*(M_2 + 1))
  {
    x <- as.numeric(freq)
    freq_table <- matrix(x,nrow=M_1 + 1,ncol=M_2 + 1,byrow = TRUE)
    rownames(freq_table) <- as.character(0:M_1)
    colnames(freq_table) <- as.character(0:M_2)
    return(freq_table)
  }else
  {
    x[observed_cat] <- as.numeric(freq)
    freq_table <- matrix(x,nrow=M_1 + 1,ncol=M_2 + 1,byrow = TRUE)
    rownames(freq_table) <- as.character(0:M_1)
    colnames(freq_table) <- as.character(0:M_2)
    return(freq_table)
  }
}
corr_indiv_log_likelihood <- function(data,mu_1,mu_2,tau_1,tau_2,
                                      phi,q_1,q_2,q_3)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  d_1 <- mu_1*(tau_1^2)*x
  d_2 <- mu_2*(tau_2^2)*x
  b_1 <- 1 + d_1 + d_2
  b_2 <- 1 + d_1
  b_3 <- 1 + d_2
  c_1 <- phi/(tau_1*tau_2)
  c_2 <- c_1 - (1/tau_1^2)
  c_3 <- c_1 - (1/tau_2^2)
  log_marg_surv_1 <- (-(1/tau_1^2))*log(b_2)
  log_marg_surv_2 <- (-(1/tau_2^2))*log(b_3)
  marg_surv_1 <- exp(log_marg_surv_1)
  marg_surv_2 <- exp(log_marg_surv_2)
  log_joint_surv <- ((-c_1)*log(b_1)) + (c_2*log(b_2)) + (c_3*log(b_3))
  joint_surv <- exp(log_joint_surv)
  c <- (1 - (marg_surv_1) - (marg_surv_2) + joint_surv)
  joint_subdist_func <- (q_3)*c
  if((j_1 != 0) && (j_2 != 0))
    {
      a <- joint_subdist_func[j_1,j_2]
      if(a > 10^(-10))
      {
        return(log(a))
      }else
      {
        return(log(10^(-10)))
      }
    }else if((j_1 != 0) && (j_2 == 0)){
      a_1 <- marg_surv_2
      a_2 <- joint_surv
      if((a_1 - a_2) > 10^(-10))
      {
        return(log(q_1[j_1]) + log(a_1 - a_2))
      }else
      {
        return(log(10^(-10)))
      }
    }else if((j_1 == 0) && (j_2 != 0)){
      a_1 <- marg_surv_1
      a_2 <- joint_surv
      if((a_1 - a_2) > 10^(-10))
      {
        return(log(q_2[j_2]) + log(a_1 - a_2))
      }else
      {
        return(log(10^(-10)))
      }
    }else{
      return(log_joint_surv)
    }
}
neg_log_likelihood_corr <- function(theta,D,M_1,M_2)
{
  lambda_vec_3 <- theta[1:M_1]
  lambda_vec_4 <- theta[(M_1 + 1):(M_1 + M_2)]
  tau_1 <- theta[M_1 + M_2 + 1]
  tau_2 <- theta[M_1 + M_2 + 2]
  phi <- theta[M_1 + M_2 + 3]
  lambda_3 <- sum(lambda_vec_3)
  lambda_4 <- sum(lambda_vec_4)
  q_4 <- (lambda_vec_3)/lambda_3
  q_5 <- (lambda_vec_4)/lambda_4
  q_6 <- outer(q_4,q_5)
  if(min((tau_1/tau_2),(tau_2/tau_1)) > phi)
  {
    return(-sum(apply(D,1,corr_indiv_log_likelihood,
                      mu_1=lambda_3,mu_2=lambda_4,q_1=q_4,
                      q_2=q_5,q_3=q_6,tau_1=tau_1,tau_2=tau_2,phi=phi)))
  }else
  {
    return(10^20)
  }
}
## fast determinant computation in R
Rcpp::cppFunction("double armaDet(const arma::mat & x) { return arma::det(x);}",depends = "RcppArmadillo")
pd_checking <- function(M,epsilon)
{
  n <- nrow(M)
  d <- replicate(n-1,0)
  for(i in 2:n)
  {
    d[i] <- armaDet(M[1:i,1:i])
  }
  d <- c(M[1,1],d)
  l <- length(which(d > epsilon))
  if(l == n)
  {
    return("the matrix is pd")
  }else
  {
    return("the matrix is not pd")
  }
}
## function to find eigen values from Rcpparmadillo
Rcpp::cppFunction("arma::vec armaEigen(arma::mat A){
                  return arma::eig_sym(A);
                  }",depends = "RcppArmadillo")
pd_checking_1 <- function(M,epsilon)
{
  eigen_values <- eigen(M,only.values = TRUE)$values
  l <- length(which(eigen_values > epsilon))
  if(l == nrow(M))
  {
    return("the matrix is pd")
  }else
  {
    return("the matrix is not pd")
  }
}
pd_checking_2 <- function(M,epsilon)
{
  eigen_values <- as.vector(armaEigen(M))
  l <- length(which(eigen_values > epsilon))
  if(l == nrow(M))
  {
    return("the matrix is pd")
  }else
  {
    return("the matrix is not pd")
  }
}
## fast matrix inverse
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x);}",depends = "RcppArmadillo")
asymptotic_sd_corr <- function(R,initial,n,mu_1,mu_2,
                               tau_1,tau_2,phi,M_1,M_2,q_1,q_2,q_3,epsilon)
{
  r <- replicate(2*(M_1+M_2+3),0)
  no_pd <- 0
  data <- corr_data_gen(n,mu_1,mu_2,tau_1,tau_2,phi,M_1,M_2,q_1,q_2,q_3)
  soln <- optim(initial,neg_log_likelihood_corr,D=data,
                M_1=M_1,M_2=M_2,method = "L-BFGS-B",lower = rep(0.01,M_1 + M_2 + 3),
                upper = c(rep(Inf,M_1 + M_2 + 2),0.99),control = list(trace=6))$par
  h <- pracma::hessian(neg_log_likelihood_corr,soln,D=data,M_1=M_1,M_2=M_2)
  if(pd_checking_2(h,epsilon) == "the matrix is pd")
  {
    Sd <- armaInv(h)
    Sd <- diag(Sd)
    Sd <- sqrt(Sd)
    no_pd <- no_pd + 1
    r <- rbind(r,c(soln,Sd))
  }
  print(no_pd)
  while(no_pd < R)
  {
    data <- corr_data_gen(n,mu_1,mu_2,tau_1,tau_2,phi,M_1,M_2,q_1,q_2,q_3)
    soln <- optim(initial,neg_log_likelihood_corr,D=data,
                  M_1=M_1,M_2=M_2,method = "L-BFGS-B",lower = rep(0.01,M_1 + M_2 + 3),
                  upper = c(rep(Inf,M_1 + M_2 + 2),0.99),control = list(trace=6))$par
    h <- pracma::hessian(neg_log_likelihood_corr,soln,D=data,M_1=M_1,M_2=M_2)
    if(pd_checking_2(h,epsilon) == "the matrix is pd")
    {
      Sd <- armaInv(h)
      Sd <- diag(Sd)
      Sd <- sqrt(Sd)
      no_pd <- no_pd + 1
      r <- rbind(r,c(soln,Sd))
    }
    print(no_pd)
  }
  return(r)
}
library("optimParallel")
asymptotic_sd_corr_par <- function(R,initial,n,tau_1,tau_2,phi,mu_1,
                                   mu_2,M_1,M_2,q_1,q_2,q_3,epsilon)
{
  r <- replicate(2*(M_1+M_2+3),0)
  no_pd <- 0
  data <- corr_data_gen(n,mu_1,mu_2,tau_1,tau_2,phi,M_1,M_2,q_1,q_2,q_3)
  soln <- optimParallel(par=initial,fn=neg_log_likelihood_corr,D=data,
                        M_1=M_1,M_2=M_2,method = "L-BFGS-B",lower = rep(0.01,M_1 + M_2 + 3),
                        upper = c(rep(Inf,M_1 + M_2 + 2),0.99),parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE))$par
  h <- pracma::hessian(neg_log_likelihood_corr,soln,D=data,M_1=M_1,M_2=M_2)
  if(pd_checking_2(h,epsilon) == "the matrix is pd")
  {
    Sd <- armaInv(h)
    Sd <- diag(Sd)
    Sd <- sqrt(Sd)
    no_pd <- no_pd + 1
    r <- rbind(r,c(soln,Sd))
  }
  print(no_pd)
  while(no_pd < R)
  {
    data <- corr_data_gen(n,mu_1,mu_2,tau_1,tau_2,phi,M_1,M_2,q_1,q_2,q_3)
    soln <- optimParallel(par=initial,fn=neg_log_likelihood_corr,D=data,
                          M_1=M_1,M_2=M_2,method = "L-BFGS-B",lower = rep(0.01,M_1 + M_2 + 3),
                          upper = c(rep(Inf,M_1 + M_2 + 2),0.99),parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE))$par
    h <- pracma::hessian(neg_log_likelihood_corr,soln,D=data,M_1=M_1,M_2=M_2)
    if(pd_checking_2(h,epsilon) == "the matrix is pd")
    {
      Sd <- armaInv(h)
      Sd <- diag(Sd)
      Sd <- sqrt(Sd)
      no_pd <- no_pd + 1
      r <- rbind(r,c(soln,Sd))
    }
    print(no_pd)
  }
  return(r)
}
neg_log_likelihood_corr_1 <- function(theta,D,M_1,M_2)
{
  lambda_vec_3 <- theta[1:M_1]
  lambda_vec_4 <- theta[(M_1 + 1):(M_1 + M_2)]
  tau_3 <- theta[M_1 + M_2 + 1]
  tau_4 <- theta[M_1 + M_2 + 2]
  phi <- theta[M_1 + M_2 + 3]
  lambda_3 <- sum(lambda_vec_3)
  lambda_4 <- sum(lambda_vec_4)
  q_4 <- (lambda_vec_3)/lambda_3
  q_5 <- (lambda_vec_4)/lambda_4
  q_6 <- outer(q_4,q_5)
  return(-sum(apply(D,1,corr_indiv_log_likelihood,
                    mu_1=lambda_3,mu_2=lambda_4,q_1=q_4,
                    q_2=q_5,q_3=q_6,tau_1=tau_3,tau_2=tau_4,phi=phi)))
}
hin <- function(theta)
{
  h <- rep(NA,1)
  for(i in 1:7)
  {
    h[i] <- theta[i]
  }
  h[8] <- min((theta[5]/theta[6]),(theta[6]/theta[5])) - theta[7]
  return(h)
}
asymptotic_sd_corr_1 <- function(R,initial,n,mu_1,mu_2,tau_1,tau_2,phi,
                                 M_1,M_2,q_1,q_2,q_3,epsilon)
{
  r <- replicate(2*(M_1+M_2+3),0)
  no_pd <- 0
  d <- corr_data_gen(n,mu_1,mu_2,tau_1,tau_2,phi,M_1,M_2,q_1,q_2,q_3)
  objfun <- function(theta)
  {
    return(neg_log_likelihood_corr_1(theta,d,2,2))
  }
  soln <- alabama::constrOptim.nl(par=initial,fn=objfun,hin=hin)$par
  h <- pracma::hessian(neg_log_likelihood_corr_1,soln,D=d,M_1=M_1,M_2=M_2)
  if(pd_checking_2(h,epsilon) == "the matrix is pd")
  {
    Sd <- armaInv(h)
    Sd <- diag(Sd)
    Sd <- sqrt(Sd)
    no_pd <- no_pd + 1
    r <- rbind(r,c(soln,Sd))
  }
  print(no_pd)
  while(no_pd < R)
  {
    d <- corr_data_gen(n,mu_1,mu_2,tau_1,tau_2,phi,M_1,M_2,q_1,q_2,q_3)
    objfun <- function(theta)
    {
      return(neg_log_likelihood_corr_1(theta,d,2,2))
    }
    soln <- alabama::constrOptim.nl(par=initial,fn=objfun,hin=hin)$par
    h <- pracma::hessian(neg_log_likelihood_corr_1,soln,D=d,M_1=M_1,M_2=M_2)
    if(pd_checking_2(h,epsilon) == "the matrix is pd")
    {
      Sd <- armaInv(h)
      Sd <- diag(Sd)
      Sd <- sqrt(Sd)
      no_pd <- no_pd + 1
      r <- rbind(r,c(soln,Sd))
    }
    print(no_pd)
  }
  return(r)
}
coverage_prob_1 <- function(theta_hat,theta_0,s)
{
  return(ifelse((theta_0 > theta_hat - (s*1.96)) &&
                  (theta_0 < theta_hat + (s*1.96)),1,0))
}
coverage_prob <- function(theta_hat_vec,theta_0_vec,sd_vec)
{
  p <- replicate(L_1+L_2+3,0)
  for(i in 1:(L_1+L_2+3))
  {
    p[i] <- coverage_prob_1(theta_hat_vec[i],
                            theta_0_vec[i],sd_vec[i])
  }
  return(p)
}
theta_org <- c(lambda_vec_1,lambda_vec_2,sigma_1,sigma_2,rho)
d <- corr_data_gen(150,lambda_1,lambda_2,sigma_1,sigma_2,rho,2,2,p_1,p_2,p_3)
soln <- system.time(optim(theta_org,neg_log_likelihood_corr,D=d,
                          M_1=2,M_2=2,method = "L-BFGS-B",lower = rep(0.01,7),upper = c(rep(Inf,6),0.99),control = list(trace=6))$par)
soln
library("doParallel")
cl <- makeCluster(15)
registerDoParallel(cl)
Rep <- 10000
n <- 50
R_4_50 <- asymptotic_sd_corr_1(Rep,theta_org,n,lambda_1,lambda_2,
                                sigma_1,sigma_2,rho,M_1=2,M_2=2,p_1,p_2,p_3,10^(-4))
R_4_50 <- R_4_50[-c(1),]
bias_50 <- apply(R_4_50[,1:(L_1+L_2+3)],2,mean) - theta_org
sample_se_50 <- apply(R_4_50[,1:(L_1+L_2+3)],2,sd)
asymptotic_se_50 <- apply(R_4_50[,(L_1+L_2+4):(2*(L_1+L_2+3))],2,mean)
no_pd_50 <- nrow(R_4_50)
prop_pd_50 <- no_pd_50/Rep
prop_pd_50
asymptotic_coverage_prob_50 <- matrix(0,nrow=no_pd_50,ncol=(L_1+L_2+3))
for(i in 1:no_pd_50)
{
  asymptotic_coverage_prob_50[i,] <- coverage_prob(R_4_50[i,1:(L_1+L_2+3)],theta_org,R_4_50[i,(L_1+L_2+4):(2*(L_1+L_2+3))])
}
coverage_prob_vec_50 <- apply(asymptotic_coverage_prob_50,2,mean)
bias_50
sample_se_50
asymptotic_se_50
coverage_prob_vec_50
Rep <- 10000
n <- 150
R_4_150 <- asymptotic_sd_corr_1(Rep,theta_org,n,lambda_1,lambda_2,
                    sigma_1,sigma_2,rho,M_1=2,M_2=2,p_1,p_2,p_3,10^(-4))
R_4_150 <- R_4_150[-c(1),]
bias_150 <- apply(R_4_150[,1:(L_1+L_2+3)],2,mean) - theta_org
sample_se_150 <- apply(R_4_150[,1:(L_1+L_2+3)],2,sd)
asymptotic_se_150 <- apply(R_4_150[,(L_1+L_2+4):(2*(L_1+L_2+3))],2,mean)
no_pd_150 <- nrow(R_4_150)
prop_pd_150 <- no_pd_150/Rep
prop_pd_150
asymptotic_coverage_prob_150 <- matrix(0,nrow=no_pd_150,ncol=(L_1+L_2+3))
for(i in 1:no_pd_150)
{
  asymptotic_coverage_prob_150[i,] <- coverage_prob(R_4_150[i,1:(L_1+L_2+3)],theta_org,R_4_150[i,(L_1+L_2+4):(2*(L_1+L_2+3))])
}
coverage_prob_vec_150 <- apply(asymptotic_coverage_prob_150,2,mean)
bias_150
sample_se_150
asymptotic_se_150
coverage_prob_vec_150
Rep <- 10000
n <- 300
R_4_300 <- asymptotic_sd_corr_1(Rep,theta_org,n,lambda_1,lambda_2,
                              sigma_1,sigma_2,rho,M_1=2,M_2=2,p_1,p_2,p_3,10^(-4))
R_4_300 <- R_4_300[-c(1),]
bias_300 <- apply(R_4_300[,1:(L_1+L_2+3)],2,mean) - theta_org
sample_se_300 <- apply(R_4_300[,1:(L_1+L_2+3)],2,sd)
asymptotic_se_300 <- apply(R_4_300[,(L_1+L_2+4):(2*(L_1+L_2+3))],2,mean)
no_pd_300 <- nrow(R_4_300)
prop_pd_300 <- no_pd_300/Rep
prop_pd_300
asymptotic_coverage_prob_300 <- matrix(0,nrow=no_pd_300,ncol=(L_1+L_2+3))
for(i in 1:no_pd_300)
{
  asymptotic_coverage_prob_300[i,] <- coverage_prob(R_4_300[i,1:(L_1+L_2+3)],theta_org,R_4_300[i,(L_1+L_2+4):(2*(L_1+L_2+3))])
}
coverage_prob_vec_300 <- apply(asymptotic_coverage_prob_300,2,mean)
bias_300
sample_se_300
asymptotic_se_300
coverage_prob_vec_300
stopCluster(cl)

create_estimate_name_corr <- function(param,i,k)
{
  a_1 <- as.character(i)
  a_2 <- as.character(k)
  if(param == "lambda")
  {
    return(paste("lambda","_","hat","_",a_1,"_",a_2,sep = ""))
  }
  if(param == "sigma")
  {
    return(paste("sigma","_","hat",sep = ""))
  }
  if(param == "rho")
  {
    return(paste("rho","_","hat",sep = ""))
  }
}
## to save it in excel file
name <- c("bias","sample_se","asymptotic_se","cov_prob")
output_matrix <- rbind(bias,sample_se,asymptotic_se,coverage_prob_vec)
rownames(output_matrix) <- name
L <- 2
colnames(output_matrix) <- c(sapply(1:L,create_estimate_name_corr,k=1,param="lambda"),
                             sapply(1:L,create_estimate_name_corr,k=2,param="lambda"),
                             "sigma_hat","rho_hat")
write.csv(output_matrix,paste("corr_gamma","sample size = ",toString(n),"&",
                              "censoring prob = ",toString(p_cen),"Rep",toString(Rep),".csv"))
