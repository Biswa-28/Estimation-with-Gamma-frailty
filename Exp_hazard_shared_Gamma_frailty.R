## true value for Exponential baseline cause-specific hazard function of the first individual
lambda_vec_1 <- c(2.4,5.8)    
## true value for Exponential baseline cause-specific hazard function of the second individual
lambda_vec_2 <- c(3.5,4.5)   
## true value of Gamma frailty's dispersion parameter
sigma <- 0.95
lambda_1 <- sum(lambda_vec_1)
lambda_2 <- sum(lambda_vec_2)
### setting the value of joint censoring probability as 0.1
p_cen <- 0.1 
## Integrand for joint survival function
integrand_1 <- function(x,mu,mu_1,mu_2,tau)
{
  return(dexp(x,rate=mu)/((1 + (x*mu_1*tau)
                           + (x*mu_2*tau))^(1/tau)))
}
## Integrand for marginal survival
integrand_2 <- function(x,mu,mu_1,tau)
{
  return(dexp(x,rate=mu)*((1 + (tau*x*mu_1))^(-(1/tau))))
}
## equation to solve for mu for a given censoring prob
f_1 <- function(mu)
{
  return((integrate(integrand_1,0,Inf,mu=mu,
                    mu_1=lambda_1,mu_2=lambda_2,tau=sigma)$value - p_cen))
}
## solution of mu_hat
mu_hat <- nleqslv::nleqslv(0,f_1)$x
### some little calculations
L_1 <- length(lambda_vec_1)
L_2 <- length(lambda_vec_2)
p_1 <- (lambda_vec_1)/lambda_1
p_2 <- (lambda_vec_2)/lambda_2
p_3 <- outer(p_1,p_2)
## function to generate failure causes for two individuals
failure_cause_gen <- function(x,tau,mu_1,mu_2,
                              L_1,L_2,q_1,q_2,q_3)
{
  a_1 <- (tau*mu_1*x)
  a_2 <- (tau*mu_2*x)
  marg_surv_1 <- (1 + a_1)^(-(1/tau))
  marg_surv_2 <- (1 + a_2)^(-(1/tau))
  joint_surv <- (1 + a_1 + a_2)^(-(1/tau))
  joint_failure_prob <- 1 - marg_surv_1 - marg_surv_2 +
    joint_surv
  prob_1 <- q_3*joint_failure_prob
  prob_1 <- rbind(q_2*(marg_surv_1 - joint_surv),prob_1)
  prob_1 <- cbind(c(joint_surv,q_1*(marg_surv_2 - joint_surv)),
                  prob_1)
  s <- which(rmultinom(1,1,prob_1)!=0)
  s_1 <- floor(s/(L_1 + 1))
  s_2 <- s - ((L_1 + 1)*s_1)
  if(s_2 != 0)
  {
    return(c(s_2 - 1,s_1))
  }else
  {
    return(c(L_1,s_1 - 1))
  }
}
## function to generate the data for a given sample size n
data_gen <- function(n,tau,mu_1,mu_2,
                     L_1,L_2,q_4,q_5,q_6)
{
  mon_time <- rexp(n,rate=mu_hat)
  F <- sapply(mon_time,failure_cause_gen,tau=sigma,
              mu_1=lambda_1,mu_2=lambda_2,L_1=L_1,L_2=L_2,
              q_1=q_4,q_2=q_5,q_3=q_6)
  return(cbind(mon_time,t(F)))
}
## Frequency table from observed data
freq_cause <- function(D,L_1,L_2)
{
  n <- nrow(D)
  M <- replicate(n,0)
  for(i in 1:n)
  {
    M[i] = ((L_1 + 1)*D[i,][2]) + D[i,][3] + 1
  }
  freq <- table(M)
  observed_cat <- as.numeric(names(freq))
  all_possible_cat <- 1:((L_1 + 1)*(L_2 + 1))
  x <- replicate((L_1 + 1)*(L_2 + 1),0)
  if(length(observed_cat) == (L_1 + 1)*(L_2 + 1))
  {
    x <- as.numeric(freq)
    freq_table <- matrix(x,nrow=L_1 + 1,ncol=L_2 + 1,byrow = TRUE)
    rownames(freq_table) <- as.character(0:L_1)
    colnames(freq_table) <- as.character(0:L_2)
    return(freq_table)
  }else
  {
    x[observed_cat] <- as.numeric(freq)
    freq_table <- matrix(x,nrow=L_1 + 1,ncol=L_2 + 1,byrow = TRUE)
    rownames(freq_table) <- as.character(0:L_1)
    colnames(freq_table) <- as.character(0:L_2)
    return(freq_table)
  }
}
## individual log-likelihood function
exp_indiv_log_likelihood <- function(data,mu_1,mu_2,
                                     q_1,q_2,q_3,tau)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  a_1 <- (tau*(mu_1)*x)
  a_2 <- (tau*(mu_2)*x)
  surv_1 <- (1 + a_1)^(-(1/tau))
  surv_2 <- (1 + a_2)^(-(1/tau))
  surv_joint <- (1 + a_1 + a_2)^(-(1/tau))
  c <- (1 - (surv_1) - (surv_2) + surv_joint)
  joint_subdist_func <- (q_3)*c
  if((j_1 != 0) && (j_2 != 0))
  {
    return(log(joint_subdist_func[j_1,j_2]))
  }else if((j_1 != 0) && (j_2 == 0)){
    return(log(q_1[j_1]) + log(surv_2 - surv_joint))
  }else if((j_1 == 0) && (j_2 != 0)){
    return(log(q_2[j_2]) + log(surv_1 - surv_joint))
  }else{
    return(log(surv_joint))
  }
}
## the overall negative log-likelihood function which is needed to minimize
neg_log_likelihood_exp <- function(D,theta,M_1,M_2)
{
  lambda_vec_3 <- theta[1:M_1]
  lambda_vec_4 <- theta[(M_1 + 1):(M_1 + M_2)]
  sigma_1 <- theta[M_1 + M_2 + 1]
  lambda_3 <- sum(lambda_vec_3)
  lambda_4 <- sum(lambda_vec_4)
  q_4 <- (lambda_vec_3)/lambda_3
  q_5 <- (lambda_vec_4)/lambda_4
  q_6 <- outer(q_4,q_5)
  return(-sum(apply(D,1,exp_indiv_log_likelihood,
                    mu_1=lambda_3,mu_2=lambda_4,q_1=q_4,
                    q_2=q_5,q_3=q_6,tau=sigma_1)))
}
library("RcppArmadillo")
## function to find eigen values from Rcpparmadillo
Rcpp::cppFunction("arma::vec armaEigen(arma::mat A){
                  return arma::eig_sym(A);
                  }",depends = "RcppArmadillo")
pd_checking <- function(M,epsilon)
{
  n <- nrow(M)
  d <- replicate(n-1,0)
  for(i in 2:n)
  {
    d[i] <- det(M[1:i,1:i])
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
asymptotic_sd_par <- function(R,initial,n,tau,mu_1,mu_2,M_1,
                              M_2,q_1,q_2,q_3,epsilon)
{
  r <- replicate(2*(M_1+M_2+1),0)
  no_pd <- 0
  data <- data_gen(n,tau,mu_1,mu_2,M_1,M_2,q_1,q_2,q_3)
  func_1 <- function(theta)
  {
    return(neg_log_likelihood_exp(data,theta,M_1,M_2))
  }
  soln <- optimParallel::optimParallel(par=initial,fn=func_1,
                                       lower = replicate(M_1+M_2+1,0.01),upper = c(replicate(M_1+M_2,9),0.9),
                                       method="L-BFGS-B",parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE))$par
  h <- pracma::hessian(neg_log_likelihood_exp,soln,D=data,M_1=M_1,M_2=M_2)
  if(pd_checking_2(h,epsilon) == "the matrix is pd")
  {
    Sd <- armaInv(h)
    Sd <- diag(Sd)
    Sd <- sqrt(Sd)
    no_pd <- no_pd + 1
    print(no_pd)
    r <- rbind(r,c(soln,Sd))
  }
  while(no_pd < R)
  {
    data <- data_gen(n,tau,mu_1,mu_2,M_1,M_2,q_1,q_2,q_3)
    func_1 <- function(theta)
    {
      return(neg_log_likelihood_exp(data,theta,M_1,M_2))
    }
    soln <- optimParallel::optimParallel(par=initial,fn=func_1,
                                         lower = replicate(M_1+M_2+1,0.01),upper = c(replicate(M_1+M_2,9),0.9),
                                         method="L-BFGS-B",parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE))$par
    h <- pracma::hessian(neg_log_likelihood_exp,soln,D=data,M_1=M_1,M_2=M_2)
    if(pd_checking_2(h,epsilon) == "the matrix is pd")
    {
      Sd <- armaInv(h)
      Sd <- diag(Sd)
      Sd <- sqrt(Sd)
      no_pd <- no_pd + 1
      print(no_pd)
      r <- rbind(r,c(soln,Sd))
    }
  }
  return(r)
}
coverage_prob_1 <- function(theta_hat,theta_0,s)
{
  return(ifelse((theta_0 > theta_hat - (s*1.96)) &&
                  (theta_0 < theta_hat + (s*1.96)),1,0))
}
## function for asymptotic coverage probability
coverage_prob <- function(theta_hat_vec,theta_0_vec,sd_vec)
{
  p <- replicate(5,0)
  for(i in 1:5)
  {
    p[i] <- coverage_prob_1(theta_hat_vec[i],
                            theta_0_vec[i],sd_vec[i])
  }
  return(p)
}
### true parameter value
theta_org <- c(lambda_vec_1,lambda_vec_2,sigma)
library("doParallel")
cl <- makeCluster(11)
setDefaultCluster(cl=cl)
Rep <- 10000
n <- 300
clusterExport(cl,c("exp_indiv_log_likelihood","neg_log_likelihood_exp"))
R_4 <- asymptotic_sd_par(Rep,theta_org,n,sigma,lambda_1,
                         lambda_2,2,2,p_1,p_2,p_3,10^(-4))
R_4 <- R_4[-c(1),]
bias <- apply(R_4[,1:5],2,mean) - theta_org
sample_se <- apply(R_4[,1:5],2,sd)
asymptotic_se <- apply(R_4[,6:10],2,mean)
no_pd <- nrow(R_4)
prop_pd <- no_pd/Rep
prop_pd
asymptotic_coverage_prob <- array(0,dim=c(no_pd,5))
for(i in 1:no_pd)
{
  asymptotic_coverage_prob[i,] <- coverage_prob(R_4[i,1:5],theta_org,R_4[i,6:10])
}
coverage_prob_vec <- apply(asymptotic_coverage_prob,2,mean)
bias
sample_se
asymptotic_se
coverage_prob_vec
stopCluster(cl)
setDefaultCluster(cl=NULL)
create_estimate_name <- function(i,k)
{
  a_1 <- as.character(i)
  a_2 <- as.character(k)
  return(paste("lambda","_","hat","_",a_1,"_",a_2,sep = ""))
}
name <- c("bias","sample_se","asymptotic_se","cov_prob")
output_matrix <- rbind(bias,sample_se,asymptotic_se,coverage_prob_vec)
rownames(output_matrix) <- name
colnames(output_matrix) <- c(sapply(1:L_1,create_estimate_name,k=1),
                             sapply(1:L_2,create_estimate_name,k=2),"sigma_hat")
write.csv(output_matrix,paste("sample size = ",toString(n),"&",
                              "censoring prob = ",toString(p_cen),"_alternate",".csv"))
