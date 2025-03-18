## Here, we will consider the hazard of exponential as 1/lambda
## and mean will be lambda.
## True values of parameter vectors
lambda_vec_1 <- c(7,6)
lambda_vec_2 <- c(8.5,10)
sigma_vec <- c(0.65,0.85)
## joint censoring probability
p_cen <- 0.1
lambda_1 <- sum(lambda_vec_1)
lambda_2 <- sum(lambda_vec_2)
integrand_1 <- function(x,mu,mu_vec_1,mu_vec_2,tau_vec)
{
  v <- -(1/tau_vec)*log(1 + tau_vec*x*((1/mu_vec_1) + (1/mu_vec_2)))
  return((exp(sum(v)))*dexp(x,rate = (1/mu)))
}
## equation to solve for mu for a given censoring prob
f_1 <- function(mu)
{
  f <- function(u)
  {
    return(integrand_1(u,mu,lambda_vec_1,lambda_vec_2,sigma_vec))
  }
  return((integrate(Vectorize(f),0,Inf)$value - p_cen))
}
## solution of mu_hat
mu_hat <- nleqslv::nleqslv(0.4,f_1)$x
f_1(mu_hat)
L_1 <- length(lambda_vec_1)
L_2 <- length(lambda_vec_2)
p_1 <- (lambda_vec_1)/lambda_1
p_2 <- (lambda_vec_2)/lambda_2
p_3 <- outer(p_1,p_2)
## integrand of marginal sub-distribution function
integrand_2 <- function(t,j,mu_vec,tau_vec)
{
  a <- 1 + tau_vec*t*(1/mu_vec)
  b <- -(1/tau_vec)*log(a)
  marg_surv <- exp(sum(b))
  return((1/mu_vec[j])*marg_surv*((a[j])^(-1)))
}
## integrand of joint sub-distribution function
integrand_3 <- function(t_1,t_2,j_1,j_2,mu_vec_1,
                        mu_vec_2,tau_vec)
{
  a_1 <- 1/mu_vec_1
  a_2 <- 1/mu_vec_2
  b_1 <- 1 + tau_vec[1]*((t_1*a_1[1]) + (t_2*a_2[1]))
  b_2 <- 1 + tau_vec[2]*((t_1*a_1[2]) + (t_2*a_2[2]))
  if(j_1 == j_2)
  {
    j = j_1
    if(j==1)
    {
      d_1 <- b_1^(-2 - (1/tau_vec[1]))
      d_2 <- b_2^(-1/tau_vec[2])
      return(a_1[1]*a_2[1]*(1 + tau_vec[1])*d_1*d_2)
    }
    if(j==2)
    {
      d_1 <- b_1^(-1/tau_vec[1])
      d_2 <- b_2^(-2 - (1/tau_vec[2]))
      return(a_1[2]*a_2[2]*(1 + tau_vec[2])*d_1*d_2)
    }
  }else
  {
    d_1 <- b_1^(-1 - (1/tau_vec[1]))
    d_2 <- b_2^(-1 - (1/tau_vec[2]))
    return(a_1[j_1]*a_2[j_2]*d_1*d_2)
  }
}
### compiling necessary C++ files for numerical integration
Rcpp::sourceCpp("numerical_integration.cpp")
Rcpp::sourceCpp("exp_hazard_shared_cause_specific_gamma_Ccode.cpp")
## marginal distribution function calculation
marg_subdis_shared_cause <- function(t,j,mu_vec,tau_vec)
{
  f_1 <- function(u)
  {
    return(integrand_2(u,j,mu_vec,tau_vec))
  }
  return((integrate(Vectorize(f_1),0,t)$value))
}
## joint sub-distribution function calculation
joint_sub_shared_cause <- function(t_1,t_2,j_1,j_2,mu_vec_1,
                          mu_vec_2,tau_vec)
{
  I <- pracma::integral2(integrand_3,0,t_1,0,t_2,j_1=j_1,j_2=j_2,
       mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,tau_vec=tau_vec)$Q
  return(I)
}

## function to generate failure causes for two individuals
failure_cause_gen <- function(x,mu_vec_1,mu_vec_2,tau_vec,M)
{
  a <- tau_vec*x
  c_1 <- a*(1/mu_vec_1)
  c_2 <- a*(1/mu_vec_2)
  d <- (1/tau_vec)
  v_3 <- -d*log(1 + c_1 + c_2)
  joint_surv <- exp(sum(v_3))
  prob_1 <- matrix(0,nrow=M,ncol=M)
  for(i in 1:M)
  {
      prob_1[i,] = sapply(1:M,g <- function(z)
      {
        h <- joint_sub_shared_cause(x,x,i,z,mu_vec_1,mu_vec_2,tau_vec)
        if(h > 10^(-9))
        {
          return(h)
        }else
        {
          return(10^(-9))
        }
      })
  }
  if(sum(prob_1) < 1 - (10^(-45)))
  {
    q_1 <- sapply(1:M,f <- function(y)
    {
      h_1 <- marg_subdis_shared_cause(x,y,mu_vec_1,tau_vec)
      h_2 <- apply(prob_1,1,sum)[y]
      if((h_1 - h_2) > 10^(-10))
      {
        return(h_1 - h_2)
      }else
      {
        return(10^(-10))
      }
    })
    q_2 <- sapply(1:M,f <- function(y)
    {
      h_3 <- marg_subdis_shared_cause(x,y,mu_vec_2,tau_vec)
      h_4 <- apply(prob_1,2,sum)[y]
      if((h_3 - h_4) > 10^(-10))
      {
        return(h_3 - h_4)
      }else
      {
        return(10^(-10))
      }
    })
    q_1 <- unlist(q_1)
    q_2 <- unlist(q_2)
    prob_1 <- rbind(q_2,prob_1)
    prob_1 <- cbind(c(joint_surv,q_1),prob_1)
    s <- which(rmultinom(1,1,prob_1)!=0)
    s_1 <- floor(s/(M + 1))
    s_2 <- s - ((M + 1)*s_1)
    if(s_2 != 0)
    {
      return(c(s_2 - 1,s_1))
    }else
    {
      return(c(M,s_1 - 1))
    }
  }else
  {
    q <- replicate(M,0)
    prob_1 <- rbind(q,prob_1)
    prob_1 <- cbind(c(0,q),prob_1)
    s <- which(rmultinom(1,1,prob_1)!=0)
    s_1 <- floor(s/(M + 1))
    s_2 <- s - ((M + 1)*s_1)
    if(s_2 != 0)
    {
      return(c(s_2 - 1,s_1))
    }else
    {
      return(c(M,s_1 - 1))
    }
  }
}
## function to generate failure causes for two individuals using C++ functions
failure_cause_gen_1 <- function(x,mu_vec_1,mu_vec_2,tau_vec,M)
{
  a <- tau_vec*x
  c_1 <- a*(1/mu_vec_1)
  c_2 <- a*(1/mu_vec_2)
  d <- (1/tau_vec)
  v_3 <- -d*log(1 + c_1 + c_2)
  joint_surv <- exp(sum(v_3))
  prob_1 <- matrix(0,nrow=M,ncol=M)
  for(i in 1:M)
  {
    prob_1[i,] = sapply(1:M,joint_sub_dis_func_cpp$approximate,t_1=x,
    t_2=x,l_1=i,mu_vec_3=mu_vec_1,mu_vec_4=mu_vec_2,tau_vec_1=tau_vec)
  }
  q_1 <- sapply(1:M,f <- function(y)
  {
    h_1 <- marg_subdis_shared_cause_cpp(t=x,l=y,mu_vec_1=mu_vec_1,tau_vec_1=tau_vec)
    h_2 <-  apply(prob_1,1,sum)[y]
    if((h_1 - h_2) > 10^(-10))
    {
      return(h_1 - h_2)
    }else
    {
      return(10^(-10))
    }
  })
  q_2 <- sapply(1:M,f <- function(y)
  {
    h_3 <- marg_subdis_shared_cause_cpp(t=x,l=y,mu_vec_1=mu_vec_2,tau_vec_1=tau_vec)
    h_4 <- apply(prob_1,2,sum)[y]
    if((h_3 - h_4) > 10^(-10))
    {
      return(h_3 - h_4)
    }else
    {
      return(10^(-10))
    }
  })
  q_1 <- unlist(q_1)
  q_2 <- unlist(q_2)
  prob_1 <- rbind(q_2,prob_1)
  prob_1 <- cbind(c(joint_surv,q_1),prob_1)
  s <- which(rmultinom(1,1,prob_1)!=0)
  s_1 <- floor(s/(M + 1))
  s_2 <- s - ((M + 1)*s_1)
  if(s_2 != 0)
  {
    return(c(s_2 - 1,s_1))
  }else
  {
    return(c(M,s_1 - 1))
  }
}
## function to generate the data for a given sample size n
shared_cause_data_gen <- function(n,mu_vec_3,mu_vec_4,tau_vec_1,M)
{
  mon_time <- rexp(n,rate=(1/mu_hat))
  F <- sapply(mon_time,failure_cause_gen,tau_vec=tau_vec_1,
  mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,M=M)
  return(cbind(mon_time,t(F)))
}
## function to generate the data for a given sample size n using C++ functions
shared_cause_data_gen_1 <- function(n,mu_vec_3,mu_vec_4,tau_vec_1,M)
{
  mon_time <- rexp(n,rate=(1/mu_hat))
  F <- sapply(mon_time,failure_cause_gen_1,tau_vec=tau_vec_1,
              mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,M=M)
  return(cbind(mon_time,t(F)))
}
## function to create frequency table
freq_cause <- function(D,M)
{
  n <- nrow(D)
  v <- replicate(n,0)
  for(i in 1:n)
  {
    v[i] = ((M + 1)*D[i,][2]) + D[i,][3] + 1
  }
  freq <- table(v)
  observed_cat <- as.numeric(names(freq))
  all_possible_cat <- 1:((M + 1)*(M + 1))
  x <- replicate((M + 1)*(M + 1),0)
  if(length(observed_cat) == (M + 1)*(M + 1))
  {
    x <- as.numeric(freq)
    freq_table <- matrix(x,nrow=M + 1,ncol=M + 1,byrow = TRUE)
    rownames(freq_table) <- as.character(0:M)
    colnames(freq_table) <- as.character(0:M)
    return(freq_table)
  }else
  {
    x[observed_cat] <- as.numeric(freq)
    freq_table <- matrix(x,nrow=M + 1,ncol=M + 1,byrow = TRUE)
    rownames(freq_table) <- as.character(0:M)
    colnames(freq_table) <- as.character(0:M)
    return(freq_table)
  }
}
## computation of theoritical cell prob
surv_1 <- integrate(integrand_2,0,Inf,mu=mu_hat,
                    mu_1=lambda_1,tau=sigma)$value
surv_2 <- integrate(integrand_2,0,Inf,mu=mu_hat,
                    mu_1=lambda_2,tau=sigma)$value
surv_joint <- integrate(integrand_1,0,Inf,mu=mu_hat,
                        mu_1=lambda_1,mu_2=lambda_2,tau=sigma)$value
subdist_func_1 <- (p_1)*(1 - surv_1)
subdist_func_2 <- (p_2)*(1 - surv_2)
c <- (1 - (surv_1) - (surv_2) + surv_joint)
joint_subdist_func <- (p_3)*c
P_cell_1 <- rbind(p_2*(surv_1 - p_cen),joint_subdist_func)
P_cell_1 <- cbind(c(p_cen,p_1*(surv_2 - p_cen)),P_cell_1)
## returning theoritical cell prob matrix
P_cell_1
d <- system.time(shared_cause_data_gen(300,lambda_vec_1,lambda_vec_2,sigma_vec,2))
d_1 <- system.time(shared_cause_data_gen_1(300,lambda_vec_1,lambda_vec_2,sigma_vec,2))
d
d_1
freq_cause(d,2)
## log-likelihood function for an individual data
shared_cause_indiv_log_lik <- function(data,mu_vec_1,mu_vec_2,tau_vec)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  a <- tau_vec*x
  c_1 <- a*(1/mu_vec_1)
  c_2 <- a*(1/mu_vec_2)
  d <- (1/tau_vec)
  v <- -d*log(1 + c_1 + c_2)
  log_joint_surv <- sum(v)
  if((j_1 != 0) && (j_2 != 0))
  {
    return(log(joint_sub_shared_cause(x,x,j_1,j_2,
              mu_vec_1,mu_vec_2,tau_vec)))
  }else if((j_1 != 0) && (j_2 == 0)){
    return(log(marg_subdis_shared_cause(x,j_1,mu_vec_1,tau_vec) -
    sum(sapply(1:2,joint_sub_shared_cause,t_1=x,t_2=x,j_1=j_1,
    mu_vec_1 = mu_vec_1,mu_vec_2 = mu_vec_2,tau_vec=tau_vec))))
  }else if((j_1 == 0) && (j_2 != 0)){
    return(log(marg_subdis_shared_cause(x,j_2,mu_vec_2,tau_vec) -
    sum(sapply(1:2,joint_sub_shared_cause,t_1=x,t_2=x,j_2=j_2,
    mu_vec_1 = mu_vec_1,mu_vec_2 = mu_vec_2,tau_vec=tau_vec))))
  }else{
    return(log_joint_surv)
  }
}
## same function as above using C++ functions instead of R 
shared_cause_indiv_log_lik_1 <- function(data,mu_vec_1,mu_vec_2,tau_vec)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  a <- tau_vec*x
  c_1 <- a*(1/mu_vec_1)
  c_2 <- a*(1/mu_vec_2)
  d <- (1/tau_vec)
  v <- -d*log(1 + c_1 + c_2)
  log_joint_surv <- sum(v)
  if((j_1 != 0) && (j_2 != 0))
  {
    return(log(joint_sub_dis_func_cpp(x,x,j_1,j_2,
                                      mu_vec_1,mu_vec_2,tau_vec)))
  }else if((j_1 != 0) && (j_2 == 0)){
    return(log(marg_subdis_shared_cause_cpp(x,j_1,mu_vec_1,tau_vec) -
                 sum(sapply(1:2,joint_sub_dis_func_cpp,t_1=x,t_2=x,l_1=j_1,
                            mu_vec_3 = mu_vec_1,mu_vec_4 = mu_vec_2,tau_vec_1=tau_vec))))
  }else if((j_1 == 0) && (j_2 != 0)){
    return(log(marg_subdis_shared_cause_cpp(x,j_2,mu_vec_2,tau_vec) -
                 sum(sapply(1:2,joint_sub_dis_func_cpp,t_1=x,t_2=x,l_2=j_2,
                            mu_vec_3 = mu_vec_1,mu_vec_4 = mu_vec_2,tau_vec_1=tau_vec))))
  }else{
    return(log_joint_surv)
  }
}
## negative log-likelihood function which we want to minimize
neg_log_lik_shared_cause <- function(theta,D,M)
{
  lambda_vec_3 <- theta[1:M]
  lambda_vec_4 <- theta[(M + 1):(2*M)]
  tau_vec_1 <- theta[((2*M) + 1):(3*M)]
  return(-sum(apply(D,1,shared_cause_indiv_log_lik,
                    mu_vec_1=lambda_vec_3,mu_vec_2=lambda_vec_4,
                    tau_vec=tau_vec_1)))
}
## negative log-likelihood using C++ functions
neg_log_lik_shared_cause_1 <- function(theta,D,M)
{
  lambda_vec_3 <- theta[1:M]
  lambda_vec_4 <- theta[(M + 1):(2*M)]
  tau_vec_1 <- theta[((2*M) + 1):(3*M)]
  return(-sum(apply(D,1,shared_cause_indiv_log_lik_1,
                    mu_vec_1=lambda_vec_3,mu_vec_2=lambda_vec_4,
                    tau_vec=tau_vec_1)))
}
## numerical gradient function of log-likelihood
num_grad_neg_log_lik <- function(theta,D,M)
{
  return(numDeriv::grad(neg_log_lik_shared_cause,theta,D=D,M=M))
}
## numerical gradient function of log-likelihood using C++ functions
num_grad_neg_log_lik_1 <- function(theta,D,M)
{
  return(numDeriv::grad(neg_log_lik_shared_cause_1,theta,D=D,M=M))
}
## fast determinant computation
Rcpp::cppFunction("double armaDet(const arma::mat & x) { return arma::det(x);}",depends = "RcppArmadillo")
##function to check whether a matrix is p.d or not
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
## to compute fast matrix inverse
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x);}",depends = "RcppArmadillo")
# function to compute asymptotic s.e for a given sample of size n
asymptotic_sd_shared_cause <- function(R,initial,n,tau_vec,mu_vec_1,
                                       mu_vec_2,M,epsilon)
{
  r <- replicate(2*(3*M),0)
  no_pd <- 0
  data <- shared_cause_data_gen(n,mu_vec_1,mu_vec_2,tau_vec,M)
  soln <- optim(initial,neg_log_lik_shared_cause,num_grad_neg_log_lik,
                D=data,M=M,method = "L-BFGS-B",lower = replicate(3*M,0.01))$par
  h <- numDeriv::hessian(neg_log_lik_shared_cause,soln,D=data,M=M)
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
    data <- shared_cause_data_gen(n,mu_vec_1,mu_vec_2,tau_vec,M)
    soln <- optim(initial,neg_log_lik_shared_cause,num_grad_neg_log_lik,
                  D=data,M=M,method = "L-BFGS-B",lower = replicate(3*M,0.01))$par
    h <- numDeriv::hessian(neg_log_lik_shared_cause,soln,D=data,M=M)
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
# function to compute asymptotic s.e for a given sample of size n parallely
asymptotic_sd_shared_cause_par <- function(R,initial,n,tau_vec,mu_vec_1,
                                       mu_vec_2,M,epsilon)
{
  r <- replicate(2*(3*M),0)
  no_pd <- 0
  data <- shared_cause_data_gen(n,mu_vec_1,mu_vec_2,tau_vec,M)
  soln <- optimParallel::optimParallel(par=initial,fn=neg_log_lik_shared_cause,
  D=data,M=2,method = "L-BFGS-B",lower=replicate(6,0.01),
  parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE))$par
  h <- numDeriv::hessian(neg_log_lik_shared_cause,soln,D=data,M=M)
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
    data <- shared_cause_data_gen(n,mu_vec_1,mu_vec_2,tau_vec,M)
    soln <- optimParallel::optimParallel(par=initial,fn=neg_log_lik_shared_cause,
    D=data,M=2,method = "L-BFGS-B",lower=replicate(6,0.01),
    parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE))$par
    h <- numDeriv::hessian(neg_log_lik_shared_cause,soln,D=data,M=M)
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
## function to compute the coverage probabilities
coverage_prob <- function(theta_hat_vec,theta_0_vec,sd_vec,M)
{
  p <- replicate(3*M,0)
  for(i in 1:(3*M))
  {
    p[i] <- coverage_prob_1(theta_hat_vec[i],
                            theta_0_vec[i],sd_vec[i])
  }
  return(p)
}
## true parameter vector
theta_org <- c(lambda_vec_1,lambda_vec_2,sigma_vec)
## parallel computation of asymptotic standard error
library("doParallel")
cl <- makeCluster(3)
setDefaultCluster(cl=cl)
## to export the functions to the cluster
clusterExport(cl,c("joint_sub_shared_cause","shared_cause_indiv_log_lik",
  "integrand_3","marg_subdis_shared_cause","integrand_2"))
Rep <- 5
n <- 50
R_4 <- asymptotic_sd_shared_cause_par(Rep,theta_org,n,sigma_vec,
lambda_vec_1,lambda_vec_2,M=2,10^(-4))
R_4
R_4 <- R_4[-c(1),]
bias <- apply(R_4[,1:6],2,mean) - theta_org
sample_se <- apply(R_4[,1:6],2,sd)
asymptotic_se <- apply(R_4[,7:(2*6)],2,mean)
no_pd <- nrow(R_4)
prop_pd <- no_pd/Rep
prop_pd
asymptotic_coverage_prob <- matrix(0,nrow=no_pd,ncol=6)
for(i in 1:no_pd)
{
  asymptotic_coverage_prob[i,] <- coverage_prob(R_4[i,1:6],theta_org,R_4[i,7:(2*6)],2)
}
coverage_prob_vec <- apply(asymptotic_coverage_prob,2,mean)
bias
sample_se
asymptotic_se
coverage_prob_vec
stopCluster(cl)
setDefaultCluster(cl=NULL)
create_estimate_name <- function(param,i,k)
{
  a_1 <- as.character(i)
  a_2 <- as.character(k)
  if(param == "lambda")
  {
    return(paste("lambda","_","hat","_",a_1,"_",a_2,sep = ""))
  }
  if(param == "sigma")
  {
    return(paste("sigma","_","hat","_",a_1,sep = ""))
  }
}
## to save it in excel file
name <- c("bias","sample_se","asymptotic_se","cov_prob")
output_matrix <- rbind(bias,sample_se,asymptotic_se,coverage_prob_vec)
rownames(output_matrix) <- name
L <- 2
colnames(output_matrix) <- c(sapply(1:L,create_estimate_name,k=1,param="lambda"),
sapply(1:L,create_estimate_name,k=2,param="lambda"),
sapply(1:L,create_estimate_name,k=2,param="sigma"))
write.csv(output_matrix,paste("shared_cause","sample size = ",toString(n),"&",
"censoring prob = ",toString(p_cen),"Rep",toString(Rep),".csv"))
