## call the data in r environment
data <- read.csv("~/SharedAndCorrCauseSpecific/R/data.csv")
## checking the structure and missingness of the data
str(data)
attach(data)
Age <- as.character(Age)
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
data_failure_cause <- array(0,dim=c(796,2))
for(i in 1:796)
{
  data_failure_cause[i,] <- failure_cause_gen(data_1[i,c(2,3)])
}
hearing_data <- cbind(data_1[,1],data_failure_cause)
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
## getting the frequency table
F <- freq_cause(hearing_data,3,3)
## getting the row and column marginals
marginal_row <- apply(F,1,sum)
marginal_col <- apply(F,2,sum)
lambda_vec_1_initial <- (marginal_row/sum(Age))[2:4]
lambda_vec_2_initial <- (marginal_col/sum(Age))[2:4]
lambda_1_initial <- sum(lambda_vec_1_initial)
lambda_2_initial <- sum(lambda_vec_2_initial)
## integrand of marginal sub-distribution function
integrand_2 <- function(t,j,mu_vec,gamma,tau_vec)
{
  a <- gamma*log(t*mu_vec)
  a <- exp(a)
  a <- 1 + tau_vec*a
  b <- -(1/tau_vec)*log(a)
  log_h <- (gamma*log(mu_vec[j])) + log(gamma) + ((gamma - 1)*log(t))
  log_val <- log_h + sum(b) - log(a[j])
  return(exp(log_val))
}
## function to compute joint survival function
wei_shared_cause_joint_surv_func <- function(t_1,t_2,mu_vec_1,mu_vec_2,
                                             gamma_1,gamma_2,tau_vec)
{
  b_1 <- 1 + tau_vec*(((t_1*mu_vec_1)^gamma_1) +
                        ((t_2*mu_vec_2)^gamma_2))
  b_2 <- -(1/tau_vec)*log(b_1)
  joint_surv <- exp(sum(b_2))
  return(joint_surv)
}
## integrand of joint sub-distribution function(number of causes is 3)
integrand_3 <- function(t_1,t_2,j_1,j_2,mu_vec_1,
                        mu_vec_2,gamma_1,gamma_2,tau_vec)
{
  joint_surv <- wei_shared_cause_joint_surv_func(t_1,t_2,mu_vec_1,mu_vec_2,
                                                 gamma_1,gamma_2,tau_vec)
  if(j_1 == j_2)
  {
    j = j_1
    b_1 <- gamma_1*log(mu_vec_1[j])
    b_2 <- gamma_2*log(mu_vec_2[j])
    log_a_1 <- b_1 + log(gamma_1) + ((gamma_1 - 1)*log(t_1))
    log_a_2 <- b_2 + log(gamma_2) + ((gamma_2 - 1)*log(t_2))
    log_c_1 <- b_1 + (gamma_1*log(t_1))
    log_c_2 <- b_2 + (gamma_2*log(t_2))
    log_d <- -2*log(1 + tau_vec[j]*(exp(log_c_1) + exp(log_c_2)))
    log_val <- log(1 + tau_vec[j]) + log_a_1 + log_a_2 + log(joint_surv) + log_d
    val <- exp(log_val)
    return(val)
  }else
  {
    b_1 <- gamma_1*log(mu_vec_1[j_1])
    b_2 <- gamma_2*log(mu_vec_2[j_1])
    b_3 <- gamma_1*log(mu_vec_1[j_2])
    b_4 <- gamma_2*log(mu_vec_2[j_2])
    log_a_1 <- b_1 + log(gamma_1) + ((gamma_1 - 1)*log(t_1))
    log_a_4 <- b_4 + log(gamma_2) + ((gamma_2 - 1)*log(t_2))
    log_c_1 <- b_1 + (gamma_1*log(t_1))
    log_c_2 <- b_2 + (gamma_2*log(t_2))
    log_c_3 <- b_3 + (gamma_1*log(t_1))
    log_c_4 <- b_4 + (gamma_2*log(t_2))
    log_d_1 <- -log(1 + tau_vec[j_1]*(exp(log_c_1) + exp(log_c_2)))
    log_d_2 <- -log(1 + tau_vec[j_2]*(exp(log_c_3) + exp(log_c_4)))
    log_val <- log_a_1 + log_a_4 + log(joint_surv) + log_d_1 + log_d_2
    val <- exp(log_val)
    return(val)
  }
}
## marginal distribution function calculation
marg_subdis_weibull_shared_cause_specific <- function(t,j,mu_vec,gamma,tau_vec)
{
  return(SharedAndCorrCauseSpecific::marg_subdis_weibull_shared_cause_specific_cpp(t,j,mu_vec,gamma,tau_vec))
}
## joint sub-distribution function calculation
joint_sub_dis_weibull_shared_cause_specific <- function(t_1,t_2,j_1,j_2,mu_vec_1,mu_vec_2,
                                                        gamma_1,gamma_2,tau_vec)
{
  return(SharedAndCorrCauseSpecific::joint_sub_dis_weibull_shared_cause_specific_cpp(t_1,t_2,
                          j_1,j_2,mu_vec_1,mu_vec_2,gamma_1,gamma_2,tau_vec,3)$approximate)
}
joint_sub_shared_cause <- function(t_1,t_2,j_1,j_2,mu_vec_1,
                                   mu_vec_2,gamma_1,gamma_2,tau_vec)
{
  integrand_new_3 <- function(u,v)
  {
    return(integrand_3(u,v,j_1,j_2,mu_vec_1,
                       mu_vec_2,gamma_1,gamma_2,tau_vec))
  }
  I <- pracma::integral2(Vectorize(integrand_new_3),0,t_1,0,t_2)$Q
  return(I)
}
wei_joint_sub_dis_shared_cause_1 <- function(t_1,t_2,j_1,j_2,mu_vec_1,
                                             mu_vec_2,gamma_1,gamma_2,tau_vec)
{
  z <- function(v)
  {
    I_1 <- sapply(v, g <- function(y)
    {
      g_1 <- function(u)
      {
        return(integrand_3(u,t_1=y,j_1=j_1,j_2=j_2,mu_vec_1=mu_vec_1,
                           mu_vec_2=mu_vec_2,gamma_1=gamma_1,gamma_2=gamma_2,tau_vec=tau_vec))
      }
      return(cubature::cubintegrate(g_1,lower=0,upper=t_1,relTol=1e-9)$integral)
    })
    return(I_1)
  }
  val <- cubature::cubintegrate(z,lower=0,upper=t_2,relTol=1e-9)$integral
  return(val)
}
## cell probability for exp hazard shared cause specific gamma frailty
prob_cell_weibull_shared_cause <- function(x,j_1,j_2,mu_vec_1,mu_vec_2,
                                       gamma_1,gamma_2,tau_vec)
{
  log_joint_surv <- -(1/tau_vec)*log(1 + tau_vec*(((x*mu_vec_1)^gamma_1) +
                                                    ((x*mu_vec_2)^gamma_2)))
  joint_surv <- exp(sum(log_joint_surv))
  if((j_1 == 1)&&(j_2 == 1))
  {
    return(joint_surv)
  }else if((j_1 > 1)&&(j_2 == 1))
  {
    return(marg_subdis_weibull_shared_cause_specific(x,j_1 - 1,mu_vec_1,gamma_1,tau_vec) -
             sum(sapply(1:3,joint_sub_dis_weibull_shared_cause_specific,t_1=x,t_2=x,j_1=j_1 - 1,
                        mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,gamma_1=gamma_1,
                        gamma_2=gamma_2,tau_vec=tau_vec)))
  }else if((j_1 == 1)&&(j_2 > 1))
  {
    return(marg_subdis_weibull_shared_cause_specific(x,j_2 - 1,mu_vec_2,gamma_2,tau_vec) -
             sum(sapply(1:3,joint_sub_dis_weibull_shared_cause_specific,t_1=x,t_2=x,j_2=j_2 - 1,
                        mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,gamma_1=gamma_1,
                        gamma_2=gamma_2,tau_vec=tau_vec)))
  }else
  {
    return(joint_sub_dis_weibull_shared_cause_specific(x,x,j_1 - 1,j_2 - 1,mu_vec_1,
                                        mu_vec_2,gamma_1,gamma_2,tau_vec))
  }
}
## initial estimate
in_est_shared_cause_1 <- function(monitoring_data,Freq_table,M,theta)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  gamma_3 <- theta[(2*M) + 1]
  gamma_4 <- theta[(2*M) + 2]
  tau_vec_1 <- theta[((2*M) + 3):((3*M) + 2)]
  c <- 0
  for(j_1 in 1:(M + 1))
  {
    for(j_2 in 1:(M + 1))
    {
      d <- monitoring_data[monitoring_data[,2] == (j_1 - 1) &
                             monitoring_data[,3] == (j_2 - 1),]
      if(is.matrix(d) == TRUE)
      {
        if(length(d[,1]) != 0)
        {
          c <- c + (Freq_table[j_1,j_2] - sum(sapply(d[,1],
                                                     prob_cell_wei_shared_cause,j_1=j_1,j_2=j_2,mu_vec_1=mu_vec_3,
                                                     mu_vec_2=mu_vec_4,gamma_1=gamma_3,gamma_2=gamma_4,tau_vec=tau_vec_1)))^2
        }
      }
      if(is.vector(d) == TRUE)
      {
        c <- c + (Freq_table[j_1,j_2] - prob_cell_wei_shared_cause(d[1],
                                                                   j_1,j_2,mu_vec_3,mu_vec_4,gamma_3,gamma_4,tau_vec_1))^2
      }
    }
  }
  return(c)
}
## log-likelihood function for an individual
wei_shared_cause_indiv_log_lik <- function(data,mu_vec_1,mu_vec_2,
                                           gamma_1,gamma_2,tau_vec)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  eps <- 10^(-25)
  if((j_1 != 0) && (j_2 != 0))
  {
    h <- joint_sub_dis_weibull_shared_cause_specific(x,x,j_1,j_2,
                                          mu_vec_1,mu_vec_2,gamma_1,gamma_2,tau_vec)
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }else if((j_1 != 0) && (j_2 == 0)){
    h <- marg_subdis_weibull_shared_cause_specific(x,j_1,mu_vec_1,gamma_1,tau_vec) -
      sum(sapply(1:3,joint_sub_dis_weibull_shared_cause_specific,t_1=x,t_2=x,j_1=j_1,
                 mu_vec_1 = mu_vec_1,mu_vec_2 = mu_vec_2,gamma_1=gamma_1,
                 gamma_2=gamma_2,tau_vec=tau_vec))
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }else if((j_1 == 0) && (j_2 != 0)){
    h <- marg_subdis_weibull_shared_cause_specific(x,j_2,mu_vec_2,gamma_2,tau_vec) -
      sum(sapply(1:3,joint_sub_dis_weibull_shared_cause_specific,t_1=x,t_2=x,j_2=j_2,
                 mu_vec_1 = mu_vec_1,mu_vec_2 = mu_vec_2,gamma_1=gamma_1,
                 gamma_2=gamma_2,tau_vec=tau_vec))
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }else{
    h <- wei_shared_cause_joint_surv_func(x,x,mu_vec_1,mu_vec_2,gamma_1,
                                          gamma_2,tau_vec)
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }
}
## negative log-likelihood function which we want to minimize
neg_log_lik_wei_shared_cause <- function(theta,D,M)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  gamma_3 <- theta[(2*M) + 1]
  gamma_4 <- theta[(2*M) + 2]
  tau_vec_1 <- theta[((2*M) + 3):((3*M) + 2)]
  b <- -apply(D,1,wei_shared_cause_indiv_log_lik,
              mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,gamma_1=gamma_3,
              gamma_2=gamma_4,tau_vec=tau_vec_1)
  return(sum(b))
}
library("doParallel")
### parallel version
neg_log_lik_wei_shared_cause_par <- function(theta,D,M,cl)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  gamma_3 <- theta[(2*M) + 1]
  gamma_4 <- theta[(2*M) + 2]
  tau_vec_1 <- theta[((2*M) + 3):((3*M) + 2)]
  b <- -parApply(cl,D,1,wei_shared_cause_indiv_log_lik,
                 mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,gamma_1=gamma_3,
                 gamma_2=gamma_4,tau_vec=tau_vec_1)
  return(sum(b))
}
library("doParallel")
cl <- makeCluster(27)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_shared_cause_specific","wei_shared_cause_joint_surv_func",
                   "joint_sub_dis_weibull_shared_cause_specific"))
neg_log_lik_wei_shared_cause_par(c(lambda_vec_1_initial,lambda_vec_2_initial,
                                   3.4,5.6,c(3.4,6.7,1.9)),hearing_data,3,cl)
stopCluster(cl)
## numerical gradient function of log-likelihood
num_grad_neg_log_lik <- function(D,theta,M)
{
  return(numDeriv::grad(neg_log_lik_shared_cause,theta,D=D,M=M))
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
library("doParallel")
cl <- makeCluster(28)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_shared_cause_specific","wei_shared_cause_joint_surv_func",
                   "joint_sub_dis_weibull_shared_cause_specific"))
neg_log_lik_wei_shared_cause_par(c(lambda_vec_1_initial,lambda_vec_2_initial,4.5,2.3,
                                   c(2.3,4.5,1.2)),hearing_data,3,cl)
neg_log_lik_wei_shared_cause_par(c(0.0235527,0.00299633,0.00416121,0.023592,0.00292502,
                                   0.00416631,3.32483,3.28406,4.06529,1.42555,1.77079),hearing_data,3,cl)
neg_log_lik_wei_shared_cause_par(c(mle_lambda_vec_1_shared_cause,
                                   mle_lambda_vec_2_shared_cause,1.42555,1.77079,mle_sigma_vec_shared_cause),hearing_data,3,cl)
stopCluster(cl)
# finding a better initial point for minimization of original
## neg-loglikelihood function
l <- seq(0.1,2,length=10)
in_est_2d_array <- array(0,dim=c(10,10))
mle_lambda_vec_1_shared_cause <- data_soln_2_par[1:3]
mle_lambda_vec_2_shared_cause <- data_soln_2_par[4:6]
mle_sigma_vec_shared_cause <- c(data_soln_2_par[7:8],7.657)
for(i in 1:length(l))
{
  for(j in 1:length(l))
  {
    in_est_2d_array[i,j] <- neg_log_lik_wei_shared_cause(c(mle_lambda_vec_1_shared_cause,
                                                           mle_lambda_vec_2_shared_cause,l[i],l[j],mle_sigma_vec_shared_cause),
                                                         hearing_data,3)
  }
}
min(in_est_2d_array)
which.min(in_est_2d_array)
### [6,6]th element is the least element and minimum value is 877.8887
### let us take c(l[6],l[6]) as initial estimate for beta_1,beta_2
library("doParallel")
cl <- makeCluster(3)
setDefaultCluster(cl=cl)
clusterExport(cl,c("prob_cell_shared_cause","marg_subdis_shared_cause",
                   "integrand_2","integrand_3","joint_sub_shared_cause","joint_surv_func"))
in_est <- optimParallel::optimParallel(par=c(lambda_vec_1_initial,lambda_vec_2_initial,replicate(5,1)),
                                       fn=in_est_shared_cause_3,monitoring_data = hearing_data,Freq_table=F,
                                       M = 3,method="L-BFGS-B",lower = replicate(11,0.01),
                                       parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
setDefaultCluster(cl=NULL)
in_est$par
in_est$convergence
in_est_shared_cause_1(in_est$par,hearing_data,F,3)
in_est_shared_cause_1(c(lambda_vec_1_initial,lambda_vec_2_initial,replicate(5,1)),
                      hearing_data,F,3)
in_est_data <- in_est$par
## min in_est_1 at sigma = 3.9,value = 4687.594,taking seq from 0.1 to 4
## Minimization of original neg-loglikelihood function for shared frailty
### this is the initial value considered using the m.l.e of exponential
### cause specific hazard with shared Gamma frailty model
initial_1 <- c(mle_lambda_vec_1_shared_cause,mle_lambda_vec_2_shared_cause,
               1.42555,1.77079,mle_sigma_vec_shared_cause)
initial_2 <- c(mle_lambda_vec_1_shared_cause,mle_lambda_vec_2_shared_cause,
               l[6],l[6],mle_sigma_vec_shared_cause)
library("doParallel")
cl <- makeCluster(23)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_shared_cause_specific","wei_shared_cause_joint_surv_func",
  "joint_sub_dis_weibull_shared_cause_specific","wei_shared_cause_indiv_log_lik"))
data_soln <- optimParallel::optimParallel(initial_2,
fn=neg_log_lik_wei_shared_cause,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(11,0.001),control = list(trace=6,ndeps=replicate(11,1e-5),maxit=200),
parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
#data_soln_1 <- c(35.9748,0.00790753,0.00864017,39.2024,0.00715227,
                 #0.00881737,1.34617,1.3219,13.3555,45.5874,18.2477)
stopCluster(cl)
data_soln_1_par <- data_soln$par
### checking hessian is p.d or not
hess_1 <- data_soln$hessian
eigen(hess_1)$values
### not all positive eigenvalues hence have to perturbate last three co-ordinates
initial_3 <- c(data_soln_1_par[1:8],28.67459,37.864,3.9675)
library("doParallel")
cl <- makeCluster(23)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_shared_cause_specific","wei_shared_cause_joint_surv_func",
                   "joint_sub_dis_weibull_shared_cause_specific","wei_shared_cause_indiv_log_lik"))
data_soln <- optimParallel::optimParallel(initial_3,
fn=neg_log_lik_wei_shared_cause,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(11,0.001),control = list(trace=6,ndeps=replicate(11,1e-5),maxit=200),
parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
### checking hessian is p.d or not
hess_2 <- data_soln$hessian
eigen(hess_2)$values
### not all positive eigenvalues hence have to perturbate last co-ordinate
data_soln_2_par <- data_soln$par
initial_4 <- c(data_soln_2_par[1:10],2.5675)
library("doParallel")
cl <- makeCluster(23)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_shared_cause_specific","wei_shared_cause_joint_surv_func",
                   "joint_sub_dis_weibull_shared_cause_specific","wei_shared_cause_indiv_log_lik"))
data_soln_3 <- optimParallel::optimParallel(initial_3,
fn=neg_log_lik_wei_shared_cause,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(11,0.001),control = list(trace=6,ndeps=replicate(11,1e-5),maxit=300),
parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
### checking hessian is p.d or not
hess_3 <- data_soln_3$hessian
eigen(hess_3)$values
initial_4 <- data_soln_3$par
cl <- makeCluster(23)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_shared_cause_specific","wei_shared_cause_joint_surv_func",
                   "joint_sub_dis_weibull_shared_cause_specific","wei_shared_cause_indiv_log_lik"))
data_soln_4 <- optimParallel::optimParallel(initial_4,
fn=neg_log_lik_wei_shared_cause,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(11,0.001),control = list(trace=6,ndeps=replicate(11,1e-5),maxit=300),
parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
### checking hessian is p.d or not
hess_4 <- data_soln_4$hessian
eigen(hess_4)$values
initial_5 <- c(data_soln_4$par[1:8],31.707321,data_soln_4$par[10:11])
cl <- makeCluster(23)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_shared_cause_specific","wei_shared_cause_joint_surv_func",
                   "joint_sub_dis_weibull_shared_cause_specific","wei_shared_cause_indiv_log_lik"))
data_soln_5 <- optimParallel::optimParallel(initial_5,
fn=neg_log_lik_wei_shared_cause,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(11,0.0001),control = list(trace=6,ndeps=replicate(11,1e-5),maxit=300),
parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
### checking hessian is p.d or not
hess_5 <- data_soln_5$hessian
eigen(hess_5)$values
### let us now apply nearPD function to compute sd
final_ml_estimate <- data_soln_5$par
var_final_ml_estimate <- diag(solve(Matrix::nearPD(hess_5)$mat))
asymptotic_sd_final_ml_estimate <- sqrt(var_final_ml_estimate)
### mle and it's estimated asymptotic standard error
final_ml_estimate
asymptotic_sd_final_ml_estimate
initial_6 <- c(data_soln_5$par[1:8],33.3591,37.9675,data_soln_5$par[11])
cl <- makeCluster(23)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_shared_cause_specific","wei_shared_cause_joint_surv_func",
                   "joint_sub_dis_weibull_shared_cause_specific","wei_shared_cause_indiv_log_lik"))
data_soln_6 <- optimParallel::optimParallel(initial_6,
fn=neg_log_lik_wei_shared_cause,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(11,0.0001),control = list(trace=6,ndeps=replicate(11,1e-5),maxit=300),
parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
### checking hessian is p.d or not
hess_6 <- data_soln_6$hessian
eigen(hess_6)$values
hess <- pracma::hessian(neg_log_lik_shared_cause,data_soln,D=hearing_data,M=3)
eigen(hess)$values
sd_mle <- diag(solve(hess))
library("doParallel")
cl <- makeCluster(30)
setDefaultCluster(cl=cl)
clusterExport(cl,c("wei_shared_cause_indiv_log_lik","wei_joint_sub_dis_shared_cause_1",
                   "integrand_3","wei_shared_cause_joint_surv_func","wei_marg_subdis_shared_cause","integrand_2"))
data_soln_1 <- c(35.9748,0.00790753,0.00864017,39.2024,0.00715227,
                 0.00881737,1.34617,1.3219,13.3555,45.5874,18.2477)
#neg_log_lik_shared_cause(hearing_data,in_est_data,3)
#neg_log_lik_shared_cause(hearing_data,data_soln_1,3)
hess <- pracma::hessian(neg_log_lik_wei_shared_cause_par,data_soln_1,D=hearing_data,M=3,cl=cl)
eigen(hess)$values
stopCluster(cl)
## Computation of AIC for shared cause specific Gamma frailty model
AIC_shared <- (2*neg_log_lik_shared_cause(hearing_data,data_soln,3)) + (2*9)
AIC_shared
