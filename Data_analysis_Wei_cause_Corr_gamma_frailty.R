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
data_failure_cause <- t(apply(data_1[,c(2,3)],1,failure_cause_gen))
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
## Marginal survival function of Weibull cause specific hazard
## with correlated frailty model
marg_surv_wei_corr <- function(x,gamma,mu_vec,tau)
{
  a <- sum(mu_vec^gamma)
  h <- (1 + ((tau^2)*a*(x^gamma)))^(-(1/(tau^2)))
  eps <- 10^(-15)
  if(h > eps)
  {
    return(h)
  }else
  {
    return(eps)
  }
}
## Joint survival function of Weibull cause specific hazard
## with correlated frailty model
joint_surv_wei_corr <- function(x,gamma_1,gamma_2,mu_vec_1,
                                mu_vec_2,tau_1,tau_2,phi)
{
  a_1 <- mu_vec_1^gamma_1
  a_2 <- mu_vec_2^gamma_2
  s_1 <- (tau_1^2)*sum(a_1)*(x^gamma_1)
  s_2 <- (tau_2^2)*sum(a_2)*(x^gamma_2)
  b_1 <- phi/(tau_1*tau_2)
  b_2 <- b_1 - 1/(tau_1^2)
  b_3 <- b_1 - 1/(tau_2^2)
  log_s <- -(b_1*log(1 + s_1 + s_2)) + (b_2*log(1 + s_1)) +
    (b_3*log(1 + s_2))
  joint_surv <- exp(log_s)
  eps <- 10^(-15)
  if(joint_surv > eps)
  {
    return(joint_surv)
  }else
  {
    return(eps)
  }
}
## Marginal sub-distribution function of Weibull cause specific hazard
## with correlated frailty model
marg_subdis_wei_corr <- function(x,j,gamma,mu_vec,tau)
{
  a <- mu_vec^gamma
  p_j <- (mu_vec[j]^gamma)/sum(a)
  h <- p_j*marg_surv_wei_corr(x,gamma,mu_vec,tau)
  eps <- 10^(-15)
  if(h > eps)
  {
    return(h)
  }else
  {
    return(eps)
  }
}
## Joint sub-distribution function of Weibull cause specific hazard
## with correlated frailty model
joint_subdis_wei_corr <- function(t_1,t_2,j_1,j_2,mu_vec_1,mu_vec_2,
                                  gamma_1,gamma_2,tau_1,tau_2,phi)
{
  a_1 <- mu_vec_1^gamma_1
  a_2 <- mu_vec_2^gamma_2
  sum_a_1 <- sum(a_1)
  sum_a_2 <- sum(a_2)
  log_constant <- log(a_1[j_1]) + log(a_2[j_2]) + log(gamma_1) +
    log(gamma_2) + (2*log(tau_1) + 2*log(tau_2))
  b_1 <- phi/(tau_1*tau_2)
  b_2 <- b_1 - 1/(tau_1^2)
  b_3 <- b_1 - 1/(tau_2^2)
  integrand_1 <- function(arg)
  {
    u <- arg[1]
    v <- arg[2]
    log_s_1 <- log(tau_1^2) + log(sum_a_1) + (gamma_1*log(u))
    log_s_2 <- log(tau_2^2) + log(sum_a_2) + (gamma_2*log(v))
    s_1 <- exp(log_s_1)
    s_2 <- exp(log_s_2)
    log_surv <- -(b_1*log(1 + s_1 + s_2)) + (b_2*log(1 + s_1)) +
      (b_3*log(1 + s_2))
    log_c_1 <- log(b_1) + log(1 + b_1) - (2*log(1 + s_1 + s_2))
    log_ingd_1 <- log_constant + log_surv + ((gamma_1 - 1)*log(u)) +
      ((gamma_2 - 1)*log(v)) + log_c_1
    ingd_1 <- exp(log_ingd_1)
    return(ingd_1)
  }
  integrand_2 <- function(arg)
  {
    u <- arg[1]
    v <- arg[2]
    log_s_1 <- log(tau_1^2) + log(sum_a_1) + (gamma_1*log(u))
    log_s_2 <- log(tau_2^2) + log(sum_a_2) + (gamma_2*log(v))
    s_1 <- exp(log_s_1)
    s_2 <- exp(log_s_2)
    log_surv <- -(b_1*log(1 + s_1 + s_2)) + (b_2*log(1 + s_1)) +
      (b_3*log(1 + s_2))
    log_c_2 <- log(b_1) + log(-b_2) - log(1 + s_1 + s_2) - log(1 + s_1)
    log_ingd_2 <- log_constant + log_surv + ((gamma_1 - 1)*log(u)) +
      ((gamma_2 - 1)*log(v)) + log_c_2
    ingd_2 <- exp(log_ingd_2)
    return(ingd_2)
  }
  integrand_3 <- function(arg)
  {
    u <- arg[1]
    v <- arg[2]
    log_s_1 <- log(tau_1^2) + log(sum_a_1) + (gamma_1*log(u))
    log_s_2 <- log(tau_2^2) + log(sum_a_2) + (gamma_2*log(v))
    s_1 <- exp(log_s_1)
    s_2 <- exp(log_s_2)
    log_surv <- -(b_1*log(1 + s_1 + s_2)) + (b_2*log(1 + s_1)) +
      (b_3*log(1 + s_2))
    log_c_3 <- log(b_1) + log(-b_3) - log(1 + s_1 + s_2) - log(1 + s_2)
    log_ingd_3 <- log_constant + log_surv + ((gamma_1 - 1)*log(u)) +
      ((gamma_2 - 1)*log(v)) + log_c_3
    ingd_3 <- exp(log_ingd_3)
    return(ingd_3)
  }
  integrand_4 <- function(arg)
  {
    u <- arg[1]
    v <- arg[2]
    log_s_1 <- log(tau_1^2) + log(sum_a_1) + (gamma_1*log(u))
    log_s_2 <- log(tau_2^2) + log(sum_a_2) + (gamma_2*log(v))
    s_1 <- exp(log_s_1)
    s_2 <- exp(log_s_2)
    log_surv <- -(b_1*log(1 + s_1 + s_2)) + (b_2*log(1 + s_1)) +
      (b_3*log(1 + s_2))
    log_c_4 <- (log(-b_2) + log(-b_3)) - log(1 + s_1) - log(1 + s_2)
    log_ingd_4 <- log_constant + log_surv + ((gamma_1 - 1)*log(u)) +
      ((gamma_2 - 1)*log(v)) + log_c_4
    ingd_4 <- exp(log_ingd_4)
    return(ingd_4)
  }
  h_1 <- cubature::cuhre(integrand_1,lowerLimit=c(0,0),upperLimit=c(t_1,t_2))$integral
  h_2 <- cubature::cuhre(integrand_2,lowerLimit=c(0,0),upperLimit=c(t_1,t_2))$integral
  h_3 <- cubature::cuhre(integrand_3,lowerLimit=c(0,0),upperLimit=c(t_1,t_2))$integral
  h_4 <- cubature::cuhre(integrand_4,lowerLimit=c(0,0),upperLimit=c(t_1,t_2))$integral
  h <- h_1 + h_2 + h_3 + h_4
  eps <- 10^(-15)
  if(h > eps)
  {
    return(h)
  }else
  {
    return(eps)
  }
}

## log-likelihood of an individual for Weibull cause specific hazard
## correlated Gamma frailty model
wei_indiv_log_lik_corr <- function(data,mu_vec_1,mu_vec_2,gamma_1,
                                   gamma_2,tau_1,tau_2,phi,q_1,q_2)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  a_1 <- mu_vec_1^gamma_1
  a_2 <- mu_vec_2^gamma_2
  s_1 <- (tau_1^2)*sum(a_1)*(x^gamma_1)
  s_2 <- (tau_2^2)*sum(a_2)*(x^gamma_2)
  b_1 <- phi/(tau_1*tau_2)
  b_2 <- b_1 - 1/(tau_1^2)
  b_3 <- b_1 - 1/(tau_2^2)
  log_s <- -(b_1*log(1 + s_1 + s_2)) + (b_2*log(1 + s_1)) +
    (b_3*log(1 + s_2))
  marg_surv_1 <- (1 + s_1)^(-(1/(tau_1^2)))
  marg_surv_2 <- (1 + s_2)^(-(1/(tau_2^2)))
  joint_surv <- exp(log_s)
  eps <- 10^(-15)
  if((j_1 != 0) && (j_2 != 0))
  {
    h <- joint_subdis_wei_corr(x,x,j_1,j_2,mu_vec_1,mu_vec_2,gamma_1,gamma_2,tau_1,tau_2,phi)
    if(h > 10^(-15))
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }else if((j_1 != 0) && (j_2 == 0)){
    h_1 <- q_1[j_1]*marg_surv_wei_corr(x,gamma_1,mu_vec_1,tau_1)
    h_2 <- sum(sapply(1:3,joint_subdis_wei_corr,t_1=x,t_2=x,j_1 = j_1,mu_vec_1=mu_vec_1,
                      mu_vec_2=mu_vec_2,gamma_1=gamma_1,gamma_2=gamma_2,phi=phi,tau_1=tau_1,tau_2=tau_2))
    if((h_1 - h_2) > eps)
    {
      return(log(h_1 - h_2))
    }else
    {
      return(log(eps))
    }
  }else if((j_1 == 0) && (j_2 != 0)){
    h_1 <- q_2[j_2]*marg_surv_wei_corr(x,gamma_2,mu_vec_2,tau_2)
    h_2 <- sum(sapply(1:3,joint_subdis_wei_corr,t_1=x,t_2=x,j_2 = j_2,mu_vec_1=mu_vec_1,
                      mu_vec_2=mu_vec_2,gamma_1=gamma_1,gamma_2=gamma_2,phi=phi,tau_1=tau_1,tau_2=tau_2))
    if((h_1 - h_2) > eps)
    {
      return(log(h_1 - h_2))
    }else
    {
      return(log(eps))
    }
  }else{
    return(log_s)
  }
}

## the overall negative log-likelihood function which is
## to minimize
neg_log_lik_wei_corr <- function(theta,D,M)
{
  lambda_vec_3 <- theta[1:M]
  lambda_vec_4 <- theta[(M + 1):(2*M)]
  beta_1 <- theta[(2*M) + 1]
  beta_2 <- theta[(2*M) + 2]
  tau_1 <- theta[(2*M) + 3]
  tau_2 <- theta[(2*M) + 4]
  phi <- theta[(2*M) + 5]
  lambda_3 <- lambda_vec_3^beta_1
  lambda_4 <- lambda_vec_4^beta_2
  q_4 <- lambda_3/sum(lambda_3)
  q_5 <- lambda_4/sum(lambda_4)
  return(-sum(apply(D,1,wei_indiv_log_lik_corr,mu_vec_1=lambda_vec_3,
                    mu_vec_2=lambda_vec_4,q_1=q_4,q_2=q_5,gamma_1=beta_1,
                    gamma_2=beta_2,tau_1=tau_1,tau_2=tau_2,phi=phi)))
}
### parallel version
neg_log_lik_wei_corr_par <- function(theta,D,M,cl)
{
  lambda_vec_3 <- theta[1:M]
  lambda_vec_4 <- theta[(M + 1):(2*M)]
  beta_1 <- theta[(2*M) + 1]
  beta_2 <- theta[(2*M) + 2]
  tau_1 <- theta[(2*M) + 3]
  tau_2 <- theta[(2*M) + 4]
  phi <- theta[(2*M) + 5]
  lambda_3 <- lambda_vec_3^beta_1
  lambda_4 <- lambda_vec_4^beta_2
  q_4 <- lambda_3/sum(lambda_3)
  q_5 <- lambda_4/sum(lambda_4)
  return(-sum(parApply(cl,D,1,wei_indiv_log_lik_corr,mu_vec_1=lambda_vec_3,
                       mu_vec_2=lambda_vec_4,q_1=q_4,q_2=q_5,gamma_1=beta_1,
                       gamma_2=beta_2,tau_1=tau_1,tau_2=tau_2,phi=phi)))
}
library("doParallel")
cl <- makeCluster(25)
setDefaultCluster(cl=cl)
clusterExport(cl,c("joint_surv_wei_corr","soln_1","soln_2",
                   "marg_surv_wei_corr","joint_subdis_wei_corr"))
neg_log_lik_wei_corr_par(c(c(0.0397,0.0327,0.0336,0.0445,0.0345,0.0361),
                           12.834,12.78,4.269,4.259,0.99),hearing_data,3,cl)
#neg_log_lik_wei_corr_par(c(initial_value_1,4.309,4.259,0.98),hearing_data,3,cl)
stopCluster(cl)
initial_value <- c(0.0397,0.0327,0.0336,0.0445,0.0345,0.0361,13.36,12.36,4.309,4.209,0.97)
initial_value_1 <- c(0.0397,0.0327,0.0336,0.0445,0.0345,0.0361,13.36,12.36)
initial_value_2 <- c(4.309,4.209,0.97)
### objective function for fixed value of variance and rho parameter
neg_log_lik_wei_corr_1 <- function(theta,D,M)
{
  lambda_vec_3 <- theta[1:M]
  lambda_vec_4 <- theta[(M + 1):(2*M)]
  beta_1 <- theta[(2*M) + 1]
  beta_2 <- theta[(2*M) + 2]
  lambda_3 <- lambda_vec_3^beta_1
  lambda_4 <- lambda_vec_4^beta_2
  q_4 <- lambda_3/sum(lambda_3)
  q_5 <- lambda_4/sum(lambda_4)
  return(-sum(apply(D,1,wei_indiv_log_lik_corr,mu_vec_1=lambda_vec_3,
                       mu_vec_2=lambda_vec_4,q_1=q_4,q_2=q_5,gamma_1=beta_1,
                       gamma_2=beta_2,tau_1=soln_2[1],tau_2=soln_2[2],phi=soln_2[3])))
}
## parallel version for fixed value of baseline hazard parameters
### parallel version
neg_log_lik_wei_corr_par_2 <- function(theta,D,M,cl)
{
  tau_1 <- theta[1]
  tau_2 <- theta[2]
  phi <- theta[3]
  lambda_3 <- lambda_vec_3^beta_1
  lambda_4 <- lambda_vec_4^beta_2
  q_4 <- lambda_3/sum(lambda_3)
  q_5 <- lambda_4/sum(lambda_4)
  return(-sum(parApply(cl,D,1,wei_indiv_log_lik_corr,mu_vec_1=lambda_vec_3,
                       mu_vec_2=lambda_vec_4,q_1=q_4,q_2=q_5,gamma_1=beta_1,
                       gamma_2=beta_2,tau_1=tau_1,tau_2=tau_2,phi=phi)))
}
cl <- makeCluster(17)
setDefaultCluster(cl=cl)
clusterExport(cl,c("wei_indiv_log_lik_corr","joint_surv_wei_corr","x",
                   "marg_surv_wei_corr","joint_subdis_wei_corr","soln_2"))
data_soln_3 <- optimParallel::optimParallel(par=soln_1,
fn=neg_log_lik_wei_corr_1,M=3,D=hearing_data,lower=replicate(8,0.01),upper=replicate(8,Inf),
              parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE),control = list(trace=6,maxit=300),
                                          method = "L-BFGS-B")$par
stopCluster(cl)
setDefaultCluster(cl=NULL)
soln_3 <- data_soln_3
cl <- makeCluster(25)
setDefaultCluster(cl=cl)
clusterExport(cl,c("wei_indiv_log_lik_corr","joint_surv_wei_corr","x",
                   "marg_surv_wei_corr","joint_subdis_wei_corr"))
library(alabama)
neg_log_lik_wei_corr_par_2 <- function(theta,D,M,cl)
{
  tau_1 <- theta[1]
  tau_2 <- theta[2]
  phi <- theta[3]
  lambda_vec_3 <- soln_1[1:3]
  lambda_vec_4 <- soln_1[4:6]
  beta_1 <- soln_1[7]
  beta_2 <- soln_1[8]
  lambda_3 <- lambda_vec_3^beta_1
  lambda_4 <- lambda_vec_4^beta_2
  q_4 <- lambda_3/sum(lambda_3)
  q_5 <- lambda_4/sum(lambda_4)
  if(min(tau_1/tau_2,tau_2/tau_1) > phi)
  {
    return(-sum(parApply(cl,D,1,wei_indiv_log_lik_corr,mu_vec_1=lambda_vec_3,
                         mu_vec_2=lambda_vec_4,q_1=q_4,q_2=q_5,gamma_1=beta_1,
                         gamma_2=beta_2,tau_1=tau_1,tau_2=tau_2,phi=phi)))
  }else
  {
    return(10^8)
  }
}
objfun=function(theta)
{
  return(neg_log_lik_wei_corr_par_2(theta,hearing_data,3,cl))
}
hin <- function(theta)
{
  return(c(theta[1],theta[2],theta[3],
           min((theta[1]/theta[2]),(theta[2]/theta[1]))-theta[3]))
}
num_gr <- function(theta)
{
  return(pracma::grad(neg_log_lik_wei_corr_par_2,theta,
                        D=hearing_data,M=3,cl=cl))
}
#soln_2 <- alabama::constrOptim.nl(soln_2,objfun,hin)$par
soln_2 <- nloptr::slsqp(x0=soln_2,fn=objfun,hin=hin)$par
stopCluster(cl)
setDefaultCluster(cl=NULL)
neg_log_lik_wei_corr_2 <- function(theta,D,M)
{
  tau_1 <- theta[1]
  tau_2 <- theta[2]
  phi <- theta[3]
  lambda_vec_3 <- soln_1[1:3]
  lambda_vec_4 <- soln_1[4:6]
  beta_1 <- soln_1[7]
  beta_2 <- soln_1[8]
  lambda_3 <- lambda_vec_3^beta_1
  lambda_4 <- lambda_vec_4^beta_2
  q_4 <- lambda_3/sum(lambda_3)
  q_5 <- lambda_4/sum(lambda_4)
  if(min(tau_1/tau_2,tau_2/tau_1) > phi)
  {
    return(-sum(apply(D,1,wei_indiv_log_lik_corr,mu_vec_1=lambda_vec_3,
                      mu_vec_2=lambda_vec_4,q_1=q_4,q_2=q_5,gamma_1=beta_1,
                      gamma_2=beta_2,tau_1=tau_1,tau_2=tau_2,phi=phi)))
  }else
  {
    return(10^9)
  }
}
cl <- makeCluster(7)
setDefaultCluster(cl=cl)
clusterExport(cl,c("wei_indiv_log_lik_corr","joint_surv_wei_corr","x",
                   "marg_surv_wei_corr","joint_subdis_wei_corr","soln_1"))
soln_2 <- optimParallel::optimParallel(par=initial_value_2,
fn=neg_log_lik_wei_corr_2,M=3,D=hearing_data,lower=replicate(3,0.01),upper=c(replicate(2,Inf),0.99),
    parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE),control = list(trace=6,maxit=300),
                                            method = "L-BFGS-B")$par
stopCluster(cl)
soln <- c(soln_1,soln_2)
library("doParallel")
cl <- makeCluster(25)
setDefaultCluster(cl=cl)
clusterExport(cl,c("joint_surv_wei_corr","soln_1","soln_2",
                   "marg_surv_wei_corr","joint_subdis_wei_corr"))
hess_mle <- pracma::hessian(neg_log_lik_wei_corr_par,c(soln_1,4.2645,4.2325,0.99),D=hearing_data,cl=cl,M=3)
stopCluster(cl)
### estimate of standard error
se_mle <- abs(diag(solve(hess_mle)))
se_mle <- sqrt(se_mle)
se_mle
