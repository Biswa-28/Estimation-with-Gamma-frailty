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
## function to compute joint survival function
## function to compute joint survival function for model with equal variance
joint_surv_weibull_hazard_corr_cause_specific <- function(t_1,t_2,mu_vec_1,
                                      mu_vec_2,gamma_1,gamma_2,tau_vec,rho_vec)
{
  log_tau_vec <- log(tau_vec)
  log_a_1 <- log(rho_vec) - (2*log_tau_vec)
  a_1 <- exp(log_a_1)
  a_2 <- a_1 - exp(-2*log_tau_vec)
  b_1 <- 1 + ((t_1*mu_vec_1)^gamma_1)*(tau_vec^2)
  b_2 <- 1 + ((t_2*mu_vec_2)^gamma_2)*(tau_vec^2)
  b_3 <- b_1 + b_2 - 1
  log_joint_surv <- -a_1*log(b_3) + (a_2*log(b_1)) + (a_2*log(b_2))
  joint_surv <- exp(sum(log_joint_surv))
  return(joint_surv)
}
## marginal distribution function calculation
marg_subdis_weibull_corr_cause_specific <- function(t,j,mu_vec,gamma,tau_vec)
{
  return(SharedAndCorrCauseSpecific::marg_subdis_weibull_corr_cause_specific_cpp(t,j,mu_vec,gamma,tau_vec))
}
## joint sub-distribution function calculation for the model with equal variance
joint_sub_dis_weibull_corr_cause_specific <- function(t_1,t_2,j_1,j_2,mu_vec_1,mu_vec_2,
                                                  gamma_1,gamma_2,tau_vec,rho_vec)
{
  return(SharedAndCorrCauseSpecific::joint_sub_dis_weibull_corr_cause_specific_cpp(t_1,t_2,j_1,j_2,
              gamma_1,gamma_2,mu_vec_1,mu_vec_2,tau_vec,rho_vec,3)$approximate)
}
## log-likelihood function for an individual for the model with equal variance
weibull_corr_cause_indiv_log_lik <- function(data,mu_vec_1,mu_vec_2,
                                         gamma_1,gamma_2,tau_vec,rho_vec)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  eps <- 10^(-25)
  if((j_1 != 0) && (j_2 != 0))
  {
    h <- joint_sub_dis_weibull_corr_cause_specific(x,x,j_1,j_2,mu_vec_1,mu_vec_2,
                                          gamma_1,gamma_2,tau_vec,rho_vec)
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }else if((j_1 != 0) && (j_2 == 0)){
    h <- marg_subdis_weibull_corr_cause_specific(x,j_1,mu_vec_1,gamma_1,tau_vec) -
      sum(sapply(1:3,joint_sub_dis_weibull_corr_cause_specific,t_1=x,t_2=x,j_1=j_1,
      mu_vec_1 = mu_vec_1,mu_vec_2 = mu_vec_2,gamma_1=gamma_1,gamma_2=gamma_2,tau_vec=tau_vec,
                 rho_vec=rho_vec))
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }else if((j_1 == 0) && (j_2 != 0)){
    h <- marg_subdis_weibull_corr_cause_specific(x,j_2,mu_vec_2,gamma_2,tau_vec) -
      sum(sapply(1:3,joint_sub_dis_weibull_corr_cause_specific,t_1=x,t_2=x,j_2=j_2,
      mu_vec_1 = mu_vec_1,mu_vec_2 = mu_vec_2,gamma_1=gamma_1,gamma_2=gamma_2,tau_vec=tau_vec,
                 rho_vec=rho_vec))
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }else{
    h <- joint_surv_weibull_hazard_corr_cause_specific(x,x,mu_vec_1,mu_vec_2,
                                             gamma_1,gamma_2,tau_vec,rho_vec)
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }
}
## negative log-likelihood function for reduced model which we want to minimize
neg_log_lik_weibull_corr_cause_specific <- function(theta,D,M)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  gamma_3 <- theta[(2*M) + 1]
  gamma_4 <- theta[(2*M) + 2]
  tau_vec_1 <- theta[((2*M) + 3):((3*M) + 2)]
  rho_vec_1 <- theta[((3*M) + 3):((4*M) + 2)]
  b <- -apply(D,1,weibull_corr_cause_indiv_log_lik,
  mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,gamma_1=gamma_3,gamma_2=gamma_4,
              tau_vec=tau_vec_1,rho_vec=rho_vec_1)
  return(sum(b))
}
library("doParallel")
### parallel version
neg_log_lik_weibull_corr_cause_specific_par <- function(theta,D,M,cl)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  gamma_3 <- theta[(2*M) + 1]
  gamma_4 <- theta[(2*M) + 2]
  tau_vec_1 <- theta[((2*M) + 3):((3*M) + 2)]
  rho_vec_1 <- theta[((3*M) + 3):((4*M) + 2)]
  b <- -parApply(cl,D,1,weibull_corr_cause_indiv_log_lik,
                 mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,gamma_1=gamma_3,gamma_2=gamma_4,
                 tau_vec=tau_vec_1,rho_vec=rho_vec_1)
  return(sum(b))
}
# finding a better initial point for minimization of original
## neg-loglikelihood function
data_soln_exp_corr_cause_specific <- c(37.599,0.02847,0.00255,38.5983,0.024056,0.00262,
                                       3.59981,10.6776,2.4805,0.999,0.999,0.999)
l <- seq(0.1,2,length=15)
in_est_2d_array <- array(0,dim=c(15,15))
c <- 0
library("doParallel")
cl <- makeCluster(35)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
"joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik",
"SharedAndCorrCauseSpecific"))
for(i in 1:length(l))
{
  for(j in 1:length(l))
  {
    in_est_2d_array[i,j] <- neg_log_lik_weibull_corr_cause_specific_par(c(data_soln_exp_corr_cause_specific[1:6],
    l[i],l[j],data_soln_exp_corr_cause_specific[7:12]),hearing_data,3,cl)
    c <- c + 1
  }
}
stopCluster(cl)
min_pos <- which.min(in_est_2d_array) ## 129 i.e c(l[8],l[9])
#### the first initial value
initial_reduced_1 <- c(data_soln_exp_corr_cause_specific[1:6],c(l[8],l[9]),
               data_soln_exp_corr_cause_specific[7:12])
### computation for reduced model
library("doParallel")
cl <- makeCluster(29)
setDefaultCluster(cl=cl)
clusterEvalQ(cl,{
  library(SharedAndCorrCauseSpecific)
  })
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
"joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik"))
initial_reduced_1 <- c(data_soln_exp_corr_cause_specific[1:6],c(l[8],l[9]),
                       data_soln_exp_corr_cause_specific[7:12])
data_soln_1 <- optimParallel::optimParallel(initial_reduced_1,
fn=neg_log_lik_weibull_corr_cause_specific,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(14,0.0001),upper = c(replicate(11,Inf),replicate(3,0.999)),
control = list(trace=6,ndeps=replicate(14,1e-5),maxit=600),parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
### checking hessian is p.d or not
data_reduced_soln_1 <- data_soln_1
data_soln_1_par <- data_soln_1$par
hess_reduced_1 <- data_reduced_soln_1$hessian
eigen(hess_reduced_1)$values  # 13 positive and 1 very large negative
cl <- makeCluster(29)
setDefaultCluster(cl=cl)
clusterEvalQ(cl,{
  library(SharedAndCorrCauseSpecific)
})
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
                   "joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik"))
#initial_3_1 <- c(data_soln_2_par[1:6],60,data_soln_2_par[8:14])
#initial_3_2 <- c(data_soln_2_par[1:6],67,data_soln_2_par[8:14]) #357.5323
initial_3_3 <- c(initial_3_2[1:7],67.5,initial_3_2[9:14]) #257.0661
initial_3_4 <- c(initial_3_3[1:10],28.561677,initial_3_3[12:14]) #248.2387
initial_3_5 <- c(initial_3_4[1:9],81.1465422,initial_3_4[11:14])
initial_3_6 <- c(initial_3_4[1:8],52.3890582,initial_3_4[10:14]) #222.7115
neg_log_lik_weibull_corr_cause_specific_par(initial_3_6,hearing_data,3,cl)
pracma::grad(neg_log_lik_weibull_corr_cause_specific_par,initial_3_6,
             D=hearing_data,M=3,cl=cl)
stopCluster(cl)
library("doParallel")
cl <- makeCluster(29)
setDefaultCluster(cl=cl)
clusterEvalQ(cl,{
  library(SharedAndCorrCauseSpecific)
})
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
                   "joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik"))
data_soln_2 <- optimParallel::optimParallel(initial_3_6,
fn=neg_log_lik_weibull_corr_cause_specific,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(14,0.0001),upper = c(replicate(11,Inf),replicate(3,0.999)),
control = list(trace=6,ndeps=replicate(14,1e-5),maxit=600),parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
### checking hessian is p.d or not
data_soln_2_par <- data_soln_2$par     ### neg-lik 149.833
hess_reduced_2 <- data_soln_2$hessian
eigen(hess_reduced_2)$values
####################_______________________######################
##other tries
cl <- makeCluster(35)
setDefaultCluster(cl=cl)
clusterEvalQ(cl,{
  library(SharedAndCorrCauseSpecific)
})
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
                   "joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik"))
initial_4_1 <- c(data_soln_2_par[1:11],c(0.9859,0.97393,0.94256)) ##149.4872
initial_4_2 <- c(data_soln_2_par[1:11],c(0.9859,0.98393,0.94256))
neg_log_lik_weibull_corr_cause_specific_par(initial_4_2,hearing_data,3,cl)
#pracma::grad(neg_log_lik_weibull_corr_cause_specific_par,initial_4_1,
 #            D=hearing_data,M=3,cl=cl)
H_4 <- pracma::hessian(neg_log_lik_weibull_corr_cause_specific_par,initial_4_1,
                D=hearing_data,M=3,cl=cl)
eigen(H_4)$values
stopCluster(cl)
cl <- makeCluster(29)
setDefaultCluster(cl=cl)
clusterEvalQ(cl,{
  library(SharedAndCorrCauseSpecific)
})
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
                   "joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik"))
data_soln_3 <- optimParallel::optimParallel(initial_4_1,
fn=neg_log_lik_weibull_corr_cause_specific,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(14,0.0001),upper = c(replicate(11,Inf),replicate(3,0.999)),
control = list(trace=6,ndeps=replicate(14,1e-6),maxit=600),parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
### checking hessian is p.d or not
data_soln_3_par <- data_soln_3$par     ### neg-lik 116.532
hess_reduced_3 <- data_soln_3$hessian
eigen(hess_reduced_3)$values
###################___________________________####################
cl <- makeCluster(29)
setDefaultCluster(cl=cl)
clusterEvalQ(cl,{
  library(SharedAndCorrCauseSpecific)
})
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
                   "joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik"))
ini_2 <- data_soln_1$par
ini_3 <- c(ini_2[1:3],38.7874,ini_2[5:14])
ini_4 <- c(ini_3[1:13],0.9918)
neg_log_lik_weibull_corr_cause_specific_par(ini_4,hearing_data,3,cl)
#pracma::grad(neg_log_lik_weibull_corr_cause_specific_par,ini_4,D=hearing_data,
             #M=3,cl=cl)
#hess_test_1 <- pracma::hessian(neg_log_lik_weibull_corr_cause_specific_par,ini_3,D=hearing_data,
                            #M=3,cl=cl)
#eigen(hess_test_1)$values
data_soln_2 <- optimParallel::optimParallel(ini_4,
fn=neg_log_lik_weibull_corr_cause_specific,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(14,0.0001),upper = c(replicate(11,Inf),replicate(3,0.9999)),
control = list(trace=6,ndeps=replicate(14,1e-5),maxit=300),parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
### checking hessian is p.d or not
data_reduced_soln_2 <- data_soln_2$par
hess_reduced_2 <- data_soln_2$hessian
eigen(hess_reduced_2)$values
initial_3 <- c(data_reduced_soln_2[1:8],60.83123,data_reduced_soln_2[10:14])
cl <- makeCluster(29)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
                   "joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik"))
neg_log_lik_weibull_corr_cause_specific_par(c(initial_3[1:7],62.8497,initial_3[9:13],0.999),hearing_data,3,cl)
data_soln_3 <- optimParallel::optimParallel(initial_3,
fn=neg_log_lik_weibull_corr_cause_specific,D = hearing_data,Freq_table=F,hessian=TRUE,
M = 3,method="L-BFGS-B",lower = replicate(14,0.0001),upper = c(replicate(11,Inf),replicate(3,0.9999)),
control = list(trace=6,ndeps=replicate(14,1e-5),maxit=300),parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
#pracma::grad(neg_log_lik_weibull_corr_cause_specific_par,c(data_reduced_soln_2[1:8],60.83123,data_reduced_soln_2[10:14]),
             #D=hearing_data,M=3,cl=cl)
#hess_original_2 <- pracma::hessian(f=neg_log_lik_weibull_corr_cause_specific_par,
#x0=c(data_reduced_soln_2[1:8],60.83123,data_reduced_soln_2[10:14]),cl=cl,D=hearing_data,M=3)
stopCluster(cl)
### checking hessian is p.d or not
data_reduced_soln_3 <- data_soln_3$par
hess_reduced_3 <- data_soln_3$hessian
eigen(hess_reduced_3)$values
cl <- makeCluster(25)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
"joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik","data_soln_2"))
clusterEvalQ(cl,{
  data_soln_2_par <- data_soln_2$par
  data_soln_2_par_1 <- data_soln_2_par[1:8]
  data_soln_2_par_2 <- data_soln_2_par[9:14]
})
### Now applying it on the parallel version of the negative log-likelihood
################# performing a profile likelihood approach ####################
############# theta_1 <- c(lambda_vec_1,lambda_vec_2,beta_1,beta_2)
############# theta_2 <- c(tau_vec,rho_1_1,rho_2_2,rho_3_3)
neg_log_lik_weibull_corr_cause_specific_1 <- function(theta,D,M)
{
  tau_vec_1 <- theta[1:M]
  rho_vec_1 <- theta[(M + 1):(2*M)]
  c <- 0
  for(i in 1:(2*M))
  {
    if(theta[i] > 0)
    {
      c <- c + 1
    }
  }
  if(c == (2*M))
  {
    b <- -apply(D,1,weibull_corr_cause_indiv_log_lik,
                   mu_vec_1=data_soln_2_par_1[1:3],mu_vec_2=data_soln_2_par_1[4:6],
                   gamma_1=data_soln_2_par_1[7],gamma_2=data_soln_2_par_1[8],
                   tau_vec=tau_vec_1,rho_vec=rho_vec_1)
    return(sum(b))
  }else
  {
    return(10^19)
  }
}
neg_log_lik_weibull_corr_cause_specific_par_1 <- function(theta,D,M,cl)
{
  tau_vec_1 <- theta[1:M]
  rho_vec_1 <- theta[(M + 1):(2*M)]
  c <- 0
  for(i in 1:M)
  {
    if(theta[i] > 0)
    {
      c <- c + 1
    }
  }
  for(i in (M + 1):(2*M))
  {
    if((theta[i] < 1)&(theta[i] > 0))
    {
      c <- c + 1
    }
  }
  if(c == (2*M))
  {
    b <- -parApply(cl,D,1,weibull_corr_cause_indiv_log_lik,
    mu_vec_1=data_soln_2_par_1[1:3],mu_vec_2=data_soln_2_par_1[4:6],
    gamma_1=data_soln_2_par_1[7],gamma_2=data_soln_2_par_1[8],
    tau_vec=tau_vec_1,rho_vec=rho_vec_1)
    return(sum(b))
  }else
  {
    return(10^19)
  }
}
objfun_Weibull_1 <- function(theta_1)
{
  return(neg_log_lik_weibull_corr_cause_specific_par_1(theta_1,hearing_data,3,cl))
}
objfun_Weibull_1(data_soln_2_par_2)
gr_objfun_Weibull_1 <- function(theta_1)
{
  return(pracma::grad(objfun_Weibull_1,theta_1))
}
#gr_objfun_Weibull_1(data_soln_2_par_2)
hin_1 <- function(theta_1)
{
  h <- rep(NA, 1)
  h[1] <- theta_1[1]
  h[2] <- theta_1[2]
  h[3] <- theta_1[3]
  h[4] <- 1 - theta_1[4]
  h[5] <- 1 - theta_1[5]
  h[6] <- 1 - theta_1[6]
  #h[7] <- theta_1[7]
  #h[8] <- theta_1[8]
  return(h)
}
hin_1.jac <- function(theta_1)
{
  j <- matrix(NA, 6, length(theta_1))
  j[1,] <- c(1,rep(0,5))
  j[2,] <- c(0,1,rep(0,4))
  j[3,] <- c(0,0,1,rep(0,3))
  j[4,] <- c(rep(0,3),-1,rep(0,2))
  j[5,] <- c(rep(0,4),-1,0)
  j[6,] <- c(rep(0,5),-1)
  j
}
hin_1(data_soln_2_par_2)
#final_mle_Weibull_1 <- NlcOptim::solnl(data_soln_2_par_2,objfun=objfun_Weibull_1,
#lb=replicate(6,0.0001),ub=c(replicate(3,Inf),replicate(3,0.999)))
final_mle_Weibull_1 <- alabama::auglag(par=data_soln_2_par_2,fn=objfun_Weibull_1,
hin=hin_1,hin.jac=hin_1.jac,control.outer = list(trace=TRUE),control.optim = list(trace=6,ndeps=rep(6,1e-5)))
#final_mle_Weibull_1 <- optimParallel::optimParallel(data_soln_2_par_1,
#neg_log_lik_weibull_corr_cause_specific_1,D=hearing_data,M=3,hessian=TRUE,lower=replicate(8,0.0001),upper=replicate(8,Inf),
#control = list(trace=6,ndeps=replicate(8,1e-6),maxit=300),parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
final_mle_Weibull_par_1 <- final_mle_Weibull_1$par
#final_mle_Weibull_par_1 <- as.vector(final_mle_Weibull_par_1)
#sub_hess_1 <- final_mle_Weibull_1$hessian
cl <- makeCluster(25)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
"joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik","data_soln_2",
"final_mle_Weibull_1"))
clusterEvalQ(cl,{
  data_soln_2_par <- data_soln_2$par
  data_soln_2_par_1 <- data_soln_2_par[1:8]
  final_mle_Weibull_par_1 <- final_mle_Weibull_1$par
  #final_mle_Weibull_par_1 <- as.vector(final_mle_Weibull_par_1)
})
neg_log_lik_weibull_corr_cause_specific_2 <- function(theta,D,M)
{
  tau_vec_1 <- theta[1:M]
  rho_vec_1 <- theta[(M + 1):(2*M)]
  b <- -apply(D,1,weibull_corr_cause_indiv_log_lik,
      mu_vec_1=final_mle_Weibull_par_1[1:3],mu_vec_2=final_mle_Weibull_par_1[4:6],
      gamma_1=final_mle_Weibull_par_1[7],gamma_2=final_mle_Weibull_par_1[8],
      tau_vec=final_mle_Weibull_par_1[1:3],rho_vec=final_mle_Weibull_par_1[4:6])
  return(sum(b))
}
neg_log_lik_weibull_corr_cause_specific_par_2 <- function(theta,D,M,cl)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  gamma_3 <- theta[(2*M) + 1]
  gamma_4 <- theta[(2*M) + 2]
  c <- 0
  for(i in 1:((2*M) + 2))
  {
    if(theta[i] > 0)
    {
      c <- c + 1
    }
  }
  if(c == ((2*M) + 2))
  {
    b <- -parApply(cl,D,1,weibull_corr_cause_indiv_log_lik,
    mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,gamma_1=gamma_3,gamma_2=gamma_4,
    tau_vec=final_mle_Weibull_par_1[1:3],rho_vec=final_mle_Weibull_par_1[4:6])
    return(sum(b))
  }else
  {
    return(10^19)
  }
}
#final_mle_Weibull_2 <- optimParallel::optimParallel(data_soln_2_par_2,
#neg_log_lik_weibull_corr_cause_specific_2,hessian=TRUE,lower=replicate(6,0.0001),upper=c(replicate(3,Inf),replicate(3,0.999)),
#control = list(trace=6,ndeps=replicate(6,1e-6),maxit=300),parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
#stopCluster(cl)
objfun_Weibull_2 <- function(theta_2)
{
  return(neg_log_lik_weibull_corr_cause_specific_par_2(theta_2,hearing_data,3,cl))
}
objfun_Weibull_2(data_soln_2_par_1)
hin_2 <- function(theta_2)
{
  h <- rep(NA, 1)
  h[1] <- theta_2[1]
  h[2] <- theta_2[2]
  h[3] <- theta_2[3]
  h[4] <- theta_2[4]
  h[5] <- theta_2[5]
  h[6] <- theta_2[6]
  h[7] <- theta_2[7]
  h[8] <- theta_2[8]
  return(h)
}
hin_2.jac <- function(theta_2)
{
  j <- matrix(NA, 8, length(theta_2))
  j <- diag(8)
  j
}
hin_2(data_soln_2_par_1)
#final_mle_Weibull_1 <- NlcOptim::solnl(data_soln_2_par_1,objfun=objfun_Weibull_1,
#lb=replicate(8,0.0001),ub=replicate(8,Inf))
#final_mle_Weibull_2 <- NlcOptim::solnl(data_soln_2_par_1,objfun=objfun_Weibull_2,
#lb=replicate(8,0.0001),ub=replicate(8,Inf))
final_mle_Weibull_2 <- alabama::auglag(par=data_soln_2_par_1,fn=objfun_Weibull_2,
hin=hin_2,hin.jac=hin_2.jac,control.outer = list(trace=TRUE),control.optim = list(trace=6,ndeps=rep(8,1e-5)))
stopCluster(cl)
final_mle_Weibull_par_2 <- final_mle_Weibull_2$par
sub_hess_2 <- final_mle_Weibull_2$hessian
##### final_mle_Weibull_par_2 is m.l.e of cause specific hazard parameters  ####
## Now to find m.l.e of frailty parameters using this
cl <- makeCluster(25)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
"joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik","final_mle_Weibull_1","final_mle_Weibull_2"))
clusterEvalQ(cl,{
  final_mle_Weibull_par_1 <- final_mle_Weibull_1$par
  final_mle_Weibull_par_2 <- final_mle_Weibull_2$par
})
### Now applying it on the parallel version of the negative log-likelihood
################# performing a profile likelihood approach ####################
############# theta_1 <- c(lambda_vec_1,lambda_vec_2,beta_1,beta_2)
############# theta_2 <- c(tau_vec,rho_1_1,rho_2_2,rho_3_3)
neg_log_lik_weibull_corr_cause_specific_1 <- function(theta,D,M)
{
  tau_vec_1 <- theta[1:M]
  rho_vec_1 <- theta[(M + 1):(2*M)]
  c <- 0
  for(i in 1:(2*M))
  {
    if(theta[i] > 0)
    {
      c <- c + 1
    }
  }
  if(c == (2*M))
  {
    b <- -apply(D,1,weibull_corr_cause_indiv_log_lik,
                mu_vec_1=final_mle_Weibull_par_2[1:3],mu_vec_2=final_mle_Weibull_par_2[4:6],
                gamma_1=final_mle_Weibull_par_2[7],gamma_2=final_mle_Weibull_par_2[8],
                tau_vec=tau_vec_1,rho_vec=rho_vec_1)
    return(sum(b))
  }else
  {
    return(10^19)
  }
}
neg_log_lik_weibull_corr_cause_specific_par_1 <- function(theta,D,M,cl)
{
  tau_vec_1 <- theta[1:M]
  rho_vec_1 <- theta[(M + 1):(2*M)]
  c <- 0
  for(i in 1:M)
  {
    if(theta[i] > 0)
    {
      c <- c + 1
    }
  }
  for(i in (M + 1):(2*M))
  {
    if((theta[i] < 1)&(theta[i] > 0))
    {
      c <- c + 1
    }
  }
  if(c == (2*M))
  {
    b <- -parApply(cl,D,1,weibull_corr_cause_indiv_log_lik,
                   mu_vec_1=final_mle_Weibull_par_2[1:3],mu_vec_2=final_mle_Weibull_par_2[4:6],
                   gamma_1=final_mle_Weibull_par_2[7],gamma_2=final_mle_Weibull_par_2[8],
                   tau_vec=tau_vec_1,rho_vec=rho_vec_1)
    return(sum(b))
  }else
  {
    return(10^19)
  }
}
objfun_Weibull_1 <- function(theta_1)
{
  return(neg_log_lik_weibull_corr_cause_specific_par_1(theta_1,hearing_data,3,cl))
}
objfun_Weibull_1(final_mle_Weibull_par_1)
gr_objfun_Weibull_1 <- function(theta_1)
{
  return(pracma::grad(objfun_Weibull_1,theta_1))
}
#gr_objfun_Weibull_1(data_soln_2_par_2)
hin_1 <- function(theta_1)
{
  h <- rep(NA, 1)
  h[1] <- theta_1[1]
  h[2] <- theta_1[2]
  h[3] <- theta_1[3]
  h[4] <- 1 - theta_1[4]
  h[5] <- 1 - theta_1[5]
  h[6] <- 1 - theta_1[6]
  #h[7] <- theta_1[7]
  #h[8] <- theta_1[8]
  return(h)
}
hin_1.jac <- function(theta_1)
{
  j <- matrix(NA, 6, length(theta_1))
  j[1,] <- c(1,rep(0,5))
  j[2,] <- c(0,1,rep(0,4))
  j[3,] <- c(0,0,1,rep(0,3))
  j[4,] <- c(rep(0,3),-1,rep(0,2))
  j[5,] <- c(rep(0,4),-1,0)
  j[6,] <- c(rep(0,5),-1)
  j
}
hin_1(final_mle_Weibull_par_1)
#final_mle_Weibull_1 <- NlcOptim::solnl(data_soln_2_par_2,objfun=objfun_Weibull_1,
#lb=replicate(6,0.0001),ub=c(replicate(3,Inf),replicate(3,0.999)))
final_mle_Weibull_1 <- alabama::auglag(par=final_mle_Weibull_par_1,fn=objfun_Weibull_1,
hin=hin_1,hin.jac=hin_1.jac,control.outer = list(trace=TRUE),control.optim = list(trace=6,ndeps=rep(6,1e-5)))
stopCluster(cl)
final_mle_Weibull_par_1 <- final_mle_Weibull_1$par
#final_mle_Weibull_par_1 <- as.vector(final_mle_Weibull_par_1)
sub_hess_1 <- final_mle_Weibull_1$hessian
#################################________________############################
### final m.l.e
final_mle_Weibull_corr_cause_specific_Gamma <- c(final_mle_Weibull_par_2,
                                                 final_mle_Weibull_par_1)
final_mle_Weibull_corr_cause_specific_Gamma
### neg log-likelihood value at m.l.e
neg_log_lik_Weibull_val <- neg_log_lik_weibull_corr_cause_specific(final_mle_Weibull_corr_cause_specific_Gamma,
                                                                   hearing_data,3)
neg_log_lik_Weibull_val
cl <- makeCluster(25)
setDefaultCluster(cl=cl)
clusterExport(cl,c("marg_subdis_weibull_corr_cause_specific","joint_surv_weibull_hazard_corr_cause_specific",
"joint_sub_dis_weibull_corr_cause_specific","weibull_corr_cause_indiv_log_lik","final_mle_Weibull_1","final_mle_Weibull_2"))
clusterEvalQ(cl,{
  final_mle_Weibull_par_1 <- final_mle_Weibull_1$par
  final_mle_Weibull_par_2 <- final_mle_Weibull_2$par
  final_mle_Weibull_corr_cause_specific_Gamma <- c(final_mle_Weibull_par_2,
                                                   final_mle_Weibull_par_1)
})
neg_log_lik_weibull_corr_cause_specific_par <- function(theta,D,M,cl)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  gamma_3 <- theta[(2*M) + 1]
  gamma_4 <- theta[(2*M) + 2]
  tau_vec_1 <- theta[((2*M) + 3):((3*M) + 2)]
  rho_vec_1 <- theta[((3*M) + 3):((4*M) + 2)]
  b <- -parApply(cl,D,1,weibull_corr_cause_indiv_log_lik,
                 mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,gamma_1=gamma_3,gamma_2=gamma_4,
                 tau_vec=tau_vec_1,rho_vec=rho_vec_1)
  return(sum(b))
}
objfun_par <- function(theta)
{
  return(neg_log_lik_weibull_corr_cause_specific_par(theta,hearing_data,3,cl))
}
final_gradient <- function(theta)
{
  return(pracma::grad(neg_log_lik_weibull_corr_cause_specific_par,
                      theta,D=hearing_data,M=3,cl=cl))
}
objfun_par(final_mle_Weibull_corr_cause_specific_Gamma)
final_hessian <- function(theta)
{
  return(pracma::jacobian(final_gradient,theta))
}
another_soln <- NlcOptim::solnl(X=data_soln_2_par,objfun=objfun_par,
lb=replicate(14,0.0001),ub=c(replicate(11,Inf),replicate(3,0.99)))
stopCluster(cl)
final_hess_1 <- another_soln$hessian
B <- final_hess_1[1:8,9:14]
A <- final_hess_1[1:8,1:8]
D <- final_hess_1[9:14,9:14]
A_inv <- as.matrix(Matrix::nearPD(armaInv(A))$mat)
D_inv <- as.matrix(Matrix::nearPD(armaInv(D))$mat)
F_1 <- as.matrix(D - t(B)%*%A_inv%*%B)  ## invertible
F_2 <- A - B%*%D_inv%*%t(B)  ## invertible
F_1 <- as.matrix(Matrix::nearPD(F_1)$mat)
F_2 <- as.matrix(Matrix::nearPD(F_2)$mat)
Z_1 <- matrix(0,nrow=6,ncol=8)
Z_2 <- matrix(0,nrow=8,ncol=6)
M_11 <- cbind(armaInv(F_1),Z_1)
M_12 <- cbind(Z_2,armaInv(F_2))
M_1 <- rbind(M_11,M_12)
M_21 <- cbind(diag(replicate(8,1),nrow=8,ncol=8),-B%*%D_inv)
M_22 <- cbind(-t(B)%*%A_inv,diag(replicate(6,1),nrow=6,ncol=6))
M_2 <- rbind(M_21,M_22)
inv_hess <- as.matrix(M_1%*%M_2)
inv_hess <- as.matrix(Matrix::nearPD(inv_hess)$mat)
### estimate of asymptotic sd of m.l.es'
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x);}",depends = "RcppArmadillo")
var_mle_Weibull <- diag(inv_hess)
asymp_sd_Weibull <- sqrt(var_mle_Weibull)
asymp_sd_Weibull
