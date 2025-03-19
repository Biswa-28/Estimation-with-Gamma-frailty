## call the data in r environment
data <- read.csv("~/SharedAndCorrCauseSpecific/R/data.csv")
## checking the structure and missingness of the data
str(data)
d <- attach(data)
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
## integrand of marginal sub-distribution function
integrand_2 <- function(t,j,mu_vec,tau_vec)
{
  a <- 1 + tau_vec*t*mu_vec
  b <- -(1/tau_vec)*log(a)
  marg_surv <- exp(sum(b))
  return((mu_vec[j])*marg_surv*((a[j])^(-1)))
}
## function to compute unconditional joint survival function
joint_surv_func <- function(t_1,t_2,mu_vec_1,mu_vec_2,tau_vec)
{
  b_1 <- log(1 + tau_vec*(t_1*mu_vec_1 + t_2*mu_vec_2))
  b_2 <- -(1/tau_vec)*b_1
  joint_surv <- exp(sum(b_2))
  return(joint_surv)
}
## integrand of joint sub-distribution function(number of causes is 3)
integrand_3 <- function(t_1,t_2,j_1,j_2,mu_vec_1,
                        mu_vec_2,tau_vec)
{
  joint_surv <- joint_surv_func(t_1,t_2,mu_vec_1,mu_vec_2,tau_vec)
  if(j_1 == j_2)
  {
    j = j_1
    a_1 <- (mu_vec_1[j])
    a_2 <- (mu_vec_2[j])
    d <- (1 + tau_vec[j]*(t_1*a_1 + t_2*a_2))^(-2)
    return((1 + tau_vec[j])*a_1*a_2*joint_surv*d)
  }else
  {
    a_1 <- (mu_vec_1[j_1])
    a_2 <- (mu_vec_1[j_2])
    a_3 <- (mu_vec_2[j_1])
    a_4 <- (mu_vec_2[j_2])
    d_1 <- (1 + tau_vec[j_1]*(t_1*a_1 + t_2*a_3))^(-1)
    d_2 <- (1 + tau_vec[j_2]*(t_1*a_2 + t_2*a_4))^(-1)
    return(a_1*a_4*joint_surv*d_1*d_2)
  }
}
## marginal distribution function calculation
#####---------------------------------------------------------- we have created a R package called SharedAndCorrCauseSpecific using C++ functions then calling it for fast computation----------------------------------------------------------######
###### its srcs' has been also uploaded here.
marg_subdis_exp_shared_cause_specific <- function(t,j,mu_vec,tau_vec)
{
  return(SharedAndCorrCauseSpecific::marg_subdis_exp_shared_cause_specific_cpp(t,j,
                                                                    mu_vec,tau_vec))
}
## joint sub-distribution function calculation
joint_sub_dis_exp_shared_cause_specific <- function(t_1,t_2,j_1,j_2,mu_vec_1,
                                   mu_vec_2,tau_vec)
{
  return(SharedAndCorrCauseSpecific::joint_sub_dis_exp_shared_cause_specific_cpp(t_1,t_2,
                                      j_1,j_2,mu_vec_1,mu_vec_2,tau_vec,3)$approximate)
}

joint_sub_dis_shared_cause_1 <- function(t_3,t_4,j_3,j_4,mu_vec_3,
                                         mu_vec_4,tau_vec_1)
{
  z <- function(v)
  {
    I_1 <- sapply(v, g <- function(y)
    {
      g_1 <- function(u)
      {
        return(integrand_3(u,t_1=y,j_1=j_3,j_2=j_4,mu_vec_1=mu_vec_3,
                           mu_vec_2=mu_vec_4,tau_vec=tau_vec_1))
      }
      return(cubature::cubintegrate(g_1,lower=0,upper=t_3,relTol=1e-9)$integral)
    })
    return(I_1)
  }
  val <- cubature::cubintegrate(z,lower=0,upper=t_4,relTol=1e-9)$integral
  return(val)
}
## cell probability for exp hazard shared cause specific gamma frailty
prob_cell_shared_cause <- function(x,j_1,j_2,mu_vec_1,mu_vec_2,tau_vec)
{
  log_joint_surv <- -(1/tau_vec)*log(1 + tau_vec*x*(mu_vec_1 + mu_vec_2))
  joint_surv <- exp(sum(log_joint_surv))
  if((j_1 == 1)&&(j_2 == 1))
  {
    return(joint_surv)
  }else if((j_1 > 1)&&(j_2 == 1))
  {
    return(marg_subdis_exp_shared_cause_specific(x,j_1 - 1,mu_vec_1,tau_vec) -
             sum(sapply(1:3,joint_sub_dis_exp_shared_cause_specific,t_1=x,t_2=x,j_1=j_1 - 1,
                        mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,tau_vec=tau_vec)))
  }else if((j_1 == 1)&&(j_2 > 1))
  {
    return(marg_subdis_exp_shared_cause_specific(x,j_2 - 1,mu_vec_2,tau_vec) -
             sum(sapply(1:3,joint_sub_dis_exp_shared_cause_specific,t_1=x,t_1=x,j_2=j_2 - 1,
                        mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,tau_vec=tau_vec)))
  }else
  {
    return(joint_sub_dis_exp_shared_cause_specific(x,x,j_1 - 1,j_2 - 1,mu_vec_1,mu_vec_2,tau_vec))
  }
}
## function to compute the initial value
in_est_shared_cause <- function(monitoring_data,Freq_table,M,theta)
{
  mu_vec_1 <- theta[1:M]
  mu_vec_2 <- theta[(M + 1):(2*M)]
  tau_vec <- theta[((2*M) + 1):(3*M)]
  sum(sapply(1:(M + 1),function(j_1)
  {
    return(sum(sapply(1:(M + 1),function(j_2)
    {
      prob_cell_new <- function(x)
      {
        return(prob_cell_shared_cause(x,j_1,j_2,mu_vec_1,
                                      mu_vec_2,tau_vec))
      }
      d <- monitoring_data[monitoring_data[,2] == (j_1 - 1) &
                             monitoring_data[,3] == (j_2 - 1),]
      if(is.matrix(d) == TRUE)
      {
        if(length(d[,1]) != 0)
        {
          c <- c + (Freq_table[j_1,j_2] - sum(sapply(d[,1],
                                                     prob_cell_new)))^2
        }
      }
      if(is.vector(d) == TRUE)
      {
        c <- c + (Freq_table[j_1,j_2] - prob_cell_new(d[1]))^2
      }
    })))
  }))
}
## initial estimate
in_est_shared_cause_1 <- function(theta,monitoring_data,Freq_table,M)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  tau_vec_1 <- theta[((2*M) + 1):(3*M)]
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
                                                     prob_cell_shared_cause,j_1=j_1,j_2=j_2,mu_vec_1=mu_vec_3,
                                                     mu_vec_2=mu_vec_4,tau_vec=tau_vec_1)))^2
        }
      }
      if(is.vector(d) == TRUE)
      {
        c <- c + (Freq_table[j_1,j_2] - prob_cell_shared_cause(d[1],
                                                               j_1,j_2,mu_vec_3,mu_vec_4,tau_vec_1))^2
      }
    }
  }
  return(c)
}
in_est_shared_cause_2 <- function(theta,monitoring_data,Freq_table,M)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  tau_vec_1 <- theta[((2*M) + 1):(3*M)]
  a <- sapply(1:(M + 1),f_1 <- function(i)
  {
    b <- sapply(1:(M + 1),f_2 <- function(j)
    {
      d <- monitoring_data[monitoring_data[,2] == (i - 1) &
                             monitoring_data[,3] == (j - 1),]
      d <- matrix(d,ncol=3)
      if(length(d[,1]) > 0)
      {
        e <- sapply(d[,1],
                    prob_cell_shared_cause,j_1=i,j_2=j,mu_vec_1=mu_vec_3,
                    mu_vec_2=mu_vec_4,tau_vec=tau_vec_1)
        return((Freq_table[i,j] - sum(e))^2)
      }else
      {
        return(c(0))
      }
    })
    return(sum(b))
  })
  return(sum(a))
}
#### parallel version
in_est_shared_cause_2_par <- function(theta,monitoring_data,Freq_table,M,cl)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  tau_vec_1 <- theta[((2*M) + 1):(3*M)]
  a <- parSapply(cl,1:(M + 1),f_1 <- function(i)
  {
    b <- parSapply(cl,1:(M + 1),f_2 <- function(j)
    {
      d <- monitoring_data[monitoring_data[,2] == (i - 1) &
                             monitoring_data[,3] == (j - 1),]
      d <- matrix(d,ncol=3)
      if(length(d[,1]) > 0)
      {
        e <- parSapply(cl,d[,1],
                       prob_cell_shared_cause,j_1=i,j_2=j,mu_vec_1=mu_vec_3,
                       mu_vec_2=mu_vec_4,tau_vec=tau_vec_1)
        return((Freq_table[i,j] - sum(e))^2)
      }else
      {
        return(c(0))
      }
    })
    return(sum(unlist(b)))
  })
  return(sum(unlist(a)))
}
## numerical deribative of above function
num_deriv_in_est_1 <- function(monitoring_data,Freq_table,M,theta)
{
  deriv_func <- pracma::grad(in_est_shared_cause_2,theta,
                             monitoring_data=monitoring_data,Freq_table=Freq_table,M=M)
  return(deriv_func)
}
## log-likelihood function for an individual
shared_cause_indiv_log_lik <- function(data,mu_vec_1,mu_vec_2,tau_vec)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  a <- tau_vec*x
  c_1 <- a*mu_vec_1
  c_2 <- a*mu_vec_2
  d <- (1/tau_vec)
  v <- -d*log(1 + c_1 + c_2)
  log_joint_surv <- sum(v)
  eps <- 10^(-25)
  if((j_1 != 0) && (j_2 != 0))
  {
    h <- joint_sub_dis_exp_shared_cause_specific(x,x,j_1,j_2,
                                      mu_vec_1,mu_vec_2,tau_vec)
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }else if((j_1 != 0) && (j_2 == 0)){
    h <- marg_subdis_exp_shared_cause_specific(x,j_1,mu_vec_1,tau_vec) - sum(sapply(1:3,joint_sub_dis_exp_shared_cause_specific,t_1=x,t_2=x,j_1=j_1,
                                                          mu_vec_1 = mu_vec_1,mu_vec_2 = mu_vec_2,tau_vec=tau_vec))
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }else if((j_1 == 0) && (j_2 != 0)){
    h <- marg_subdis_exp_shared_cause_specific(x,j_2,mu_vec_2,tau_vec) -
      sum(sapply(1:3,joint_sub_dis_exp_shared_cause_specific,t_1=x,t_2=x,j_2=j_2,
                 mu_vec_1 = mu_vec_1,mu_vec_2 = mu_vec_2,tau_vec=tau_vec))
    if(h > eps)
    {
      return(log(h))
    }else
    {
      return(log(eps))
    }
  }else{
    return(log_joint_surv)
  }
}
## negative log-likelihood function which we want to minimize
## the argument about which we want to minimize shoud be first
neg_log_lik_shared_cause <- function(theta,D,M)
{
  lambda_vec_3 <- theta[1:M]
  lambda_vec_4 <- theta[(M + 1):(2*M)]
  tau_vec_1 <- theta[((2*M) + 1):(3*M)]
  return(-sum(apply(D,1,shared_cause_indiv_log_lik,
                    mu_vec_1=lambda_vec_3,mu_vec_2=lambda_vec_4,
                    tau_vec=tau_vec_1)))
}
### parallel version
neg_log_lik_shared_cause_par <- function(theta,D,M,cl)
{
  lambda_vec_3 <- theta[1:M]
  lambda_vec_4 <- theta[(M + 1):(2*M)]
  tau_vec_1 <- theta[((2*M) + 1):(3*M)]
  return(-sum(parApply(cl,D,1,shared_cause_indiv_log_lik,
                       mu_vec_1=lambda_vec_3,mu_vec_2=lambda_vec_4,
                       tau_vec=tau_vec_1)))
}
cl <- makeCluster(28)
setDefaultCluster(cl=cl)
clusterExport(cl,c("prob_cell_shared_cause","marg_subdis_shared_cause","parSapply",
                   "integrand_2","integrand_3","joint_sub_dis_shared_cause_1","joint_surv_func"))
#neg_log_lik_shared_cause_par(c(43.6802,0.015,0.019,13.2426,0.032,0.024,12.7383,15.96,13.4317),hearing_data,3,cl)
neg_log_lik_shared_cause_par(soln_2,hearing_data,3,cl)
neg_log_lik_shared_cause_par(c(35.9582,0.0619684,0.0785706,39.2107,0.0568886,
0.0840602,12.9245,45.5892,18.2817),hearing_data,3,cl)
#hess <- pracma::hessian(neg_log_lik_shared_cause_par,c(35.9582,0.0619684,0.0785706,39.2107,0.0568886,0.0840602,12.9245,45.5892,18.2817),
                        #D=hearing_data,M=3,cl=cl)
#in_est_shared_cause_2_par(c(lambda_vec_1_initial,lambda_vec_2_initial,c(0.6,4.5,2.3)),hearing_data,F,3,cl)
stopCluster(cl)
## numerical gradient function of log-likelihood
num_grad_neg_log_lik <- function(D,theta,M)
{
  return(pracma::grad(neg_log_lik_shared_cause,theta,D=D,M=M))
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
ui <- diag(replicate(9,1),nrow=9,ncol=9)
ci <- replicate(9,0)
# finding a better initial point for minimization of original
## neg-loglikelihood function
l <- seq(0.1,2,length=5)
in_val_3d <- array(0,dim=c(5,5,5))
c <- 0
for(i in 1:length(l))
{
  for(j in 1:length(l))
  {
    for(k in 1:length(l))
    {
      in_val_3d[i,j,k] <- in_est_shared_cause_1(c(lambda_vec_1_initial,lambda_vec_2_initial,
                                                  c(l[i],l[j],l[k])),hearing_data,F,3)
      c <- c + 1
      print(c)
    }
  }
}
which.min(in_val_3d)
library("doParallel")
cl <- makeCluster(19)
setDefaultCluster(cl=cl)
clusterExport(cl,c("prob_cell_shared_cause","marg_subdis_shared_cause",
                   "integrand_2","integrand_3","joint_sub_dis_shared_cause_1","joint_surv_func"))
in_est <- optimParallel::optimParallel(par=c(lambda_vec_1_initial,lambda_vec_2_initial,c(l[1],l[5],l[5])),
                                       fn=in_est_shared_cause_1,monitoring_data = hearing_data,Freq_table=F,
                                       M = 3,lower=replicate(9,0.01),method="L-BFGS-B",control = list(trace=6,ndeps=replicate(9,1e-5)),
                                       parallel = list(cl=NULL,forward=FALSE,loginfo=TRUE))
stopCluster(cl)
#cl <- makeCluster(19)
#setDefaultCluster(cl=cl)
#clusterExport(cl,c("prob_cell_shared_cause","marg_subdis_shared_cause",
#"integrand_2","integrand_3","joint_sub_shared_cause","joint_surv_func"))
#in_est <- optimParallel::optimParallel(par=c(lambda_vec_1_initial,lambda_vec_2_initial,c(1,1,1)),
#fn=in_est_shared_cause_2,monitoring_data = hearing_data,Freq_table=F,
#M = 3,lower=replicate(9,0.01),method="L-BFGS-B",
#parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE))
#stopCluster(cl)
#setDefaultCluster(cl=NULL)
in_est$par
in_est$convergence
in_est_shared_cause_1(hearing_data,F,3,in_est$par)
in_est_shared_cause_1(hearing_data,F,3,c(lambda_vec_1_initial,lambda_vec_2_initial,c(l[1],l[5],l[5])))
in_est_data <- in_est$par
## min in_est_1 at sigma = 3.9,value = 4687.594,taking seq from 0.1 to 4
## Minimization of original neg-loglikelihood function for shared frailty
library("doParallel")
#cl <- makePSOCKcluster(7)
#registerDoParallel(cl)
#data_soln <- constrOptim(in_est_data,neg_log_lik_shared_cause,num_grad_neg_log_lik,
#M=3,D=hearing_data,ui=ui,ci=ci,method = "BFGS")$par
#hess <- numDeriv::hessian(neg_log_lik_shared_cause,data_soln,D=hearing_data,M=3)
#stopCluster(cl)
#data_soln
cl <- makeCluster(19)
setDefaultCluster(cl=cl)
clusterExport(cl,c("prob_cell_shared_cause","marg_subdis_shared_cause",
                   "integrand_2","integrand_3","joint_sub_dis_shared_cause_1","joint_surv_func",
                   "shared_cause_indiv_log_lik"))
data_soln <- optimParallel::optimParallel(par=c(lambda_vec_1_initial,lambda_vec_2_initial,c(l[1],l[5],l[5])),
                                          fn=neg_log_lik_shared_cause,M=3,D=hearing_data,lower=replicate(9,0.01),control=list(trace=6),
                                          parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE),method = "L-BFGS-B")$par
stopCluster(cl)
setDefaultCluster(cl=NULL)
library(doParallel)
cl <- makeCluster(19)
setDefaultCluster(cl=cl)
clusterExport(cl,c("prob_cell_shared_cause","marg_subdis_exp_shared_cause_specific",
                   "joint_sub_dis_exp_shared_cause_specific","joint_surv_func",
                   "shared_cause_indiv_log_lik"))
data_soln <- optimParallel::optimParallel(par=c(35.9582,0.0051,0.0067,39.2107,0.0053,0.0061,12.9245,45.5892,18.2817),hessian = TRUE,
                                          fn=neg_log_lik_shared_cause,M=3,D=hearing_data,lower=replicate(9,0.001),control=list(trace=6),
                                          parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE),method = "L-BFGS-B")
stopCluster(cl)
setDefaultCluster(cl=NULL)
library(doParallel)
cl <- makeCluster(19)
setDefaultCluster(cl=cl)
clusterExport(cl,c("prob_cell_shared_cause","marg_subdis_exp_shared_cause_specific",
                   "joint_sub_dis_exp_shared_cause_specific","joint_surv_func",
                   "shared_cause_indiv_log_lik"))
data_soln <- optimParallel::optimParallel(par=soln_1,hessian = TRUE,
                                          fn=neg_log_lik_shared_cause,M=3,D=hearing_data,lower=replicate(9,0.001),control=list(trace=6),
                                          parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE),method = "L-BFGS-B")
stopCluster(cl)
setDefaultCluster(cl=NULL)
library(doParallel)
cl <- makeCluster(19)
setDefaultCluster(cl=cl)
clusterExport(cl,c("prob_cell_shared_cause","marg_subdis_exp_shared_cause_specific",
                   "joint_sub_dis_exp_shared_cause_specific","joint_surv_func",
                   "shared_cause_indiv_log_lik"))
data_soln_2_par <- c(35.958200536,0.006048146,0.003377293,39.210699763,0.005502268,
                     0.003451682,12.924503787,45.589199923,9.474994132)
data_soln_3 <- optimParallel::optimParallel(par=c(data_soln_2_par[1:8],7.657),hessian = TRUE,
                                            fn=neg_log_lik_shared_cause,M=3,D=hearing_data,lower=replicate(9,0.001),control=list(trace=6),
                                            parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE),method = "L-BFGS-B")
stopCluster(cl)
setDefaultCluster(cl=NULL)
data_soln_3_par <- data_soln_3$par
hess_3 <- data_soln_3$hessian
eigen(hess_3)$values
soln_1 <- c(35.9582,0.00619684,0.00785706,39.2107,0.00568886,0.00840602,12.9245,45.5892,18.2817)
soln_2 <- c(35.9582,0.19684,0.285706,39.2107,0.468886,0.840602,12.9245,45.5892,18.2817)
hess <- data_soln$hessian
eigen(hess)$values ## one eigen value negative and very close to zero
pd_hess <- Matrix::nearPD(hess)$mat
### inverse of nearest pd matrix
inv_pd_hess <- solve(pd_hess)
### estimate of asymptotic s.d of m.l.e
sd_mle_1 <- sqrt(diag(inv_pd_hess))
cl <- makeCluster(28)
setDefaultCluster(cl=cl)
clusterExport(cl,c("prob_cell_shared_cause","marg_subdis_shared_cause",
                   "integrand_2","integrand_3","joint_sub_dis_shared_cause_1","joint_surv_func",
                   "shared_cause_indiv_log_lik","parApply"))
hess_1 <- pracma::hessian(neg_log_lik_shared_cause_par,soln_1,D=hearing_data,M=3,cl=cl)
eigen(hess_1)$values
sd_mle_1 <- sqrt(diag(solve(hess_1)))
sd_mle_1
stopCluster(cl)
## checking hessian at data_soln_3(i.e at m.l.e) is p.d or not
eigen(hess_3)$values
## estimating asymptotic standard error of the m.l.e
sd_mle <- sqrt(diag(solve(hess_3)))
sd_mle
## Computation of AIC for shared frailty model
AIC_shared_cause_specific <- (2*neg_log_lik_shared_cause(data_soln_3_par,hearing_data,3)) + (2*9)
AIC_shared_cause_specific        #1838.964
