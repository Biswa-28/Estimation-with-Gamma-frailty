## checking the structure and missingness of the data
str(data)
d <- attach(data)
Age <- Age[-c(380)]
## transforming the monitoring time in months
chng_year_to_mths <- function(x)## checking the structure and missingness of the data
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
delete_rows <- as.numeric(rownames(subset(data_1,Loss == "Y" & is.na(Type) == TRUE)))
data_1 <- data_1[-c(delete_rows),]
failure_cause_gen <- function(x)
{
  loss <- x[1]
  type <- x[2]
  if((loss == "N") & (is.na(type) == TRUE))     ## 1
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
  if((loss == "Y") & (type == "R-Conductive"))  ## 11
  {
    return(c(0,2))
  }
  if((loss == "Y") & (type == "L-Conductive, R-Mixed"))  ## 12
  {
    return(c(2,3))
  }
  if((loss == "Y") & (type == "L-Conductive"))  ## 13
  {
    return(c(2,0))
  }
}
data_failure_cause <- t(apply(data_1[,c(2,3)],1,failure_cause_gen))
data_failure_cause <- array(0,dim = c(796,2))
for(i in 1:796)
{
  data_failure_cause[i,] <- t(failure_cause_gen(data_1[i,c(2,3)]))
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
library("devtools")
devtools::load_all()
Rcpp::sourceCpp("corr_cause_specific.cpp")
corrspecific::marg_subdis_corr_cause_cpp(13.4,lambda_vec_1_initial,1,c(0.8,0.4,0.3))
marg_subdis_corr_cause_cpp(19.40,lambda_vec_2_initial,2,c(0.8,0.4,0.3))
## integrand of marginal sub-distribution function
integrand_2 <- function(t,j,mu_vec,tau_vec)
{
  a <- 1 + (tau_vec^2)*t*(1/mu_vec)
  b <- -((1/tau_vec)^2)*log(a)
  marg_surv <- exp(sum(b))
  return((1/mu_vec[j])*marg_surv*((a[j])^(-1)))
}
## integrand of joint sub-distribution function
integrand_3 <- function(t_1,t_2,j_1,j_2,mu_vec_1,
                        mu_vec_2,tau_vec_1,tau_vec_2,phi_vec)
{
  a_1 <- 1/mu_vec_1
  a_2 <- 1/mu_vec_2
  c_1 <- 1/tau_vec_1
  c_2 <- 1/tau_vec_2
  b <- (phi_vec*c_1*c_2)
  b_1 <- 1 + ((tau_vec_1^2)*t_1*(a_1))
  b_2 <- 1 + ((tau_vec_2^2)*t_2*(a_2))
  b_3 <- b_1 + b_2 - 1
  v_1 <- -(b*log(b_3))
  v_2 <- ((b - (c_1^2))*log(b_1))
  v_3 <- ((b - (c_2^2))*log(b_2))
  joint_surv <- exp(sum(v_1 + v_2 + v_3))
  if(j_1 == j_2)
  {
    j = j_1
    e_1 <- b[j]
    d_1 <- e_1*(1 + e_1)*((1/b_3[j])^2)
    d_2 <- (e_1*((c_1[j]^2) - e_1))*(1/(b_1[j]*b_3[j]))
    d_3 <- (e_1*((c_2[j]^2) - e_1))*(1/(b_2[j]*b_3[j]))
    d_4 <- ((c_1[j]^2) - e_1)*((c_2[j]^2) - e_1)*(1/(b_1[j]*b_2[j]))
    return((a_1[j])*(a_2[j])*((tau_vec_1[j])^2)*
             ((tau_vec_2[j])^2)*joint_surv*(d_1 + d_2 + d_3 + d_4))
  }else
  {
    e_1 <- b[j_1]
    e_2 <- b[j_2]
    d_1 <- (e_1*e_2)/(b_3[j_1]*b_3[j_2])
    d_2 <- (e_1*((c_2[j_2]^2) - e_2))/(b_3[j_1]*b_2[j_2])
    d_3 <- (e_2*((c_1[j_1]^2) - e_1))/(b_3[j_2]*b_1[j_1])
    d_4 <- (((c_1[j_1]^2) - e_1)*((c_2[j_2]^2) - e_2))/(b_1[j_1]*b_2[j_2])
    return(a_1[j_1]*a_2[j_2]*((tau_vec_1[j_1])^2)*((tau_vec_2[j_2])^2)*
             joint_surv*(d_1 + d_2 + d_3 + d_4))
  }
}
marg_subdis_corr_cause_cpp <- function(t,mu_vec,j,tau_vec)
{
  return(corrspecific::marg_subdis_corr_cause_cpp_1(t,mu_vec,j,tau_vec))
}
## consider common sigma vector for the cause sets for simplicity
joint_sub_dis_corr_cause_cpp <- function(t_1,t_2,j_3,j_4,mu_vec_3,
                                           mu_vec_4,tau_vec,phi_vec_1)
{
  return(corrspecific::joint_sub_dis_corr_cause_cpp_1(t_1,t_2,j_3,j_4,
                                      mu_vec_3,mu_vec_4,tau_vec,phi_vec_1,3)$approximate)
}
joint_sub_dis_corr_cause_cpp(150,150,1,3,lambda_vec_1_initial,lambda_vec_2_initial,
                               c(0.7,0.4,0.1),c(0.78,0.45,0.88))
## consider common sigma vector for the cause sets for simplicity
## log-likelihood function for an individual
corr_cause_indiv_log_lik <- function(data,mu_vec_1,mu_vec_2,
                                     tau_vec,phi_vec)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  a_1 <- (1/mu_vec_1)
  a_2 <- (1/mu_vec_2)
  t <- (tau_vec^2)
  b_1 <- 1 + (t*x*a_1)
  b_2 <- 1 + (t*x*a_2)
  b_3 <- b_1 + b_2 - 1
  c <- (1/tau_vec)
  d <- phi_vec*(c^2)
  v_1 <- -d*log(b_3)
  v_2 <- ((d - (c^2))*log(b_1))
  v_3 <- ((d - (c^2))*log(b_2))
  log_joint_surv_corr <- sum(v_1 + v_2 + v_3)
  eps <- 10^(-50)
  if((j_1 != 0) && (j_2 != 0))
  {
    h <- joint_sub_dis_corr_cause_cpp(x,x,j_1,j_2,
                                        mu_vec_1,mu_vec_2,tau_vec,phi_vec)
    if(h > eps)
    {
      return(log(h))
    }else{
      return(log(eps))
    }
  }else if((j_1 != 0) && (j_2 == 0)){
    h <- marg_subdis_corr_cause_cpp(x,mu_vec_1,j_1,tau_vec) -
      sum(sapply(1:3,joint_sub_dis_corr_cause_cpp,t_1=x,t_2=x,j_3=j_1,
                 mu_vec_3 = mu_vec_1,mu_vec_4 = mu_vec_2,tau_vec=tau_vec,phi_vec_1=phi_vec))
    if(h > eps)
    {
      return(log(h))
    }else{
      return(log(eps))
    }
  }else if((j_1 == 0) && (j_2 != 0)){
    h <- marg_subdis_corr_cause_cpp(x,mu_vec_2,j_2,tau_vec) -
      sum(sapply(1:3,joint_sub_dis_corr_cause_cpp,t_1=x,t_2=x,j_4=j_2,
                 mu_vec_3 = mu_vec_1,mu_vec_4 = mu_vec_2,tau_vec=tau_vec,phi_vec_1=phi_vec))
    if(h > eps)
    {
      return(log(h))
    }else{
      return(log(eps))
    }
  }else{
    return(log_joint_surv_corr)
  }
}
## negative log-likelihood function which we want to minimize
neg_log_lik_corr_cause <- function(theta,D,M)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  tau_vec <- theta[((2*M) + 1):(3*M)]
  phi_vec_1 <- theta[((3*M) + 1):(4*M)]
  return(-sum(apply(D,1,corr_cause_indiv_log_lik,
                    mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,
                    tau_vec=tau_vec,phi_vec=phi_vec_1)))
}
## parallel version of the above
neg_log_lik_corr_cause_par <- function(theta,D,M,cl)
{
  mu_vec_3 <- theta[1:M]
  mu_vec_4 <- theta[(M + 1):(2*M)]
  tau_vec <- theta[((2*M) + 1):(3*M)]
  phi_vec_1 <- theta[((3*M) + 1):(4*M)]
  return(-sum(parApply(cl,D,1,corr_cause_indiv_log_lik,
                    mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,
                    tau_vec=tau_vec,phi_vec=phi_vec_1)))
}
## numerical gradient function of log-likelihood
num_grad_neg_log_lik <- function(theta,D,M)
{
  return(pracma::grad(neg_log_lik_corr_cause,theta,D=D,M=M))
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
cl <- makeCluster(25)
setDefaultCluster(cl=cl)
clusterExport(cl,c("corr_cause_indiv_log_lik","joint_sub_dis_corr_cause_cpp_1",
                   "joint_sub_dis_corr_cause_cpp","marg_subdis_corr_cause_cpp"))
soln <- optimParallel::optimParallel(par=c(lambda_vec_1_initial,lambda_vec_2_initial,
        c(10.1,7.1,6.2),c(0.9,0.9,0.9)),fn=neg_log_lik_corr_cause,
D=hearing_data,M=3,method = "L-BFGS-B",lower=replicate(4*3,0.01),control=list(trace=6),
upper=c(replicate(3*3,Inf),0.99,0.99,0.99),parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE))$par
stopCluster(cl)
cl <- makeCluster(25)
setDefaultCluster(cl=cl)
clusterExport(cl,c("corr_cause_indiv_log_lik","joint_sub_dis_corr_cause_cpp_1",
                   "joint_sub_dis_corr_cause_cpp","marg_subdis_corr_cause_cpp"))
system.time(neg_log_lik_corr_cause(soln,hearing_data,3))
system.time(neg_log_lik_corr_cause_par(soln,hearing_data,3,cl))
stopCluster(cl)
cl <- makeCluster(25)
setDefaultCluster(cl=cl)
clusterExport(cl,c("corr_cause_indiv_log_lik","joint_sub_dis_corr_cause_cpp_1",
                   "joint_sub_dis_corr_cause_cpp","marg_subdis_corr_cause_cpp"))
hess <- pracma::hessian(neg_log_lik_corr_cause_par,soln,D=hearing_data,M=3,cl=cl)
stopCluster(cl)
