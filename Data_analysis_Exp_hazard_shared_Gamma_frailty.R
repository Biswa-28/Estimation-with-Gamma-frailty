## call the data in r environment
data <- read.csv("~/data.csv")
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
marg_surv <- function(x,tau,mu)
{
  return((1 + tau*mu*x)^(-(1/tau)))
}
joint_surv <- function(x,tau,mu_1,mu_2)
{
  return((1 + tau*x*(mu_1+mu_2))^(-(1/tau)))
}
## to find estimate of sigma
s_initial <- function(sigma)
{
  return(mean(sapply(Age,joint_surv,tau=sigma,mu_1=lambda_1_initial,
                     mu_2=lambda_2_initial)) - (F[1,1]/sum(F)))
}
sigma <- seq(0.01,0.1,length=100)
p <- sapply(sigma,s_initial)
plot(sigma,p,type="l",ylim=c(0,0.3))
## function to generate failure causes for two individuals
failure_cause_gen <- function(x,tau,mu_1,mu_2,
                              M_1,M_2,q_1,q_2,q_3)
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
data_gen <- function(n,tau,mu_1,mu_2,
                     M_1,M_2,q_4,q_5,q_6)
{
  mon_time <- rexp(n,rate=mu_hat)
  F <- sapply(mon_time,failure_cause_gen,tau=sigma,
              mu_1=lambda_1,mu_2=lambda_2,M_1=M_1,M_2=M_2,
              q_1=q_4,q_2=q_5,q_3=q_6)
  return(cbind(mon_time,t(F)))
}
## calculation of cell probabilities for exponential hazard shared frailty model
prob_cell <- function(x,mu_1,mu_2,tau,
                      q_1,q_2,q_3,j_1,j_2)
{
  if((j_1 == 1)&&(j_2 == 1))
  {
    return(joint_surv(x,tau,mu_1,mu_2))
  }else if((j_1 > 1)&&(j_2 == 1))
  {
    return(q_1[j_1 - 1]*(marg_surv(x,tau,mu_2) - joint_surv(x,tau,mu_1,mu_2)))
  }else if((j_1 == 1)&&(j_2 > 1))
  {
    return(q_2[j_2 - 1]*(marg_surv(x,tau,mu_1) - joint_surv(x,tau,mu_1,mu_2)))
  }else
  {
    return(q_3[j_1 - 1,j_2 - 1]*(1 - marg_surv(x,tau,mu_1) - marg_surv(x,tau,mu_2) + joint_surv(x,tau,mu_1,mu_2)))
  }
}
## function to minimize for initial value for shared frailty model
in_est_1 <- function(monitoring_data,Freq_table,M_1,M_2,theta)
{
  mu_vec_1 <- theta[1:M_1]
  mu_vec_2 <- theta[(M_1 + 1):(M_1 + M_2)]
  tau <- theta[M_1+M_2+1]
  mu_1 <- sum(mu_vec_1)
  mu_2 <- sum(mu_vec_2)
  q_4 <- mu_vec_1/mu_1
  q_5 <- mu_vec_2/mu_2
  q_6 <- outer(q_4,q_5)
  c <- 0
  for(j_1 in 1:(M_1 + 1))
  {
    for(j_2 in 1:(M_2 + 1))
    {
      c <- c + (Freq_table[j_1,j_2] -  sum(sapply(monitoring_data[,1],prob_cell,
           mu_1=mu_1,mu_2=mu_2,tau=tau,q_1=q_4,q_2=q_5,q_3=q_6,j_1=j_1,j_2=j_2)))^2
    }
  }
  return(c)
}
## log-likelihood for a particular sample point
## for shared frailty model
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
## the overall negative log-likelihood function for shared frailty model which is
## to minimize
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
