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
data_failure_cause <- t(apply(data_1[,c(2,3)],1,failure_cause_gen))
hearing_data <- cbind(data_1[,1],data_failure_cause)
## getting the frequency table
F <- freq_cause(hearing_data,3,3)
## getting the row and column marginals
marginal_row <- apply(F,1,sum)
marginal_col <- apply(F,2,sum)
lambda_vec_1_initial <- (marginal_row/sum(Age))[2:4] 
lambda_vec_2_initial <- (marginal_col/sum(Age))[2:4]
lambda_1_initial <- sum(lambda_vec_1_initial)
lambda_2_initial <- sum(lambda_vec_2_initial)
marg_surv_wei <- function(x,beta,mu_vec,tau)
{
  a <- sum(mu_vec^beta)
  return((1 + (tau*a*(x^beta)))^(-(1/tau)))
}
joint_surv_wei <- function(x,beta_1,beta_2,mu_vec_1,mu_vec_2,tau)
{
  a_1 <- sum(mu_vec_1^beta_1)
  a_2 <- sum(mu_vec_2^beta_2)
  return((1 + tau*(a_1*(x^beta_1) + a_2*(x^beta_2)))^(-(1/tau)))
}
## function to generate failure causes for two individuals
weibull_failure_cause_gen <- function(x,tau,beta_1,beta_2,
                                      a_1,a_2,M_1,M_2,q_1,q_2,q_3)
{
  d_1 <- (tau*a_1*(x^beta_1))
  d_2 <- (tau*a_2*(x^beta_2))
  marg_surv_1 <- (1 + d_1)^(-(1/tau))
  marg_surv_2 <- (1 + d_2)^(-(1/tau))
  joint_surv <- (1 + d_1 + d_2)^(-(1/tau))
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
weibull_data_gen <- function(n,tau,a_1,a_2,beta_1,beta_2,
                             M_1,M_2,q_1,q_2,q_3)
{
  mon_time <- rexp(n,rate=mu_hat)
  F <- sapply(mon_time,weibull_failure_cause_gen,tau=tau,
              a_1=a_1,a_2=a_2,M_1=M_1,M_2=M_2,beta_1=beta_1,
              beta_2=beta_2,q_1=q_1,q_2=q_2,q_3=q_3)
  return(cbind(mon_time,t(F)))
}

## calculation of cell probabilities for exponential hazard shared frailty model
prob_cell_wei_sh <- function(x,mu_vec_1,mu_vec_2,beta_1,beta_2,
                      tau,q_1,q_2,q_3,j_1,j_2)
{
  if((j_1 == 1)&&(j_2 == 1))
  {
    return(joint_surv_wei(x,beta_1,beta_2,mu_vec_1,mu_vec_2,tau))
  }else if((j_1 > 1)&&(j_2 == 1))
  {
    return(q_1[j_1 - 1]*(marg_surv_wei(x,beta_2,mu_vec_2,tau) - 
           joint_surv_wei(x,beta_1,beta_2,mu_vec_1,mu_vec_2,tau)))
  }else if((j_1 == 1)&&(j_2 > 1))
  {
    return(q_2[j_2 - 1]*(marg_surv_wei(x,beta_1,mu_vec_1,tau) - 
          joint_surv_wei(x,beta_1,beta_2,mu_vec_1,mu_vec_2,tau)))
  }else
  {
    return(q_3[j_1 - 1,j_2 - 1]*(1 - marg_surv_wei(x,beta_1,mu_vec_1,tau) - 
                  marg_surv_wei(x,beta_2,mu_vec_2,tau) + 
        joint_surv_wei(x,beta_1,beta_2,mu_vec_1,mu_vec_2,tau)))
  }
}
## function to minimize for initial value for 
## weibull hazard shared frailty model
in_est_wei_sh <- function(theta,monitoring_data,Freq_table,M_1,M_2)
{
  mu_vec_1 <- theta[1:M_1]
  mu_vec_2 <- theta[(M_1 + 1):(M_1 + M_2)]
  beta_1 <- theta[M_1+M_2+1]
  beta_2 <- theta[M_1+M_2+2]
  tau <- theta[M_1+M_2+3]
  a_1 <- sum(mu_vec_1^beta_1)
  a_2 <- sum(mu_vec_2^beta_2)
  q_4 <- (mu_vec_1^beta_1)/a_1
  q_5 <- (mu_vec_2^beta_2)/a_2
  q_6 <- outer(q_4,q_5)
  c <- 0
  for(j_1 in 1:(M_1 + 1))
  {
    for(j_2 in 1:(M_2 + 1))
    {
      c <- c + (Freq_table[j_1,j_2] -  sum(sapply(monitoring_data[,1],prob_cell_wei_sh,
      mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,beta_1=beta_1,beta_2=beta_2,tau=tau,
      q_1=q_4,q_2=q_5,q_3=q_6,j_1=j_1,j_2=j_2)))^2
    }
  }
  return(c)
}
in_est_wei_sh_1 <- function(theta,monitoring_data,Freq_table,M_1,M_2)
{
  mu_vec_1 <- theta[1:M_1]
  mu_vec_2 <- theta[(M_1 + 1):(M_1 + M_2)]
  beta_1 <- theta[M_1+M_2+1]
  beta_2 <- theta[M_1+M_2+2]
  tau <- theta[M_1+M_2+3]
  a_1 <- sum(mu_vec_1^beta_1)
  a_2 <- sum(mu_vec_2^beta_2)
  q_4 <- (mu_vec_1^beta_1)/a_1
  q_5 <- (mu_vec_2^beta_2)/a_2
  q_6 <- outer(q_4,q_5)
  foreach(j_1 = 1:(M_1 + 1),.combine = sum) %do% {
    foreach(j_2 = 1:(M_2 + 1),.combine = sum) %do% {
      (Freq_table[j_1,j_2] -  sum(sapply(monitoring_data[,1],prob_cell_wei_sh,
      mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,beta_1=beta_1,beta_2=beta_2,tau=tau,
      q_1=q_4,q_2=q_5,q_3=q_6,j_1=j_1,j_2=j_2)))^2
    }
  }
  return(c)
}
weibull_indiv_log_lik <- function(data,mu_1,mu_2,q_1,q_2,
                                  q_3,gamma_1,gamma_2,tau)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  surv_1 <- (1 + (tau*(mu_1)*(x^gamma_1)))^(-(1/tau))
  surv_2 <- (1 + (tau*(mu_2)*(x^gamma_2)))^(-(1/tau))
  surv_joint <- (1 + tau*((mu_1*(x^gamma_1)) + 
                            (mu_2*(x^gamma_2))))^(-(1/tau))
  c <- 1 - (surv_1) - (surv_2) + surv_joint
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
## the overall negative log-likelihood function which is 
## to minimize
neg_log_likelihood_weibull <- function(theta,D,M_1,M_2)
{
  lambda_vec_3 <- theta[1:M_1]
  lambda_vec_4 <- theta[(M_1 + 1):(M_1 + M_2)]
  beta_1 <- theta[M_1 + M_2 + 1]
  beta_2 <- theta[M_1 + M_2 + 2]
  sigma_1 <- theta[M_1 + M_2 + 3]
  lambda_3 <- sum(lambda_vec_3)
  lambda_4 <- sum(lambda_vec_4)
  q_4 <- (lambda_vec_3)/lambda_3
  q_5 <- (lambda_vec_4)/lambda_4
  q_6 <- outer(q_4,q_5)
  return(-sum(apply(D,1,weibull_indiv_log_lik,
                    mu_1=lambda_3,mu_2=lambda_4,q_1=q_4,q_2=q_5,q_3=q_6,
                    gamma_1=beta_1,gamma_2=beta_2,tau=sigma_1)))
}
wei_sh_indiv_gr <- function(data,mu_vec_1,mu_vec_2,beta_1,beta_2,
                            mu_1,mu_2,a_1,a_2,tau,M_1,M_2)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  if((j_1 == 0) && (j_2 == 0))
  {
    Q_1 <- -M_1*M_2
    Q_2 <- (x^beta_1)
    Q_3 <- (x^beta_2)
    d_1 <- tau*Q_2*a_1
    d_2 <- tau*Q_3*a_2
    e_1 <- sum((mu_vec_1^beta_1)*log(mu_vec_1))
    e_2 <- sum((mu_vec_2^beta_2)*log(mu_vec_2))
    H_4 <- 1 + d_1 + d_2
    G_1 <- Q_1*beta_1*Q_2*(mu_vec_1^(beta_1 - 1))
    G_2 <- Q_1*beta_2*Q_3*(mu_vec_2^(beta_2 - 1))
    G_3 <- Q_1*Q_2*(a_1*log(x) + e_1)
    G_4 <- Q_1*Q_3*(a_2*log(x) + e_2)
    G_5 <- Q_1*(Q_2*a_1 + Q_3*a_2)*(1/tau)
    return(c(G_1,G_2,G_3,G_4,G_5)/H_4)
  }else if((j_1 == 0) && (j_2 != 0))
  {
    Q_1 <- (x^beta_1)
    Q_2 <- (x^beta_2)
    d_1 <- tau*Q_1*a_1
    d_2 <- tau*Q_2*a_2
    d_3 <- d_1 + d_2
    e_1 <- sum((mu_vec_1^beta_1)*log(mu_vec_1))
    e_2 <- sum((mu_vec_2^beta_2)*log(mu_vec_2))
    Q_3 <- (1 + d_1)^(-(1/tau))
    Q_4 <- (1 + d_3)^(-(1/tau))
    H_2 <- Q_3 - Q_4
    Q_5 <- (1 + d_1)^(-(1/tau) - 1)
    Q_6 <- (1 + d_3)^(-(1/tau) - 1)
    Q_7 <- M_1/mu_2
    Q_8 <- M_1/H_2
    G_1 <- (-Q_5 + Q_6)*beta_1*Q_1*(mu_vec_1^(beta_1 - 1))*Q_8
    G_2 <- Q_6*beta_2*Q_2*(mu_vec_2^(beta_2 - 1))*Q_8
    G_3 <- (-Q_5 + Q_6)*Q_1*(a_1*log(x) + e_1)*Q_8
    G_4 <- Q_6*Q_2*(a_2*log(x) + e_2)*Q_8
    G_5 <- ((Q_3)*((1/tau^2)*log(1 + d_1) - ((Q_1*a_1)/(tau*(1 + d_1)))) -
           (Q_4)*((log(1 + d_3)/(tau^2)) - 
           ((Q_1*a_1 + Q_2*a_2)/(tau*(1 + d_3)))))*Q_8
    if(j_2 == 1)
    {
      return(c(G_1,G_2[j_2] + (M_1/mu_vec_2[j_2]) - Q_7,
               G_2[2:M_2] - Q_7,G_3,G_4,G_5))
    }else if(j_2 == M_2)
    {
      return(c(G_1,G_2[1:(M_2 - 1)] - Q_7,
               G_2[j_2] + (M_1/mu_vec_2[j_2]) - Q_7,G_3,G_4,G_5))
    }else
    {
      return(c(G_1,G_2[1:(j_2 - 1)] - Q_7,
               G_2[j_2] + (M_1/mu_vec_2[j_2]) - Q_7,
               G_2[(j_2 + 1):M_2] - Q_7,G_3,G_4,G_5))
    }
  }else if((j_1 != 0) && (j_2 == 0))
  {
    Q_1 <- (x^beta_1)
    Q_2 <- (x^beta_2)
    d_1 <- tau*Q_1*a_1
    d_2 <- tau*Q_2*a_2
    d_3 <- d_1 + d_2
    e_1 <- sum((mu_vec_1^beta_1)*log(mu_vec_1))
    e_2 <- sum((mu_vec_2^beta_2)*log(mu_vec_2))
    Q_3 <- (1 + d_2)^(-(1/tau))
    Q_4 <- (1 + d_3)^(-(1/tau))
    H_3 <- Q_3 - Q_4
    Q_5 <- (1 + d_2)^(-(1/tau) - 1)
    Q_6 <- (1 + d_3)^(-(1/tau) - 1)
    Q_7 <- M_2/mu_1
    Q_8 <- M_2/H_3
    G_1 <- Q_6*beta_1*Q_1*(mu_vec_1^(beta_1 - 1))*Q_8
    G_2 <- (-Q_5 + Q_6)*beta_2*Q_2*(mu_vec_2^(beta_2 - 1))*Q_8
    G_3 <- Q_6*Q_1*(a_1*log(x) + e_1)*Q_8
    G_4 <- (-Q_5 + Q_6)*Q_2*(a_2*log(x) + e_2)*Q_8
    G_5 <- ((Q_3)*((1/tau^2)*log(1 + d_2) - ((Q_2*a_2)/(tau*(1 + d_2)))) -
              (Q_4)*((log(1 + d_3)/(tau^2)) - 
              ((Q_1*a_1 + Q_2*a_2)/(tau*(1 + d_3)))))*Q_8
    if(j_1 == 1)
    {
      return(c(G_1[j_1] + (M_2/mu_vec_1[j_1]) - Q_7,
               G_1[2:M_1] - Q_7,G_2,G_3,G_4,G_5))
    }else if(j_1 == M_1)
    {
      return(c(G_1[1:(M_1 - 1)] - Q_7,G_1[j_1] + (M_2/mu_vec_1[j_1]) - Q_7,
               G_2,G_3,G_4,G_5))
    }else
    {
      return(c(G_1[1:(j_1 - 1)] - Q_7,
             G_1[j_1] + (M_2/mu_vec_1[j_1]) - Q_7,
             G_1[(j_1 + 1):M_1] - Q_7,G_2,G_3,G_4,G_5))
    }
  }else
  {
    Q_1 <- (x^beta_1)
    Q_2 <- (x^beta_2)
    d_1 <- tau*Q_1*a_1
    d_2 <- tau*Q_2*a_2
    d_3 <- d_1 + d_2
    e_1 <- sum((mu_vec_1^beta_1)*log(mu_vec_1))
    e_2 <- sum((mu_vec_2^beta_2)*log(mu_vec_2))
    Q_3 <- (1 + d_1)^(-(1/tau))
    Q_4 <- (1 + d_2)^(-(1/tau))
    Q_5 <- (1 + d_3)^(-(1/tau))
    H_1 <- 1 - Q_3 - Q_4 + Q_5
    Q_6 <- (1 + d_1)^(-(1/tau) - 1)
    Q_7 <- (1 + d_2)^(-(1/tau) - 1)
    Q_8 <- (1 + d_3)^(-(1/tau) - 1)
    Q_9 <- 1/mu_1
    Q_10 <- 1/mu_2
    Q_11 <- 1/H_1
    G_1 <- (Q_6 - Q_8)*beta_1*Q_1*(mu_vec_1^(beta_1 - 1))*Q_11
    G_2 <- (Q_7 - Q_8)*beta_2*Q_2*(mu_vec_2^(beta_2 - 1))*Q_11
    G_3 <- (Q_6 - Q_8)*Q_1*(a_1*log(x) + e_1)*Q_11
    G_4 <- (Q_7 - Q_8)*Q_2*(a_2*log(x) + e_2)*Q_11
    G_5 <- (-Q_3*((1/tau^2)*log(1 + d_1) - ((Q_1*a_1)/(tau*(1 + d_1)))) -
            Q_4*((1/tau^2)*log(1 + d_2) - ((Q_2*a_2)/(tau*(1 + d_2)))) +
            Q_5*((log(1 + d_3)/(tau^2)) - 
           ((Q_1*a_1 + Q_2*a_2)/(tau*(1 + d_3)))))*Q_11
    if((j_1 == 1) & (j_2 == 1))
    {
      return(c(G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,G_1[2:M_1] - Q_9,
             G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,G_2[2:M_2] - Q_10,
             G_3,G_4,G_5))
    }else if((j_1 == M_1) & (j_2 == 1))
    {
      return(c(G_1[1:(M_1 - 1)] - Q_9,G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,
               G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,G_2[2:M_2] - Q_10,
               G_3,G_4,G_5))
    }else if((j_1 == 1) & (j_2 == M_2))
    {
      return(c(G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,G_1[2:M_1] - Q_9,
               G_2[1:(M_2 - 1)] - Q_10,G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,
               G_3,G_4,G_5))
    }else if((j_1 == M_1) & (j_2 == M_2))
    {
      return(c(G_1[1:(M_1 - 1)] - Q_9,G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,
               G_2[1:(M_2 - 1)] - Q_10,G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,
               G_3,G_4,G_5))
    }else
    {
      return(c(G_1[1:(j_1 - 1)] - Q_9,G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,
               G_1[(j_1 + 1):M_1] - Q_9,
               G_2[1:(j_2 - 1)] - Q_10,G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,
               G_2[(j_2 + 1):M_2] - Q_10,
               G_3,G_4,G_5))
    }
  }
}
wei_sh_indiv_gr_1 <- function(data,mu_vec_1,mu_vec_2,beta_1,beta_2,
                            mu_1,mu_2,a_1,a_2,tau,M_1,M_2)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  if((j_1 == 0) && (j_2 == 0))
  {
    Q_1 <- -M_1*M_2
    Q_2 <- (x^beta_1)
    Q_3 <- (x^beta_2)
    d_1 <- tau*Q_2*a_1
    d_2 <- tau*Q_3*a_2
    e_1 <- sum((mu_vec_1^beta_1)*log(mu_vec_1))
    e_2 <- sum((mu_vec_2^beta_2)*log(mu_vec_2))
    H_4 <- 1 + d_1 + d_2
    G_1 <- Q_1*beta_1*Q_2*(mu_vec_1^(beta_1 - 1))
    G_2 <- Q_1*beta_2*Q_3*(mu_vec_2^(beta_2 - 1))
    G_3 <- Q_1*Q_2*(a_1*log(x) + e_1)
    G_4 <- Q_1*Q_3*(a_2*log(x) + e_2)
    G_5 <- Q_1*(Q_2*a_1 + Q_3*a_2)*(1/tau)
    return(c(G_1,G_2,G_3,G_4,G_5)/H_4)
  }else if((j_1 == 0) && (j_2 != 0))
  {
    Q_1 <- (x^beta_1)
    Q_2 <- (x^beta_2)
    d_1 <- tau*Q_1*a_1
    d_2 <- tau*Q_2*a_2
    d_3 <- d_1 + d_2
    e_1 <- sum((mu_vec_1^beta_1)*log(mu_vec_1))
    e_2 <- sum((mu_vec_2^beta_2)*log(mu_vec_2))
    Q_3 <- (1 + d_1)^(-(1/tau))
    Q_4 <- (1 + d_3)^(-(1/tau))
    H_2 <- Q_3 - Q_4
    Q_5 <- (1 + d_1)^(-(1/tau) - 1)
    Q_6 <- (1 + d_3)^(-(1/tau) - 1)
    Q_7 <- M_1/mu_2
    Q_8 <- M_1/H_2
    G_1 <- (-Q_5 + Q_6)*beta_1*Q_1*(mu_vec_1^(beta_1 - 1))*Q_8
    G_2 <- Q_6*beta_2*Q_2*(mu_vec_2^(beta_2 - 1))*Q_8
    G_3 <- (-Q_5 + Q_6)*Q_1*(a_1*log(x) + e_1)*Q_8
    G_4 <- Q_6*Q_2*(a_2*log(x) + e_2)*Q_8
    G_5 <- ((Q_3)*((1/tau^2)*log(1 + d_1) - ((Q_1*a_1)/(tau*(1 + d_1)))) -
              (Q_4)*((log(1 + d_3)/(tau^2)) - 
                       ((Q_1*a_1 + Q_2*a_2)/(tau*(1 + d_3)))))*Q_8
    if(j_2 == 1)
    {
      return(c(G_1,G_2[j_2] + (M_1/mu_vec_2[j_2]) - Q_7,
               G_2[2:M_2] - Q_7,G_3,G_4,G_5))
    }else if(j_2 == M_2)
    {
      return(c(G_1,G_2[1:(M_2 - 1)] - Q_7,
               G_2[j_2] + (M_1/mu_vec_2[j_2]) - Q_7,G_3,G_4,G_5))
    }else
    {
      return(c(G_1,G_2[1:(j_2 - 1)] - Q_7,
               G_2[j_2] + (M_1/mu_vec_2[j_2]) - Q_7,
               G_2[(j_2 + 1):M_2] - Q_7,G_3,G_4,G_5))
    }
  }else if((j_1 != 0) && (j_2 == 0))
  {
    Q_1 <- (x^beta_1)
    Q_2 <- (x^beta_2)
    d_1 <- tau*Q_1*a_1
    d_2 <- tau*Q_2*a_2
    d_3 <- d_1 + d_2
    e_1 <- sum((mu_vec_1^beta_1)*log(mu_vec_1))
    e_2 <- sum((mu_vec_2^beta_2)*log(mu_vec_2))
    Q_3 <- (1 + d_2)^(-(1/tau))
    Q_4 <- (1 + d_3)^(-(1/tau))
    H_3 <- Q_3 - Q_4
    Q_5 <- (1 + d_2)^(-(1/tau) - 1)
    Q_6 <- (1 + d_3)^(-(1/tau) - 1)
    Q_7 <- M_2/mu_1
    Q_8 <- M_2/H_3
    G_1 <- Q_6*beta_1*Q_1*(mu_vec_1^(beta_1 - 1))*Q_8
    G_2 <- (-Q_5 + Q_6)*beta_2*Q_2*(mu_vec_2^(beta_2 - 1))*Q_8
    G_3 <- Q_6*Q_1*(a_1*log(x) + e_1)*Q_8
    G_4 <- (-Q_5 + Q_6)*Q_2*(a_2*log(x) + e_2)*Q_8
    G_5 <- ((Q_3)*((1/tau^2)*log(1 + d_2) - ((Q_2*a_2)/(tau*(1 + d_2)))) -
              (Q_4)*((log(1 + d_3)/(tau^2)) - 
                       ((Q_1*a_1 + Q_2*a_2)/(tau*(1 + d_3)))))*Q_8
    if(j_1 == 1)
    {
      return(c(G_1[j_1] + (M_2/mu_vec_1[j_1]) - Q_7,
               G_1[2:M_1] - Q_7,G_2,G_3,G_4,G_5))
    }else if(j_1 == M_1)
    {
      return(c(G_1[1:(M_1 - 1)] - Q_7,G_1[j_1] + (M_2/mu_vec_1[j_1]) - Q_7,
               G_2,G_3,G_4,G_5))
    }else
    {
      return(c(G_1[1:(j_1 - 1)] - Q_7,
               G_1[j_1] + (M_2/mu_vec_1[j_1]) - Q_7,
               G_1[(j_1 + 1):M_1] - Q_7,G_2,G_3,G_4,G_5))
    }
  }else
  {
    Q_1 <- (x^beta_1)
    Q_2 <- (x^beta_2)
    d_1 <- tau*Q_1*a_1
    d_2 <- tau*Q_2*a_2
    d_3 <- d_1 + d_2
    e_1 <- sum((mu_vec_1^beta_1)*log(mu_vec_1))
    e_2 <- sum((mu_vec_2^beta_2)*log(mu_vec_2))
    Q_3 <- (1 + d_1)^(-(1/tau))
    Q_4 <- (1 + d_2)^(-(1/tau))
    Q_5 <- (1 + d_3)^(-(1/tau))
    H_1 <- 1 - Q_3 - Q_4 + Q_5
    Q_6 <- (1 + d_1)^(-(1/tau) - 1)
    Q_7 <- (1 + d_2)^(-(1/tau) - 1)
    Q_8 <- (1 + d_3)^(-(1/tau) - 1)
    Q_9 <- 1/mu_1
    Q_10 <- 1/mu_2
    Q_11 <- 1/H_1
    G_1 <- (Q_6 - Q_8)*beta_1*Q_1*(mu_vec_1^(beta_1 - 1))*Q_11
    G_2 <- (Q_7 - Q_8)*beta_2*Q_2*(mu_vec_2^(beta_2 - 1))*Q_11
    G_3 <- (Q_6 - Q_8)*Q_1*(a_1*log(x) + e_1)*Q_11
    G_4 <- (Q_7 - Q_8)*Q_2*(a_2*log(x) + e_2)*Q_11
    G_5 <- (-Q_3*((1/tau^2)*log(1 + d_1) - ((Q_1*a_1)/(tau*(1 + d_1)))) -
              Q_4*((1/tau^2)*log(1 + d_2) - ((Q_2*a_2)/(tau*(1 + d_2)))) +
              Q_5*((log(1 + d_3)/(tau^2)) - 
                     ((Q_1*a_1 + Q_2*a_2)/(tau*(1 + d_3)))))*Q_11
    if(j_1 == 1)
    {
      if(j_2 == 1)
      {
        return(c(G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,G_1[2:M_1] - Q_9,
                 G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,G_2[2:M_2] - Q_10,
                 G_3,G_4,G_5))
      }else if(j_2 == M_2)
      {
        return(c(G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,G_1[2:M_1] - Q_9,
                 G_2[1:(M_2 - 1)] - Q_10,G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,
                 G_3,G_4,G_5))
      }else
      {
        return(c(G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,G_1[2:M_1] - Q_9,
                 G_2[1:(j_2 - 1)] - Q_10,G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,
                 G_2[(j_2 + 1):M_2] - Q_10,
                 G_3,G_4,G_5))
      }
    }else if(j_1 == M_1)
    {
      if(j_2 == 1)
      {
        return(c(G_1[1:(M_1 - 1)] - Q_9,G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,
                 G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,G_2[2:M_2] - Q_10,
                 G_3,G_4,G_5))
      }else if(j_2 == M_2)
      {
        return(c(G_1[1:(M_1 - 1)] - Q_9,G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,
                 G_2[1:(M_2 - 1)] - Q_10,G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,
                 G_3,G_4,G_5))
      }else
      {
        return(c(G_1[1:(M_1 - 1)] - Q_9,G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,
                 G_2[1:(j_2 - 1)] - Q_10,G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,
                 G_2[(j_2 + 1):M_2] - Q_10,
                 G_3,G_4,G_5))
      }
    }else
    {
      if(j_2 == 1)
      {
        return(c(G_1[1:(j_1 - 1)] - Q_9,G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,
                 G_1[(j_1 + 1):M_1] - Q_9,
                 G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,G_2[2:M_2] - Q_10,
                 G_3,G_4,G_5))
      }else if(j_2 == M_2)
      {
        return(c(G_1[1:(j_1 - 1)] - Q_9,G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,
                 G_1[(j_1 + 1):M_1] - Q_9,
                 G_2[1:(M_2 - 1)] - Q_10,G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,
                 G_3,G_4,G_5))
      }else
      {
        return(c(G_1[1:(j_1 - 1)] - Q_9,G_1[j_1] + 1/mu_vec_1[j_1] - Q_9,
                 G_1[(j_1 + 1):M_1] - Q_9,
                 G_2[1:(j_2 - 1)] - Q_10,G_2[j_2] + 1/mu_vec_2[j_2] - Q_10,
                 G_2[(j_2 + 1):M_2] - Q_10,
                 G_3,G_4,G_5))
      }
    }
  }  
}
## negative gradient for weibull hazard shared gamma frailty model
neg_gr_wei_sh <- function(D,theta,M_1,M_2)
{
  mu_vec_3 <- theta[1:M_1]
  mu_vec_4 <- theta[(M_1 + 1):(M_1 + M_2)]
  beta_3 <- theta[M_1 + M_2 + 1]
  beta_4 <- theta[M_1 + M_2 + 2]
  sigma_1 <- theta[M_1 + M_2 + 3]
  mu_3 <- sum(mu_vec_3)
  mu_4 <- sum(mu_vec_4)
  a_3 <- sum(mu_vec_3^beta_3)
  a_4 <- sum(mu_vec_4^beta_4)
  m <- -apply(D,1,wei_sh_indiv_gr_1,tau=sigma_1,mu_vec_1=mu_vec_3,
       mu_vec_2=mu_vec_4,mu_1=mu_3,mu_2=mu_4,a_1=a_3,a_2=a_4,
       beta_1=beta_3,beta_2=beta_4,M_1=M_1,M_2=M_2)
  if(typeof(m) != "list")
  {
    return(apply(m,1,sum))
  }else
  {
    m_1 <- matrix(unlist(m),ncol=M_1+M_2+3,byrow=TRUE)
    return(apply(m_1,2,sum))
  }
}
## numerical gradient
num_gr_wei_shared <- function(D,theta,M_1,M_2)
{
  return(pracma::grad(neg_log_likelihood_weibull,theta,D=D,M_1=M_1,
                        M_2=M_2))
}
## gradient of prob_cell function for weibull hazard and shared gamma frailty
gr_prob_cell_wei_sh <- function(x,j_1,j_2,mu_vec_1,mu_vec_2,beta_1,beta_2,
                         mu_1,mu_2,a_1,a_2,tau,M_1,M_2)
{
  q_1 <- (mu_vec_1^beta_1)/a_1
  q_2 <- (mu_vec_2^beta_2)/a_2
  q_3 <- outer(q_1,q_2)
  p <- prob_cell_wei_sh(x,mu_vec_1,mu_vec_2,beta_1,beta_2,tau,q_1,q_2,q_3,j_1,j_2)
  g <- -wei_sh_indiv_gr_1(c(x,j_1 - 1,j_2 - 1),mu_vec_1,mu_vec_2,beta_1,
                       beta_2,mu_1,mu_2,a_1,a_2,tau,M_1,M_2)
  if(typeof(g) != "list")
  {
    return(p*g)
  }else
  {
    return(p*(matrix(unlist(g),nrow=1,ncol = length(mu_vec_1) + 
                       length(mu_vec_2) + 3)))
  }
}
## gradient of function in_est_wei_sh for shared frailty model
gr_in_est_wei_sh <- function(monitoring_data,Freq_table,M_1,M_2,theta)
{
  mu_vec_1 <- theta[1:M_1]
  mu_vec_2 <- theta[(M_1 + 1):(M_1 + M_2)]
  beta_1 <- theta[M_1 + M_2 + 1]
  beta_2 <- theta[M_1 + M_2 + 2]
  tau <- theta[M_1+M_2+3]
  mu_1 <- sum(mu_vec_1)
  mu_2 <- sum(mu_vec_2)
  a_1 <- sum(mu_vec_1^beta_1)
  a_2 <- sum(mu_vec_2^beta_2) 
  q_4 <- (mu_vec_1^beta_1)/a_1
  q_5 <- (mu_vec_2^beta_2)/a_2
  q_6 <- outer(q_4,q_5)
  c_1 <- replicate(M_1,0)
  c_2 <- replicate(M_2,0)
  c_3 <- 0
  c_4 <- 0
  c_5 <- 0
  d <- matrix(0,nrow = M_1 + 1,ncol = M_2 + 1)
  for(j_1 in 1:(M_1 + 1))
  {
    for(j_2 in 1:(M_2 + 1))
    {
      d[j_1,j_2] <- 2*(Freq_table[j_1,j_2] - sum(sapply(monitoring_data[,1],
                    prob_cell_wei_sh,mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,tau=tau,
                    beta_1=beta_1,beta_2=beta_2,q_1=q_4,q_2=q_5,q_3=q_6,j_1=j_1,j_2=j_2)))
      m <- sapply(monitoring_data[,1],gr_prob_cell_wei_sh,mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,
           beta_1=beta_1,beta_2=beta_2,mu_1=mu_1,mu_2=mu_2,a_1=a_1,a_2=a_2,tau=tau,M_1=M_1,M_2=M_2,j_1=j_1,j_2=j_2)
      m_1 <- matrix(m[1:M_1,],nrow=M_1,ncol=length(monitoring_data[,1]))
      c_1 <- c_1 + d[j_1,j_2]*apply(m_1,1,sum)
      c_3 <- c_3 + d[j_1,j_2]*sum(m[M_1 + M_2 + 1,])
      c_4 <- c_4 + d[j_1,j_2]*sum(m[M_1 + M_2 + 2,])
      c_5 <- c_5 + d[j_1,j_2]*sum(m[M_1 + M_2 + 3,])
    }
  }
  for(j_2 in 1:(M_2 + 1))
  {
    for(j_1 in 1:(M_1 + 1))
    {
      m <- sapply(monitoring_data[,1],gr_prob_cell_wei_sh,mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,
           beta_1=beta_1,beta_2=beta_2,mu_1=mu_1,mu_2=mu_2,a_1=a_1,a_2=a_2,tau=tau,M_1=M_1,M_2=M_2,j_1=j_1,j_2=j_2)
      m_2 <- matrix(m[(M_1 + 1):(M_1 + M_2),],nrow=M_2,ncol=length(monitoring_data[,1]))
      c_2 <- c_2 + d[j_1,j_2]*apply(m_2,1,sum)
    }
  }
  return(c(c_1,c_2,c_3,c_4,c_5))
}
library("doParallel")
cl <- makeCluster(3)
setDefaultCluster(cl=cl)
clusterExport(cl,c("prob_cell_wei_sh","joint_surv_wei",
"marg_surv_wei"))
in_est <- optimParallel::optimParallel(par=c(lambda_vec_1_initial,lambda_vec_2_initial,c(0.1,0.7,1)),
fn=in_est_wei_sh,monitoring_data = hearing_data,Freq_table=F,
M_1 = 3,M_2 = 3,lower=replicate(9,0.01),method="L-BFGS-B",control=list(trace=6),
parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE))
stopCluster(cl)
setDefaultCluster(cl=NULL)
in_est$par
in_est$convergence
in_est_wei_sh(in_est$par,hearing_data,F,3,3)
in_est_wei_sh(c(lambda_vec_1_initial,lambda_vec_2_initial,c(0.1,0.7,1)),
              hearing_data,F,3,3)
in_est_data <- in_est$par
cl <- makeCluster(3)
setDefaultCluster(cl=cl)
clusterExport(cl,"weibull_indiv_log_lik")
data_soln <- optimParallel::optimParallel(par=in_est_data,
fn=neg_log_likelihood_weibull,M_1=3,M_2=3,D=hearing_data,lower=replicate(9,0.01),
parallel = list(cl=NULL,forward=FALSE,loginfo = TRUE),control = list(trace=6),
method = "L-BFGS-B")$par
stopCluster(cl)
setDefaultCluster(cl=NULL)
hess <- pracma::hessian(neg_log_likelihood_weibull,data_soln,D=hearing_data,M_1=3,M_2=3)
neg_log_likelihood_weibull(in_est_data,hearing_data,3,3)
neg_log_likelihood_weibull(data_soln,hearing_data,3,3)
eigen(hess)$values
sd_mle <- diag(solve(hess))
