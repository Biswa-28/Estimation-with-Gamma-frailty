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

  g_1 <- function(x,tau,mu_1,mu_2)
{
  return(x/(1 + tau*((mu_1*x) + (mu_2*x))))
}
g_2 <- function(x,tau,mu_1,mu_2)
{
  a_1 <- (mu_1*x)
  a_2 <- (mu_2*x)
  return((a_1 + a_2)/((1 + tau*(a_1 + a_2))))
}
marg_surv <- function(x,tau,mu)
{
  return((1 + tau*mu*x)^(-(1/tau)))
}
joint_surv <- function(x,tau,mu_1,mu_2)
{
  return((1 + tau*x*(mu_1+mu_2))^(-(1/tau)))
}
h_1 <- function(x,tau,mu_1,mu_2)
{
  return(marg_surv(x,tau,mu_1) -
           joint_surv(x,tau,mu_1,mu_2))
}
h_2 <- function(x,tau,mu_1,mu_2)
{
  return(marg_surv(x,tau,mu_2) -
           joint_surv(x,tau,mu_1,mu_2))
}
h_3 <- function(x,tau,mu_1,mu_2)
{
  return(1 - marg_surv(x,tau,mu_1) - marg_surv(x,tau,mu_2)
         + joint_surv(x,tau,mu_1,mu_2))
}
g_3 <- function(x,tau,mu)
{
  return(((1 + tau*x*mu)^(-1 - (1/tau))))
}
g_4 <- function(x,tau,mu_1,mu_2)
{
  return((1 + tau*x*(mu_1+mu_2))^(-1-(1/tau)))
}
g_5 <- function(x,tau,mu)
{
  a <- x*tau*mu
  return((log(1 + a)/(tau^2)) -
           ((x*mu)/(tau*(1 + a))))
}
g_6 <- function(x,tau,mu_1,mu_2)
{
  a <- x*tau*(mu_1 + mu_2)
  return(log(1 + a)/(tau^2) -
           (x*(mu_1 + mu_2))/(tau*(1 + a)))
}
g_7 <- function(x,tau,mu)
{
  a <- x*mu
  return(((1/(tau^2))*log(1 + tau*a)) - (1 + (1/tau))*(a/(1 + tau*a)))
}
g_8 <- function(x,tau,mu_1,mu_2)
{
  a <- x*(mu_1 + mu_2)
  return((1/(tau^2))*log(1 + tau*a) - (1 + (1/tau))*(a/(1 + tau*a)))
}
g_9 <- function(x,tau,mu)
{
  a <- x*mu
  return(-(2/(tau^3))*log(1 + tau*a) + 2*(1/(tau^2))*(a/(1 + tau*a))
         + (1/tau)*((a/(1 + tau*a))^2))
}
g_10 <- function(x,tau,mu_1,mu_2)
{
  a <- x*(mu_1 + mu_2)
  return(-(2/(tau^3))*log(1 + tau*a) + 2*(1/(tau^2))*(a/(1 + tau*a))
         + (1/tau)*((a/(1 + tau*a))^2))
}
  ## individual gradient function for shared frailty model
indiv_gr <- function(data,tau,mu_vec_1,mu_vec_2,
                     mu_1,mu_2,L_1,L_2)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  if((j_1 == 0) && (j_2 == 0))
  {
    G_1 <- g_1(x,tau,mu_1,mu_2)
    return(L_1*L_2*c(replicate(L_1 + L_2,G_1),
                     (1/tau)*g_2(x,tau,mu_1,mu_2) - 
                       (log(1 + tau*(mu_1 + mu_2)*x))/(tau^2)))
  }else if((j_1 == 0) && (j_2 != 0))
  {
    H_1 <- h_1(x,tau,mu_1,mu_2)
    Q_1 <- (L_1*(g_3(x,tau,mu_1)*x - g_4(x,tau,mu_1,mu_2)*x))/H_1
    Q_2 <- -(L_1*(g_4(x,tau,mu_1,mu_2)*x))/H_1
    Q_3 <- (L_1/mu_2)
    Q_4 <- L_1*(-g_5(x,tau,mu_1)*marg_surv(x,tau,mu_1) + 
                  g_6(x,tau,mu_1,mu_2)*joint_surv(x,tau,mu_1,mu_2))/H_1
    if(j_2 == 1)
    {
      return(c(replicate(L_1,Q_1),
               Q_2 - (L_1/mu_vec_2[j_2]) + Q_3,
               replicate(L_2 - 1,Q_2 + Q_3),Q_4))
    }else if(j_2 == L_2)
    {
      return(c(replicate(L_1,Q_1),
               replicate(L_2 - 1,Q_2 + Q_3),
               Q_2 - (L_1/mu_vec_2[j_2]) + Q_3,Q_4))
    }else
    {
      return(c(replicate(L_1,Q_1),
               replicate(j_2 - 1,Q_2 + Q_3),
               Q_2 - (L_1/mu_vec_2[j_2]) + Q_3,
               replicate(L_2 - j_2,Q_2 + Q_3),Q_4))
    }
  }else if((j_1 != 0) && (j_2 == 0))
  {
    H_2 <- h_2(x,tau,mu_1,mu_2)
    Q_1 <- (-L_2*g_4(x,tau,mu_1,mu_2)*x)/H_2
    Q_2 <- (L_2*(g_3(x,tau,mu_2)*x - g_4(x,tau,mu_1,mu_2)*x))/H_2
    Q_3 <- (L_2/mu_1)
    Q_4 <- L_2*(-g_5(x,tau,mu_2)*marg_surv(x,tau,mu_2) + 
                  g_6(x,tau,mu_1,mu_2)*joint_surv(x,tau,mu_1,mu_2))/H_2
    if(j_1 == 1)
    {
      return(c(Q_1 - (L_2/mu_vec_1[j_1]) + Q_3,
               replicate(L_1 - 1,Q_1 + Q_3),
               replicate(L_2,Q_2),Q_4))
    }else if(j_1 == L_1)
    {
      return(c(replicate(L_1 - 1,Q_1 + Q_3),
               Q_1 - (L_2/mu_vec_1[j_1]) + Q_3,
               replicate(L_2,Q_2),Q_4))
    }else
    {
      return(c(replicate(j_1 - 1,Q_1 + Q_3),
               Q_1 - (L_2/mu_vec_1[j_1]) + Q_3,
               replicate(L_1 - j_1,Q_1 + Q_3),
               replicate(L_2,Q_2),Q_4))
    }
  }else
  {
    H_3 <- h_3(x,tau,mu_1,mu_2)
    Q_1 <- (-g_3(x,tau,mu_1)*x + g_4(x,tau,mu_1,mu_2)*x)/H_3
    Q_2 <- (-g_3(x,tau,mu_2)*x + g_4(x,tau,mu_1,mu_2)*x)/H_3
    Q_3 <- (1/mu_1)
    Q_4 <- (1/mu_2)
    Q_5 <- (g_5(x,tau,mu_1)*marg_surv(x,tau,mu_1) 
            + g_5(x,tau,mu_2)*marg_surv(x,tau,mu_2)
            - g_6(x,tau,mu_1,mu_2)*joint_surv(x,tau,mu_1,mu_2))/H_3
    if((j_1 == 1) && (j_2 == 1))
    {
      return(c(Q_1 - (1/mu_vec_1[j_1]) + Q_3,
               replicate(L_1 - 1,Q_1 + Q_3),
               Q_2 - (1/mu_vec_2[j_2]) + Q_4,
               replicate(L_2 - 1,Q_2 + Q_4),Q_5))
    }else if((j_1 == 1) && (j_2 == L_2))
    {
      return(c(Q_1 - (1/mu_vec_1[j_1]) + Q_3,
               replicate(L_1 - 1,Q_1 + Q_3),
               replicate(L_2 - 1,Q_2 + Q_4),
               Q_2 - (1/mu_vec_2[j_2]) + Q_4,Q_5))
    }else if((j_1 == L_1) && (j_2 == 1))
    {
      return(c(replicate(L_1 - 1,Q_1 + Q_3),
               Q_1 - (1/mu_vec_1[j_1]) + Q_3,
               Q_2 - (1/mu_vec_2[j_2]) + Q_4,
               replicate(L_2 - 1,Q_2 + Q_4),Q_5))
    }else if((j_1 == L_1) && (j_2 == L_2))
    {
      return(c(replicate(L_1 - 1,Q_1 + Q_3),
               Q_1 - (1/mu_vec_1[j_1]) + Q_3,
               replicate(L_2 - 1,Q_2 + Q_4),
               Q_2 - (1/mu_vec_2[j_2]) + Q_4,Q_5))
    }else
    {
      return(c(replicate(j_1 - 1,Q_1 + Q_3),
               Q_1 - (1/mu_vec_1[j_1]) + Q_3,
               replicate(L_1 - j_1,Q_1 + Q_3),
               replicate(j_2 - 1,Q_2 + Q_4),
               Q_2 - (1/mu_vec_1[j_2]) + Q_4,
               replicate(L_2 - j_2,Q_2 + Q_4),Q_5))
    }
  }
}
## numerical gradient of the negative log-likelihood function using numDeriv package
  num_gradient <- function(D,theta,M_1,M_2)
{
  return(numDeriv::grad(neg_log_likelihood_exp,theta,D=D,M_1=M_1,
                        M_2=M_2))
}
## gradient of prob_cell function for shared frailty
gr_prob_cell <- function(x,j_1,j_2,mu_vec_1,mu_vec_2,
                         mu_1,mu_2,tau,M_1,M_2)
{
  q_1 <- mu_vec_1/mu_1
  q_2 <- mu_vec_2/mu_2
  q_3 <- outer(q_1,q_2)
  p <- prob_cell(x,mu_1,mu_2,tau,q_1,q_2,q_3,j_1,j_2)
  g <- indiv_gr(c(x,j_1 - 1,j_2 - 1),tau,mu_vec_1,mu_vec_2,mu_1,mu_2,M_1,M_2)
  if(typeof(g) != "list")
  {
    return(p*g)
  }else
  {
    return(p*(matrix(unlist(g),nrow=1,ncol = length(mu_vec_1) +
                       length(mu_vec_2) + 1)))
  }
}
## gradient of function in_est_1 for shared frailty model
gr_in_est_1 <- function(monitoring_data,Freq_table,M_1,M_2,theta)
{
  mu_vec_1 <- theta[1:M_1]
  mu_vec_2 <- theta[(M_1 + 1):(M_1 + M_2)]
  tau <- theta[M_1+M_2+1]
  mu_1 <- sum(mu_vec_1)
  mu_2 <- sum(mu_vec_2)
  q_4 <- mu_vec_1/mu_1
  q_5 <- mu_vec_2/mu_2
  q_6 <- outer(q_4,q_5)
  c_1 <- replicate(M_1,0)
  c_2 <- replicate(M_2,0)
  c_3 <- 0
  d <- matrix(0,nrow = M_1 + 1,ncol = M_2 + 1)
  for(j_1 in 1:(M_1 + 1))
  {
    for(j_2 in 1:(M_2 + 1))
    {
      d[j_1,j_2] <- 2*(Freq_table[j_1,j_2] - sum(sapply(monitoring_data[,1],
                    prob_cell,mu_1=mu_1,mu_2=mu_2,tau=tau,q_1=q_4,q_2=q_5,q_3=q_6,
                    j_1=j_1,j_2=j_2)))
      m <- sapply(monitoring_data[,1],gr_prob_cell,mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,
                  mu_1=mu_1,mu_2=mu_2,tau=tau,M_1=M_1,M_2=M_2,j_1=j_1,j_2=j_2)
      m_1 <- matrix(m[1:M_1,],nrow=M_1,ncol=length(monitoring_data[,1]))
      c_1 <- c_1 + d[j_1,j_2]*apply(m_1,1,sum)
      c_3 <- c_3 + d[j_1,j_2]*sum(m[M_1 + M_2 + 1,])
    }
  }
  for(j_2 in 1:(M_2 + 1))
  {
    for(j_1 in 1:(M_1 + 1))
    {
      m <- sapply(monitoring_data[,1],gr_prob_cell,mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,
                  mu_1=mu_1,mu_2=mu_2,tau=tau,M_1=M_1,M_2=M_2,j_1=j_1,j_2=j_2)
      m_2 <- matrix(m[(M_1 + 1):(M_1 + M_2),],nrow=M_2,ncol=length(monitoring_data[,1]))
      c_2 <- c_2 + d[j_1,j_2]*apply(m_2,1,sum)
    }
  }
  return(c(c_1,c_2,c_3))
}
ui <- rbind(c(1,replicate(6,0)),c(0,1,replicate(5,0)),c(0,0,1,0,0,0,0),
            c(0,0,0,1,0,0,0),c(0,0,0,0,1,0,0),c(replicate(5,0),1,0),
            c(replicate(6,0),1))
ci <- replicate(7,0)
 # finding a better initial point for minimization of original
## neg-loglikelihood function
library("doParallel")
cl <- makePSOCKcluster(7)
registerDoParallel(cl)
in_est <- constrOptim(c(lambda_vec_1_initial,lambda_vec_2_initial,3.8),
                      in_est_1,monitoring_data = hearing_data,Freq_table=F,
                      M_1 = 3,M_2 = 3,ui=ui,ci=ci,grad=gr_in_est_1,method="BFGS")
stopCluster(cl)
in_est$par
in_est$convergence
in_est_1(hearing_data,F,3,3,in_est$par)
in_est_1(hearing_data,F,3,3,c(lambda_vec_1_initial,lambda_vec_2_initial,3.8))
in_est_data <- in_est$par
## min in_est_1 at sigma = 3.9,value = 4687.594,taking seq from 0.1 to 4
## Minimization of original neg-loglikelihood function for shared frailty
  library("doParallel")
cl <- makePSOCKcluster(7)
registerDoParallel(cl)
data_soln <- constrOptim(in_est_data,neg_log_likelihood_exp,gradient,
                         M_1=3,M_2=3,D=hearing_data,ui=ui,ci=ci,method = "BFGS")$par
hess <- numDeriv::hessian(neg_log_likelihood_exp,data_soln,D=hearing_data,M_1=3,M_2=3)
stopCluster(cl)
data_soln
neg_log_likelihood_exp(hearing_data,in_est_data,3,3)
neg_log_likelihood_exp(hearing_data,data_soln,3,3)
H <- neg_hessian(hearing_data,data_soln,3,3)
eigen(H)$values
eigen(hess)$values
sd_mle <- diag(solve(hess))
## Computation of AIC for shared frailty model
AIC_shared <- 2*neg_log_likelihood_exp(hearing_data,data_soln,3,3) + (2*7)
AIC_shared

  
