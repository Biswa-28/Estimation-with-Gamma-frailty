## call the data in r environment
data <- read.csv("data.csv")
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
## cell probability for correlated frailty model
prob_cell_corr <- function(x,mu_1,mu_2,tau_1,tau_2,phi,
                           q_1,q_2,q_3,j_1,j_2)
{
  a_1 <- mu_1*(tau_1^2)*x
  a_2 <- mu_2*(tau_2^2)*x
  b_1 <- 1 + a_1 + a_2
  b_2 <- 1 + a_1
  b_3 <- 1 + a_2
  c_1 <- phi/(tau_1*tau_2)
  c_2 <- c_1 - (1/(tau_1^2))
  c_3 <- c_1 - (1/(tau_2^2))
  marg_surv_1 <- (b_2)^(-(1/(tau_1^2)))
  marg_surv_2 <- (b_3)^(-(1/(tau_2^2)))
  joint_surv <-  (b_1^(-c_1))*(b_2^c_2)*(b_3^c_3)
  if((j_1 == 1)&&(j_2 == 1))
  {
    return(joint_surv)
  }else if((j_1 > 1)&&(j_2 == 1))
  {
    return(q_1[j_1 - 1]*(marg_surv_2 - joint_surv))
  }else if((j_1 == 1)&&(j_2 > 1))
  {
    return(q_2[j_2 - 1]*(marg_surv_1 - joint_surv))
  }else
  {
    return(q_3[j_1 - 1,j_2 - 1]*(1 - marg_surv_1 - marg_surv_2 + joint_surv))
  }
}
## for correlated gamma model
in_est_corr <- function(monitoring_data,Freq_table,M_1,M_2,theta)
{
  mu_vec_1 <- theta[1:M_1]
  mu_vec_2 <- theta[(M_1 + 1):(M_1 + M_2)]
  tau_1 <- theta[M_1+M_2+1]
  tau_2 <- theta[M_1+M_2+2]
  phi <- theta[M_1+M_2+3]
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
      c <- c + (Freq_table[j_1,j_2] -  sum(sapply(monitoring_data[,1],prob_cell_corr,phi=phi,
      mu_1=mu_1,mu_2=mu_2,tau_1=tau_1,tau_2=tau_2,q_1=q_4,q_2=q_5,q_3=q_6,j_1=j_1,j_2=j_2)))^2
    }
  }
  return(c)
}
## log-likelihood for a particular sample point 
## for correlated frailty model
corr_indiv_log_likelihood <- function(data,tau_1,tau_2,mu_1,mu_2,
                                      phi,q_1,q_2,q_3)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  a_1 <- mu_1*(tau_1^2)*x
  a_2 <- mu_2*(tau_2^2)*x
  b_1 <- 1 + a_1 + a_2
  b_2 <- 1 + a_1
  b_3 <- 1 + a_2
  c_1 <- phi/(tau_1*tau_2)
  c_2 <- c_1 - (1/(tau_1^2))
  c_3 <- c_1 - (1/(tau_2^2))
  marg_surv_1 <- (b_2)^(-(1/(tau_1^2)))
  marg_surv_2 <- (b_3)^(-(1/(tau_2^2)))
  if(phi < min((tau_1/tau_2),(tau_2/tau_1)))
  {
    log_joint_surv <- ((-c_1)*log(b_1)) + (c_2*log(b_2)) + (c_3*log(b_3)) 
    joint_surv <- exp(log_joint_surv)
    c <- (1 - (marg_surv_1) - (marg_surv_2) + joint_surv)
    joint_subdist_func <- (q_3)*c
    if((j_1 != 0) && (j_2 != 0))
    {
      return(log(joint_subdist_func[j_1,j_2]))
    }else if((j_1 != 0) && (j_2 == 0)){
      return(log(q_1[j_1]) + log(marg_surv_2 - joint_surv))
    }else if((j_1 == 0) && (j_2 != 0)){
      return(log(q_2[j_2]) + log(marg_surv_1 - joint_surv))
    }else{
      return(log_joint_surv)
    }
  }
}
## negative log-likelihood of correlated frailty model
## The argument about which it is to be maximized should be first
## argument for optimparallel
neg_log_likelihood_corr <- function(theta,D,M_1,M_2)
{
  lambda_vec_3 <- theta[1:M_1]
  lambda_vec_4 <- theta[(M_1 + 1):(M_1 + M_2)]
  tau_3 <- theta[M_1 + M_2 + 1]
  tau_4 <- theta[M_1 + M_2 + 2]
  phi <- theta[M_1 + M_2 + 3]
  lambda_3 <- sum(lambda_vec_3)
  lambda_4 <- sum(lambda_vec_4)
  q_4 <- (lambda_vec_3)/lambda_3
  q_5 <- (lambda_vec_4)/lambda_4
  q_6 <- outer(q_4,q_5)
  return(-sum(apply(D,1,corr_indiv_log_likelihood,
                    mu_1=lambda_3,mu_2=lambda_4,q_1=q_4,
                    q_2=q_5,q_3=q_6,tau_1=tau_3,tau_2=tau_4,phi=phi)))
}
## individual gradient for correlated frailty model
corr_indiv_neg_gr <- function(data,mu_vec_1,mu_vec_2,mu_1,mu_2,
                              phi,tau_1,tau_2,M_1,M_2)
{
  x <- data[1]
  j_1 <- data[2]
  j_2 <- data[3]
  Q_1 <- 1 + (tau_1^2)*mu_1*x
  Q_2 <- 1 + (tau_2^2)*mu_2*x
  Q_3 <- Q_1 + Q_2 - 1 
  Q_4 <- phi/(tau_1*tau_2)
  Q_5 <- Q_4 - 1/(tau_1^2)
  Q_6 <- Q_4 - 1/(tau_2^2)
  if((j_1 == 0)&(j_2 == 0))
  {
    G_1 <- -(Q_4*(tau_1^2)*x)/(Q_3) + ((Q_5*(tau_1^2)*x)/Q_1)
    G_2 <- -(Q_4*(tau_2^2)*x)/(Q_3) + ((Q_6*(tau_2^2)*x)/Q_2)
    G_3 <- -(Q_4/phi)*(log(Q_3) - log(Q_1) - log(Q_2)) 
    G_4 <- (Q_4*(1/tau_1)*log(Q_3)) - (2*tau_1*Q_4*x*mu_1)/(Q_3) +
      (Q_5)*(2*x*tau_1*mu_1)*(1/Q_1) +
      (2/(tau_1^3) - Q_4*(1/tau_1))*log(Q_1) - (Q_4*log(Q_2)*(1/tau_1))
    G_5 <- (Q_4*(1/tau_2)*log(Q_3)) - (2*tau_2*Q_4*x*mu_2)/(Q_3) +
      (Q_6)*(2*x*tau_2*mu_2)*(1/Q_2) +
      (2/(tau_2^3) - Q_4*(1/tau_2))*log(Q_2) - (Q_4*log(Q_1)*(1/tau_2)) 
    v <- c(replicate(M_1,G_1),replicate(M_2,G_2),G_3,G_4,G_5)
    if(typeof(v) != "list")
    {
      return(-M_1*M_2*v)
    }else
    {
      return(-M_1*M_2*unlist(v))
    }
  }else if((j_1 == 0) && (j_2 != 0))
  {
    a <- (Q_1)^(-(1/(tau_1^2)))
    d_1 <- (Q_3)^(-Q_4)
    d_2 <- (Q_1)^(Q_5)
    d_3 <- (Q_2)^(Q_6)
    H_2 <- a - (d_1*d_2*d_3)
    G_1 <- (-((Q_1)^(-(1/(tau_1^2)) - 1)) + 
              (phi*tau_1/tau_2)*(Q_3^(-Q_4 - 1))*(Q_1^Q_5)*(Q_2^Q_6) -
              (Q_3^(-Q_4))*((phi*tau_1/tau_2) - 1)*(Q_1^(Q_5 - 1))*(Q_2^(Q_6)))*x
    G_2 <- ((phi*tau_2/tau_1)*(Q_3^(-Q_4 - 1))*(Q_1^Q_5)*(Q_2^Q_6) -
              (Q_3^(-Q_4))*(Q_1^Q_5)*Q_6*(tau_2^2)*(Q_2^(Q_6 - 1)) - 
              (Q_3^(-Q_4))*(Q_2^Q_6)*Q_5*(tau_1^2)*(Q_1^(Q_5 - 1)))*x
    G_3 <- (Q_4/phi)*((Q_3^(-Q_4))*log(Q_3)*(Q_1^Q_5)*(Q_2^Q_6) - 
                        (Q_3^(-Q_4))*(Q_1^Q_5)*(Q_2^Q_6)*log(Q_1) - 
                        (Q_3^(-Q_4))*(Q_1^Q_5)*(Q_2^Q_6)*log(Q_2))
    G_4 <- (Q_1^(-1/(tau_1^2)))*((2/(tau_1^3))*log(Q_1) - (1/(tau_1^2))*(2*tau_1*x*mu_1)*(1/Q_1)) - 
      (Q_3^(-Q_4))*((1/tau_1)*Q_4*log(Q_3) - (2*Q_4*tau_1*x*mu_1/Q_3))*(Q_1^Q_5)*(Q_2^Q_6) - 
      (Q_3^(-Q_4))*(Q_1^Q_5)*((2/(tau_1^3) - (Q_4*1/tau_1))*log(Q_1) + Q_5*(2*x*mu_1*tau_1/Q_1))*(Q_2^Q_6) +
      (Q_3^(-Q_4))*(Q_1^Q_5)*log(Q_2)*(Q_4/tau_1)*(Q_2^Q_6)
    G_5 <- -(Q_3^(-Q_4))*((Q_4/tau_2)*log(Q_3) - (Q_4*(2*mu_2*tau_2*x/Q_3)))*
      (Q_1^(Q_5))*(Q_2^(Q_6)) + (Q_3^(-Q_4))*(Q_1^(Q_5))*(Q_4/tau_2)*log(Q_1)*(Q_2^(Q_6)) -
      (Q_3^(-Q_4))*(Q_1^(Q_5))*(Q_2^(Q_6))*((2/(tau_2^3) - Q_4/tau_2)*log(Q_2) + Q_6*(2*mu_2*x*tau_2/Q_2))
    G_6 <- (M_1/mu_2)
    if(j_2 == 1)
    {
      v <- c(replicate(M_1,G_1),c(G_2 - (M_1/mu_vec_2[j_2]) + G_6,
                                  replicate(M_2 - 1,G_2 + G_6)),G_3,G_4,G_5)
    }else if(j_2 == M_2)
    {
      v <- c(replicate(M_1,G_1),c(replicate(M_2 - 1,G_2 + G_6),
                                  G_2 - (M_1/mu_vec_2[j_2]) + G_6),G_3,G_4,G_5)
    }else
    {
      v <- c(replicate(M_1,G_1),c(replicate(j_2 - 1,G_2 + G_6),
                                  G_2 - (M_1/mu_vec_2[j_2]) + G_6,replicate(M_2 - j_2,G_2 + G_6)),G_3,G_4,G_5)
    }
    if(typeof(v) != "list")
    {
      return(-(M_1/H_2)*v)
    }else
    {
      return(-(M_1/H_2)*unlist(v))
    } 
  }else if((j_1 != 0) && (j_2 == 0))
  {
    b <- (Q_2)^(-(1/(tau_2^2)))
    d_1 <- (Q_3)^(-Q_4)
    d_2 <- (Q_1)^(Q_5)
    d_3 <- (Q_2)^(Q_6)
    H_3 <- b - (d_1*d_2*d_3)
    G_1 <- (Q_4*(Q_3^(-Q_4-1))*(tau_1^2)*(Q_1^(Q_5))*(Q_2^(Q_6)) - 
              (Q_3^(-Q_4))*Q_5*(Q_1^(Q_5 - 1))*(tau_1^2)*(Q_2^(Q_6)))*x 
    G_2 <- (-(Q_2^(-1 - (1/(tau_2^2)))) + 
              Q_4*(Q_3^(-Q_4-1))*(tau_2^2)*(Q_1^(Q_5))*(Q_2^(Q_6)) - 
              (Q_3^(-Q_4))*(Q_1^(Q_5))*Q_6*(Q_2^(Q_6 - 1))*(tau_2^2))*x
    G_3 <- ((Q_3^(-Q_4))*log(Q_3)*(Q_1^(Q_5))*(Q_2^(Q_6)) - 
              (Q_3^(-Q_4))*(Q_1^(Q_5))*log(Q_1)*(Q_2^(Q_6)) - 
              (Q_3^(-Q_4))*(Q_1^(Q_5))*(Q_2^(Q_6))*log(Q_2))*(Q_4/phi)
    G_4 <- -(Q_3^(-Q_4))*(log(Q_3)*(Q_4/tau_1) - (Q_4*2*tau_1*mu_1*x/Q_3))*(Q_1^(Q_5))*(Q_2^(Q_6)) - 
      (Q_3^(-Q_4))*(Q_1^(Q_5))*((Q_5*2*x*mu_1*tau_1/Q_1) + log(Q_1)*((2/(tau_1^3)) - (Q_4/tau_1)))*(Q_2^(Q_6)) + 
      (Q_3^(-Q_4))*(Q_1^Q_5)*log(Q_2)*(Q_4/tau_1)*(Q_2^Q_6)
    G_5 <- (Q_2^(-1/(tau_2^2)))*((2/(tau_2^3))*log(Q_2) - (1/tau_2)*(2*x*mu_2/Q_2))
    -(Q_3^(-Q_4))*((Q_4/tau_2)*log(Q_3) - (Q_4*(2*mu_2*tau_2*x/Q_3)))*
      (Q_1^(Q_5))*(Q_2^(Q_6)) + (Q_3^(-Q_4))*(Q_1^(Q_5))*(Q_4/tau_2)*log(Q_1)*(Q_2^(Q_6)) -
      (Q_3^(-Q_4))*(Q_1^(Q_5))*(Q_2^(Q_6))*((2/(tau_2^3) - Q_4/tau_2)*log(Q_2) + Q_6*(2*mu_2*x*tau_2/Q_2))
    G_6 <- M_2/mu_1
    if(j_1 == 1)
    {
      v <- c(c(G_1 - (M_2/mu_vec_1[j_1]) + G_6,
               replicate(M_1 - 1,G_1 + G_6)),replicate(M_2,G_2),G_3,G_4,G_5)
    }else if(j_1 == M_1)
    {
      v <- c(c(replicate(M_1 - 1,G_1 + G_6),G_1 - (M_2/mu_vec_1[j_1]) + G_6),
             replicate(M_2,G_2),G_3,G_4,G_5)
    }else
    {
      v <- c(c(replicate(j_1 - 1,G_1 + G_6),G_1 - (M_2/mu_vec_1[j_1]) + G_6,
               replicate(M_1 - j_1,G_1 + G_6)),replicate(M_2,G_2),G_3,G_4,G_5)
    }
    if(typeof(v) != "list")
    {
      return(-(M_2/H_3)*v)
    }else
    {
      return(-(M_2/H_3)*unlist(v))
    } 
  }else
  {
    a <- (Q_1)^(-(1/(tau_1^2)))
    b <- (Q_2)^(-(1/(tau_2^2)))
    d_1 <- (Q_3)^(-Q_4)
    d_2 <- (Q_1)^(Q_5)
    d_3 <- (Q_2)^(Q_6)
    H_4 <- 1 - a - b + (d_1*d_2*d_3) 
    G_1 <- ((Q_1)^(-(1/(tau_1^2)) - 1) - 
              (phi*tau_1/tau_2)*(Q_3^(-Q_4 - 1))*(Q_1^Q_5)*(Q_2^Q_6) +
              (Q_3^(-Q_4))*((phi*tau_1/tau_2) - 1)*(Q_1^(Q_5 - 1))*(Q_2^(Q_6)))*x
    G_2 <- ((Q_2^(-1 - (1/(tau_2^2)))) - 
              Q_4*(Q_3^(-Q_4-1))*(tau_2^2)*(Q_1^(Q_5))*(Q_2^(Q_6)) + 
              (Q_3^(-Q_4))*(Q_1^(Q_5))*Q_6*(Q_2^(Q_6 - 1))*(tau_2^2))*x
    G_3 <- (-(Q_3^(-Q_4))*log(Q_3)*(Q_1^(Q_5))*(Q_2^(Q_6)) + 
              (Q_3^(-Q_4))*(Q_1^(Q_5))*log(Q_1)*(Q_2^(Q_6)) + 
              (Q_3^(-Q_4))*(Q_1^(Q_5))*(Q_2^(Q_6))*log(Q_2))*(Q_4/phi)
    G_4 <- -(Q_1^(-1/(tau_1^2)))*((2/(tau_1^3))*log(Q_1) - (1/(tau_1^2))*(2*tau_1*x*mu_1)*(1/Q_1)) + 
      (Q_3^(-Q_4))*((1/tau_1)*Q_4*log(Q_3) - (2*Q_4*tau_1*x*mu_1/Q_3))*(Q_1^Q_5)*(Q_2^Q_6) + 
      (Q_3^(-Q_4))*(Q_1^Q_5)*((2/(tau_1^3) - (Q_4*1/tau_1))*log(Q_1) + Q_5*(2*x*mu_1*tau_1/Q_1))*(Q_2^Q_6) -
      (Q_3^(-Q_4))*(Q_1^Q_5)*log(Q_1)*(Q_4/tau_1)*(Q_2^Q_6)
    G_5 <- -(Q_2^(-1/(tau_2^2)))*((2/(tau_2^3))*log(Q_2) - (1/tau_2)*(2*x*mu_2/Q_2)) + 
      (Q_3^(-Q_4))*((Q_4/tau_2)*log(Q_3) - (Q_4*(2*mu_2*tau_2*x/Q_3)))*
      (Q_1^(Q_5))*(Q_2^(Q_6)) - (Q_3^(-Q_4))*(Q_1^(Q_5))*(Q_4/tau_2)*log(Q_1)*(Q_2^(Q_6)) +
      (Q_3^(-Q_4))*(Q_1^(Q_5))*(Q_2^(Q_6))*((2/(tau_2^3) - Q_4/tau_2)*log(Q_2) + Q_6*(2*mu_2*x*tau_2/Q_2))
    G_6 <- 1/mu_1
    G_7 <- 1/mu_2
    if((j_1 == 1)&(j_2 == 1))
    {
      v <- c(c(G_1 - (1/mu_vec_1[j_1]) + G_6,replicate(M_1 - 1,G_1 + G_6)),
             c(G_2 - (1/mu_vec_2[j_2]) + G_7,replicate(M_2 - 1,G_2 + G_7)),G_3,G_4,G_5)
    }else if((j_1 == 1)&(j_2 == M_2))
    {
      v <- c(c(G_1 - (1/mu_vec_1[j_1]) + G_6,replicate(M_1 - 1,G_1 + G_6)),
             c(replicate(M_2 - 1,G_2 + G_7),G_2 - (1/mu_vec_2[j_2]) + G_7),G_3,G_4,G_5)
    }else if((j_1 == M_1)&(j_2 == 1))
    {
      v <- c(c(replicate(M_1 - 1,G_1 + G_6),G_1 - (1/mu_vec_1[j_1]) + G_6),
             c(G_2 - (1/mu_vec_2[j_2]) + G_7,replicate(M_2 - 1,G_2 + G_7)),G_3,G_4,G_5)
    }else if((j_1 == M_1)&(j_2 == M_2))
    {
      v <- c(c(replicate(M_1 - 1,G_1 + G_6),G_1 - (1/mu_vec_1[j_1]) + G_6),
             c(replicate(M_2 - 1,G_2 + G_7),G_2 - (1/mu_vec_2[j_2]) + G_7),G_3,G_4,G_5)
    }else
    {
      v <- c(c(replicate(j_1 - 1,G_1 + G_6),G_1 - (1/mu_vec_1[j_1]) + G_6,replicate(M_1 - j_1,G_1 + G_6)),
             c(replicate(j_2 - 1,G_2 + G_7),G_2 - (1/mu_vec_2[j_2]) + G_7,replicate(M_2 - j_2,G_2 + G_7)),G_3,G_4,G_5)
    }
    if(typeof(v) != "list")
    {
      return(-(1/H_4)*v)
    }else
    {
      return(-(1/H_4)*unlist(v))
    } 
  }
}
## negative gradient for correlated frailty model
corr_neg_gradient <- function(D,theta,M_1,M_2)
{
  mu_vec_3 <- theta[1:M_1]
  mu_vec_4 <- theta[(M_1 + 1):(M_1 + M_2)]
  tau_1 <- theta[M_1 + M_2 + 1]
  tau_2 <- theta[M_1 + M_2 + 2]
  phi <- theta[M_1 + M_2 + 3]
  mu_3 <- sum(mu_vec_3)
  mu_4 <- sum(mu_vec_4)
  m <- apply(D,1,corr_indiv_neg_gr,tau_1=tau_1,tau_2=tau_2,
             mu_vec_1=mu_vec_3,mu_vec_2=mu_vec_4,phi=phi,
             mu_1=mu_3,mu_2=mu_4,M_1=M_1,M_2=M_2)
  if(typeof(m) != "list")
  {
    return(apply(m,1,sum))
  }else
  {
    m_1 <- matrix(unlist(m),ncol=M_1+M_2+1,byrow=TRUE)
    return(apply(m_1,2,sum))
  }
}
## numerical gradient for correlated frailty
corr_num_gradient <- function(D,theta,M_1,M_2)
{
  return(numDeriv::grad(neg_log_likelihood_corr,D=D,M_1=M_1,M_2=M_2))
}
## gradient of prob_cell_corr function
gr_prob_cell_corr <- function(x,j_1,j_2,mu_vec_1,mu_vec_2,tau_1,
                              tau_2,phi,mu_1,mu_2,M_1,M_2)
{
  q_1 <- mu_vec_1/mu_1
  q_2 <- mu_vec_2/mu_2
  q_3 <- outer(q_1,q_2)
  p <- prob_cell_corr(x,mu_1,mu_2,tau_1,tau_2,phi,q_1,q_2,q_3,j_1,j_2)
  g <- corr_indiv_neg_gr(c(x,j_1 - 1,j_2 - 1),mu_vec_1,mu_vec_2,mu_1,mu_2,phi,tau_1,tau_2,M_1,M_2)
  if(typeof(g) != "list")
  {
    return(p*g)
  }else
  {
    return(p*(matrix(unlist(g),nrow=1,ncol = length(mu_vec_1) + 
                       length(mu_vec_2) + 3)))
  }
}
## gradient of function in_est_corr for correlated frailty model
gr_in_est_corr <- function(monitoring_data,Freq_table,M_1,M_2,theta)
{
  mu_vec_1 <- theta[1:M_1]
  mu_vec_2 <- theta[(M_1 + 1):(M_1 + M_2)]
  tau_1 <- theta[M_1+M_2+1]
  tau_2 <- theta[M_1+M_2+2]
  phi <- theta[M_1+M_2+3]
  mu_1 <- sum(mu_vec_1)
  mu_2 <- sum(mu_vec_2)
  q_4 <- mu_vec_1/mu_1
  q_5 <- mu_vec_2/mu_2
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
      prob_cell_corr,mu_1=mu_1,mu_2=mu_2,tau_1=tau_1,tau_2=tau_2,
      phi=phi,q_1=q_4,q_2=q_5,q_3=q_6,j_1=j_1,j_2=j_2)))
      m <- sapply(monitoring_data[,1],gr_prob_cell_corr,mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,
                  mu_1=mu_1,mu_2=mu_2,tau_1=tau_1,tau_2=tau_2,phi=phi,M_1=M_1,M_2=M_2,j_1=j_1,j_2=j_2)
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
      m <- sapply(monitoring_data[,1],gr_prob_cell_corr,mu_vec_1=mu_vec_1,mu_vec_2=mu_vec_2,
                  mu_1=mu_1,mu_2=mu_2,tau_1=tau_1,tau_2=tau_2,phi=phi,M_1=M_1,M_2=M_2,j_1=j_1,j_2=j_2)
      m_2 <- matrix(m[(M_1 + 1):(M_1 + M_2),],nrow=M_2,ncol=length(monitoring_data[,1]))
      c_2 <- c_2 + d[j_1,j_2]*apply(m_2,1,sum)
    }
  }
  return(c(c_1,c_2,c_3,c_4,c_5))
}
gr_in_est_corr(hearing_data,F,3,3,replicate(9,44531))
#### correlated frailty model fitting ####
library("doParallel")
data_soln <- c(0.65008783,0.04032277,0.08070887,0.62244563,0.05502849,0.08690448)
## sum function to get initial estimate for final parameter estimation
## in correlated frailty model
S_corr <- function(tau_1,tau_2,phi)
{
  return(in_est_corr(hearing_data,F,3,3,
                     c(data_soln[1:6],tau_1,tau_2,phi)))
}
cl <- makePSOCKcluster(27)
registerDoParallel(cl)
tau <- seq(2,3,by = 0.1)
r <- c(seq(-0.9,-0.1,by=0.1),seq(0.1,0.9,by=0.1))
mat <- matrix(0,nrow=(11^2),ncol=4)
n <- 0
for(i in 1:11)
{
  for(j in 1:11)
  {
    n <- n + 1
    mat[n,] <- c(tau[i],tau[j],0.99,S_corr(tau[i],tau[j],0.99))
  }
}
min_pos <- which.min(mat[,4]) # 37th row
mat[min_pos,] #tau_1=2.3,tau_2=2.3,r=0.99,value = 17549.37(which is min)
stopCluster(cl)
pho_vec <- seq(0.89,0.99,by=0.001)
f_value <- sapply(pho_vec,S_corr,tau_1=2.3,tau_2=2.3)
f_value
plot(pho_vec,f_value,type="l")
#initial_est_corr <- BB::spg(c(1,1.2,0.6),
#                    in_est_corr,monitoring_data = hearing_data[1:20,],Freq_table=F,M_1 = 3,M_2 = 3,method = 2,
#                    lower = c(replicate(2,0.01),0.1),upper = c(replicate(2,10^2),0.99),control = list(checkGrad=FALSE))
#initial_est_corr$par
#initial_est_corr$convergence
#in_est_corr(hearing_data,F,3,3,initial_est_corr$par)
#in_est_corr(hearing_data,F,3,3,c(lambda_vec_1_initial,lambda_vec_2_initial,1,1.2,0.6))
#in_est_corr_data <- initial_est_corr$par
data_soln <- c(0.65008783,0.04032277,0.08070887,0.62244563,0.05502849,0.08690448)
in_est_corr_data_1 <- c(data_soln,2.3,2.3,0.3)
library("doParallel")
cl <- makeCluster(34)
setDefaultCluster(cl=cl)
clusterExport(cl,"corr_indiv_log_likelihood")
data_soln_corr <- optimParallel::optimParallel(par=in_est_corr_data_1,fn=neg_log_likelihood_corr,
D = hearing_data,M_1 = 3,M_2 = 3,lower = replicate(9,0.01),
upper = c(replicate(8,Inf),1),method = "L-BFGS-B",control=list(trace=6),
parallel=list(cl=NULL,forward=FALSE,loginfo=TRUE))$par#control=list(checkGrad=FALSE))
hess_corr <- numDeriv::hessian(neg_log_likelihood_corr,data_soln_corr,D=hearing_data,M_1=3,M_2=3)
stopCluster(cl)
data_soln_corr
neg_log_likelihood_corr(in_est_corr_data_1,hearing_data,3,3)
neg_log_likelihood_corr(as.vector(data_soln_corr),hearing_data,3,3)
con=function(theta)
{
  f=NULL
  tau_1=theta[7]
  tau_2=theta[8]
  rho=theta[9]
  f=rho-min((tau_1/tau_2),(tau_2/tau_1))
  return(list(ceq=NULL,c=f))
}
objfun=function(theta)
{
  return(neg_log_likelihood_corr(hearing_data,theta,3,3))
}
in_est_corr_data_1 <- c(data_soln,2.3,2.3,0.9)
data_soln_corr_1 <- NlcOptim::solnl(in_est_corr_data_1,objfun,con,tolFun = 1e-4,
                    lb = c(replicate(8,0),0),ub = c(replicate(8,Inf),1),tolCon = 1e-4)
data_soln_corr_1
eigen(data_soln_corr_1$hessian)
neg_log_likelihood_corr(in_est_corr_data_1,hearing_data,3,3)
#H <- neg_hessian(hearing_data,data_soln,3,3) 
#eigen(H)$values
eigen(hess_corr)$values
library(alabama)
objfun=function(theta)
{
  return(neg_log_likelihood_corr(theta,hearing_data,3,3))
}
hin <- function(theta)
{
  return(c(theta[1],theta[2],theta[3],theta[4],
           theta[5],theta[6],theta[7],theta[8],
           theta[9],
           min((theta[7]/theta[8]),(theta[8]/theta[7]))-theta[9]))
}
hin_1 <- function(theta)
{
  h <- numeric(10)
  h[1] <- -theta[1]
  h[2] <- -theta[2]
  h[3] <- -theta[3]
  h[4] <- -theta[4]
  h[5] <- -theta[5]
  h[6] <- -theta[6]
  h[7] <- -theta[7]
  h[8] <- -theta[8]
  h[9] <- -theta[9]
  h[10] <- -min((theta[7]/theta[8]),(theta[8]/theta[7])) + theta[9]
  return(h)
}
gr <- function(theta)
{
  return(corr_neg_gradient(hearing_data,theta,3,3))
}
num_gr <- function(theta)
{
  return(numDeriv::grad(neg_log_likelihood_corr,theta,
         D=hearing_data,M_1=3,M_2=3))
}
jac_hin_1 <- function(theta)
{
  rbind(diag(replicate(9,1),nrow=9,ncol=9),
        numDeriv::grad(hin_1[10]))
}
in_est_corr_data_1 <- c(data_soln,2.4,2.4,0.99)
data_soln_corr_1 <- constrOptim.nl(in_est_corr_data_1,objfun,num_gr,hin)$par
data_soln_corr_2 <- nloptr::nloptr(x0=in_est_corr_data_1,eval_f=objfun,eval_grad_f = num_gr,
                    eval_g_ineq = hin_1,opts = list("algorithm"="NLOPT_LD_MMA"))$par
neg_log_likelihood_corr(data_soln_corr_1,hearing_data,3,3)
gr(data_soln_corr_1)
hess <- pracma::hessian(neg_log_likelihood_corr,data_soln_corr_1,
                          D=hearing_data,M_1=3,M_2=3)
eigen(hess)$values
V <- as.matrix(bdsmatrix::gchol(hess))
sd_mle_corr <- diag(solve(t(V)%*%V))
