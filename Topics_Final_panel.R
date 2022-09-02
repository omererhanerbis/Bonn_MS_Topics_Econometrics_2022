set.seed(2016300192)
library(MASS)
library(caret)
library(DRDID)

# Monte_Carlo_DRDID
n_monte_carlo <- 10
n_sample <- 1000

DGP1_vec_dr_imp <- c()
DGP1_vec_ipw <- c()
DGP1_vec_or <- c()
DGP1_vec_dr <- c()
DGP1_vec_twfe <- c()
DGP2_vec_dr_imp <- c()
DGP2_vec_ipw <- c()
DGP2_vec_or <- c()
DGP2_vec_dr <- c()
DGP2_vec_twfe <- c()
DGP3_vec_dr_imp <- c()
DGP3_vec_ipw <- c()
DGP3_vec_or <- c()
DGP3_vec_dr <- c()
DGP3_vec_twfe <- c()
DGP4_vec_dr_imp <- c()
DGP4_vec_ipw <- c()
DGP4_vec_or <- c()
DGP4_vec_dr <- c()
DGP4_vec_twfe <- c()

DGP1_vec_dr_imp_cover <- c()
DGP1_vec_ipw_cover    <- c()
DGP1_vec_or_cover     <- c()
DGP1_vec_dr_cover     <- c()
DGP1_vec_twfe_cover   <- c()
DGP2_vec_dr_imp_cover <- c()
DGP2_vec_ipw_cover    <- c()
DGP2_vec_or_cover     <- c()
DGP2_vec_dr_cover     <- c()
DGP2_vec_twfe_cover   <- c()
DGP3_vec_dr_imp_cover <- c()
DGP3_vec_ipw_cover    <- c()
DGP3_vec_or_cover     <- c()
DGP3_vec_dr_cover     <- c()
DGP3_vec_twfe_cover   <- c()
DGP4_vec_dr_imp_cover <- c()
DGP4_vec_ipw_cover    <- c()
DGP4_vec_or_cover     <- c()
DGP4_vec_dr_cover     <- c()
DGP4_vec_twfe_cover   <- c()


DGP1_vec_dr_imp_inf <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP1_vec_ipw_inf    <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP1_vec_or_inf     <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP1_vec_dr_inf     <- matrix(ncol = n_monte_carlo, nrow = n_sample)
# DGP1_twfe_inf       <- matrix(ncol = n_sample, nrow = 1000)
DGP2_vec_dr_imp_inf <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP2_vec_ipw_inf    <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP2_vec_or_inf     <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP2_vec_dr_inf     <- matrix(ncol = n_monte_carlo, nrow = n_sample)
# DGP2_twfe_inf       <- matrix(ncol = n_sample, nrow = 1000)
DGP3_vec_dr_imp_inf <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP3_vec_ipw_inf    <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP3_vec_or_inf     <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP3_vec_dr_inf     <- matrix(ncol = n_monte_carlo, nrow = n_sample)
# DGP3_twfe_inf       <- matrix(ncol = n_sample, nrow = 1000)
DGP4_vec_dr_imp_inf <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP4_vec_ipw_inf    <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP4_vec_or_inf     <- matrix(ncol = n_monte_carlo, nrow = n_sample)
DGP4_vec_dr_inf     <- matrix(ncol = n_monte_carlo, nrow = n_sample)
# DGP4_twfe_inf       <- matrix(ncol = n_sample, nrow = 1000)


# Functions
f_reg <- function(vector_W)
{
  return(210 + 27.4*vector_W[1] + 13.7*(sum(vector_W[2:4])))
}


f_ps <- function(vector_W)
{
  return(0.75*(-vector_W[1] + vector_W[2]/2 - vector_W[3]/4 - vector_W[4]/10))
}

nu <- function(vector_W, treatment)
{
  holder <- c()
  for(i in 1:nrow(vector_W))
  {
    selection <- rnorm(1, mean = as.numeric(treatment[i] * f_reg(vector_W[i,])), sd = 1)
    holder    <- c(holder, selection)
  }
  return(holder)
}
# # treatment choice
p <- function(vector_W)
{
  return(exp(f_ps(vector_W)) / (1 + exp(f_ps(vector_W))))
}





start_time <- Sys.time()


for(m in 1:n_monte_carlo)
{
  # Correct(X) and Observed(Z) Covariates creation
  X <- mvrnorm(n = n_sample, mu = rep(0, 4), Sigma = matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1), nrow = 4))
  X <- as.data.frame(X)
  
  Z_tilde <- data.frame(Z_tilde_1 = as.numeric(), Z_tilde_2 = as.numeric(), Z_tilde_3 = as.numeric(), Z_tilde_4 = as.numeric())
  for(i in 1:nrow(X))
  {
    Z_tilde[i,] <- c(exp(X[i,1]/2), 10+X[i,2]/(1+exp(X[i,1])), (0.6+X[i,1]*X[i,3]/25)^3, (20+X[i,2]+X[i,4])^2)
  }
  Z <- as.data.frame(scale(Z_tilde, center = T, scale = T))
  
  
  # threshold for treatment choice
  U <- runif(n_sample, 0, 1)
  
  
  # noise epsilon_indexD
  epsilon_0  <- rnorm(n_sample)
  epsilon_10 <- rnorm(n_sample/2)
  epsilon_11 <- rnorm(n_sample/2)
  
  # FASTENING THE CODE
  p_Z <- p(Z)
  p_X <- p(X)
  f_reg_Z <- f_reg(Z)
  f_reg_X <- f_reg(X)
  nu_Z_0_1000 <- nu(Z, rep(0, 1000))
  nu_Z_1_1000 <- nu(Z, rep(1, 1000))
  nu_X_0_1000 <- nu(X, rep(0, 1000))
  nu_X_1_1000 <- nu(X, rep(1, 1000))
  
  #################### DGP1 #################### 
  D <- as.numeric(p_Z >= U)
  
  
  # Y_tD
  Y_00 <- f_reg_Z + nu(Z, D) + epsilon_0
  Y_10 <- 2 * f_reg_Z + nu_Z_0_1000 + epsilon_10
  Y_11 <- 2 * f_reg_Z + nu_Z_1_1000 + epsilon_11
  
  # complete data-frame creations
  time_0_correct <- data.frame(re = Y_00, t = 0, id = 1:1000, d = D,
                               Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  
  time_1_treat_0_all <- data.frame(re = Y_10, t = 1, id = 1:1000, d = D,
                                   Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  time_1_treat_1_all <- data.frame(re = Y_11, t = 1, id = 1:1000, d = D,
                                   Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  
  time_1_correct <- rbind(time_1_treat_0_all[time_1_treat_0_all$d==0,],
                          (time_1_treat_1_all[time_1_treat_1_all$d==1,]))
  
  
  # finalizing the correct complete data-frame
  df <- rbind(time_0_correct, time_1_correct)
  colnames(df)[1] <- "re"
  
  
  dr_imp_1 <- drdid(yname="re", tname = "t", idname = "id", dname = "d",
                  xformla= ~ Z1+Z2+Z3+Z4,
                  data = df, panel = TRUE, inffunc = T)
  ipw_1 <- ipwdid(yname="re", tname = "t", idname = "id", dname = "d",
                xformla= ~ Z1+Z2+Z3+Z4,
                data = df, panel = TRUE, inffunc = T)
  or_1 <- ordid(yname="re", tname = "t", idname = "id", dname = "d",
              xformla= ~ Z1+Z2+Z3+Z4,
              data = df, panel = TRUE, inffunc = T)
  dr_1 <- drdid(yname="re", tname = "t", idname = "id", dname = "d",
                    xformla= ~ Z1+Z2+Z3+Z4, estMethod = "trad",
                    data = df, panel = TRUE, inffunc = T)
  twfe_1 <- twfe_did_panel(y1=time_1_correct$Z_tilde_1, y0=time_0_correct$Z_tilde_1,
                           D=D, covariates = as.matrix(cbind(df$Z1[1:n_sample], df$Z2[1:n_sample], df$Z3[1:n_sample], df$Z4[1:n_sample])), inffunc = F)
  
  DGP1_vec_dr_imp <- c(DGP1_vec_dr_imp, dr_imp_1$ATT)
  DGP1_vec_ipw <-    c(DGP1_vec_ipw, ipw_1$ATT)
  DGP1_vec_or <-     c(DGP1_vec_or, or_1$ATT)
  DGP1_vec_dr <-     c(DGP1_vec_dr, dr_1$ATT)
  DGP1_vec_twfe <-   c(DGP1_vec_twfe, twfe_1$ATT)
  DGP1_vec_dr_imp_cover <- c(DGP1_vec_dr_imp_cover, as.numeric(dr_imp_1$lci < 0) * as.numeric(dr_imp_1$uci > 0))
  DGP1_vec_ipw_cover    <- c(DGP1_vec_ipw_cover   , as.numeric(ipw_1$lci < 0) * as.numeric(ipw_1$uci > 0))
  DGP1_vec_or_cover     <- c(DGP1_vec_or_cover    , as.numeric(or_1$lci < 0) * as.numeric(or_1$uci > 0))
  DGP1_vec_dr_cover     <- c(DGP1_vec_dr_cover    , as.numeric(dr_1$lci < 0) * as.numeric(dr_1$uci > 0))
  DGP1_vec_twfe_cover   <- c(DGP1_vec_twfe_cover  , as.numeric(twfe_1$lci < 0) * as.numeric(twfe_1$uci > 0))
  
  DGP1_vec_dr_imp_inf[,m] <- dr_imp_1$att.inf.func
  DGP1_vec_ipw_inf[,m]       <- ipw_1$att.inf.func
  DGP1_vec_or_inf[,m]         <- or_1$att.inf.func
  DGP1_vec_dr_inf[,m]         <- dr_1$att.inf.func
  
  
  
  
  #################### DGP2 #################### 
  D <- as.numeric(p_X >= U)
  
  
  # Y_tD
  Y_00 <- f_reg_Z + nu(Z, D) + epsilon_0
  Y_10 <- 2 * f_reg_Z + nu_Z_0_1000 + epsilon_10
  Y_11 <- 2 * f_reg_Z + nu_Z_1_1000 + epsilon_11
  
  # complete data-frame creations
  time_0_correct <- data.frame(re = Y_00, t = 0, id = 1:1000, d = D,
                               Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  
  time_1_treat_0_all <- data.frame(re = Y_10, t = 1, id = 1:1000, d = D,
                                   Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  time_1_treat_1_all <- data.frame(re = Y_11, t = 1, id = 1:1000, d = D,
                                   Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  
  time_1_correct <- rbind(time_1_treat_0_all[time_1_treat_0_all$d==0,],
                          (time_1_treat_1_all[time_1_treat_1_all$d==1,]))
  
  
  # finalizing the correct complete data-frame
  df <- rbind(time_0_correct, time_1_correct)
  colnames(df)[1] <- "re"
  
  
  dr_imp_2 <- drdid(yname="re", tname = "t", idname = "id", dname = "d",
                  xformla= ~ Z1+Z2+Z3+Z4,
                  data = df, panel = TRUE, inffunc = T)
  ipw_2 <- ipwdid(yname="re", tname = "t", idname = "id", dname = "d",
                xformla= ~ Z1+Z2+Z3+Z4,
                data = df, panel = TRUE, inffunc = T)
  or_2 <- ordid(yname="re", tname = "t", idname = "id", dname = "d",
              xformla= ~ Z1+Z2+Z3+Z4,
              data = df, panel = TRUE, inffunc = T)
  dr_2 <- drdid(yname="re", tname = "t", idname = "id", dname = "d",
                xformla= ~ Z1+Z2+Z3+Z4, estMethod = "trad",
                data = df, panel = TRUE, inffunc = T)
  twfe_2 <- twfe_did_panel(y1=time_1_correct$Z_tilde_1, y0=time_0_correct$Z_tilde_1,
                           D=D, covariates = as.matrix(cbind(df$Z1[1:n_sample], df$Z2[1:n_sample], df$Z3[1:n_sample], df$Z4[1:n_sample])), inffunc = F)
  
  
  DGP2_vec_dr_imp <- c(DGP2_vec_dr_imp, dr_imp_2$ATT)
  DGP2_vec_ipw <-    c(DGP2_vec_ipw, ipw_2$ATT)
  DGP2_vec_or <-     c(DGP2_vec_or, or_2$ATT)
  DGP2_vec_dr <-     c(DGP2_vec_dr, dr_2$ATT)
  DGP2_vec_twfe <-   c(DGP2_vec_twfe, twfe_2$ATT)
  DGP2_vec_dr_imp_cover <- c(DGP2_vec_dr_imp_cover, as.numeric(dr_imp_2$lci < 0) * as.numeric(dr_imp_2$uci > 0))
  DGP2_vec_ipw_cover    <- c(DGP2_vec_ipw_cover   , as.numeric(ipw_2$lci < 0) * as.numeric(ipw_2$uci > 0))
  DGP2_vec_or_cover     <- c(DGP2_vec_or_cover    , as.numeric(or_2$lci < 0) * as.numeric(or_2$uci > 0))
  DGP2_vec_dr_cover     <- c(DGP2_vec_dr_cover    , as.numeric(dr_2$lci < 0) * as.numeric(dr_2$uci > 0))
  DGP2_vec_twfe_cover   <- c(DGP2_vec_twfe_cover  , as.numeric(twfe_2$lci < 0) * as.numeric(twfe_2$uci > 0))
  
  
  DGP2_vec_dr_imp_inf[,m] <- dr_imp_2$att.inf.func
  DGP2_vec_ipw_inf[,m]       <- ipw_2$att.inf.func
  DGP2_vec_or_inf[,m]         <- or_2$att.inf.func
  DGP2_vec_dr_inf[,m]         <- dr_2$att.inf.func
  
  
  #################### DGP3 #################### 
  D <- as.numeric(p_Z >= U)
  
  
  # Y_tD
  Y_00 <- f_reg_X + nu(X, D) + epsilon_0
  Y_10 <- 2 * f_reg_X + nu_X_0_1000 + epsilon_10
  Y_11 <- 2 * f_reg_X + nu_X_1_1000 + epsilon_11
  
  # complete data-frame creations
  time_0_correct <- data.frame(re = Y_00, t = 0, id = 1:1000, d = D,
                               Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  
  time_1_treat_0_all <- data.frame(re = Y_10, t = 1, id = 1:1000, d = D,
                                   Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  time_1_treat_1_all <- data.frame(re = Y_11, t = 1, id = 1:1000, d = D,
                                   Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  
  time_1_correct <- rbind(time_1_treat_0_all[time_1_treat_0_all$d==0,],
                          (time_1_treat_1_all[time_1_treat_1_all$d==1,]))
  
  
  # finalizing the correct complete data-frame
  df <- rbind(time_0_correct, time_1_correct)
  colnames(df)[1] <- "re"
  
  
  dr_imp_3 <- drdid(yname="re", tname = "t", idname = "id", dname = "d",
                  xformla= ~ Z1+Z2+Z3+Z4,
                  data = df, panel = TRUE, inffunc = T)
  ipw_3 <- ipwdid(yname="re", tname = "t", idname = "id", dname = "d",
                xformla= ~ Z1+Z2+Z3+Z4,
                data = df, panel = TRUE, inffunc = T)
  or_3 <- ordid(yname="re", tname = "t", idname = "id", dname = "d",
              xformla= ~ Z1+Z2+Z3+Z4,
              data = df, panel = TRUE, inffunc = T)
  dr_3 <- drdid(yname="re", tname = "t", idname = "id", dname = "d",
                xformla= ~ Z1+Z2+Z3+Z4, estMethod = "trad",
                data = df, panel = TRUE, inffunc = T)
  twfe_3 <- twfe_did_panel(y1=time_1_correct$V1, y0=time_0_correct$V1,
                           D=D, covariates = as.matrix(cbind(df$Z1[1:n_sample], df$Z2[1:n_sample], df$Z3[1:n_sample], df$Z4[1:n_sample])), inffunc = F)
  
  
  DGP3_vec_dr_imp <- c(DGP3_vec_dr_imp, dr_imp_3$ATT)
  DGP3_vec_ipw <-    c(DGP3_vec_ipw, ipw_3$ATT)
  DGP3_vec_or <-     c(DGP3_vec_or, or_3$ATT)
  DGP3_vec_dr <-     c(DGP3_vec_dr, dr_3$ATT)
  DGP3_vec_twfe <-   c(DGP3_vec_twfe, twfe_3$ATT)
  DGP3_vec_dr_imp_cover <- c(DGP3_vec_dr_imp_cover, as.numeric(dr_imp_3$lci < 0) * as.numeric(dr_imp_3$uci > 0))
  DGP3_vec_ipw_cover    <- c(DGP3_vec_ipw_cover   , as.numeric(ipw_3$lci < 0) * as.numeric(ipw_3$uci > 0))
  DGP3_vec_or_cover     <- c(DGP3_vec_or_cover    , as.numeric(or_3$lci < 0) * as.numeric(or_3$uci > 0))
  DGP3_vec_dr_cover     <- c(DGP3_vec_dr_cover    , as.numeric(dr_3$lci < 0) * as.numeric(dr_3$uci > 0))
  DGP3_vec_twfe_cover   <- c(DGP3_vec_twfe_cover  , as.numeric(twfe_3$lci < 0) * as.numeric(twfe_3$uci > 0))
  
  
  DGP3_vec_dr_imp_inf[,m] <- dr_imp_3$att.inf.func
  DGP3_vec_ipw_inf[,m]       <- ipw_3$att.inf.func
  DGP3_vec_or_inf[,m]         <- or_3$att.inf.func
  DGP3_vec_dr_inf[,m]         <- dr_3$att.inf.func
  
  
  
  
  
  #################### DGP4 #################### 
  D <- as.numeric(p_X >= U)
  
  
  # Y_tD
  Y_00 <- f_reg_X + nu(X, D) + epsilon_0
  Y_10 <- 2 * f_reg_X + nu_X_0_1000  + epsilon_10
  Y_11 <- 2 * f_reg_X + nu_X_1_1000  + epsilon_11
  
  # complete data-frame creations
  time_0_correct <- data.frame(re = Y_00, t = 0, id = 1:1000, d = D,
                               Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  
  time_1_treat_0_all <- data.frame(re = Y_10, t = 1, id = 1:1000, d = D,
                                   Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  time_1_treat_1_all <- data.frame(re = Y_11, t = 1, id = 1:1000, d = D,
                                   Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4])
  
  time_1_correct <- rbind(time_1_treat_0_all[time_1_treat_0_all$d==0,],
                          (time_1_treat_1_all[time_1_treat_1_all$d==1,]))
  
  
  # finalizing the correct complete data-frame
  df <- rbind(time_0_correct, time_1_correct)
  colnames(df)[1] <- "re"
  
  
  dr_imp_4 <- drdid(yname="re", tname = "t", idname = "id", dname = "d",
                  xformla= ~ Z1+Z2+Z3+Z4,
                  data = df, panel = TRUE, inffunc = T)
  ipw_4 <- ipwdid(yname="re", tname = "t", idname = "id", dname = "d",
                xformla= ~ Z1+Z2+Z3+Z4,
                data = df, panel = TRUE, inffunc = T)
  or_4 <- ordid(yname="re", tname = "t", idname = "id", dname = "d",
              xformla= ~ Z1+Z2+Z3+Z4,
              data = df, panel = TRUE, inffunc = T)
  dr_4 <- drdid(yname="re", tname = "t", idname = "id", dname = "d",
                xformla= ~ Z1+Z2+Z3+Z4, estMethod = "trad",
                data = df, panel = TRUE, inffunc = T)
  twfe_4 <- twfe_did_panel(y1=time_1_correct$V1, y0=time_0_correct$V1,
                           D=D, covariates = as.matrix(cbind(df$Z1[1:n_sample], df$Z2[1:n_sample], df$Z3[1:n_sample], df$Z4[1:n_sample])), inffunc = F)
  
  
  DGP4_vec_dr_imp <- c(DGP4_vec_dr_imp, dr_imp_4$ATT)
  DGP4_vec_ipw    <- c(DGP4_vec_ipw, ipw_4$ATT)
  DGP4_vec_or     <- c(DGP4_vec_or, or_4$ATT)
  DGP4_vec_dr <-     c(DGP4_vec_dr, dr_4$ATT)
  DGP4_vec_twfe <-   c(DGP4_vec_twfe, twfe_4$ATT)
  DGP4_vec_dr_imp_cover <- c(DGP4_vec_dr_imp_cover, as.numeric(dr_imp_4$lci < 0) * as.numeric(dr_imp_4$uci > 0))
  DGP4_vec_ipw_cover    <- c(DGP4_vec_ipw_cover   , as.numeric(ipw_4$lci < 0) * as.numeric(ipw_4$uci > 0))
  DGP4_vec_or_cover     <- c(DGP4_vec_or_cover    , as.numeric(or_4$lci < 0) * as.numeric(or_4$uci > 0))
  DGP4_vec_dr_cover     <- c(DGP4_vec_dr_cover    , as.numeric(dr_4$lci < 0) * as.numeric(dr_4$uci > 0))
  DGP4_vec_twfe_cover   <- c(DGP4_vec_twfe_cover  , as.numeric(twfe_4$lci < 0) * as.numeric(twfe_4$uci > 0))
  
  
  DGP4_vec_dr_imp_inf[,m] <- dr_imp_4$att.inf.func
  DGP4_vec_ipw_inf[,m]       <- ipw_4$att.inf.func
  DGP4_vec_or_inf[,m]         <- or_4$att.inf.func
  DGP4_vec_dr_inf[,m]         <- dr_4$att.inf.func
  
}

end_time <- Sys.time()
end_time - start_time

# DR model treatment estimates

sapply(data.frame(DGP1_vec_dr_imp, DGP1_vec_ipw, DGP1_vec_or, DGP1_vec_dr), mean)
sapply(data.frame(DGP2_vec_dr_imp, DGP2_vec_ipw, DGP2_vec_or, DGP2_vec_dr), mean)
sapply(data.frame(DGP3_vec_dr_imp, DGP3_vec_ipw, DGP3_vec_or, DGP3_vec_dr), mean)
sapply(data.frame(DGP4_vec_dr_imp, DGP4_vec_ipw, DGP4_vec_or, DGP4_vec_dr), mean)

# Coverage Probability

sapply(data.frame(DGP1_vec_dr_imp_cover, DGP1_vec_ipw_cover, DGP1_vec_or_cover, DGP1_vec_dr_cover), mean)
sapply(data.frame(DGP2_vec_dr_imp_cover, DGP2_vec_ipw_cover, DGP2_vec_or_cover, DGP2_vec_dr_cover), mean)
sapply(data.frame(DGP3_vec_dr_imp_cover, DGP3_vec_ipw_cover, DGP3_vec_or_cover, DGP3_vec_dr_cover), mean)
sapply(data.frame(DGP4_vec_dr_imp_cover, DGP4_vec_ipw_cover, DGP4_vec_or_cover, DGP4_vec_dr_cover), mean)

# Variance
c(mean(colMeans(DGP1_vec_dr_imp_inf^2)), mean(colMeans(DGP1_vec_ipw_inf^2)), mean(colMeans(DGP1_vec_or_inf^2)), mean(colMeans(DGP1_vec_dr_inf^2)))
c(mean(colMeans(DGP2_vec_dr_imp_inf^2)), mean(colMeans(DGP2_vec_ipw_inf^2)), mean(colMeans(DGP2_vec_or_inf^2)), mean(colMeans(DGP2_vec_dr_inf^2)))
c(mean(colMeans(DGP3_vec_dr_imp_inf^2)), mean(colMeans(DGP3_vec_ipw_inf^2)), mean(colMeans(DGP3_vec_or_inf^2)), mean(colMeans(DGP3_vec_dr_inf^2)))
c(mean(colMeans(DGP4_vec_dr_imp_inf^2)), mean(colMeans(DGP4_vec_ipw_inf^2)), mean(colMeans(DGP4_vec_or_inf^2)), mean(colMeans(DGP4_vec_dr_inf^2)))

# c(mean(DGP1_vec_twfe),mean(DGP2_vec_twfe),mean(DGP3_vec_twfe),mean(DGP4_vec_twfe))
# c(mean(DGP1_vec_twfe_cover),mean(DGP2_vec_twfe_cover),mean(DGP3_vec_twfe_cover),mean(DGP4_vec_twfe_cover))
# c(mean(colMeans(DGP1_vec_twfe_inf^2)), mean(colMeans(DGP1_vec_twfe_inf^2)), mean(colMeans(DGP1_vec_twfe_inf^2)), mean(colMeans(DGP1_vec_twfe_inf^2)))

# hold_EF <- c()
# hold_EF_trial <- c()
# for(sira in 1:10)
# {
#   epsilon_0  <- rnorm(n_sample)
#   epsilon_10 <- rnorm(n_sample/2)
#   epsilon_11 <- rnorm(n_sample/2)
#   dft <- df
#   dft$Z1 <- rep(rnorm(1000, 0,1),2)
#   dft$Z2 <- rep(rnorm(1000, 0,1),2)
#   dft$Z3 <- rep(rnorm(1000, 0,1),2)
#   dft$Z4 <- rep(rnorm(1000, 0,1),2)
#   dft$d <- c(rep(1,500), rep(0,500), rep(1,500), rep(0,500))
#   dft$re <- c(rep(210,1000)+epsilon_0+rnorm(1000,dft$d*210, 1), rep(420,1000)+c(epsilon_11, epsilon_10)+rnorm(1000,dft$d*210, 1))
#   dft$id <- rep(1:1000, 2)
#   true_trial <- drdid(yname="re", tname = "t", idname = "id", dname = "d", data = dft, xformla= ~ Z1+Z2+Z3+Z4, panel = TRUE, inffunc = T)
#   true <-  drdid(yname="re", tname = "t", idname = "id", dname = "d", data = dft, panel = TRUE, inffunc = T)
#   hold_EF<- c(hold_EF, mean(true$att.inf.func^2))
#   hold_EF_trial <- c(hold_EF_trial, mean(true_trial$att.inf.func^2))
# }
# 
# mean(hold_EF)
# mean(hold_EF_trial)


