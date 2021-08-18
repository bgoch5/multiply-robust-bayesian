########################
# functions.R #
# Multiply Robust Bayesian Procedures 
# Gochanour et al. (2021+) 

########## PROJ ##########

fPROJ <- function(dat, w) {
  return(fPROJ1(dat, w) - fPROJ0(dat, w))
}

fgammaPROJ <- function(dat, C) { # Change C back to 500 later
  w0 <- rexp(dim(dat)[1] * C, 1)
  my_w0 <- matrix(w0, dim(dat)[1], C, byrow = T)
  fp <- function(a0) {
    return(fPROJ(dat, a0))
  }
  store <- apply(my_w0, 2, fp)
  my_mean <- apply(store, 1, mean) # Take mean accross each row
  low_quant <- apply(store, 1, lowquant_fun)
  hi_quant <- apply(store, 1, hiquant_fun)
  # cov_ind=as.numeric(low_quant<theta0 & hi_quant>theta0)
  # ci_len=hi_quant-low_quant
  my_result <- c(my_mean, low_quant, hi_quant)
  return(my_result)
}

fPROJ1 <- function(dat, w, freg1 = reg1, freg2 = reg2,
                   fprop1 = prop1,
                   fprop2 = prop2,
                   fdf_reg=df_reg1,
                   fdf_prop=df_prop) {
  
  r=dat$r
  y1=dat$y1
  ry=y1[r==1]
  rw <- w[r == 1]
  n <- length(y1)
  nr <- length(ry)
  nm <- n - nr
  
  fdf_reg[,"rw"]=rw
  fdf_prop[,"w"]=w
  
  
  ##################### END MODIFY #############################
  
  ### estimator under case (1010):
  glm11 <- glm(formula = fprop1, family = binomial(logit), data = fdf_prop, weights = w)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  etheta_MRP1010 <- 1 / sum(w) * (sum(r * y1 * w / esp1) + sum((1 - r / esp1) * esa1 * w))
  
  ### estimator under case 0110:
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  etheta_MRP0110 <- 1 / sum(w) * (sum(r * y1 * w / esp1) + sum((1 - r / esp1) * esa1 * w))
  
  ### estimator under case (1001):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  etheta_MRP1001 <- 1 / sum(w) * (sum(r * y1 * w / esp1) + sum((1 - r / esp1) * esa1 * w))
  
  ### estimator under case (0101):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  etheta_MRP0101 <- 1 / sum(w) * (sum(r * y1 * w / esp1) + sum((1 - r / esp1) * esa1 * w))
  
  ### estimator under case (1011):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  esp2 <- replace(esp2, esp2 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  lm_p <- lm(formula = r ~ 0 + esp1 + esp2)
  ebeta_p <- lm_p$coefficient
  ebeta_ps <- ebeta_p^{2} / sum(ebeta_p^{2})
  ep_1011 <- ebeta_ps[1] * esp1 + ebeta_ps[2] * esp2
  em_1011 <- esa1
  etheta_MRP1011 <- 1 / sum(w) * (sum(r * y1 * w / ep_1011) + sum((1 - r / ep_1011) * em_1011 * w))
  
  ### estimator under case (0111):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  esp2 <- replace(esp2, esp2 < .005, .005)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  lm_p <- lm(formula = r ~ 0 + esp1 + esp2)
  ebeta_p <- lm_p$coefficient
  ebeta_ps <- ebeta_p^{2} / sum(ebeta_p^{2})
  ep_0111 <- ebeta_ps[1] * esp1 + ebeta_ps[2] * esp2
  em_0111 <- esa1
  etheta_MRP0111 <- 1 / sum(w) * (sum(r * y1 * w / ep_0111) + sum((1 - r / ep_0111) * em_0111 * w))
  
  ### estimator under case (1110):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  resa1 <- esa1[r == 1]
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  resa2 <- esa2[r == 1]
  lm_m <- lm(formula = ry ~ 0 + resa1 + resa2)
  ebeta_m <- lm_m$coefficient
  ebeta_ms <- ebeta_m^{2} / sum(ebeta_m^{2})
  ep_1110 <- esp1
  em_1110 <- ebeta_ms[1] * esa1 + ebeta_ms[2] * esa2
  etheta_MRP1110 <- 1 / sum(w) * (sum(r * y1 * w / ep_1110) + sum((1 - r / ep_1011) * em_1110 * w))
  
  ### estimator under case (1101):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  resa1 <- esa1[r == 1]
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  resa2 <- esa2[r == 1]
  lm_m <- lm(formula = ry ~ 0 + resa1 + resa2)
  ebeta_m <- lm_m$coefficient
  ebeta_ms <- ebeta_m^{2} / sum(ebeta_m^{2})
  ep_1101 <- esp1
  em_1101 <- ebeta_ms[1] * esa1 + ebeta_ms[2] * esa2
  etheta_MRP1101 <- 1 / sum(w) * (sum(r * y1 * w / ep_1101) + sum((1 - r / ep_1101) * em_1101 * w))
  
  ### estimator under case (1111):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  esp2 <- replace(esp2, esp2 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  resa1 <- esa1[r == 1]
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  resa2 <- esa2[r == 1]
  lm_p <- lm(formula = r ~ 0 + esp1 + esp2)
  ebeta_p <- lm_p$coefficient
  ebeta_ps <- ebeta_p^{2} / sum(ebeta_p^{2})
  ep_1111 <- ebeta_ps[1] * esp1 + ebeta_ps[2] * esp2
  lm_m <- lm(formula = ry ~ 0 + resa1 + resa2)
  ebeta_m <- lm_m$coefficient
  ebeta_ms <- ebeta_m^{2} / sum(ebeta_m^{2})
  em_1111 <- ebeta_ms[1] * esa1 + ebeta_ms[2] * esa2
  etheta_MRP1111 <- 1 / sum(w) * (sum(r * y1 * w / ep_1111) + sum((1 - r / ep_1111) * em_1111 * w))
  
  etheta_MRP <- c(
    etheta_MRP1010, etheta_MRP0110, etheta_MRP1001, etheta_MRP0101, etheta_MRP1011,
    etheta_MRP0111, etheta_MRP1110, etheta_MRP1101, etheta_MRP1111
  )
  return(etheta_MRP)
}

fPROJ0 <- function(dat, w, freg1 = reg1, freg2 = reg2,
                   fprop1 = prop0_1,
                   fprop2 = prop0_2,
                   fdf_reg=df_reg0,
                   fdf_prop=df_prop) {
  
  r=dat$r
  y0=dat$y0
  ry=y0[r==0]
  rw <- w[r == 0]
  n <- length(y0)
  nr <- length(ry)
  nm <- n - nr
  fdf_reg[,"rw"]=rw
  fdf_prop[,"w"]=w
  rw=w[r==0]
  ### estimator under case (1010):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  etheta_MRP1010 <- 1 / sum(w) * (sum((1 - r) * y0 * w / esp1) + sum((1 - (1 - r) / esp1) * esa1 * w))
  
  ### estimator under case (0110):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  etheta_MRP0110 <- 1 / sum(w) * (sum((1 - r) * y0 * w / esp1) + sum((1 - (1 - r) / esp1) * esa1 * w))
  
  ### estimator under case (1001):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  etheta_MRP1001 <- 1 / sum(w) * (sum((1 - r) * y0 * w / esp1) + sum((1 - (1 - r) / esp1) * esa1 * w))
  
  ### estimator under case (0101):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  etheta_MRP0101 <- 1 / sum(w) * (sum((1 - r) * y0 * w / esp1) + sum((1 - (1 - r) / esp1) * esa1 * w))
  
  ### estimator under case (1011):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  esp2 <- replace(esp2, esp2 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  lm_p <- lm(formula = (1 - r) ~ 0 + esp1 + esp2)
  ebeta_p <- lm_p$coefficient
  ebeta_ps <- ebeta_p^{2} / sum(ebeta_p^{2})
  ep_1011 <- ebeta_ps[1] * esp1 + ebeta_ps[2] * esp2
  em_1011 <- esa1
  etheta_MRP1011 <- 1 / sum(w) * (sum((1 - r) * y0 * w / ep_1011) + sum((1 - (1 - r) / ep_1011) * em_1011 * w))
  
  ### estimator under case (0111):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  esp2 <- replace(esp2, esp2 < .005, .005)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  lm_p <- lm(formula = (1 - r) ~ 0 + esp1 + esp2)
  ebeta_p <- lm_p$coefficient
  ebeta_ps <- ebeta_p^{2} / sum(ebeta_p^{2})
  ep_0111 <- ebeta_ps[1] * esp1 + ebeta_ps[2] * esp2
  em_0111 <- esa1
  etheta_MRP0111 <- 1 / sum(w) * (sum((1 - r) * y0 * w / ep_0111) + sum((1 - (1 - r) / ep_0111) * em_0111 * w))
  
  ### estimator under case (1110):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  resa1 <- esa1[r == 0]
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  resa2 <- esa2[r == 0]
  lm_m <- lm(formula = ry ~ 0 + resa1 + resa2)
  ebeta_m <- lm_m$coefficient
  ebeta_ms <- ebeta_m^{2} / sum(ebeta_m^{2})
  ep_1110 <- esp1
  em_1110 <- ebeta_ms[1] * esa1 + ebeta_ms[2] * esa2
  etheta_MRP1110 <- 1 / sum(w) * (sum((1 - r) * y0 * w / ep_1110) + sum((1 - (1 - r) / ep_1110) * em_1110 * w))
  
  ### estimator under case (1101):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  resa1 <- esa1[r == 0]
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  resa2 <- esa2[r == 0]
  lm_m <- lm(formula = ry ~ 0 + resa1 + resa2)
  ebeta_m <- lm_m$coefficient
  ebeta_ms <- ebeta_m^{2} / sum(ebeta_m^{2})
  ep_1101 <- esp1
  em_1101 <- ebeta_ms[1] * esa1 + ebeta_ms[2] * esa2
  etheta_MRP1101 <- 1 / sum(w) * (sum((1 - r) * y0 * w / ep_1101) + sum((1 - (1 - r) / ep_1101) * em_1101 * w))
  
  ### estimator under case (1111):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  esp1 <- replace(esp1, esp1 < .005, .005)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  esp2 <- replace(esp2, esp2 < .005, .005)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  resa1 <- esa1[r == 0]
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  resa2 <- esa2[r == 0]
  lm_p <- lm(formula = (1 - r) ~ 0 + esp1 + esp2)
  ebeta_p <- lm_p$coefficient
  ebeta_ps <- ebeta_p^{ 2} / sum(ebeta_p^{2})
  ep_1111 <- ebeta_ps[1] * esp1 + ebeta_ps[2] * esp2
  lm_m <- lm(formula = ry ~ 0 + resa1 + resa2)
  ebeta_m <- lm_m$coefficient
  ebeta_ms <- ebeta_m^
    {
      2
    } / sum(ebeta_m^{
      2
    })
  em_1111 <- ebeta_ms[1] * esa1 + ebeta_ms[2] * esa2
  etheta_MRP1111 <- 1 / sum(w) * (sum((1 - r) * y0 * w / ep_1111) + sum((1 - (1 - r) / ep_1111) * em_1111 * w))
  
  etheta_MRP <- c(
    etheta_MRP1010, etheta_MRP0110, etheta_MRP1001, etheta_MRP0101, etheta_MRP1011,
    etheta_MRP0111, etheta_MRP1110, etheta_MRP1101, etheta_MRP1111
  )
  return(etheta_MRP)
}


########## 6. CAL ##########

fCAL <- function(dat, w) {
  return(fCAL1(dat, w) - fCAL0(dat, w))
}

fgammaCAL <- function(dat, C) { 
  w0 <- rexp(dim(dat)[1] * C, 1)
  my_w0 <- matrix(w0, dim(dat)[1], C, byrow = T)
  fp <- function(a0) {
    return(fCAL(dat, a0))
  }
  store <- apply(my_w0, 2, fp)
  my_mean <- apply(store, 1, mean) # Take mean accross each row
  low_quant <- apply(store, 1, lowquant_fun)
  hi_quant <- apply(store, 1, hiquant_fun)
  # cov_ind=as.numeric(low_quant<theta0 & hi_quant>theta0)
  # ci_len=hi_quant-low_quant
  my_result <- c(my_mean, low_quant, hi_quant)
  return(my_result)
}

fCAL1 <- function(dat, w, freg1 = reg1, freg2 = reg2,
                  fprop1 = prop1,
                  fprop2 = prop2,
                  fdf_reg=df_reg1,
                  fdf_prop=df_prop) {
  
  r=dat$r
  y1=dat$y1
  ry=y1[r==1]
  rw <- w[r == 1]
  n <- length(y1)
  nr <- length(ry)
  nm <- n - nr
  rw=w[r==1]
  fdf_reg[,"rw"]=rw
  fdf_prop[,"w"]=w
  
  ### estimator under case (1010):
  # Model updates after this one
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg) # rw is ws for those in this treatment group
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 1]
  ruu2 <- esa1[r == 1]
  
  uu1 <- esp1
  uu2 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2))
  my_u2 <- as.matrix(cbind(uu1, uu2))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP_1010 <- eomegaP_1010 <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1))
  etheta_MRC1010 <- sum(eomegaP_1010 * ry)
  
  
  ### estimator under case (0110):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 1]
  ruu2 <- esa1[r == 1]
  
  uu1 <- esp1
  uu2 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2))
  my_u2 <- as.matrix(cbind(uu1, uu2))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP_0110 <- eomegaP_1010 <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1))
  etheta_MRC0110 <- sum(eomegaP_0110 * ry)
  
  
  ### estimator under case (1001):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 1]
  ruu2 <- esa1[r == 1]
  
  uu1 <- esp1
  uu2 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2))
  my_u2 <- as.matrix(cbind(uu1, uu2))
  
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- eomegaP_1010 <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1))
  etheta_MRC1001 <- sum(eomegaP * ry)
  
  ### estimator under case (0101):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 1]
  ruu2 <- esa1[r == 1]
  
  uu1 <- esp1
  uu2 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2))
  my_u2 <- as.matrix(cbind(uu1, uu2))
  
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- eomegaP_1010 <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1))
  etheta_MRC0101 <- sum(eomegaP * ry)
  
  ### estimator under case (1011):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  # esp2<-replace(esp2,esp2<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 1]
  ruu2 <- esp2[r == 1]
  ruu3 <- esa1[r == 1]
  
  uu1 <- esp1
  uu2 <- esp2
  uu3 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2, ruu3))
  my_u2 <- as.matrix(cbind(uu1, uu2, uu3))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1))
  
  
  
  etheta_MRC1011 <- sum(eomegaP * ry)
  
  
  
  ### estimator under case (0111):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  # esp2<-replace(esp2,esp2<.01,.01)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 1]
  ruu2 <- esp2[r == 1]
  ruu3 <- esa1[r == 1]
  
  uu1 <- esp1
  uu2 <- esp2
  uu3 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2, ruu3))
  my_u2 <- as.matrix(cbind(uu1, uu2, uu3))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1))
  etheta_MRC0111 <- sum(eomegaP * ry)
  
  ### estimator under case (1110):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  
  ruu1 <- esp1[r == 1]
  ruu2 <- esa1[r == 1]
  ruu3 <- esa2[r == 1]
  
  uu1 <- esp1
  uu2 <- esa1
  uu3 <- esa2
  
  my_u <- as.matrix(cbind(ruu1, ruu2, ruu3))
  my_u2 <- as.matrix(cbind(uu1, uu2, uu3))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1))
  etheta_MRC1110 <- sum(eomegaP * ry)
  
  ### estimator under case (1101):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  # esa1=predict(glm21,df_prop)
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  
  ruu1 <- esp1[r == 1]
  ruu2 <- esa1[r == 1]
  ruu3 <- esa2[r == 1]
  
  uu1 <- esp1
  uu2 <- esa1
  uu3 <- esa2
  
  my_u <- as.matrix(cbind(ruu1, ruu2, ruu3))
  my_u2 <- as.matrix(cbind(uu1, uu2, uu3))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1))
  etheta_MRC1101 <- sum(eomegaP * ry)
  
  ### estimator under case (1111):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  # esp2<-replace(esp2,esp2<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  
  ruu1 <- esp1[r == 1]
  ruu2 <- esp2[r == 1]
  ruu3 <- esa1[r == 1]
  ruu4 <- esa2[r == 1]
  
  uu1 <- esp1
  uu2 <- esp2
  uu3 <- esa1
  uu4 <- esa2
  
  my_u <- as.matrix(cbind(ruu1, ruu2, ruu3, ruu4))
  my_u2 <- as.matrix(cbind(uu1, uu2, uu3, uu4))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]) + elambdaP[4] * (ruu4 - my_mu[4]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]) + elambdaP[4] * (ruu4 - my_mu[4]))^(-1))
  
  etheta_MRC1111 <- sum(eomegaP * ry)
  
  etheta_MRC <- c(
    etheta_MRC1010, etheta_MRC0110, etheta_MRC1001, etheta_MRC0101,
    etheta_MRC1011, etheta_MRC0111, etheta_MRC1110, etheta_MRC1101, etheta_MRC1111
  )
  return(etheta_MRC)
}


fCAL0 <- function(dat, w, freg1 = reg1, freg2 = reg2,
                  fprop1 = prop0_1,
                  fprop2 = prop0_2,
                  fdf_reg=df_reg0,
                  fdf_prop=df_prop) {
  
  r=dat$r
  y0=dat$y0
  ry=y0[r==0]
  rw <- w[r == 0]
  n <- length(y0)
  nr <- length(ry)
  nm <- n - nr
  rw=w[r==0]
  fdf_reg[,"rw"]=rw
  fdf_prop[,"w"]=w
  
  ### estimator under case (1010):
  # Model updates after this one
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg) # rw is ws for those in this treatment group
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 0]
  ruu2 <- esa1[r == 0]
  
  uu1 <- esp1
  uu2 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2))
  my_u2 <- as.matrix(cbind(uu1, uu2))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP_1010 <- eomegaP_1010 <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1))
  etheta_MRC1010 <- sum(eomegaP_1010 * ry)
  
  
  ### estimator under case (0110):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 0]
  ruu2 <- esa1[r == 0]
  
  uu1 <- esp1
  uu2 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2))
  my_u2 <- as.matrix(cbind(uu1, uu2))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  
  eomegaP_0110 <- eomegaP_1010 <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1))
  etheta_MRC0110 <- sum(eomegaP_0110 * ry)
  
  
  ### estimator under case (1001):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 0]
  ruu2 <- esa1[r == 0]
  
  uu1 <- esp1
  uu2 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2))
  my_u2 <- as.matrix(cbind(uu1, uu2))
  
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- eomegaP_1010 <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1))
  etheta_MRC1001 <- sum(eomegaP * ry)
  
  ### estimator under case (0101):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 0]
  ruu2 <- esa1[r == 0]
  
  uu1 <- esp1
  uu2 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2))
  my_u2 <- as.matrix(cbind(uu1, uu2))
  
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- eomegaP_1010 <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]))^(-1))
  etheta_MRC0101 <- sum(eomegaP * ry)
  
  ### estimator under case (1011):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  # esp2<-replace(esp2,esp2<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 0]
  ruu2 <- esp2[r == 0]
  ruu3 <- esa1[r == 0]
  
  uu1 <- esp1
  uu2 <- esp2
  uu3 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2, ruu3))
  my_u2 <- as.matrix(cbind(uu1, uu2, uu3))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1))
  etheta_MRC1011 <- sum(eomegaP * ry)
  
  ### estimator under case (0111):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  # esp2<-replace(esp2,esp2<.01,.01)
  glm21 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  
  ruu1 <- esp1[r == 0]
  ruu2 <- esp2[r == 0]
  ruu3 <- esa1[r == 0]
  
  uu1 <- esp1
  uu2 <- esp2
  uu3 <- esa1
  
  my_u <- as.matrix(cbind(ruu1, ruu2, ruu3))
  my_u2 <- as.matrix(cbind(uu1, uu2, uu3))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1))
  etheta_MRC0111 <- sum(eomegaP * ry)
  
  ### estimator under case (1110):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  
  ruu1 <- esp1[r == 0]
  ruu2 <- esa1[r == 0]
  ruu3 <- esa2[r == 0]
  
  uu1 <- esp1
  uu2 <- esa1
  uu3 <- esa2
  
  my_u <- as.matrix(cbind(ruu1, ruu2, ruu3))
  my_u2 <- as.matrix(cbind(uu1, uu2, uu3))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1))
  etheta_MRC1110 <- sum(eomegaP * ry)
  
  ### estimator under case (1101):
  glm11 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  
  ruu1 <- esp1[r == 0]
  ruu2 <- esa1[r == 0]
  ruu3 <- esa2[r == 0]
  
  uu1 <- esp1
  uu2 <- esa1
  uu3 <- esa2
  
  my_u <- as.matrix(cbind(ruu1, ruu2, ruu3))
  my_u2 <- as.matrix(cbind(uu1, uu2, uu3))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]))^(-1))
  etheta_MRC1101 <- sum(eomegaP * ry)
  
  ### estimator under case (1111):
  glm11 <- glm(formula = fprop1, family = binomial(logit), weights = w, data = fdf_prop)
  esp1 <- glm11$fitted.values
  # esp1<-replace(esp1,esp1<.01,.01)
  glm12 <- glm(formula = fprop2, family = binomial(logit), weights = w, data = fdf_prop)
  esp2 <- glm12$fitted.values
  # esp2<-replace(esp2,esp2<.01,.01)
  glm21 <- glm(formula = freg1, weights = rw, data = fdf_reg)
  ebeta1 <- glm21$coefficient
  esa1 <- predict(glm21, df_prop)
  glm22 <- glm(formula = freg2, weights = rw, data = fdf_reg)
  ebeta2 <- glm22$coefficient
  esa2 <- predict(glm22, df_prop)
  
  ruu1 <- esp1[r == 0]
  ruu2 <- esp2[r == 0]
  ruu3 <- esa1[r == 0]
  ruu4 <- esa2[r == 0]
  
  uu1 <- esp1
  uu2 <- esp2
  uu3 <- esa1
  uu4 <- esa2
  
  my_u <- as.matrix(cbind(ruu1, ruu2, ruu3, ruu4))
  my_u2 <- as.matrix(cbind(uu1, uu2, uu3, uu4))
  
  my_mu <- as.matrix(apply(w %*% my_u2, 2, sum) / sum(w))
  elambdaP <- Lag2(u = my_u, ds = rw, mu = my_mu)
  
  eomegaP <- rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]) + elambdaP[4] * (ruu4 - my_mu[4]))^(-1) / sum(rw * (1 + elambdaP[1] * (ruu1 - my_mu[1]) + elambdaP[2] * (ruu2 - my_mu[2]) + elambdaP[3] * (ruu3 - my_mu[3]) + elambdaP[4] * (ruu4 - my_mu[4]))^(-1))
  etheta_MRC1111 <- sum(eomegaP * ry)
  
  etheta_MRC <- c(
    etheta_MRC1010, etheta_MRC0110, etheta_MRC1001, etheta_MRC0101,
    etheta_MRC1011, etheta_MRC0111, etheta_MRC1110, etheta_MRC1101, etheta_MRC1111
  )
  return(etheta_MRC)
}

######### Misc functions ##########
lowquant_fun <- function(x) {
  return(quantile(x, probs = 0.025))
}
hiquant_fun <- function(x) {
  return(quantile(x, probs = 0.975))
}


fsamp <- function(tt) {
  ress <- sample(tt, 1)
  return(ress)
}

######################################################
#  R function for solving g_1(lambda)=0 in Wu (2005) #
######################################################

#########################################
# R Function for EL Lagrange Multiplier #
#     x is of dimension 2 or higher     #
# Input: u=(x_1,x_2,...,x_n)            #
#        ds=(d_1,d_2,...,d_n)           #
#        (Design Wights: d_i=1/pi_i)    #
#        ds=(1,1,...1) for iid data     #
#        mu: benchmark means for x      #
# Output: lambda=M                      #
#                                       #
# Written by Changbao Wu, March, 2000   #
# Modified by Randy Sitter, June, 2006  #
#########################################

Lag2=function(u,ds,mu){ 
  n=length(ds)
  u=u-rep(1,n)%*%t(mu)
  M=0*mu
  dif=1
  tol=1e-8
  k=0
  while(dif>tol & k<=50){
  D1=t(u)%*%((ds/(1+u%*%M))*rep(1,n))
  DD=-t(u)%*%(c((ds/(1+u%*%M)^2))*u)
  D2=solve(DD,D1,tol=1e-40)
  dif=max(abs(D2))
  rule=1
  while(rule>0){
    rule=0
    if(min(1+t(M-D2)%*%t(u))<=0) rule=rule+1
    if(rule>0) D2=D2/2
               }
    M=M-D2
    k=k+1
                   }
  if(k>=50) M=0*mu
  return(as.vector(M))
          }