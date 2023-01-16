# Load the necessary packages
library(forecast)
library(dplyr)
library(rugarch)

# load the data in Merck.csv
mydata <- read.delim("Merck.csv", header = TRUE,  sep = ",")

date <- as.Date(mydata$Date, format = "%d/%m/%Y")
y <- mydata$y
r <- diff(log(y))

# plot the returns
plot(date[-1], r, type = "l", xlab = "", ylab = "returns")

# 1
# (a)
submods <- c("sGARCH", "gjrGARCH")
ARMA_TGARCH_est <- list()
ic_arma_tgarch <- matrix(nrow = 4 ^ 2 * 2 ^ 3,
                          ncol = 7)
colnames(ic_arma_tgarch) <- c("pm", "qm", "ph", "qh",
                              "thresh", "aic", "bic")
i <- 0; t0 <- proc.time()
for (pm in 0:3)
{
  for (qm in 0:3)
  {
    for (ph in 0:1)
    {
      for (qh in 0:1)
      {
        for (thresh in 0:1)
        {
          i <- i + 1
          ic_arma_tgarch[i, 1:5] <- c(pm, qm, ph, qh, thresh)
          
          if (ph == 0 && qh == 0)
          {
            if (!thresh)
            {
              # For models without ARCH/GARCH effects use arfimaspec and arfimafit
              try(silent = T, expr =
              {
                ARMA_TGARCH_mod <- arfimaspec(
                  mean.model = list(armaOrder = c(pm, qm)))
                
                ARMA_TGARCH_est[[i]] <- arfimafit(ARMA_TGARCH_mod, r)
                
                ic_arma_tgarch[i,6:7] <- infocriteria(
                  ARMA_TGARCH_est[[i]])[1:2]
              })
            }
          }
          else
          {
            # For models with ARCH/GARCH effects use ugarchspec and ugarchfit
            # Specify model sGARCH or gjrGARCH, if the specification has (or not)
            # a Threshold term.
            try(silent = T, expr =
            {
              ARMA_TGARCH_mod <- ugarchspec(
                mean.model = list(armaOrder = c(pm, qm)),
                variance.model = list(model = submods[1 + thresh],
                                      garchOrder = c(ph, qh)))
              
              ARMA_TGARCH_est[[i]] <- ugarchfit(ARMA_TGARCH_mod, r,
                                                solver = 'hybrid')
              
              ic_arma_tgarch[i,6:7] <- infocriteria(ARMA_TGARCH_est[[i]])[1:2]
            })
          }
          cat("\r", "Processed ARMA(",
              pm, ",", qm, ")-",
              c("", "T")[1 + thresh], "GARCH(",
              ph, ",", qh, ")",
              c(" ", "")[1 + thresh], sep = "")
          
        }
      }
    }
  }
}
cat("\n", proc.time() - t0, "\n")

# Capture the first 20 models in terms of AIC/BIC
ic_aic_arma_tgarch <- ic_arma_tgarch[order(ic_arma_tgarch[,6]),][1:20,]
ic_bic_arma_tgarch <- ic_arma_tgarch[order(ic_arma_tgarch[,7]),][1:20,]

# find the intersection of AIC and BIC preferred sets
ic_int_arma_tgarch <- intersect(as.data.frame(ic_aic_arma_tgarch),
                                as.data.frame(ic_bic_arma_tgarch))

# We select the entire intersection set.
adq_set_arma_tgarch <- as.matrix(arrange(as.data.frame(
                                  ic_int_arma_tgarch), pm, qm, ph, qh, thresh))
adq_idx_arma_tgarch <- match(data.frame(t(adq_set_arma_tgarch[, 1:5])),
                                  data.frame(t(ic_arma_tgarch[, 1:5])))

# Check the residuals
nmods <- length(adq_idx_arma_tgarch)
sacf_tgarch <- matrix(nrow = nmods, ncol = 15)
colnames(sacf_tgarch) <- c("pm", "qm", "ph", "qh", "thresh", 1:10)
for (i in 1:nmods)
{
  sacf_tgarch[i,1:5] <- adq_set_arma_tgarch[i,1:5]
  sacf_tgarch[i,6:15] <-
                  acf(ARMA_TGARCH_est[[adq_idx_arma_tgarch[i]]]@fit$z,
                                       lag = 10, plot = F)$acf[2:11]
}

# (b)
# Plot estimated volatilities
title_tgarch <- rep(NA, times = nmods)
for (i in 1:nmods)
{
  title_tgarch[i] <- paste("ARMA(",
                     as.character(adq_set_arma_tgarch[i, 1]), ",",
                     as.character(adq_set_arma_tgarch[i, 2]),
                     ")-",
                     c("", "T")[1 + adq_set_arma_tgarch[i,5]],
                     "GARCH(",
                     as.character(adq_set_arma_tgarch[i, 3]), ",",
                     as.character(adq_set_arma_tgarch[i, 4]), ")",
                     sep = "")
  plot(date[-1], sqrt(
        ARMA_TGARCH_est[[adq_idx_arma_tgarch[i]]]@fit$var),
        type = "l", xlab = "", ylab = "volatilities",
        ylim = c(0, 0.08), main = title_tgarch[i])
}


# (c)
# To test for leverage effects
for (i in 1:nmods)
{
  if (adq_set_arma_tgarch[i, 5] == 1)
  {
    # this is a specification with a threshold
    lambda_est <- ARMA_TGARCH_est[[
              adq_idx_arma_tgarch[i]]]@fit$coef["gamma1"]
    lambda_tvl <- ARMA_TGARCH_est[[
              adq_idx_arma_tgarch[i]]]@fit$tval["gamma1"]
    cat(paste0(title_tgarch[i], ": lambda_hat = ",
                                   round(lambda_est, 2),
                                ", t-value = ",
                                   round(lambda_tvl, 2)),
                                "\n")
  }
}

# (d)
# Plot forecasted volatilities
for (i in 1:nmods)
{
  plot(ugarchboot(ARMA_TGARCH_est[[adq_idx_arma_tgarch[i]]],
                  method = "Partial"), which = 3)
}

# 2
# (a)
# Estimate ARMA(p,q)-GARCH(1,1) models with and without GARCH-M.
ARMA_GARCHM_est <- list()
ic_arma_garchm <- matrix( nrow = 4 ^ 2 * 2, ncol = 5 )
colnames(ic_arma_garchm) <- c("p", "q", "m", "aic", "bic")
i <- 0; t0 <- proc.time()
for (p in 0:3)
{
  for (q in 0:3)
  {
    for (m in 0:1)
    {
      i <- i + 1
      ic_arma_garchm[i, 1:3] <- c(p, q, m)
        
      try(silent = T, expr =
      {
        ARMA_GARCHM_mod <- ugarchspec(
                            mean.model = list(armaOrder = c(p, q),
                                              archm = m),
                            variance.model = list(garchOrder = c(1, 1)))
                    
        ARMA_GARCHM_est[[i]] <- ugarchfit(ARMA_GARCHM_mod, r,
                                          solver = 'hybrid')
                  
        ic_arma_garchm[i,4:5] <- infocriteria(ARMA_GARCHM_est[[i]])[1:2]
      })
      cat("\r", "Processed ARMA(",p, ",", q, ")",
              c("-GARCH(1,1) ", "-GARCHM(1,1)")[1 + m], sep = "")
          
    }
  }
}
cat("\n", proc.time() - t0, "\n")

ic_aic_arma_garchm <- ic_arma_garchm[order(ic_arma_garchm[,4]),][1:20,]
ic_bic_arma_garchm <- ic_arma_garchm[order(ic_arma_garchm[,5]),][1:20,]

# find the intersection of AIC and BIC preferred sets
ic_int_arma_garchm <- intersect(as.data.frame(ic_aic_arma_garchm),
                                as.data.frame(ic_bic_arma_garchm))

# We select the entire intersection set.
adq_set_arma_garchm <- as.matrix(arrange(as.data.frame(
                                ic_int_arma_garchm), p, q, m))
adq_idx_arma_garchm <- match(data.frame(t(adq_set_arma_garchm[, 1:3])),
                                  data.frame(t(ic_arma_garchm[, 1:3])))

# Check the residuals
nmods <- length(adq_idx_arma_garchm)
sacf_garchm <- matrix(nrow = nmods, ncol = 13)
colnames(sacf_garchm) <- c("p", "q", "m", 1:10)
for (i in 1:nmods)
{
  sacf_garchm[i,1:3] <- adq_set_arma_garchm[i,1:3]
  sacf_garchm[i,4:13] <-
                  acf(ARMA_GARCHM_est[[adq_idx_arma_garchm[i]]]@fit$z,
                                       lag = 10, plot = F)$acf[2:11]
}

# (b)
title_garchm <- rep(NA, times = nmods)
for (i in 1:nmods)
{
  title_garchm[i] <- paste("ARMA(",
        as.character(adq_set_arma_garchm[i, 1]), ",",
        as.character(adq_set_arma_garchm[i, 2]),
        c(")-GARCH(1,1)", ")-GARCHM(1,1)")[1 + adq_set_arma_garchm[i,3]],
        sep = "")
  plot(date[-1], sqrt(
        ARMA_GARCHM_est[[adq_idx_arma_garchm[i]]]@fit$var),
        type = "l", xlab = "", ylab = "volatilities",
        ylim = c(0, 0.08), main = title_garchm[i])
}


# (c)
# To test for time-varying risk premium
for (i in 1:nmods)
{
  if (adq_set_arma_garchm[i, 3] == 1)
  {
    # this is a specification with GARCH-in-mean
    archm_est <- ARMA_GARCHM_est[[
                          adq_idx_arma_garchm[i]]]@fit$coef["archm"]
    archm_tvl <- ARMA_GARCHM_est[[
                          adq_idx_arma_garchm[i]]]@fit$tval["archm"]
    cat(paste0(title_garchm[i], ": archm_hat = ",
                                   round(archm_est, 2),
                                ", t-value = ",
                                   round(archm_tvl, 2)),
                                "\n")
  }
}


# (d)
#
for (i in 1:nmods)
{
  plot(ugarchboot(ARMA_GARCHM_est[[adq_idx_arma_garchm[i]]],
                                   method = "Partial"), which = 3)
}