install.packages("forecast")
install.packages("dplyr")
install.packages("zoo")
install.packages("aTSA")

library(forecast)
library(dplyr)
library(zoo)
library(aTSA)

# 1
#
# (a)
#
# load the data in usdata.csv
mydata <- read.delim("usdata.csv", header = TRUE,  sep = ",")

dates <- as.yearqtr(mydata$obs)
y <- mydata$GDP
r <- mydata$FFR

plot(dates, y, type = "l", xlab = "Time (Quarters)",
     main = "Log Real US GDP per Capita")

# (b)
#
# Estimate ADF regressions for different values of p
# plus combinations of constant and/or trend.
TT <- length(y)
ADF_est <- list()
ic <- matrix(nrow = 30, ncol = 5)
colnames(ic) <- c("cons", "trend", "p", "aic", "bic")
i <- 0
for (const in 0:1)
{
  for (p in 0:9)
  {
    i <- i + 1
    ADF_est[[i]] <- Arima(diff(y), xreg = y[-TT],
                          order = c(p, 0, 0),
                          include.mean = as.logical(const),
                          include.drift = F)
    ic[i,] <- c(const, 0, p, ADF_est[[i]]$aic,
                             ADF_est[[i]]$bic)
  }
  
  if (const)
  {
    # only add a specification with trend if there is a
    # constant (i.e., exclude no constant with trend)
    for (p in 0:9)
    {
      i <- i + 1
      ADF_est[[i]] <- Arima(diff(y), xreg = y[-TT],
                            order = c(p, 0, 0),
                            include.mean = as.logical(const),
                            include.drift = T)
      ic[i,] <- c(const, 1, p, ADF_est[[i]]$aic,
                               ADF_est[[i]]$bic)
    }
  }
}

ic_aic <- ic[order(ic[,4]),][1:10,]
ic_bic <- ic[order(ic[,5]),][1:10,]

# Select top 5 BIC, add to an object and obtain index
adq_set <- as.matrix(arrange(as.data.frame(ic_bic[1:5,]),
                   const, trend, p))
adq_idx <- match(data.frame(t(adq_set[, 1:3])),
                 data.frame(t(ic[, 1:3])))

# Check residuals
for (i in 1:length(adq_idx))
{
  checkresiduals(ADF_est[[adq_idx[i]]])
}

# All residuals look OK.

# (c)
#
adf.test(y)

# (d)
#
# repeat for differenced y series
TT <- length(diff(y))
ADF_est_diff <- list()
ic_diff <- matrix( nrow = 30, ncol = 5 )
colnames(ic_diff) <- c("cons", "trend", "p", "aic", "bic")
i <- 0
for (const in 0:1)
{
  for (p in 0:9)
  {
    i <- i + 1
    ADF_est_diff[[i]] <- Arima(diff(diff(y)),
                               xreg = diff(y)[-TT],
                               order = c(p, 0, 0),
                               include.mean = as.logical(const),
                               include.drift = F)
    ic_diff[i,] <- c(const, 0, p, ADF_est_diff[[i]]$aic,
                                  ADF_est_diff[[i]]$bic)
  }
  
  if (const)
  {
    for (p in 0:9)
    {
      i <- i + 1
      ADF_est_diff[[i]] <- Arima(diff(diff(y)),
                                 xreg = diff(y)[-TT],
                                 order = c(p, 0, 0),
                                 include.mean = as.logical(const),
                                 include.drift = T)
      ic_diff[i,] <- c(const, 1, p, ADF_est_diff[[i]]$aic,
                                    ADF_est_diff[[i]]$bic)
    }
  }
}

ic_aic_diff <- ic_diff[order(ic_diff[,4]),][1:10,]
ic_bic_diff <- ic_diff[order(ic_diff[,5]),][1:10,]

# AIC and BIC agree on the top 5

adq_set_diff <- as.matrix(arrange(as.data.frame(
                      ic_bic_diff[1:5,]), const, trend, p))
adq_idx_diff <- match(data.frame(t(adq_set_diff[, 1:3])),
                      data.frame(t(ic_diff[, 1:3])))

for (i in 1:length(adq_idx_diff))
{
  checkresiduals(ADF_est_diff[[adq_idx_diff[i]]])
}

# all residuals look ok

adf.test(diff(y))

# (f)
TT <- length(y)
ARIMA_est <- list()
ic_arima <- matrix(nrow = (2 * 2 + 2) * 4 ^ 2, ncol = 7)
colnames(ic_arima) <- c("d", "cons", "trend", "p", "q", "aic", "bic")
i <- 0
for (d in 0:1)
{
  for (const in 0:1)
  {
    for (p in 0:3)
    {
      for (q in 0:3)
      {
        i <- i + 1
        d1 <- as.logical(d)
        c1 <- as.logical(const)
        
        try(silent = T, expr =
        {
        ARIMA_est[[i]] <- Arima(y, order = c(p, d, q),
                                include.constant = c1)
        
        ic_arima[i,] <- c(d, const, 0, p, q,
                          ARIMA_est[[i]]$aic,
                          ARIMA_est[[i]]$bic)
        })

        if (const)
        {
          # only add a specification with trend if there is a
          # constant (i.e., exclude no constant with trend)
          i <- i + 1
          
          if (d1)
          {
            x <- c(0,cumsum(1:(TT - 1)))
          }
          else
          {
            x <- NULL
          }

          try(silent = T, expr =
          {
          ARIMA_est[[i]] <- Arima(y, order = c(p, d, q),
                                  xreg = x,
                                  include.constant = c1,
                                  include.drift = T)
          
          ic_arima[i,] <- c(d, const, 1, p, q,
                            ARIMA_est[[i]]$aic,
                            ARIMA_est[[i]]$bic)
          })
        }
      }
    }
  }
}

ic_aic_arima <- ic_arima[order(ic_arima[,6]),][1:10,]
ic_bic_arima <- ic_arima[order(ic_arima[,7]),][1:10,]

# find the intersection of AIC and BIC preferred sets
ic_int_arima <- intersect(as.data.frame(ic_aic_arima),
                          as.data.frame(ic_bic_arima))

# Add ARIMA(1,1,0) and ARIMA(2,1,0) to the four in the intersecting set.

adq_set_arima <- as.matrix(arrange(as.data.frame(
                           rbind(ic_int_arima,
                                 ic_bic_arima[c(1, 3),])),
                                      d, const, trend, p))
adq_idx_arima <- match(data.frame(t(adq_set_arima[, 1:5])),
                       data.frame(t(ic_arima[, 1:5])))

# Check the residuals for specifications in the adequate set.
for (i in 1:length(adq_idx_arima))
{
  checkresiduals(ARIMA_est[[adq_idx_arima[i]]])
}

# All residuals look OK.


# 2

# (a)

plot(dates, r, type = "l", xlab = "Time (Quarters)",
     main = "Federal Funds Rate")

# (b)
TT <- length(r)
ADF_est <- list()
ic <- matrix( nrow = 30, ncol = 5 )
colnames(ic) <- c("cons", "trend", "p", "aic", "bic")
i <- 0
for (const in 0:1)
{
  for (p in 0:9)
  {
    i <- i + 1
    ADF_est[[i]] <- Arima(diff(r), xreg = r[-TT],
                          order = c(p, 0, 0),
                          include.mean = as.logical(const),
                          include.drift = F)
    ic[i,] <- c(const, 0, p, ADF_est[[i]]$aic,
                ADF_est[[i]]$bic)
  }
  
  if (const)
  {
    # only add a specification with trend if there is a
    # constant (i.e., exclude no constant with trend)
    for (p in 0:9)
    {
      i <- i + 1
      ADF_est[[i]] <- Arima(diff(r), xreg = r[-TT],
                            order = c(p, 0, 0),
                            include.mean = as.logical(const),
                            include.drift = T)
      ic[i,] <- c(const, 1, p, ADF_est[[i]]$aic,
                  ADF_est[[i]]$bic)
    }
  }
}

ic_aic <- ic[order(ic[,4]),][1:10,]
ic_bic <- ic[order(ic[,5]),][1:10,]

# find the intersection of AIC and BIC preferred sets
ic_int <- intersect(as.data.frame(ic_aic),
                    as.data.frame(ic_bic))

# Adding intersecting set to an object and obtain index
adq_set <- as.matrix(arrange(as.data.frame(ic_int),
                             const, trend, p))
adq_idx <- match(data.frame(t(adq_set[, 1:3])),
                 data.frame(t(ic[, 1:3])))

for (i in 1:length(adq_idx))
{
  checkresiduals(ADF_est[[adq_idx[i]]])
}

# ADF test.
adf.test(r, nlag = 10)

#Repeat for diff(r)
TT <- length(diff(r))
ADF_est_diff <- list()
ic_diff <- matrix( nrow = 30, ncol = 5 )
colnames(ic_diff) <- c("cons", "trend", "p", "aic", "bic")
i <- 0
for (const in 0:1)
{
  for (p in 0:9)
  {
    i <- i + 1
    ADF_est_diff[[i]] <- Arima(diff(diff(r)),
                               xreg = diff(r)[-TT],
                               order = c(p, 0, 0),
                               include.mean = as.logical(const),
                               include.drift = F)
    ic_diff[i,] <- c(const, 0, p, ADF_est_diff[[i]]$aic,
                     ADF_est_diff[[i]]$bic)
  }
  
  if (const)
  {
    # only add a specification with trend if there is a
    # constant (i.e., exclude no constant with trend)
    for (p in 0:9)
    {
      i <- i + 1
      ADF_est_diff[[i]] <- Arima(diff(diff(r)),
                                 xreg = diff(r)[-TT],
                                 order = c(p, 0, 0),
                                 include.mean = as.logical(const),
                                 include.drift = T)
      ic_diff[i,] <- c(const, 1, p, ADF_est_diff[[i]]$aic,
                       ADF_est_diff[[i]]$bic)
    }
  }
}

ic_aic_diff <- ic_diff[order(ic_diff[,4]),][1:10,]
ic_bic_diff <- ic_diff[order(ic_diff[,5]),][1:10,]

# find the intersection of AIC and BIC preferred sets
ic_int_diff <- intersect(as.data.frame(ic_aic_diff),
                         as.data.frame(ic_bic_diff))

# Capture the selected set
adq_set_diff <- as.matrix(arrange(
                          as.data.frame(ic_int_diff),
                                       const, trend, p))
adq_idx_diff <- match(data.frame(t(adq_set_diff[, 1:3])),
                      data.frame(t(ic_diff[, 1:3])))

# Check residuals
for (i in 1:length(adq_idx_diff))
{
  checkresiduals(ADF_est_diff[[adq_idx_diff[i]]])
}

# ADF test

adf.test(diff(r), nlag = 10)
