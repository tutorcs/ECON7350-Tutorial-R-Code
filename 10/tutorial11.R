# Load necessary packages
library(zoo)
library(vars)
library(pracma)

# load the data
mydata <- read.delim("money_dem.csv", header = TRUE,  sep = ",")

date <- as.yearqtr(mydata$DATE)
lrgdp <- log(mydata$RGDP)
price <- mydata$GDP/mydata$RGDP
lrm2 <- log(mydata$M2) - log(price)
rs <- mydata$TB3mo
x <- cbind(lrgdp, lrm2, rs)

# 1
# (a)
VAR_est <- list()
ic_var <- matrix(nrow = 20, ncol = 3)
colnames(ic_var) <- c("p", "aic", "bic")
for (p in 1:20){
  VAR_est[[p]] <- VAR(x, p)
  ic_var[p,] <- c(p, AIC(VAR_est[[p]]),
                     BIC(VAR_est[[p]]))
}

ic_aic_var <- ic_var[order(ic_var[,2]),]
ic_bic_var <- ic_var[order(ic_var[,3]),]

# Create adequate set and obtain indices
adq_set_var <- as.matrix(ic_var[2:4,])
adq_idx_var <- c(2:4)

# Check the residuals for adequate set
for (i in 2:4){
  print(paste0("Checking VAR(", i, ")"))
  print(serial.test(VAR_est[[i]], type = "BG"))
}

# Check the residuals for all 20 VAR(p) models
LMtest <- matrix(nrow = 20, ncol = 2)
colnames(LMtest) <- c("p", "p-value")
for (p in 1:20){
  LMtest[p,1] <- p
  LMtest[p,2] <- serial.test(VAR_est[[p]], type = "BG")[["serial"]][["p.value"]]
  print(paste0("Checking VAR(", p, ")"))
  print(serial.test(VAR_est[[p]], type = "BG"))
}

# Proceed with p = 8, 9, 10 as the adequate set.
adq_set_var <- as.matrix(ic_var[8:10,])
adq_idx_var <- c(8:10)

# (b)
# intercept and slope coefficients
nmods <- length(adq_idx_var)
n <- ncol(x)
for (i in 1:nmods){
  p <- adq_idx_var[i]
  print(paste0("VAR(", p, ") has ",
               n * (1 + n * p),
               " coefficients."))
}

# (c)
nmods <- length(adq_idx_var)
for (i in 1:nmods){
  p <- adq_idx_var[i]
  print(paste0("VAR(", p,
        "): Maximum absolute eigenvalue is ",
        max(vars::roots(VAR_est[[p]]))))
}


# (d)
hrz = 12
VAR_fcst <- list()
xlim <- c(length(date) - 3 * hrz,
              length(date) + hrz)
ylim <- c(lrgdp[xlim[1]],
          max(lrgdp) + 0.2)
for (i in 1:nmods){
  p <- adq_idx_var[i]
  VAR_fcst[[i]] <- predict(VAR_est[[p]],
                           n.ahead = hrz)
  plot(VAR_fcst[[i]], names = "lrgdp",
           xlim = xlim, ylim = ylim,
           main = paste0("Forecast for Log Real GDP of VAR(",p,")"),
           xlab = "Horizon",
           ylab = "RRP")
}

# 2
# (a)
orders <- perms(1:3)
vnames <- c("lrgdp", "lrm2", "rs")
for (i in 1:3){
  for (j in 1:3){
    for (k in 1:nrow(orders)){
      title_i_j_k <- paste0("Response of ",
                            vnames[i],
                            " to a shock in ",
                            vnames[j],
                            "; x = (",
                            paste0(vnames[orders[k,]],
                                  collapse = ", "),
                            ")'")

      irf_i_j_k <- irf(VAR(x[,orders[k,]], 8),
                          n.ahead = 40,
                          response = vnames[i],
                          impulse = vnames[j],
                          boot = TRUE)
        
      plot(irf_i_j_k, main = title_i_j_k)
      cat("\r", title_i_j_k, "  ", sep = "")
    }
  }
}

# (b)
FEVD_est <- fevd(VAR_est[[8]], n.ahead = 40)
plot(FEVD_est, mar = c(2,1,2,1),
               oma = c(0,1,0,1))

# (c)
for (i in 1:3){
  ctest_i <- causality(VAR_est[[8]],
                       cause = vnames[i])
  print(ctest_i$Granger)
}