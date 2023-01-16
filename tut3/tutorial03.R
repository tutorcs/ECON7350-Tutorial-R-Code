# Clear objects from environment
rm(list = ls())

# Install and load package "forecast"
install.packages("forecast")
library(forecast)

# 1
# (a)
# load the data in Merck.csv
mydata <- read.delim("Merck.csv", header = TRUE,  sep = ",")
sel_sample <- mydata$Date >= as.Date("2011-01-01") &
              mydata$Date <= as.Date("2012-01-31")
y <- as.matrix(mydata$Adj_Close[sel_sample])

# (b)
# Construct variables:
# Change in price (first difference of y)
Dy <- diff(y)
# Log returns of y
r <- as.matrix(log(y[2:nrow(y)]) - log(lag(y[1:nrow(y) - 1])))

# Change variable names
colnames(y) <- "Stock Prices of Merck (MRK)"
colnames(Dy) <- "Changes in Stock Prices"
colnames(r) <- "Log Returns"

# (c)
# Time series plot of y and Dy
sel2011 <- mydata$Date[sel_sample] <= as.Date("2011-12-31")
dates <- as.Date(mydata$Date[sel_sample])
plot(dates[sel2011], y[sel2011], type = "l", xlab = "Time (2011)",
     ylab = colnames(y))
plot(dates[sel2011], Dy[sel2011], type = "l", xlab = "Time (2011)",
     ylab = colnames(Dy))

# (d)
# ACF and PACF of y
acf(y[sel2011], main = colnames(y))
pacf(y[sel2011], main = colnames(y))
# ACF and PACF of Dy
acf(Dy[sel2011], main = colnames(Dy))
pacf(Dy[sel2011], main = colnames(Dy))

# (e)
# Create an empty matrix to store the results
ic <- matrix(nrow = 25, ncol = 4)
colnames(ic) <- c("p", "q", "aic", "bic")
# Add AIC and BIC information for the first model (ARMA(0,0), to the first row)
fit_p_q <- Arima(Dy, order = c(0, 0, 0))
ic[1,] = c(0, 0, fit_p_q[["aic"]], fit_p_q[["bic"]])

# Simplify the process with a loop
for (p in 0:4){
  for (q in 0:4){
    fit_p_q <- Arima(Dy, order = c(p, 0, q))
    ic[p * 5 + q + 1,] = c(p, q, fit_p_q[["aic"]], fit_p_q[["bic"]])
  }
}

# (f)
# Adding top 10 specifications preferred by AIC and BIC
ic_aic <- ic[order(ic[,3], decreasing = FALSE),][1:10,]
ic_bic <- ic[order(ic[,4], decreasing = FALSE),][1:10,]

# Select the models that make top 10 in both lists
adq_set = list(c(1, 0, 1), c(1, 0, 2), c(2, 0, 1), c(3, 0, 0))

# (g)
# Perform residual analysis
checkresiduals(Arima(Dy[sel2011], order = c(1, 0, 1)))

# Simplify with a loop
for (i in 1:length(adq_set)){
  checkresiduals(Arima(Dy[sel2011], order = adq_set[[i]]))
}

# (h)
# Forecast MRK in January 2012
# Count how many days in the test sample (forecast horizon)
hrz = sum(sel_sample) - sum(sel2011)
# Set three days as reference points
xticks <- c(sum(sel_sample) - 3 * hrz + c(1, 2 * hrz, 3 * hrz))
# Subset by observations in the test set
actual_Dy <- as.matrix(Dy[!sel2011])
# Create an empty vector to store forecasts
fcst_Dy <- vector(mode = "list", length(adq_set))
# Generate forecasts with a loop
for (i in 1:length(adq_set)){
  model_p_q <- adq_set[[i]]
  fcst_Dy[[i]] <- forecast(Arima(Dy[sel2011], order = model_p_q),
                           h = hrz, level = c(68, 95))

  title_p_q <- paste("ARMA(", as.character(model_p_q[1]), ", ",
                              as.character(model_p_q[3]), ")", sep = "")
  
  plot(fcst_Dy[[i]], include = hrz * 2, ylab = colnames(Dy),
       main = title_p_q, xaxt = "n")
  lines(sum(sel2011) + 1:hrz, actual_Dy)
  axis(1, at = xticks, labels = dates[xticks])
}

# (i)
# Repeat the forecasting steps with the "y" variable using an ARMA(2,1)
actual_y <- as.matrix(y[!sel2011])
fcst_y_lev = forecast(Arima(y[sel2011], order = c(2, 0, 1)),
                      h = hrz, level = c(68, 95))

plot(fcst_y_lev, include = hrz * 2, ylab = colnames(y),
     main = "ARMA(2, 1)", xaxt = "n", ylim = c(26.1, 33.4))
lines(sum(sel2011) + 1:hrz, actual_y)
axis(1, at = xticks, labels = dates[xticks])

# Now repeat (g) using ARIMA on "y" instead of ARMA on "Dy"
# Capture the value of Adj_Close on 2010-12-31 to be our y0
y0 <- mydata$Adj_Close[sum(mydata$Date < as.Date("2011-01-01") - 1)]
# Adding it to the training set
y_ext = as.matrix(c(y0, y[sel2011]))
# Produce forecasts
fcst_y <- vector(mode = "list", length(adq_set))
for (i in 1:length(adq_set)){
  model_p_q <- adq_set[[i]]
  model_p_q[2] = 1
  fcst_y[[i]] <- forecast(Arima(y_ext, order = model_p_q, include.constant = T),
                          h = hrz, level = c(68, 95))
  
  title_p_q <- paste("ARIMA(", as.character(model_p_q[1]), ", ", 
                     as.character(model_p_q[2]), ", ",
                     as.character(model_p_q[3]), ")", sep = "")
  
  plot(fcst_y[[i]], include = hrz * 2, ylab = colnames(y),
       main = title_p_q, xaxt = "n")
  lines(1 + sum(sel2011) + 1:hrz, actual_y)
  axis(1, at = 1 + xticks, labels = dates[xticks])
}

# (j)
# Repeat the process for r
# Time series plot of r
plot(dates[sel2011], r[sel2011], type = "l", xlab = "Time (2011)",
     ylab = colnames(r))

# ACF and PACF of r
acf(r[sel2011], main = colnames(y))
pacf(r[sel2011], main = colnames(y))

# Model selection using AIC/BIC
ic <- matrix(nrow = 25, ncol = 4)
colnames(ic) <- c("p", "q", "aic", "bic")
for (p in 0:4)
{
  for (q in 0:4)
  {
    fit_p_q <- Arima(r, order = c(p, 0, q))
    c(p * 5 + q + 1, p, q)
    ic[p * 5 + q + 1,] = c(p, q, fit_p_q[["aic"]], fit_p_q[["bic"]])
  }
}

ic_aic <- ic[order(ic[,3], decreasing = FALSE),][1:10,]
ic_bic <- ic[order(ic[,4], decreasing = FALSE),][1:10,]

adq_set = list(c(1, 0, 1), c(1, 0, 2), c(2, 0, 1), c(3, 0, 0))

# Check the residuals
for (i in 1:length(adq_set))
{
  checkresiduals(Arima(r[sel2011], order = adq_set[[i]]))
}

# Produce the forecasts
hrz <- sum(sel_sample) - sum(sel2011)
xticks <- c(sum(sel_sample) - 3 * hrz + c(1, 2 * hrz, 3 * hrz))
actual_r <- as.matrix(r[!sel2011])
fcst_r <- vector(mode = "list", length(adq_set))
for (i in 1:length(adq_set))
{
  model_p_q <- adq_set[[i]]
  fcst_r[[i]] <- forecast(Arima(r[sel2011], order = model_p_q),
                          h = hrz, level = c(68, 95))
  
  title_p_q <- paste("ARMA(", as.character(model_p_q[1]), ", ",
                     as.character(model_p_q[3]), ")", sep = "")
  
  plot(fcst_r[[i]], include = hrz * 2, ylab = colnames(r),
       main = title_p_q, xaxt = "n")
  lines(sum(sel2011) + 1:hrz, actual_r)
  axis(1, at = xticks, labels = dates[xticks])
}