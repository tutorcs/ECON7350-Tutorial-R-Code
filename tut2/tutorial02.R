# Clean up the environment
rm(list = ls())

# Load the data
mydata <- read.delim("arma.csv", header = TRUE,  sep = ",")

# 3
# Generating acf and pacf graphs
# DGP1
acf(mydata$DGP1, main = colnames(mydata[2]))
pacf(mydata$DGP1, main = colnames(mydata[2]))

# DGP2 onwards
# Saving time with loops
for (i in 2:8){
  acf(mydata[1 + i], main = colnames(mydata[1 + i]))
  pacf(mydata[1 + i], main = colnames(mydata[1 + i]))
}