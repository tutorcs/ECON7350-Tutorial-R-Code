#Problem 1
#(a)
mydata <- read.delim("consumption.txt", header = TRUE, sep = "")
#(b) this code is used to generate graphs
plot(
  mydata$INC,
  mydata$CONS,
  main = "Consumption Data",
  xlab = "Income",
  ylab = "Consumption",
  pch = 17,
  col = "lightblue3"
)
#(c)
mydata[8, 1] <- 90
detach(mydata)
attach(mydata)
plot(
  INC,
  CONS,
  main = "Consumption Data",
  xlab = "Income",
  ylab = "Consumption",
  pch = 17,
  col = "cadetblue4"
)
#(d)
summary(mydata)
#(e)
cor(mydata)
cor(mydata$INC, mydata$CONS)
#(f)
DCONS <- 0.5 * CONS
LCONS <- log(CONS)
INC2 <- INC ^ 2
SQRTINC <- sqrt(INC)
#To attach to a data frame
mydata$DCONS <- 0.5 * CONS
mydata$LCONS <- log(CONS)
mydata$INC2 <- INC ^ 2
mydata$SQRTINC <- sqrt(INC)
#(G)
#Deleting variables from the environment directly
rm(DCONS, SQRTINC)
#Deleting variables from the data frame
mydata <- subset(mydata, select = -c(DCONS, SQRTINC))
mydata <- subset(mydata, select = c(CONS, INC, INC2, LCONS))
#(h)
rm(list = ls())

#Problem 2
#(a)
fultonfish <- read.delim("fultonfish.dat", header = F, sep = "")
colnames(fultonfish)[1:4] <- c("date", "lprice", "quan", "lquan")
#(b)
mean(fultonfish$quan)
sd(fultonfish$quan)
#(c) and (d)
t.test(fultonfish$quan, mu = 7200)
#we reject H0 if P-value < 0.05
#(e)
attach(fultonfish)
plot(lquan, lprice)
plot(
  lquan,
  lprice,
  main = "Log Price and Log Quantity",
  xlab = "log(Quantity)",
  ylab = "log(Price) of Whiting per pound",
  pch = 16,
  col = "darkred"
)
#to tidy up your code, highlight and then ctrl+shift+a
cor(lquan, lprice)
#(f)
save(list = ls(all = TRUE), file = "tutorial01.Rdata")
