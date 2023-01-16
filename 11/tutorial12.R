library(dplyr)
library(zoo)
library(vars)
library(urca)

mydata <- read.delim("term_structure.csv", header = TRUE,  sep = ",")

dates <- as.yearmon(mydata$obs, format = "%YM%m")
T <- length(dates)

i90d <- mydata$I90D[-T]
i180d <- mydata$I180D[-T]
i3y <- mydata$I3Y[-T]
i5y <- mydata$I5Y[-T]
x <- cbind(i90d, i180d, i3y, i5y)
n <- ncol(x)

# 1.
VAR_est <- list()
ic_var <- matrix(nrow = 20,  ncol = 3)
colnames(ic_var) <- c("p", "aic", "bic")
for (p in 1:20){
  VAR_est[[p]] <- VAR(x, p)
  ic_var[p,] <- c(p, AIC(VAR_est[[p]]),
                     BIC(VAR_est[[p]]))
}

ic_aic_var <- ic_var[order(ic_var[,2]),]
ic_bic_var <- ic_var[order(ic_var[,3]),]

adq_set_var <- as.matrix(ic_var[2:5,])
adq_idx_var <- c(2:5)

nmods <- length(adq_idx_var)
for (i in 1:nmods){
  p <- adq_idx_var[i]
  print(paste0("Checking VAR(", p, ")"))
  print(serial.test(VAR_est[[p]], lags.bg = 1,
                                  type = "BG"))
}

adq_set_var <- as.matrix(ic_var[c(3:5),])
adq_idx_var <- c(3:5)

# 2.
nmods <- length(adq_idx_var)
for (i in 1:nmods){
  p <- adq_idx_var[i]
  print(paste0("VAR(", p,
        "): Maximum absolute eigenvalue is ",
        max(vars::roots(VAR_est[[p]]))))
}

nmods <- length(adq_idx_var)
for (i in 1:nmods){
  p <- adq_idx_var[i]
  print(paste0("VAR(", p, ")"))
  print(summary(ca.jo(x, type = "trace", K = p)))
}

# 3.
VECM_est <- list()
ic_vecm <- matrix(nrow = 4 * (1 + n),  ncol = 4)
colnames(ic_vecm) <- c("p", "r", "aic", "bic")
i <- 0
for (p in 3:6){
  for (r in 0:n){
    i <- i + 1
    if (r == n){
      VECM_est[[i]] <- VAR(x, p)
    }
    else if (r == 0){
      VECM_est[[i]] <- VAR(diff(x), p - 1)
    }
    else {
      VECM_est[[i]] <- vec2var(ca.jo(x, K = p), r)
    }
    ic_vecm[i,] <- c(p, r, AIC(VECM_est[[i]]),
                           BIC(VECM_est[[i]]))
  }
}

ic_aic_vecm <- ic_vecm[order(ic_vecm[,3]),][1:10,]
ic_bic_vecm <- ic_vecm[order(ic_vecm[,4]),][1:10,]

ic_int_vecm <- intersect(as.data.frame(ic_aic_vecm),
                         as.data.frame(ic_bic_vecm))

adq_set_vecm <- as.matrix(arrange(as.data.frame(
                                  ic_int_vecm), p, r))
adq_idx_vecm <- match(data.frame(t(adq_set_vecm[, 1:2])),
                           data.frame(t(ic_vecm[, 1:2])))

nmods <- length(adq_idx_vecm)
for (i in 1:nmods){
  p <- adq_set_vecm[i, 1]
  r <- adq_set_vecm[i, 2]
  print(paste0("Checking VECM(", p, ", ", r, ")"))
  print(serial.test(VECM_est[[adq_idx_vecm[i]]],
                    lags.bg = 1,
                    type = "BG"))
}

# 4.
p <- 3
for (r in 2:4){
  if (r < n){
    vec_pr <- ca.jo(x, type = "trace", spec = "transitory")
    beta_pr <- cajorls(vec_pr, r)$beta
    
    bpvals_pr <- beta_pr
    bpvals_pr[,] <- NA
    
    for (i in (r + 1):n){
      for (j in 1:r){
        H <- beta_pr
        H[i, j] <- 0
        bpvals_pr[i, j] <- blrtest(vec_pr, H, r)@pval[1]
      }
    }
    
    print(paste0("Results for VECM(", p, ", ", r, "):"))
    print(beta_pr)
    print(bpvals_pr)
  }
}

# 5.
hrz <- 60
vnames <- c("i90d", "i180d", "i3y", "i5y")
nmods <- length(adq_idx_vecm)
for (i in 1:n){
  for (j in 1:n){
    for (imod in 1:nmods){
      p <- adq_set_vecm[imod, 1]
      r <- adq_set_vecm[imod, 2]
      title_i_j <- paste0("VECM(", p, ", ", r,
                          "): Response of ", vnames[i],
                          " to a shock in ", vnames[j])
  
      irf_i_j <- irf(VECM_est[[adq_idx_vecm[imod]]],
                            n.ahead = hrz,
                            response = vnames[i],
                            impulse = vnames[j],
                            boot = TRUE)
          
      plot(irf_i_j, main = title_i_j)
      cat("\r", title_i_j, "  ", sep = "")
    }
  }
}