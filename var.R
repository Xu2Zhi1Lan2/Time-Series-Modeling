### Example: estimate money demand function
### Theory: 
### 1.交易性需求（Baumol模型）：price level，income level，interest rate
### 2.预防性需求（Whalen模型）
### 3.投机性需求（Tobin模型）
### Data: Johansen and Juselius (1990) Denmark 1974Q1 – 1987Q3, logarithm of real money M2
### logarithm of real income, logarithm of price deflator, bond rate, bank deposit rate
# LRM : logarithm of real money M2 (LRM)
# LRY : logarithm of real income (LRY)
# LPY : logarithm of price deflator (LPY)
# IBO : bond rate (IBO)
# IDE : bank deposit rate (IDE)
# the period 1974:Q1 – 1987:Q3

library(urca)
library(vars) 
library(lmtest)

############################################################################
### data preperation
data(denmark)
df.lev <- denmark[,c("LRM","LRY","LPY","IDE")]
m.lev  <- as.matrix(df.lev)
nr_lev <- nrow(df.lev)

### plot of time series data: LPY is definitely integrated data
str.main <- c(
  "LRM=ln(real money M2)", "LRY=ln(real income)", 
  "LPY=ln(price deflator)", "IDE=bank deposit rate")

x11(width=12, height = 6); 
par(mfrow=c(2,2), mar=c(5,3,3,3))
for(i in 1:4) {
  matplot(m.lev[,i], axes=FALSE,
          type="l", col = "blue", 
          main = str.main[i])
  
  axis(2) # show y axis
  
  # show x axis and replace it with 
  # an user defined sting vector
  axis(1, at=seq_along(1:nrow(df.lev)),
       labels=denmark$ENTRY, las=2)
}

### 1-st order difference the data
df.diff <- diff(as.matrix(df.lev), lag = 1)
colnames(df.diff) <- c("dLRM", "dLRY", "dLPY", "dIDE")
m.diff <- as.matrix(df.diff)

############################################################################
### Model specification and estimation
VARselect(df.diff, lag.max = 4,
          type = "const", season = 4)

vare_diff <- VAR(df.diff, p = 1, 
                 type = "const", season = 4)

summary(vare_diff)

### Model diagnostics
serial.test(vare_diff, lags.pt = 4, type = "PT.asymptotic") #H0: no serial correlation
arch.test(vare_diff, lags.multi = 4, multivariate.only = TRUE) #H0: no heteroskadasticity
normality.test(vare_diff, multivariate.only = TRUE) #H0: normal
plot(stability(vare_diff, type = "OLS-CUSUM")) 

############################################################################
### Prediction
varf_diff <- predict(vare_diff, n.ahead = 8)
x11(); par(mai=rep(0.4,4)); plot(varf_diff)
x11(); par(mai=rep(0.4,4)); fanchart(varf_diff)

### recover level data prediction
m.varf_lev_ft <- rbind(m.lev, matrix(NA, 8, 4))
m.ft_df <- do.call(cbind,lapply(varf_diff$fcst, 
                                function(x) x[,"fcst"]))

for(h in (nr_lev+1):(nr_lev+8)) {
  hf <- h - nr_lev
  m.varf_lev_ft[h,] <- m.varf_lev_ft[h-1,] + m.ft_df[hf,]
}

x11(width=8, height = 8); 
par(mfrow=c(4,1), mar=c(2,2,2,2))

for(i in 1:4) {
  df <- m.varf_lev_ft[,i]
  matplot(df, type="l", col = "blue", 
          main = str.main[i]) 
  abline(v=nr_lev, col="blue")
}

############################################################################
### Granger Causality, H0: Time series X does not Granger-cause time series Y
grangertest(LRM ~ LRY, order = 1, data = denmark) # p = 0.1379
grangertest(LRM ~ LPY, order = 1, data = denmark) # p = 0.1465

grangertest(LRM ~ LRY, order = 4, data = denmark) # p = 0.5091
grangertest(LRM ~ LPY, order = 4, data = denmark) # p = 0.01426

causality(vare_diff, cause = "dLPY")
causality(vare_diff, cause = "dLRM")
causality(vare_diff, cause = "dLRY")
causality(vare_diff, cause = "dIDE")

############################################################################
### Impulse Response Analysis
LPYirf <- irf(vare_diff, impulse = "dLPY", response = "dLRM", n.ahead = 20, boot = TRUE)
plot(LPYirf, ylab = "dLRM", main = "LPY's shock to LRM") # no significant response

LRMirf <- irf(vare_diff, impulse = "dLRM", response = "dLPY", n.ahead = 20, boot = TRUE)
plot(LRMirf, ylab = "dLPY", main = "LRM's shock to LPY") # negative response

LRYirf <- irf(vare_diff, impulse = "dLRM", response = "dLRY", n.ahead = 20, boot = TRUE)
plot(LRYirf, ylab = "dLRY", main = "LRM's shock to LRY") # positive response

############################################################################
### Counterfactual simulation
library(svars)

View(cf)

v1 <- vars::VAR(USA, lag.max = 10, ic = "AIC" )
x1 <- id.dc(v1)
x2 <- cf(x1, series = 2)
plot(x2)



