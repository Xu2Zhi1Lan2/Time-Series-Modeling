###############################################################
### Example 1: Engle-Granger two step method
set.seed (2023) 
e1 <- rnorm (100) 
e2 <- rnorm (100) 
y1 <- cumsum(e1) # I(1) data
y2 <- 0.575*y1 + e2 # I(1) data, linearly related to y1

## step 1: cointegrating regression
model1 <- lm(y2 ~ -1+y1) 
summary(model1) # R^2 = 0.9619, hat alpha = (1, -0.59322)

# check stationarity of error term
library(fUnitRoots)
error <- residuals(model1) 
adfTest(error,type = "nc")
unitrootTest(error, type = "nc") # no unit root

# check spurious regression and serial correlation
library(lmtest)
dwtest(model1) # CRDW = 2.1452 > R^2, not spurious regression

## step 2: error-correction model
error.lagged <- error[-c(99,100)] # first order diff and first order lag, so minus 2 degree of freedom
dy1 <- diff(y1) 
dy2 <- diff(y2) 
diff.dat <- data.frame(embed(cbind(dy1,dy2),2)) # x1 is the lag of dy1, x2 is the lag of dy2
colnames(diff.dat) <- c('dy1.1','dy2.1','dy1','dy2') 
ecm <- lm(dy2 ~ error.lagged + dy1.1 + dy2.1, data= diff.dat)
summary(ecm)

###############################################################
### Example 2: VECM
library (urca) 
set.seed (2023) 
e1 <- rnorm (250,0,0.5) 
e2 <- rnorm (250,0,0.5) 
e3 <- rnorm (250,0,0.5) 
u1.ar1 <- arima.sim(model = list(ar= 0.75), 
                    innov = e1 , n = 250) # simulate an AR(1) error
u2.ar1 <- arima.sim(model = list(ar = 0.3), 
                    innov = e2 ,n = 250) 
y3 <- cumsum(e3) # I(1) data
y1 <- 0.8*y3 + u1.ar1 # I(1) data
y2 <- -0.2*y3 + u2.ar1 # I(1) data
y.mat <- data.frame(y1,y2,y3) 

# fit VECM
vecm <- ca.jo(y.mat)  # fit VECM with Johansen method
summary(vecm) # fail to reject r<=2, reject r<=1, so r = 2
              # also estimation of beta and alpha
vecm.r2 <- cajorls(vecm , r = 2) 
summary(vecm.r2$rlm)

# recover VECM as VAR 
library(vars)
vecm.level <- vec2var(vecm, r = 2) # ca.jo object to VAR object
arch.test(vecm.level) 
normality.test(vecm.level) 
serial.test(vecm.level) 
predict(vecm.level) 
irf(vecm.level, boot = T) 

