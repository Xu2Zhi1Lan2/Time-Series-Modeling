#################################################################################
### Part 1: descriptive statistics
### prepare data
copper <- scan(file="copper1800t1996.txt")  # load data 
year <- 1800:1996  # years

par(mfcol=c(2,1))
plot(year,copper,xlab='year',ylab='price',type='l')
acf(copper) # acf decreases under threshold after finite years (stationary)

### use loess to check whether there's regime change in lag 1 relation (graphically)
y <- copper[2:197]
x <- copper[1:196]
m1 <- loess(y~x)  ## local smoothing
sx <- sort(x,index=T)
par(mfcol=c(1,1))
ix <- sx$ix
plot(x,y,xlab='x(t-1)',ylab='x(t)')
lines(x[ix],m1$fitted[ix],col="red") # indicates two regimes

### Chow test
library(strucchange)
sctest(y~x, type = "Chow") # no significant structural change

#################################################################################
### Part 2: fit ARIMA
library(forecast)
copper_arima = auto.arima(copper, 
                          trace = T, 
                          seasonal= T, 
                          stepwise=F,    
                          approximation=F) # ARMA(3,2) is preferred

summary(copper_arima) # AIC = 465.48, BIC = 485.18
autoplot(copper_arima) # check if stationary
checkresiduals(copper_arima) # check linearity

#################################################################################
### Part 3: fit TAR
library(TSA)

m3 <- tar(copper,3,3,1,method = "MAIC", # use AIC
          estimate.thd = T,
          order.select = T) # d and r should be estimated by min AIC or other objective func

m3$thd 

m3$AIC # 427, even if there's no significant structural change, 
       # TAR provides a model with stronger predictability

m3$p1
m3$qr1$coefficients

m3$p2
m3$qr2$coefficients




