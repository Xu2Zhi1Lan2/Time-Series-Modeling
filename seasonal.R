### Goal: Investigate the seasonal effect in Time series data 
data(AirPassengers)
AP <- AirPassengers
plot(AP, ylab = "Passengers (1000's)") # time series plot

### Decomposition
season = cycle(AP) # cycle (12 months in a year)
Time = seq(1,length(AP))

AP.lm <- lm(AP ~ Time + factor(season) - 1) # additive seasonal effect
AP.fitted = predict(AP.lm)
plot(AP.fitted,type = "l")

AP.lm2 <- lm(log(AP) ~ Time + factor(season) - 1) # multiplicative seasonal effect
AP.fitted2 = exp(predict(AP.lm2))
plot(AP.fitted2,type = "l")

### Harmonic Seasonal Model
### use Fourier predictor: sin(2*pi*time/12), sin(4*pi*time/12),.....
time = 1:length(AP)
sin.var = matrix(0,length(AP),6)
cos.var = matrix(0,length(AP),6)
for (k in 1:6) {
  sin.var[,k] = sin(2*pi*k*time/12)
  cos.var[,k] = cos(2*pi*k*time/12)
}

HR1 <- lm(AP ~ time+sin.var+cos.var) # additive seasonal effect
summary(HR1)
plot(predict(HR1),type = "l",xlab = "",ylab = "",col = "red") 

HR2 <- lm(log(AP) ~ time+sin.var+cos.var) # multiplicative seasonal effect
summary(HR2)
plot(exp(predict(HR2)),type = "l",xlab = "",ylab = "",col = "red") 

### seasonal difference
AP.sdiff = diff(AP,lag = 12)
plot(AP.sdiff,type = "l")
acf(AP.sdiff, lag.max = 48)

AP.slogdiff = diff(log(AP),lag = 12)
plot(AP.slogdiff,type = "l")
acf(AP.slogdiff, lag.max = 48)

AP.sdiffar <- ar(AP.sdiff, aic = T, method = "ols")
plot(AP.sdiffar$aic)

AP.slogdiffar <- ar(AP.slogdiff, aic = T, method = "ols")
plot(AP.slogdiffar$aic)

### seasonal ARIMA
AP.SARIMA <- arima(AP, order = c(0,1,1), seasonal = list(order = c(0,1,1),period = 12))
plot(AP + AP.SARIMA$residuals,type = "l") # fitting

### X-12-ARIMA (read paper "Seasonal Adjustment with the R Packages x12 and x12GUI" 2014 Journal of Statistical Software for detailed tutorial)
library(x12)
s <- new("x12Single", ts = AirPassengers, tsName = "air")
s <- x12(s)
plot(s, trend = TRUE, sa = TRUE, forecast = TRUE)
plotRsdAcf(s, which = "acf2")
