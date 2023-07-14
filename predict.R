### Goal: back-test for predictability and predict the future
library(glmnet)
library(forecast)
data(AirPassengers)
AP <- AirPassengers
plot(AP, ylab = "Passengers (1000's)") # time series plot

###########################################################
### Review: Harmonic Seasonal Model and Seasonal ARIMA
#############################################################

### Harmonic Seasonal Model
### use Fourier predictor: sin(2*pi*time/12), sin(4*pi*time/12),.....
time = 1:length(AP)
sin.var = matrix(0,length(AP),6)
cos.var = matrix(0,length(AP),6)
for (k in 1:6) {
  sin.var[,k] = sin(2*pi*k*time/12)
  cos.var[,k] = cos(2*pi*k*time/12)
}
HR2 <- lm(log(AP) ~ time+sin.var+cos.var) # multiplicative seasonal effect
summary(HR2)
plot(exp(predict(HR2)),type = "l",xlab = "",ylab = "",col = "red") 

### seasonal ARIMA
AP.SARIMA <- arima(AP, order = c(0,1,1), seasonal = list(order = c(0,1,1),period = 12))
plot(AP + AP.SARIMA$residuals,type = "l") # fitting


################################################################
### Variable Selection in Harmonic Seasonal Model
#################################################################

### LASSO and Elastic Net
Y = AP[2:144]
X = cbind(time[2:144],sin.var[2:144,],cos.var[2:144,])
H1pred = cv.glmnet(X,log(Y),family = "gaussian",nfolds=5,type.measure="mse")
coef(H1pred, s = 0.1)
Y.fit = exp(predict(H1pred, newx = X, s="lambda.1se", family="gaussian", type="response"))
plot(Y.fit~time[2:144],type = "l")

X0 = cbind(time[2:144],sin.var[2:144,],cos.var[2:144,])
H0pred = cv.glmnet(X0,log(Y),family = "gaussian",nfolds=5,type.measure="mse")
coef(H0pred, s = 0.1)
Y.fit0 = exp(predict(H0pred, newx = X0, s="lambda.1se", family="gaussian", type="response"))
plot(Y.fit0~time[2:144],type = "l")
# drawbacks: no standard error for CIs and classical Bootstrap fails to provide consistent CI

### Forward Stepwise Selection
Y = AP[2:144]
intercept_only <- lm(Y ~ 1)
all <- lm(Y ~ time[2:144] + AP[1:143] + sin.var[2:144,] + cos.var[2:144,])
forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)
forward$anova # tends to overfitting
Y.pred <- Y + forward$residuals
plot(Y.pred ~ time[2:144],type = "l")

### Naive Selection (low efficient, and only works for )
summary(all)
all1 <- lm(log(Y) ~ time[2:144] + sin.var[2:144,1:5] + cos.var[2:144,1:5])
Y.pred <- exp(predict(all1))
plot(Y.pred ~ time[2:144],type = "l")

##################################################################
### Back-test
#################################################################

### Method 1: recursive framework
### 1949-1955 (84 months) as in-sample and 1956-1960 (60 months) as out-of-sample
### consider one-period prediction
Y = as.vector(AP)

Time = 1:144
In_sample = 84
Out_of_sample = 60

sin1 = sin(2*pi*Time/12)
sin2 = sin(4*pi*Time/12)
sin3 = sin(6*pi*Time/12)
sin4 = sin(8*pi*Time/12)
sin5 = sin(10*pi*Time/12)
sin6 = sin(12*pi*Time/12)

cos1 = cos(2*pi*Time/12)
cos2 = cos(4*pi*Time/12)
cos3 = cos(6*pi*Time/12)
cos4 = cos(8*pi*Time/12)
cos5 = cos(10*pi*Time/12)
cos6 = cos(12*pi*Time/12)

X = data.frame(Time = Time,
               sin1 = sin1,
               sin2 = sin2,
               sin3 = sin3,
               sin4 = sin4,
               sin5 = sin5,
               sin6 = sin6,
               cos1 = cos1,
               cos2 = cos2,
               cos3 = cos3,
               cos4 = cos4,
               cos5 = cos5,
               cos6 = cos6)

X = cbind(Time,sin.var,cos.var)

resid.HSM = rep(0,Out_of_sample)
resid.ARIMA = rep(0,Out_of_sample)

for (t in 1:Out_of_sample) {
  ## Harmonic seasonal model
  y = Y[1:(In_sample+t-1)]
  train = data.frame(X[1:(In_sample+t-1),], y = y)
  newx = X[(In_sample+t),]
  HSM = lm(log(y)~Time+sin1+sin2+sin3+sin4+sin5+sin6+cos1+cos2+cos3+cos4+cos5+cos6,data = train)
  resid.HSM[t] = exp(predict(HSM, newx, se.fit = T)$fit) - Y[In_sample+t]
  
  ## Seasonal ARIMA
  AP.SARIMA <- arima(Y[1:(In_sample+t-1)], order = c(0,1,1), seasonal = list(order = c(0,1,1),period = 12))
  resid.ARIMA[t] = as.numeric(predict(AP.SARIMA, n.ahead = 1, newxreg = NULL,se.fit = TRUE)$pred - Y[In_sample+t])
}

dm.test(resid.HSM,resid.ARIMA,alternative = "greater")

### Method 2: rolling window framework
### window size = 7 years = 84 months

resid.HSM2 = rep(0,Out_of_sample)
resid.ARIMA2 = rep(0,Out_of_sample)

for (t in 1:Out_of_sample) {
  ## Harmonic seasonal model
  y = Y[t:(In_sample+t-1)]
  train = data.frame(X[t:(In_sample+t-1),], y = y)
  newx = X[(In_sample+t),]
  HSM = lm(log(y)~Time+sin1+sin2+sin3+sin4+sin5+sin6+cos1+cos2+cos3+cos4+cos5+cos6,data = train)
  resid.HSM2[t] = exp(predict(HSM, newx, se.fit = T)$fit) - Y[In_sample+t]
  
  ## Seasonal ARIMA
  AP.SARIMA <- arima(Y[t:(In_sample+t-1)], order = c(0,1,1), seasonal = list(order = c(0,1,1),period = 12))
  resid.ARIMA2[t] = as.numeric(predict(AP.SARIMA, n.ahead = 1, newxreg = NULL,se.fit = TRUE)$pred - Y[In_sample+t])
}

dm.test(resid.HSM2,resid.ARIMA2,alternative = "greater")

##############################################################
### prediction with ARIMA
AP.SARIMA <- arima(Y, order = c(0,1,1), seasonal = list(order = c(0,1,1),period = 12))
AP.pred <- predict(AP.SARIMA, n.ahead = 12, newxreg = NULL,se.fit = TRUE)
Pred = 1:156
YPred = rep(0,156); YPred[1:144] = Y; YPred[145:156] = AP.pred$pred
LowerCI = AP.pred$pred - 1.96*AP.pred$se
UpperCI = AP.pred$pred + 1.96*AP.pred$se
plot(YPred~Pred, type = "l",col = "red",xlim = c(0,160),ylim = c(0,700),xlab = "Time",ylab = "Airpassengers")
par(new = T)
plot(LowerCI,col = "blue",lty = 2,xlim = c(0,160),ylim = c(0,700),xlab = "",ylab = "")
par(new = T)
plot(UpperCI,col = "blue",lty = 2,xlim = c(0,160),ylim = c(0,700),xlab = "",ylab = "")
abline(v=145)

