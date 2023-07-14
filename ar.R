### Goal: Time series data pre-analysis and model suggestion
data(AirPassengers)
AP <- AirPassengers
AP
class(AP)
start(AP); end(AP); frequency(AP)
plot(AP, ylab = "Passengers (1000's)") # time series plot
plot(aggregate(AP)) # monthly aggregation
plot(cycle(AP)) # cycle (12 months in a year)
boxplot(AP ~ cycle(AP)) # comparison among different months

####################################################################################
### Pre-analysis I: Trend and Cycle Decomposition
AP.decom <- decompose(AP, type = "additive")
plot(AP.decom$seasonal)
plot(AP.decom$trend)

### Fourier series
### use predictor: 1, time, sin(2*pi*time/12), sin(4*pi*time/12),.....
time = 1:length(AP)
sin.var = matrix(0,length(AP),6)
cos.var = matrix(0,length(AP),6)
for (k in 1:6) {
  sin.var[,k] = sin(2*pi*k*time/12)
  cos.var[,k] = cos(2*pi*k*time/12)
}
HR1 <- lm(AP ~ time+sin.var+cos.var)
summary(HR1)

plot(AP, ylab = "Passengers (1000's)", ylim = c(100,600))
par(new = T)
plot(predict(HR1),type = "l",xlab = "",ylab = "",ylim = c(100,600),col = "red") # failed to capture heteroscedasticity

########################################################################################
### Pre-analysis II: Stationary vs Non-stationary 

### simulated example to illustrate how ACFs of stationary and non-stationary TS differ
t = 0:400
y_stationary <- rnorm(length(t),mean=1,sd=1) 
y_integrated <- cumsum(rnorm(length(t),mean=1,sd=1))+t/400 
y_stationary <- y_stationary/max(y_stationary) 
y_integrated <- y_integrated/max(y_integrated) 
plot.new()
frame()
par(mfcol=c(2,2))
plot(t,y_stationary,
     type='l',col='red',
     xlab = "time (t)",
     ylab = "Y(t)",
     main = "Stationary Series")
acf(y_stationary,lag.max = length(y_stationary),
    xlab = "lag #", ylab = 'ACF',main=' ')
plot(t,y_integrated,
     type='l',col='red',
     xlab = "time (t)",
     ylab = "Y(t)",
     main = "Integrated Series")
acf(y_integrated,lag.max = length(y_integrated),
    xlab = "lag #", ylab = 'ACF', main=' ')

### plot the acf of AP
plot.new()
par(mfcol=c(1,1))
acf(AP,lag.max = length(AP),
    xlab = "lag #", ylab = 'ACF', main=' ') # probably non-stationary time series

### Augmented Dickey–Fuller (ADF) t-test for unit root
library(tseries)
adf.test(AP) # no unit root

### If not stationary, try log, diff or log diff
# log
LAP <- log(AP)
plot(LAP)
# difference
DAP <- diff(AP)
plot(DAP)
# log difference
LDAP <- diff(log(AP))
plot(LDAP)

#####################################################################################
### Step I: Identification

### partial acf plot (threshold +/- 1.96/sqrt(T))
pacf(AP,lag.max = 36)
pacf(AP.decom$trend[!is.na(AP.decom$trend)],lag.max = 36) # p = 1

#####################################################################################
### Step E: Estimation

### Fit an AR(p) model
AP.AR <- ar(AP, aic = T, method = "ols") # use OLS to estimate and AIC to determine the order
plot(AP.AR$aic~seq(0,AP.AR$order.max,1),type = "l")

### Fit an ARMA(p,q) model
AP.ARMA <- arma(AP,order = c(1,1))

### Fit an ARIMA(p,d,q) model
AP.ARIMA <- arima(AP,order = c(1,1,1))

### Fit a seasonal ARIMA(p,d,q) model
AP.SARIMA <- arima(AP, order = c(1,1,1), seasonal = list(order = c(1,1,1),period = 12))

#####################################################################################
### Step D: Diagnostics

### Ljung-Box test: null hypothesis rho = 0
Box.test(AP, lag=36, type="Ljung-Box") # small p-value

### residual plot
plot(AP.AR$resid[!is.na(AP.AR$resid)]) # potential heteroscedasticity
acf(AP.AR$resid[!is.na(AP.AR$resid)]) # actually, residuals have no sequential correlation

plot(AP.ARMA$resid[!is.na(AP.ARMA$resid)]) # obvious heteroscedasticity

plot(AP.ARIMA$resid[!is.na(AP.ARIMA$resid)]) # obvious heteroscedasticity

plot(AP.SARIMA$resid[!is.na(AP.SARIMA$resid)]) # potential heteroscedasticity

### ARCH test
library(aTSA)
arch.test(AP.ARIMA) # heteroscedastic

### remedy: GARCH
library(rugarch)
AP_GARCH_spec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(14,0)))
AP_GARCH_spec
AP_GARCH <- ugarchfit(spec = AP_GARCH_spec, data = AP)
AP_GARCH

plot(AP_GARCH@fit$resid[!is.na(AP_GARCH@fit$resid)],type = "l") # potential heteroscedasticity

