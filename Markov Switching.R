###################################################################################
library(MSwM)
### Part 1: A simulated example
## Generate data: T = 300, k = 2, p = 1
set.seed(2023)
n <- 300; np1 <- n+1

St <- c(rep(1,51),rep(2,20),rep(1,30),rep(2,30),rep(1,70),rep(2,10),
        rep(1,50),rep(2,20),rep(1,20)) # specify a sequence of state variables
zt <- runif(np1) ## auto-regressor
at <- rnorm(np1)  ## noise

x <- at[1]  ### set initial value
for (i in 2:np1){
  tmp = 0
  if(St[i]==1){
    tmp = 2.0 + 0.8*x[i-1] + 1.5*at[i]
  }
  else{
    tmp = -2.0+0.6*zt[i]-0.8*x[i-1]+at[i]
  }
  x <- c(x,tmp)
} # generate data w.r.t. state variables

xt <- x[-1]; zt <- zt[-1]; St <- St[-1]  ## Remove the initial value
x <- cbind(xt,rep(1,300),zt) # Start model fitting
colnames(x) <- c("xt","cnst","zt")
X <- data.frame(x)

## estimate an MSM
sim.lm <- lm(xt~-1+cnst+zt,data=X)
sim.MSM <- msmFit(sim.lm,k=2,p=1,sw=c(T,T,T,T)) # sw is a logical vector indicating which coefficients have switching.
summary(sim.MSM)

plotDiag(sim.MSM) # residual serial correlation check

par(mfcol=c(2,1))
plotDiag(sim.MSM,which=1) # Residual plot
plotDiag(sim.MSM,which=2) # Q-Q plot of residuals

par(mfrow=c(2,1))
plot(sim.MSM@Fit@filtProb[,1],type='l',main='Filtered Probability Regime 1',
       ylab='prob',xlab='time')
plot(sim.MSM@Fit@smoProb[,1],type='l',main='Smoothed Probability Regime 1',
       ylab='prob',xlab='time')

### Part 2: A real data example
## prepare the data
data=read.table("GDPC1.txt",header=T) # US quarterly GDP data
gdp=diff(log(data[,4])) # take diff log to make it stationary (gdp growth rate)
par(mfrow=c(1,1))
plot(gdp,type = "l") # 150 = year 1980, there is an obvious change 
                     # named the Great Moderation stopped by 2008 Financial crisis
                     # hence, we consider k = 3
x <- cbind(gdp,rep(1,length(gdp)))
colnames(x) <- c("gdp","cnst")
X <- data.frame(x)

## estimate an MSM
library(forecast)
GDP.arma <- auto.arima(gdp,seasonal = F,stepwise = F, trace = F)
summary(GDP.arma) # sugests p = 3

GDP.lm <- lm(gdp~-1+cnst,data=X) 
GDP.MSM <- msmFit(GDP.lm,k=3,p=3,sw=c(T,T,T,T,T))
summary(GDP.MSM)

plotDiag(GDP.MSM) 

par(mfcol=c(2,1))
plotDiag(GDP.MSM,which=1) 
plotDiag(GDP.MSM,which=2) 

par(mfrow=c(3,1))
plot(GDP.MSM@Fit@filtProb[,1],type='l',main='Filtered Probability Regime 1',
       ylab='prob',xlab='time')
plot(GDP.MSM@Fit@filtProb[,2],type='l',main='Filtered Probability Regime 2',
       ylab='prob',xlab='time')
plot(GDP.MSM@Fit@filtProb[,3],type='l',main='Filtered Probability Regime 3',
       ylab='prob',xlab='time')

plot(GDP.MSM@Fit@smoProb[,1],type='l',main='Smoothed Probability Regime 1',
       ylab='prob',xlab='time')
plot(GDP.MSM@Fit@smoProb[,2],type='l',main='Smoothed Probability Regime 2',
       ylab='prob',xlab='time')
plot(GDP.MSM@Fit@smoProb[,3],type='l',main='Smoothed Probability Regime 3',
       ylab='prob',xlab='time')


