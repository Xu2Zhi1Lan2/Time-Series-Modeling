### Example: multi-regime stock return distribution
### S&P500 monthly data 1987Jan-2016Dec
library(HiddenMarkov)

### read data
data <- read.csv(file="SP500.csv",header=T)

sp <- data[,3]
time <- as.Date(data[,1],format="%d-%b-%Y")
time <- as.Date(ifelse(time > "0016-12-30", 
                       format(time, "19%y-%m-%d"), 
                       format(time, "20%y-%m-%d")))

plot(time,sp,type='l',xlab='Time',ylab='SP500')
abline(h=0)

hist(sp,nclass=40,main="")

### 3-regime Gaussian mixture
# initial values
Pi <- diag(rep(0.8,3))+matrix(rep(1,9)/3*0.2,ncol=3,nrow=3)
delta <- rep(1,3)/3
dist <- 'norm'      # mixture normal
pm <- list(mean=c(-0.3,0.05,0.2),sd=c(0.1,0.1,0.5))

# estimation
temp <- dthmm(sp,Pi,delta,dist,pm)   #setting model
out3 <- BaumWelch(temp,control=bwcontrol(prt=FALSE)) #Baum Welch estimation
summary(out3)

# smoothed State  probability
par(mfrow=c(3,1))
plot(time,out3$u[,1],type='l',xlab='time',ylim=c(0,1),ylab='prob',main='State 1 Prob')
plot(time,out3$u[,2],type='l',xlab='time',ylim=c(0,1),ylab='prob',main='State 2 Prob')
plot(time,out3$u[,3],type='l',xlab='time',ylim=c(0,1),ylab='prob',main='State 3 Prob')

# most likely state path
par(mfrow=c(1,1))
v3 <- Viterbi(out3)   # most likely path by Viterbi
plot(time,v3,cex=0.5,xlab='time',ylab='state') #plot most likely path

# diagnostics
res3 <- residuals(out3) 
pacf(res3**2,main='PACF residual squared')    # pacf residual squared
pacf(sp**2, main='PACF SP return squared')    # pacf sp return squared
qqnorm(res3, main='residual')    # qq plot residual
qqnorm(sp,main='SP return')      # qq plot sp return
