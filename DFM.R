##########################################
### use dynamic factor model
library(dlm)

## step 1: buildDFM3
buildDFM3=function(parm){
  hh <- c(parm[1],parm[2],0,parm[3])
  gg <- c(1,1,parm[4:11])
  qq <- c(parm[12:13],0,parm[14])
  rr <- parm[15:19]
  H <- matrix(hh,ncol=2,nrow=2,byrow=T)	
  G <- matrix(gg,ncol=2,nrow=5,byrow=T)
  W <- matrix(qq,ncol=2,nrow=2,byrow=T)
  V <- diag(rr)
  return(list(
    m0=rep(0,2),
    C0=5*diag(2),
    FF=G,
    GG=H,
    V=V%*%t(V),
    W=W%*%t(W)
  ))
}

## step 2: out-of-sample prediction 
# DFM out-of-sample prediction
outsample.prediction.DFM=function(x,np,n.Start,n.End,par.init){
  xpred <- x
  SS <- rep(0,5)
  for(i in (n.Start-np):(n.End-np)){
    MLE <- dlmMLE(x[1:i,],par.init,buildDFM3)
    myMode <- buildDFM3(MLE$par)
    filter <- dlmFilter(x[1:i,],myMode)
    pp <- dlmForecast(filter,nAhead=np)
    xpred[i+np,] <- pp$f[np,]
    SS <- SS+(x[i+np,]-xpred[i+np,])**2
  }
  return(list(xpred=xpred,SS=SS,meanSS=mean(SS)))
}

# PCA out-of-sample prediction
outsample.prediction.PCA=function(x,np,n.Start,n.End){
  xpred <- x
  SS <- rep(0,5)
  for(i in (n.Start-np):(n.End-np)){
    PCA <- princomp(x[1:i,], cor=F)
    ff <- PCA$score[,1:2]
    out <- VAR(ff)
    pp <- predict(out,n.ahead=np)
    pp0 <- c(pp$fcst$Comp.1[np,1],pp$fcst$Comp.2[np,1])
    xpred[i+np,] <- PCA$loading[,1:2]%*%pp0
    SS <- SS+(x[i+np,]-xpred[i+np,])**2
  }
  return(list(xpred=xpred,SS=SS,meanSS=mean(SS)))
}

# AR prediction
outsample.prediction.uAR=function(x,np,n.Start,n.End){
  xpred <- x
  SS <- rep(0,5)
  for(i in (n.Start-np):(n.End-np)){
    for(j in 1:5){
      out <- arima(x[1:i,j],c(1,0,0))
      xpred[i+np,j] <- predict(out,n.ahead=np)$pred[np]
    }
    SS <- SS+(x[i+np,]-xpred[i+np,])**2
  }
  return(list(xpred=xpred,SS=SS,meanSS=mean(SS)))
}

## step 3: model estimation and comparison
xx.nm <- matrix(scan("GDP_panel_nm.dat"),ncol=5)

# model 1: PCA (use first two PCs to build a VAR)
PCAout <- princomp(xx.nm, cor=F)
ff <- PCAout$scores[,1:2]
library(vars)
out <- VAR(ff)

# model 2: dynamic factor model
hh3 <- c(0.7,0.2,0.8)
gg3 <- rep(c(0.7,0.7), 4)
qq3 <- c(1,-1,1)
rr3 <- rep(0.3,5)
parm3 <- c(hh3,gg3,qq3,rr3)  # initial parameter values
outMLE3 <- dlmMLE(xx.nm,parm3,buildDFM3)
myModel3 <- buildDFM3(outMLE3$par)
outSmooth3 <- dlmSmooth(xx.nm,myModel3)
fit3 <- xx.nm
for(i in 1:5){
  fit3[,i] <- outSmooth3$s[-1,]%*%myModel3$F[i,]
}
print((sum(xx.nm-fit3)**2)/sum(xx.nm**2)) # relative fit error

# compare the predictions
for(np in 1:3){
  pred.d <- outsample.prediction.DFM(xx.nm,np,88,107,outMLE3$par)
  pred.p <- outsample.prediction.PCA(xx.nm,np,88,107)
  pred.ar <- outsample.prediction.uAR(xx.nm,np,88,107)
}

