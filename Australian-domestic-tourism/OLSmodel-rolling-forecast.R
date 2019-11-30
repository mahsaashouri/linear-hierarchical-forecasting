
lapply( c("quantmod","forecast","dplyr","plyr","stringr","hts","data.table", "Matrix"),require,character.only = T)
## reading data
TourismData <-ts(read.csv("TourismData_v3.csv", header = TRUE)[-c(1,2)],start=1,frequency =12)
ausgts <- gts(TourismData, characters = list(c(1, 1, 1), 3),
              gnames = c("State", "Zone", "Region", "Purpose","State x Purpose", "Zone x Purpose"))
########hierarchy+ARIMA-hierarchy+ets
k<-24
n<-nrow(TourismData)
train_tourist <-window(ausgts,start = c(1, 1),end = c(1, (n-k)))
validation_tourist <-window(ausgts,start = c(1, ((n-k)+1)),end = c(1, n))
ally <- aggts(ausgts)

## set the reconceliation matrix using smatrix
gmat<-GmatrixG(ausgts$groups)
smatrix<-SmatrixM(gmat)
wvec<- InvS4g(ausgts$groups)
lambda <- diag(wvec)
rec.adj <- as.matrix(smatrix%*%solve(t(smatrix)%*%solve(lambda)%*%smatrix)%*%t(smatrix)%*%solve(lambda))

##### ARIMA & ARIMAlog & ETS & ETSlog

# Set up array for forecasts 
h <- NROW(validation_tourist$bts)
fc.arima.ets <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=2))
dimnames(fc.arima.ets) <- list(
  Horizon = paste0("h=",seq(h)),
  Series = colnames(ally),
  Method = c("ETS","ARIMA")
)
for(i in seq(NCOL(ally))){
  for(j in 1:h)
  { 
    austrain.1 <- window(ausgts, start=c(1,1),end = c(1, (n - h) + (j - 1)))
    ally.1<-aggts(austrain.1)
    # ETS forecast
    fc.arima.ets[j,i,"ETS"] <- forecast(ets(ally.1[,i]), h=1)$mean
    # ARIMA forecast
    fc.arima.ets[j,i,"ARIMA"] <- forecast(auto.arima(ally.1[,i]), h=1)$mean
  }
}
### setting negative base forecasts zero  
fc.arima.ets[fc.arima.ets<0]<-0

### reconciled the results
fc.arima.ets.rec <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=2,Reconciled=2))
dimnames(fc.arima.ets.rec) <- list(
  Horizon = paste0("h=",seq(h)),
  Series = colnames(ally),
  Method = c("ETS", "ARIMA"),
  Reconciled = c("Reconciled","Unreconciled")
)

## Compute unreconciled forecasts
for(i in seq(dim(fc.arima.ets.rec)[3]))
  fc.arima.ets.rec[,,i,"Unreconciled"] <- fc.arima.ets[,,i]
## Compute reconciled forecasts
for(i in seq(dim(fc.arima.ets.rec)[3])){
  for(j in seq(dim(fc.arima.ets.rec)[1]))
    fc.arima.ets.rec[j,,i,"Reconciled"] <- rec.adj %*% fc.arima.ets[j,,i] 
}

errors <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=dim(fc.arima.ets)[3], Reconciled=2))
dimnames(errors) <- list(
  Horizon = dimnames(fc.arima.ets)[[1]],
  Series = colnames(ally),
  Method = dimnames(fc.arima.ets)[[3]],
  Reconciled = c("Reconciled","Unreconciled")
)

# Compute errors for unreconciled forecasts
nbts <- NCOL(aus)
nseries <- NCOL(ally)
for(i in seq(dim(errors)[3]))
  errors[,,i,"Unreconciled"] <- aggts(validation_tourist) - fc.arima.ets[,,i]

# Compute errors for reconciled forecasts
for(i in seq(dim(errors)[3]))
{
  errors[,,i,"Reconciled"] <-  aggts(validation_tourist) - fc.arima.ets.rec[,,i,"Reconciled"]
}


#### creating forecast: OLS and OLSlog
#### OLS rolling function
OLSmodel<-function(X,freq,maxlag,h){
  X<-as.vector(X)
  trend<-seq(NROW(X))
  season<-forecast::seasonaldummy(ts(X,frequency = freq))
  Xlag<-quantmod::Lag(X,k=1:maxlag)
  X_mat<-cbind.data.frame(X,trend,season,Xlag)
  n <- nrow(X_mat)
  fore_base_OLS<-matrix(NA,nrow = h,ncol=1)
  for (i in 1:h) {
    train.1 <- X_mat[1:((n - h) + (i - 1)), ]
    valid.1 <- X_mat[(n - h) + i, ]
    fit <- lm(X ~. , data = train.1)
    fore <- predict.lm( fit , newdata = valid.1)
    fore_base_OLS[i,]<-fore
  }
  return(fore_base_OLS)
}

### computing base forecasts
fc.OLS.base <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=1))
dimnames(fc.OLS.base) <- list(
  Horizon = paste0("h=",seq(h)),
  Series = colnames(ally),
  Method = c("OLS")
)
for(i in seq(NCOL(ally)))
{
  # OLS forecasts
  fc.OLS.base[,i,"OLS"] <- OLSmodel(ally[,i],12,12,24)
}
### setting negative base forecasts zero
fc.OLS.base[fc.OLS.base<0]<-0

#### OLS reconcile
fc.OLS.rec <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=1,Reconciled=2))
dimnames(fc.OLS.rec) <- list(
  Horizon = paste0("h=",seq(h)),
  Series = colnames(ally),
  Method = c("OLS"),
  Reconciled = c("Reconciled","Unreconciled")
)

## Compute unreconciled forecasts
for(i in seq(dim(fc.OLS.rec)[3]))
  fc.OLS.rec[,,i,"Unreconciled"] <- fc.OLS.base[,,i]
## Compute reconciled forecasts
for(i in seq(dim(fc.OLS.rec)[3])){
  for(j in seq(dim(fc.OLS.rec)[1]))
    fc.OLS.rec[j,,i,"Reconciled"] <- rec.adj %*% fc.OLS.base[j,,i] 
}

errors <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=dim(fc.OLS.base)[3], Reconciled=2))
dimnames(errors) <- list(
  Horizon = dimnames(fc)[[1]],
  Series = colnames(ally),
  Method = dimnames(fc.OLS.base)[[3]],
  Reconciled = c("Reconciled","Unreconciled")
)

# Compute errors for unreconciled forecasts
nbts <- NCOL(aus)
nseries <- NCOL(ally)
for(i in seq(dim(errors)[3]))
  errors[,,i,"Unreconciled"] <- aggts(validation_tourist) - fc.OLS.base[,,i]

# Compute errors for reconciled forecasts
for(i in seq(dim(errors)[3]))
{
  errors[,,i,"Reconciled"] <-  aggts(validation_tourist) - fc.OLS.rec[,,i,"Reconciled"]
}


