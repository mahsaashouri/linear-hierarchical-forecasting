library(tidyverse)
library(hts)
library(Martix)
source("olsfc.R")

# Reading data and adding separate category columns
aus <- ts(readr::read_csv("TourismData_v3.csv")[, -(1:2)],
  start = 1998, frequency = 12)
ausgts <- gts(aus, characters = list(c(1, 1, 1), 3),
              gnames = c("State", "Zone", "Region", "Purpose",
                         "State x Purpose", "Zone x Purpose"))

# Splitting data into training and validation sets
austrain <- window(ausgts, end=c(2014,12))
austest  <- window(ausgts, start=c(2015,1))

easter.info<-as.data.frame(easter(aus,easter.mon = TRUE))
easter.info.train<-easter.info[(1:204),]
easter.info.test<- easter.info[(205:228),]
# Construct matrix of all time series including aggregates
ally <- aggts(austrain)

# Set up array for forecasts
h <- NROW(austest$bts)
fc <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=5))
dimnames(fc) <- list(
  Horizon = paste0("h=",seq(h)),
  Series = colnames(ally),
  Method = c("ETS","ARIMA","ARIMAX","OLS","OLSX")
)

# Create forecasts for all methods
for(i in seq(NCOL(ally)))
{
  # ETS forecasts
  fc[,i,"ETS"] <- pmax(forecast(ets(ally[,i]), h=h)$mean,0)
  # ARIMA forecasts
  fc[,i,"ARIMA"] <- pmax(forecast(auto.arima(ally[,i]), h=h)$mean,0)
  # ARIMAX forecasts 
  fc[,i,"ARIMAX"] <- pmax(forecast(auto.arima(ally[,i],xreg = easter.info.train), xreg=easter.info.test, h=h)$mean,0)
  # OLS forecasts
  fc[,i,"OLS"] <- pmax(olsfc(ally[,i], h=h, maxlag = 12, nolag = c(1,12)),0)
  # OLSX forecasts
  fc[,i,"OLSX"] <- pmax(olsfc.external(ally[,i], externaldata=easter.info, h=h, maxlag = 12, nolag = c(1,12)),0)
}

## set the reconceliation matrix using smatrix
gmat<-GmatrixG(ausgts$groups)
smatrix<-SmatrixM(gmat)
wvec<- InvS4g(ausgts$groups)
lambda <- diag(wvec)
rec.adj <- as.matrix(smatrix%*%solve(t(smatrix)%*%solve(lambda)%*%smatrix)%*%t(smatrix)%*%solve(lambda))


## Set up array for unreconciled and reconciled forecasts
fc.rec <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=dim(fc)[3], Reconciled=2))
dimnames(fc.rec) <- list(
  Horizon = dimnames(fc)[[1]],
  Series = colnames(ally),
  Method = dimnames(fc)[[3]],
  Reconciled = c("Reconciled","Unreconciled")
)

## Compute unreconciled forecasts
for(i in seq(dim(fc.rec)[3]))
  fc.rec[,,i,"Unreconciled"] <- fc[,,i]

# Compute reconciled forecasts
for(i in seq(dim(fc.rec)[3])){
  for(j in seq(dim(fc.rec)[1]))
    fc.rec[j,,i,"Reconciled"] <- rec.adj %*% fc[j,,i] 
}


errors <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=dim(fc)[3], Reconciled=2))
dimnames(errors) <- list(
  Horizon = dimnames(fc)[[1]],
  Series = colnames(ally),
  Method = dimnames(fc)[[3]],
  Reconciled = c("Reconciled","Unreconciled")
)

# Compute errors for unreconciled forecasts
nbts <- NCOL(aus)
nseries <- NCOL(ally)
for(i in seq(dim(errors)[3]))
  errors[,,i,"Unreconciled"] <- aggts(austest) - fc[,,i]

# Compute errors for reconciled forecasts
for(i in seq(dim(errors)[3]))
{
  errors[,,i,"Reconciled"] <-  aggts(austest) - fc.rec[,,i,"Reconciled"]
}
