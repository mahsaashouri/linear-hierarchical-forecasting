
lapply( c("quantmod","forecast","dplyr","plyr","stringr","hts","data.table", "Matrix"),require,character.only = T)
wikipedia_data <-read.csv("wikipedia_data.csv", header = TRUE)
wikipedia_wide <- wikipedia_data$views %>%
  matrix(nrow = 394, ncol = 913) %>%
  as.data.frame() %>%
  ts(frequency = 7)
colnames(wikipedia_wide) <- unique(wikipedia_data$cat_column) %>% substr(1,14)

##################
## using hierarchies and groupings up to 2-way combinations
##################
wikigts <- gts(wikipedia_wide, character=c(7,2,2,3),
               gnames = c("Access",
                          "Agent",
                          "Language",
                          "Purpose",
                          "Access x Agent",
                          "Access x Language",
                          "Access x Purpose",
                          "Agent x Language",
                          "Agent x Purpose",
                          "Language x Purpose"
               ))
# Splitting data into training and validation sets
wikitrain<-window(wikigts,end=c(1,366))
wikitest<-window(wikigts,start=c(1,367))
# Construct matrix of all time series including aggregates
ally <- aggts(wikigts)

# Set the reconciled matrix
gmat<-GmatrixG(wikigts$groups)
smatrix<-SmatrixM(gmat)
wvec<- InvS4g(wikigts$groups)
lambda <- diag(wvec)
rec.adj.test <- as.matrix(smatrix%*%solve(t(smatrix)%*%solve(lambda)%*%smatrix)%*%t(smatrix)%*%solve(lambda))

## computing ETS and ARIMA base forecasts
n<-nrow(wikipedia_wide)
h <- NROW(wikitest$bts)
fc.arima.ets <- array(NA, c(Horizon=h, Series=NCOL(ally), Method=2))
dimnames(fc.arima.ets) <- list(
  Horizon = paste0("h=",seq(h)),
  Series = colnames(ally),
  Method = c("ETS","ARIMA")
)
for(i in seq(NCOL(ally))){
  for(j in 1:h)
  { 
    wikitrain.1 <- window(wikigts, start=c(1,1),end = c(1, (n - h) + (j - 1)))
    ally.1<-aggts(wikitrain.1)
    # ETS forecast
    fc.arima.ets[j,i,"ETS"] <- forecast(ets(ally.1[,i]), h=1)$mean
    # ARIMA forecast
    fc.arima.ets[j,i,"ARIMA"] <- forecast(auto.arima(ally.1[,i]), h=1)$mean
  }
}
### setting negative base forecasts zero  
fc.arima.ets[fc.arima.ets<0]<-0

### reconciled the base forecasts
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
for(i in seq(dim(errors)[3]))
  errors[,,i,"Unreconciled"] <- aggts(wikitest) - fc.arima.ets[,,i]

# Compute errors for reconciled forecasts
for(i in seq(dim(errors)[3]))
{
  errors[,,i,"Reconciled"] <-  aggts(wikitest) - fc.arima.ets.rec[,,i,"Reconciled"]
}

## OLS
#### OLS rolling function
OLSmodel<-function(X,freq,maxlag,h){
  X<-as.vector(X)
  trend1 <-seq(NROW(X))
  trend2 <- (seq(NROW(X)))^2
  season <- forecast::seasonaldummy(ts(X,frequency = freq))
  Xlag <- quantmod::Lag(X,k=1:maxlag)
  X_mat <- cbind.data.frame(X, trend1, trend2, season, Xlag)
  n <- nrow(X_mat)
  fore_base_OLS <- matrix(NA,nrow = h,ncol=1)
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
  fc.OLS.base[,i,"OLS"] <- OLSmodel(ally[,i],7,7,28)
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
  errors[,,i,"Unreconciled"] <- aggts(wikitest) - fc.OLS.base[,,i]

# Compute errors for reconciled forecasts
for(i in seq(dim(errors)[3]))
{
  errors[,,i,"Reconciled"] <-  aggts(wikitest) - fc.OLS.rec[,,i,"Reconciled"]
}


