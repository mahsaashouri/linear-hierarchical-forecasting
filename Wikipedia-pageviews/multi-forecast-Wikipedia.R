### Wikipedia dataset
##########################
### 28-step-ahead
##########################

lapply(c("quantmod", "forecast", "dplyr", "plyr", "stringr", "hts", "data.table", "Matrix"),
  require,
  character.only = TRUE
)
wikipedia_data <- read.csv("wikipedia_data.csv", header = TRUE)

## Data reshaping
# 394: length of each series and 913: number of series
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
                          "Language x Purpose"))
                          
# Splitting data into training and test sets
wikitrain <- window(wikigts, end = c(1, 366))
wikitest <- window(wikigts, start = c(1, 367))

# Construct matrix of all time series including aggregates
ally <- aggts(wikitrain)

# Set the reconciled matrix
gmat<-GmatrixG(wikigts$groups)
smatrix<-SmatrixM(gmat)
wvec<- InvS4g(wikigts$groups)
lambda <- diag(wvec)
rec.adj.test <- as.matrix(smatrix%*%solve(t(smatrix)%*%solve(lambda)%*%smatrix)%*%t(smatrix)%*%solve(lambda))


# Set up array for forecasts
h <- NROW(wikitest$bts)
fc <- array(NA, c(Horizon = h, Series = NCOL(ally), Method = 3))
dimnames(fc) <- list(
  Horizon = paste0("h=", seq(h)),
  Series = colnames(ally),
  Method = c("ETS", "ARIMA", "OLS")
)

# Create forecasts for all methods
for (i in seq(NCOL(ally))) {
  fc[, i, "ETS"] <- pmax(forecast(ets(ally[, i]), h = h)$mean, 0)
  fc[, i, "ARIMA"] <- pmax(forecast(auto.arima(ally[, i]), h = h)$mean, 0)
  fc[, i, "OLS"] <- pmax(olsfc(ally[, i], h = h, maxlag = 7), 0)
}


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

