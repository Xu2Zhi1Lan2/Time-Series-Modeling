### data: percentage changes in quarterly personal consumption expenditure, 
### personal disposable income, production, savings, 
## and the unemployment rate for the US, 1970 to 2016
library(fpp2)
library(ggplot2)
library(forecast)
data(uschange)
summary(uschange)

########################################################################
#### descriptive statistics

# plot the data
autoplot(uschange, facets=TRUE) +
  xlab("Year") + ylab("") +
  ggtitle("Quarterly percentage changes")

# plot the ACF
ggAcf(uschange)

# plot the PACF
ggPacf(uschange)
ggPacf(uschange[,-4])

# decomposition of "Consumption"
uschange_freq_4 = uschange[,"Consumption"] %>% ts(., frequency=4) # freq = 4 quarters
uschange_stl = stl(uschange_freq_4, s.window = "periodic")
autoplot(uschange_stl) # looks seasonal effect is not that significant, SARIMAX shouldn't be a lot better

########################################################################
#### ARIMAX
uschange_arimax = auto.arima(uschange[,1], # specify main trend
                             xreg = uschange[,2:5], # specify exogenous variables
                             trace = TRUE, 
                             seasonal= FALSE, # if F, restrict to non-seasonal model
                             stepwise=FALSE,    # if yes, fo stepwise selection; otherwise searches all models
                             approximation=FALSE) # if yes, approximate SS in information criteria like AIC, BIC
summary(uschange_arimax) # AIC = 125.5, BIC = 144.85
autoplot(uschange_arimax) # check if stationary
checkresiduals(uschange_arimax) # check linearity

#### SARIMAX
uschange_sarimax = auto.arima(uschange[,1], 
                              xreg = uschange[,2:5], 
                              trace = TRUE, 
                              seasonal= TRUE, # allow a SARIMAX model
                              stepwise=FALSE,
                              approximation=FALSE) # seasonal model is not preferred
summary(uschange_sarimax) # keep income and saving
autoplot(uschange_sarimax) # nearly non-stationary
checkresiduals(uschange_sarimax)

# retrain model
uschange_sarimax2 = auto.arima(uschange[,1],
                              xreg = uschange[,c(2,4)], 
                              trace = FALSE, 
                              seasonal= TRUE, 
                              stepwise= FALSE,
                              approximation=FALSE)
summary(uschange_sarimax2) # getting worse, can't drop the two predictors

# compare with SARIMA
uschange_sarimax = auto.arima(uschange[,1], 
                              trace = TRUE, 
                              seasonal= TRUE, 
                              stepwise=FALSE,
                              approximation=FALSE) 
summary(uschange_sarimax) # checking AIC, BIC, exogenous predictors do provide additional information
