library(forecast)
library(zoo)
library(dplyr)
library(MLmetrics)

#####
#####Question 1
#####

nasa_data <- read.csv("C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/GLB.Ts+dSST_Processed.csv") #please select the data file

nasa_data <- mutate(nasa_data, temperature_degrees = nasa_data$Mean_Temp + 14, temperature_fahrenheit = nasa_data$Mean_Temp + 57.2) # converting anomalies to actual temperature

nasa_data$DateMon <- as.yearmon(paste(nasa_data$Year, nasa_data$Month), "%Y %m") # creating a date column

write.csv(nasa_data,"C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/nasa_data_1880.csv")

##Time Series Models

# Create exponential smoothing models: additive vs multiplicative noise (first A vs M), additive vs multiplicative trend (second A vs M) and no seasonality vs automatic detection (third N vs Z) trend and no seasonlity (AAN), multiplicative (MMN)
nasa_ts <- ts(nasa_data$temperature_degrees,start=1880, frequency=12) 
fit <- decompose(nasa_ts) #decompose using "classical" method, multiplicative form
plot(fit)

nasa_ts_AAN <- ets(nasa_ts, model="AAN")
nasa_ts_AAZ <- ets(nasa_ts, model="AAZ", damped=FALSE)
nasa_ts_MMZ <- ets(nasa_ts, model="MMZ", damped=FALSE)# seasonality change to automatic results in "N" for this data

# Create their prediction
nasa_ts_AAN_pred <- forecast(nasa_ts_AAN, h=972, level=c(0.8, 0.90))
nasa_ts_AAZ_pred <- forecast(nasa_ts_AAZ, h=972, level=c(0.8, 0.90))
nasa_ts_MMZ_pred <- forecast(nasa_ts_MMZ, h=972, level=c(0.8, 0.90))

par(mfrow=c(1,3)) # This command sets the plot window to show 1 row of 4 plots
plot(nasa_ts_AAN_pred, xlab="Year", ylab="Predicted Climate") # does not captures trend
plot(nasa_ts_AAZ_pred, xlab="Year", ylab="Predicted Climate")
plot(nasa_ts_MMZ_pred, xlab="Year", ylab="Predicted Climate")

#Create a trigonometric box-cox autoregressive trend seasonality (TBATS) model
nasa_tbats <- tbats(nasa_ts)
nasa_tbats_pred <-forecast(nasa_tbats, h=972, level=c(0.8, 0.90))
par(mfrow=c(1,1))
plot(nasa_tbats_pred, xlab="Year", ylab="Predicted Climate")

par(mfrow=c(1,3)) # Lets look at the three models with seasonality on one graph on the same scale
plot(nasa_ts_AAZ_pred, xlab="Year", ylab="Predicted Climate")
plot(nasa_ts_MMZ_pred, xlab="Year", ylab="Predicted Climate")
plot(nasa_tbats_pred, xlab="Year", ylab="Predicted Climate") # does not capture the trend

###Arima models

#Auto Arima Model #(2,0,2)(2,1,0)
par(mfrow=c(1,1))
nasa_arima <- auto.arima(nasa_ts,seasonal=TRUE, trace = TRUE, D = 1)
nasa_arima_pred <- forecast(nasa_arima, h=972, level=c(0.8, 0.90))
plot(nasa_arima_pred, xlab="Year", ylab="Predicted Climate")
plot(nasa_arima)
###Introducing multiple seasonality
nasa_msts <- msts(nasa_data$temperature_degrees, seasonal.periods = c(12,120))
fit <- decompose(nasa_msts) 
plot(fit)

#MSTS Auto Arima (3,1,1)(0,1,0)
nasa_arima2 <- auto.arima(nasa_msts, seasonal = TRUE, trace = TRUE, D = 1)
plot(nasa_arima2)
nasa_arima_pred2 <- forecast(nasa_arima2, h=972, level=c(0.8,0.90))
plot(nasa_arima_pred2, xlab = 'Year', ylab = "Predicted Climate")

# ARIMA with regressors
monthMatrix <- cbind(month = model.matrix(~as.factor(nasa_data$Month)))
monthMatrix <- monthMatrix[,-1]
colnames(monthMatrix) <- c("Feb", "Mar", "Apr", "May", "Jun","Jul","Aug","Sep","Oct","Nov", "Dec")
matrix_of_regressors <- monthMatrix
nasa_arima3 <- auto.arima(nasa_data$temperature_degrees, xreg=matrix_of_regressors)
nasa_arima3
xreg.pred <- matrix_of_regressors[-c(973:1683),]
nasa_arima_pred3 <- forecast(nasa_arima3, h = 972, xreg = xreg.pred, level = c(0.8, 0.90))
plot(nasa_arima_pred3, xlab = "Year", ylab = "Predicted Climate")

# ARIMA on residuals
nasa_msts2 <- tslm(nasa_msts ~ trend + season)
summary(nasa_msts2)
residarima1 <- auto.arima(nasa_msts2$residuals)
residarima1
residualsArimaForecast <- forecast(residarima1, h=972)
residualsF <- as.numeric(residualsArimaForecast$mean)
regressionForecast <- forecast(nasa_msts2, h=972)
regressionF <- as.numeric(regressionForecast$mean)
forecastR <- regressionF + residualsF
print(forecastR)
for (i in 1:972) {points(i+1680, forecastR[i], col='red', pch=19, cex=0.5)}

### Comparing models -- Time series Cross Validation (Rolling Horizon Holdout)

f_AAN  <- function(y, h) forecast(ets(y, model="AAN"), h = h) 
errors_AAN <- tsCV(nasa_ts, f_AAN, h=1, window=60)

f_MMN  <- function(y, h) forecast(ets(y, model="MMN"), h = h)
errors_MMN <- tsCV(nasa_ts, f_MMN, h=1, window=60)

f_auto_arima2 <- function(y, h) forecast(Arima(y), order = c(3,1,1), seasonal=c(0,1,0), h = h) #MSTS Auto Arima (3,1,1)(0,1,0)
errors_auto.arima2 <- tsCV(nasa_msts, f_auto_arima2, h=1, window=60)

mape_AAN <-  mean(abs(errors_AAN/nasa_ts), na.rm=TRUE)*100 #0.6423788
mape_MMN <- mean(abs(errors_MMN/nasa_ts), na.rm=TRUE)*100 #0.6425854
mape_auto.arima2 <-mean(abs(errors_auto.arima2/nasa_msts), na.rm=TRUE)*100 #0.8414727

errors <- data.frame("AAN" = mape_AAN, "MMN" = mape_MMN,"auto_arima" = mape_auto.arima2)
write.csv(errors,"C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/nasa_errors.csv")

nasa_arima_pred2

nasa_date <- data.frame("Date" =seq(as.Date("2020/4/1"), by = "month", length.out =972))
nasa_prediction <- cbind(nasa_date,nasa_arima_pred2)

write.csv(nasa_prediction,"C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/nasa_prediction.csv", row.names = FALSE)

#### Performing Arima Models on the Met data

met_data <- read.csv("C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/HadCRUT4_median_TS_monthly.csv")

names(met_data)[names(met_data) == "ï..Year"] <- "Year"
names(met_data)[names(met_data) == "Temp"] <- "Mean_Temp"

met_data <- mutate(met_data, temperature_degrees = met_data$Mean_Temp + 14, temperature_fahrenheit = met_data$Mean_Temp + 57.2) # converting anomalies to actual temperature

met_ts <- ts(met_data$temperature_degrees,start=1850, frequency=12) 
fit <- decompose(met_ts) #decompose using "classical" method, multiplicative form
plot(fit)

met_msts <- msts(met_data$temperature_degrees, seasonal.periods = c(12,120))
fit <- decompose(met_msts) 
plot(fit) 

#MSTS Auto Arima (3,1,1)(0,1,0)
met_arima <- auto.arima(met_msts, seasonal = TRUE, trace = TRUE, D = 1)
plot(met_arima)
met_arima_pred <- forecast(met_arima, h=972, level=c(0.8,0.90))
plot(met_arima_pred, xlab = 'Year', ylab = "Predicted Climate")

f_auto_arima_met <- function(y, h) forecast(Arima(y), order = c(3,1,1), seasonal=c(0,1,0), h = h) #MSTS Auto Arima (3,1,1)(0,1,0)
errors_auto.arima_met <- tsCV(met_msts, f_auto_arima_met, h=1, window=60)
mape_auto.arima_met <-mean(abs(errors_auto.arima_met/met_msts), na.rm=TRUE)*100 #0.8898724

met_date <- data.frame("Date" =seq(as.Date("2020/3/1"), by = "month", length.out =972))
met_prediction <- cbind(met_date,met_arima_pred)

write.csv(met_prediction,"C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/met_prediction.csv", row.names = FALSE)

####
###Question 2, point predictions January and July 2030 2050 and 2100
####

answer_2_nasa_prep <- data.frame("Date" = nasa_prediction$Date, "point_forecast" = nasa_prediction$`Point Forecast`, "lo_90" = nasa_prediction$`Lo 90`, "hi_90" = nasa_prediction$`Hi 90`)

Jan_2030 <- answer_2_nasa_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2030-01-01")
Jul_2030 <- answer_2_nasa_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2030-07-01")
Jan_2050 <- answer_2_nasa_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2050-01-01")
Jul_2050 <- answer_2_nasa_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2050-07-01")
Jan_2100 <- answer_2_nasa_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2100-01-01")
Jul_2100 <- answer_2_nasa_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2100-07-01")

answer_2_nasa <- rbind(Jan_2030,Jul_2030,Jan_2050,Jul_2050,Jan_2100,Jul_2100)

write.csv(answer_2_nasa,"C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/answer_2_nasa.csv", row.names = FALSE)

answer_2_met_prep <- data.frame("Date" = met_prediction$Date, "point_forecast" = met_prediction$`Point Forecast`, "lo_90" = met_prediction$`Lo 90`, "hi_90" = met_prediction$`Hi 90`)

Jan_2030 <- answer_2_met_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2030-01-01")
Jul_2030 <- answer_2_met_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2030-07-01")
Jan_2050 <- answer_2_met_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2050-01-01")
Jul_2050 <- answer_2_met_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2050-07-01")
Jan_2100 <- answer_2_met_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2100-01-01")
Jul_2100 <- answer_2_met_prep %>% select(Date, point_forecast, lo_90, hi_90) %>% filter(Date == "2100-07-01")

answer_2_met <- rbind(Jan_2030,Jul_2030,Jan_2050,Jul_2050,Jan_2100,Jul_2100)

write.csv(answer_2_met,"C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/answer_2_met.csv", row.names = FALSE)

###
#### Conducting tests on the final models 1.Auto.arima (3,1,1) & 2.Arima with Regressors 
###
nasa_data2 <- read.csv("C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/GLB.Ts+dSST_Processed.csv") #please select the data file
nasa_data_train <- mutate(nasa_data2, temperature_degrees = nasa_data2$Mean_Temp + 14, temperature_fahrenheit = nasa_data2$Mean_Temp + 57.2) # converting anomalies to actual temperature
nasa_subset1_train <- subset(nasa_data_train, Year >=1880) # Baseline is 1951 - 1980
nasa_data_train <- subset(nasa_subset1_train, Year <= 2018) #limiting the data to 2019 
nasa_data_test <- subset(nasa_subset1_train, Year > 2018)
nasa_ts_train <- ts(nasa_data_train$temperature_degrees,start=1880, frequency=12) 
nasa_msts_train <- msts(nasa_data_train$temperature_degrees, seasonal.periods = c(12,120))

##auto.arima(3,1,1) our model - Test
nasa_arima2_train <- auto.arima(nasa_msts_train, seasonal = TRUE, trace = TRUE, D = 1)
plot(nasa_arima2_train)
nasa_arima_pred2_train <- forecast(nasa_arima2_train, h=15, level=c(0.8,0.90))
plot(nasa_arima_pred2_train, xlab = 'Year', ylab = "Predicted Climate")
nasa_arima_pred2_train_date <- data.frame("Date" =seq(as.Date("2019/1/1"), by = "month", length.out =15))
nasa_arima_pred2_train_final <- cbind(nasa_arima_pred2_train_date,nasa_arima_pred2_train)
nasa_arima_pred2_final_data <- data.frame("Date" =nasa_arima_pred2_train_final$Date,"point_forecast" = nasa_arima_pred2_train_final$`Point Forecast`, "lo_90" = nasa_arima_pred2_train_final$`Lo 90`, "hi_90" = nasa_arima_pred2_train_final$`Hi 90`) 
nasa_arima_pred2_final_data <- mutate(nasa_arima_pred2_final_data, "Actual" = nasa_data_test$temperature_degrees)
mep_nasa_arima_pred2 <- mean(((nasa_arima_pred2_final_data$Actual - nasa_arima_pred2_final_data$point_forecast)/nasa_arima_pred2_final_data$Actual)*100)
##0.3514449
mape_nasa_arima_pred2 <- MLmetrics::MAPE(nasa_arima_pred2_final_data$point_forecast,nasa_arima_pred2_final_data$Actual)
## 0.006679404
write.csv(nasa_arima_pred2_final_data,"C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/auto_arima_pred.csv")

# ARIMA with regressors - Test
monthMatrix <- cbind(month = model.matrix(~as.factor(nasa_data_train$Month)))
monthMatrix <- monthMatrix[,-1]
colnames(monthMatrix) <- c("Feb", "Mar", "Apr", "May", "Jun","Jul","Aug","Sep","Oct","Nov", "Dec")
matrix_of_regressors <- monthMatrix
nasa_arima3_train <- auto.arima(nasa_data_train$temperature_degrees, xreg=matrix_of_regressors)
xreg.pred <- matrix_of_regressors[-c(16:31),]
nasa_arima_pred3_train <- forecast(nasa_arima3_train, h = 15, xreg = xreg.pred, level = c(0.8, 0.90))
plot(nasa_arima_pred3_train , xlab = "Year", ylab = "Predicted Climate")
nasa_arima_pred3_train_date <- data.frame("Date" =seq(as.Date("2019/1/1"), by = "month", length.out =1652))
nasa_arima_pred3_train_final <- cbind(nasa_arima_pred3_train_date,nasa_arima_pred3_train)
nasa_arima_pred3_train_final <- head(nasa_arima_pred3_train_final,15)
nasa_arima_pred3_train
nasa_arima_pred3_final_data <- data.frame("Date" =nasa_arima_pred3_train_final$Date,"point_forecast" = nasa_arima_pred3_train_final$`Point Forecast`, "lo_90" = nasa_arima_pred3_train_final$`Lo 90`, "hi_90" = nasa_arima_pred3_train_final$`Hi 90`) 
nasa_arima_pred3_final_data <- mutate(nasa_arima_pred3_final_data, "Actual" = nasa_data_test$temperature_degrees)
mep_nasa_arima_pred3 <- mean(((nasa_arima_pred3_final_data$Actual - nasa_arima_pred3_final_data$point_forecast)/nasa_arima_pred3_final_data$Actual)*100)
## 0.9566396%
mape_nasa_arima_pred3 <- MLmetrics::MAPE(nasa_arima_pred3_final_data$point_forecast,nasa_arima_pred3_final_data$Actual)
## 0.009613902
write.csv(nasa_arima_pred3_final_data,"C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/arima_regressors_pred.csv")
mape.errors <- data.frame("auto.arima(3,1,1)" = mape_nasa_arima_pred2 , "arima.regressors" = mape_nasa_arima_pred3)
write.csv(mape.errors,"C:/Users/mhais/Desktop/MMA/867/Team Assignment 1/mape_values.csv", row.names = FALSE)


####
####Additional tests
####
nasa_pred <- read.csv("C://Users//du_xf//Downloads//nasa_prediction.csv", header = TRUE)
met_pred <- read.csv("C://Users//du_xf//Downloads//met_prediction.csv", header = TRUE)
tail(nasa_pred)
tail(met_pred)

nasa_pred_2100 <- nasa_pred[958:969,]
nasadata_2019 <- nasadata_2019 + 14
nasadata_diff <- as.data.frame(cbind(nasadata_2019$Temperature, nasa_pred_2100[,2]))
names(nasadata_diff) <- c("2019", "2100")
nasadata_diff$diff <- nasadata_diff$`2100` - nasadata_diff$`2019`
mean(nasadata_diff$`2100`)
sd(nasadata_diff$`2100`)
avg_diff <- mean(nasadata_diff$diff)
std_dev_diff <- sd(nasadata_diff$diff)
CI90_U <- avg_diff + qnorm(0.95)*std_dev_diff
CI90_L <- avg_diff - qnorm(0.95)*std_dev_diff
nasadata_temp_increase <- as.data.frame(cbind(avg_diff, std_dev_diff, CI90_L, CI90_U))

met_pred_2100 <- met_pred[959:970,]
UKMet_2019 <- UKMetData_copy[(nrow(UKMetData_copy)-14):(nrow(UKMetData_copy)-3),2]
UKMet_2019$Temperature <- UKMet_2019$Temperature + 14
UKMet_diff <- as.data.frame(cbind(UKMet_2019[,1], met_pred_2100[,2])) 
names(UKMet_diff) <- c("2019", "2100")
UKMet_diff$Diff <- UKMet_diff$`2100` - UKMet_diff$`2019`
mean(UKMet_diff$`2100`)
sd(UKMet_diff$`2100`)
avg_diff <- mean(UKMet_diff$Diff)
std_dev_diff <- sd(UKMet_diff$Diff)
hist(UKMet_diff$Diff)
CI90_U <- avg_diff + qnorm(0.95)*std_dev_diff
CI90_L <- avg_diff - qnorm(0.95)*std_dev_diff
ukmet_temp_increase <- as.data.frame(cbind(avg_diff, std_dev_diff, CI90_L, CI90_U))

