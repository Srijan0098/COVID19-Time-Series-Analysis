library(tseries)
library(ggplot2)
library(zoo)
library(TTR)
library(itsmr)
library(randtests)
library(astsa)
library(lmtest)
library(forecast)
library(urca)
library(dplyr)

data <- read.csv('/home/srijan/Downloads/WHO-COVID-19-global-data.csv')
data <- data %>% filter(Country == 'India')
View(data)
death <- data[,c('Date_reported','New_deaths')]
View(death)

# Plotting the data

par(mfrow=c(1,1))
daily_deaths <- zoo(death$New_deaths)
plot(daily_deaths, type='l', xlab = "Days", ylab="New Deaths", main=" Daily Covid-19 Deaths in India")

# Smoothing by 7-day moving average

d = 7
death_ext7 <- c(rep(daily_deaths[1],3),
               daily_deaths,rep(daily_deaths[length(daily_deaths)],3))  # padding both sides
death_smooth7 <- stats::filter(death_ext7, sides=2, c(rep(1,d))/d)      # moving average
death_smooth7 <- na.omit(death_smooth7)

plot.ts(death_smooth7,xlab='Days',ylab="Daily Deaths",main="Smoothed daily death (using 7-day moving average)")

# Separating the 2 waves

wave1 <- death_smooth7[70:420]
wave2 <- death_smooth7[420:length(death_smooth7)]

plot(wave1, type='l', xlab = "Days", ylab="Daily New Deaths (Smoothed)",
     main="Smoothed (7-days) Daily New Covid-19 Deaths in India - First Wave")
plot(wave2, type='l', xlab = "Days", ylab="Daily New Deaths (Smoothed)",
     main="Smoothed (7-days) Daily New Covid-19 Deaths in India - Second Wave")


##########   ANALYSIS OF FIRST WAVE   ##########

# Checking for stationarity
acf2(wave1)

# Performing Unit root test to verify non-stationarity

summary(ur.df(wave1,type="drift"))  # CONCLUSION: Null Hypothesis is accepted i.e. 
wave1                               #             there is a unit root => the data is non-stationary

log_wave1 <- log(wave1)
plot(log_wave1, type='l', xlab = "Days", ylab="New Deaths", main=" Daily Covid-19 Deaths in India (log scaled) during First Wave")

acf2(log_wave1)
summary(ur.df(log_wave1,type="drift"))

# Performing differencing to try to make the data stationary
par(mfrow=c(1,1))
diff_log_wave1 <- diff(log_wave1,1)
plot(diff_log_wave1, type='l', xlab = "Days", ylab="Growth in Death", main=" Rate of Growth in daily deaths in India during First Wave of Covid")
abline(h = mean(diff_log_wave1), col="purple", lwd=2,lty=3)
plot(diff_log_wave1**2, type='l', xlab = "Days", ylab="Squared Growth in Deaths", main=" Squared Rate of Growth in daily deaths in India during First Wave of Covid")
acf2(diff_log_wave1)
summary(ur.df(diff_log_wave1,type="drift"))  # No unit root

par(mfrow=c(3,1))
chart.RollingMean(ts(diff_log_wave1),width=30)
chart.RollingPerformance(ts(diff_log_wave1),width=30,FUN = 'sd.annualized',main='Rolling 1-month moving average on Variance of Residuals')
chart.RollingPerformance(ts(diff_log_wave1**2),width=30,FUN = 'sd.annualized',main='Rolling 1-month moving average on square Variance of Residuals')
par(mfrow=c(1,1))


################ ARMA GARCH MODELLING ###################

chart.Histogram(diff_log_wave1,
                methods = c('add.density', 'add.normal'),
                colorset = c("steelblue", 'dark green', 'red'),main='Histogram of growth rate in death during first wave')
legend("topright",inset=0.1,legend = c("Kernel","Normal Distribution"),col=c("dark green", "red"),lty = 1:1,cex=1.1)
skewness(diff_log_wave1) # Negatively skewed
qqnorm(growth_rate2_diff1)
qqline(growth_rate2_diff1)

AICC <- function(n,k,AIC){
  return(AIC+ (2*((k)**2+k)/(n-k-1)))
}

w1_AIC_norm<-c()
w1_BIC_norm<-c()
#NIF_Big_norm <-c()
name_vec_w1<- c()

for (a in seq(0,3)){
  for (b in seq(0,3)){
    for (i in seq(1,3)){
      for (j in seq(1,3)){
        try(k <- ugarchspec(mean.model = list(armaOrder=c(a,b)),variance.model = list(model="sGARCH",garchOrder=c(i,j)),distribution.model="sstd"))
        try(l <-ugarchfit(data=diff_log_wave1, spec = k, out.sample = 20))
        try({w1_AIC_norm <- c(w1_AIC_norm,infocriteria(l)[1])
        w1_BIC_norm<- c(w1_BIC_norm,infocriteria(l)[2])
        name_vec_w1<- c(name_vec_w1,paste(a,b,i,j))})
      }
    }
    
  }
}

WAVE1_deaths<-data.frame(order=name_vec_w1, AIC=w1_AIC_norm,BIC=w1_BIC_norm)
WAVE1_deaths

for (i in seq(1,length(WAVE1_deaths$AIC))){
  WAVE1_deaths$AICC[i] <- AICC(length(diff_log_wave1),sum(as.numeric(str_split(WAVE2_deaths$order[i]," ")[[1]])),WAVE2_deaths$AIC[i])
}

write.csv(WAVE1_deaths,"/home/srijan/srijan_ds/Sem3/BDA_TSA/Covid_project/wave1_deaths.csv", row.names = FALSE)

which(WAVE1_deaths==min(WAVE1_deaths$AICC),arr.ind=TRUE)
WAVE1_deaths

# AICc----
# BEST MODEL (1,1,1,1)

k <- ugarchspec(mean.model = list(armaOrder=c(1,1)),
                 variance.model = list(model="sGARCH",garchOrder=c(1,1)),
                 distribution.model="sstd")
l <- ugarchfit(data=diff_log_wave1, spec = k, out.sample = 20)
l

forecast <-ugarchforecast(l,data=diff_log_wave1,n.ahead=20)
forecast

##################################
######### Rolling Forecast #######
##################################


fit_roll <- ugarchfit(k, data = diff_log_wave1, out.sample =20)
fore_roll <- ugarchforecast(fit_roll, n.ahead=20, n.roll=20)
fore_roll
f4 <- fore_roll@forecast
pred4 <- f4$seriesFor[,1]
sum((diff_log_wave1[(length(diff_log_wave1)-19):length(diff_log_wave1)] - pred4)^2)/20
par(mfrow=c(1,2))
plot(fore_roll,which=1)
plot(fore_roll,which=2)



# BIC----
# BEST MODEL (2,2,1,1)

k2 <- ugarchspec(mean.model = list(armaOrder=c(2,2)),
                variance.model = list(model="sGARCH",garchOrder=c(1,1)),
                distribution.model="sstd")
l2 <- ugarchfit(data=diff_log_wave1, spec = k2, out.sample = 20)
l2

forecast <-ugarchforecast(l,data=diff_log_wave1,n.ahead=20)
forecast

##################################
######### Rolling Forecast #######
##################################


fit_roll <- ugarchfit(k2, data = diff_log_wave1, out.sample =20)
fore_roll <- ugarchforecast(fit_roll, n.ahead=20, n.roll=20)
fore_roll
f5 <- fore_roll@forecast
pred5 <- f5$seriesFor[,1]
sum((diff_log_wave1[(length(diff_log_wave1)-19):length(diff_log_wave1)] - pred5)^2)/20
par(mfrow=c(1,2))
plot(fore_roll,which=1)
plot(fore_roll,which=2)

acf2(daily_deaths)



########## ANALYSIS of 2nd WAVE ##########

wave2

# Checking for stationarity
acf2(wave2)
wave2
# Performing Unit root test to verify non-stationarity

summary(ur.df(wave2,type="drift"))  # CONCLUSION:  Unit root

plot(wave2,type='l')
                         
par(mfrow=c(1,1))
log_wave2 <- log(wave2)
plot(log_wave2, type='l', xlab = "Days", ylab="New Deaths", main=" Daily Covid-19 Deaths in India (log scaled) during Second Wave")

acf2(log_wave2)
summary(ur.df(log_wave2,type="drift")) #There is unit root

# Performing differencing to try to make the data stationary
par(mfrow=c(1,1))
diff_log_wave2 <- diff(log_wave2,1)
plot(diff_log_wave2, type='l', xlab = "Days", ylab="New Deaths", main=" Rate of Growth in daily deaths in India during Second Wave of Covid")
abline(h = mean(diff_log_wave2), col="purple", lwd=2,lty=3)
plot(diff_log_wave2**2, type='l', xlab = "Days", ylab="New Deaths", main=" Squared Rate of Growth in daily deaths in India during Second Wave of Covid")

diff_log_wave2
acf2(diff_log_wave2)
summary(ur.df(diff_log_wave2,type="drift")) #No unit root

par(mfrow=c(3,1))
chart.RollingMean(ts(diff_log_wave2),width=30)
chart.RollingPerformance(ts(diff_log_wave2),width=30,main='Rolling 1-month moving average on Variance of Residuals')
chart.RollingPerformance(ts(diff_log_wave2**2),width=30,main='Rolling 1-month moving average on squared Variance of Residuals')
par(mfrow=c(1,1))

acf2(diff_log_wave2)
summary(ur.df(diff_log_wave2,type='drift') )# No unit root

########### ARMA GARCH MODELLING ###########
chart.Histogram(diff_log_wave2,
                methods = c('add.density', 'add.normal'),
                colorset = c("steelblue", 'dark green', 'red'),main='Histogram of growth rate in death during second wave')
legend("topright",inset=0.07,legend = c("Kernel","Normal Distribution"),col=c("dark green", "red"),lty = 1:1,cex=1.1)
skewness(diff_log_wave2) # Negatively skewed
qqnorm(growth_rate2_diff1)
qqline(growth_rate2_diff1)

AICC <- function(n,k,AIC){
  return(AIC+ (2*((k)**2+k)/(n-k-1)))
}

w2_AIC_norm <- c()
w2_BIC_norm <- c()
#NIF_Big_norm <-c()
name_vec_w2 <- c()

for (a in seq(0,3)){
  for (b in seq(0,3)){
    for (i in seq(1,3)){
      for (j in seq(1,3)){
        try(k2 <- ugarchspec(mean.model = list(armaOrder=c(a,b)),variance.model = list(model="sGARCH",garchOrder=c(i,j)),distribution.model="sstd"))
        try(l2 <-ugarchfit(data=diff_log_wave2, spec = k2, out.sample = 20))
        try({w2_AIC_norm <- c(w2_AIC_norm,infocriteria(l2)[1])
        w2_BIC_norm<- c(w2_BIC_norm,infocriteria(l2)[2])
        name_vec_w2<- c(name_vec_w2,paste(a,b,i,j))})
      }
    }
    
  }
}

WAVE2_deaths<-data.frame(order=name_vec_w2, AIC=w2_AIC_norm,BIC=w2_BIC_norm)
WAVE2_deaths

for (i in seq(1,length(WAVE2_deaths$AIC))){
  WAVE2_deaths$AICC[i] <- AICC(length(diff_log_wave2),sum(as.numeric(str_split(WAVE2_deaths$order[i]," ")[[1]])),WAVE2_deaths$AIC[i])
}

write.csv(WAVE2_deaths,"/home/srijan/srijan_ds/Sem3/BDA_TSA/Covid_project/wave2_deaths.csv", row.names = FALSE)

which(WAVE2_deaths==min(WAVE2_deaths$AICC),arr.ind=TRUE)
WAVE2_deaths

#Best : (1,1,1,1)

k3 <- ugarchspec(mean.model = list(armaOrder=c(1,1)),
                variance.model = list(model="sGARCH",garchOrder=c(1,1)),
                distribution.model="sstd")
l3 <- ugarchfit(data=diff_log_wave2, spec = k2, out.sample = 20)
l3

forecast <-ugarchforecast(l3,data=diff_log_wave2,n.ahead=20)
forecast

##################################
######### Rolling Forecast #######
##################################


fit_roll3 <- ugarchfit(k3, data = diff_log_wave2, out.sample =20)
fore_roll3 <- ugarchforecast(fit_roll3, n.ahead=20, n.roll=20)
fore_roll3
f6 <- fore_roll3@forecast
pred6 <- f6$seriesFor[,1]
sum((diff_log_wave2[(length(diff_log_wave2)-19):length(diff_log_wave2)] - pred6)^2)/20
par(mfrow=c(1,2))
plot(fore_roll3,which=1)
plot(fore_roll3,which=2)

