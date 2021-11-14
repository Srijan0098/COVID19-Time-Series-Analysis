library(tseries)
library(ggplot2)
library(zoo)
library(itsmr)
library(TTR)
library(randtests)
library(astsa)
library(lmtest)
library(forecast)
library(urca)
library(dplyr)
library(quantmod)
library(parallel)
library(rugarch)
library(PerformanceAnalytics)
library(MuMIn)

par(mfrow=c(1,1))

data <- read.csv("/home/srijan/Downloads/covid-india.csv")
View(data)
new_cases <- zoo(data$new_cases,data$date)
plot(new_cases, type='l', xlab = "Days", ylab="Daily New Cases",
     main=" Daily New Covid-19 Cases in India")

# Smoothing the data using 7-day moving average

d = 7
cases_ext7 <- c(rep(new_cases[1],3),
                new_cases,rep(new_cases[length(new_cases)],3))    # padding both sides
cases_smooth7 <- stats::filter(cases_ext7, sides=2, c(rep(1,d))/d)      # moving average
cases_smooth7 <- cases_smooth7[-c(1,2,3,length(cases_smooth7)-2,length(cases_smooth7)-1,length(cases_smooth7))]
plot(cases_smooth7, type='l', xlab = "Days", ylab="Daily New Cases (Smoothed)",
     main="Smoothed (7-days) Daily New Covid-19 Cases in India ")

# Separating the 2 waves

wave1 <- cases_smooth7[1:300]
wave2 <- cases_smooth7[300:length(cases_smooth7)]

plot(wave1, type='l', xlab = "Days", ylab="Daily New Cases (Smoothed)",
     main="Smoothed (7-days) Daily New Covid-19 Cases in India - First Wave")
plot(wave2, type='l', xlab = "Days", ylab="Daily New Cases (Smoothed)",
     main="Smoothed (7-days) Daily New Covid-19 Cases in India - Second Wave")

##########   ANALYSIS OF FIRST WAVE   ##########

# Defining growth rate (log-difference of smoothed covid cases)

acf2(wave1)
summary(ur.df(wave1,type="drift"))

log_wave1 <- log(wave1)
acf2(log_wave1)
summary(ur.df(log_wave1,type="drift"))
growth_rate1 <- diff(log_wave1)
plot(growth_rate1, type='l', xlab = "Days", ylab="Growth Rate",
     main="Growth Rate of Daily Covid Cases in India during First Wave")
abline(h = mean(growth_rate1), col="purple", lwd=2,lty=3)
plot(growth_rate1**2, type='l', xlab = "Days", ylab="Growth Rate Squared",
     main="Squared Growth Rate of Daily Covid Cases in India during First Wave")

par(mfrow=c(3,1))
chart.RollingMean(ts(growth_rate1),width = 30)
chart.RollingPerformance(ts(growth_rate1),width = 30,main = '1-month rolling average on Variance of growth rate')
chart.RollingPerformance(ts(growth_rate1**2),width = 30,main = '1-month rolling average on squared Variance of growth rate')
acf2(growth_rate1)
summary(ur.df(growth_rate1,type="drift"))
par(mfrow=c(1,1))

chart.Histogram(growth_rate1,
                xlab = "Growth rate",
                methods = c('add.density', 'add.normal'),
                colorset = c("steelblue", 'dark green', 'red'),
                main='Histogram of growth rate')
legend("topright",inset=0.1,legend = c("Kernel","Normal Distribution"),col=c("dark green", "red"),lty = 1:1,cex=1.1)

skewness(growth_rate1)
qqnorm(growth_rate1)
qqline(growth_rate1)

shapiro.test(growth_rate1)
Box.test(diff(growth_rate1),type =  "Ljung-Box")
auto.arima(growth_rate1)

# GARCH model

AICC <- function(n,k,AIC){
  return(AIC+(2*((k)**2+k)/(n-k-1)))
}

NIF_AIC_norm2<-c()
NIF_BIC_norm2<-c()
NIF_Big_norm2 <-c()
name_vec2<- c()

for (a in seq(0,3)){
  for (b in seq(0,3)){
    for (i in seq(1,4)){
      for (j in seq(1,4)){
        try(k2 <- ugarchspec(mean.model = list(armaOrder=c(a,b)),variance.model = list(model="sGARCH",garchOrder=c(i,j)),distribution.model="sstd"))
        try(l2 <-ugarchfit(data=growth_rate1, spec = k2, out.sample = 20))
        try({NIF_AIC_norm2 <- c(NIF_AIC_norm2,infocriteria(l2)[1])
        NIF_BIC_norm2<- c(NIF_BIC_norm2,infocriteria(l2)[2])
        name_vec2<- c(name_vec2,paste(a,b,i,j))})
      }
    }
    
  }
}
NIF_AIC_norm2
WAVE1_results<-data.frame(order=name_vec2, AIC=NIF_AIC_norm2,BIC=NIF_BIC_norm2)
WAVE1_results
library(stringr)
for (i in seq(1,length(WAVE1_results$AIC))){
  WAVE1_results$AICC[i] <- AICC(length(growth_rate1),sum(as.numeric(str_split(WAVE1_results$order[i]," ")[[1]])),WAVE1_results$AIC[i])
  
}
WAVE1_results
which(WAVE1_results==min(WAVE1_results$AICC),arr.ind=TRUE)

write.csv(WAVE1_results,"/home/srijan/srijan_ds/Sem3/BDA_TSA/Covid_project/wave1_results.csv", row.names = FALSE)

# AICc----
# BEST MODEL (1,1,1,1)

k <- ugarchspec(mean.model = list(armaOrder=c(1,1)),
                variance.model = list(model="sGARCH",garchOrder=c(1,1)),
                distribution.model="sstd")
l <- ugarchfit(data=growth_rate1, spec = k, out.sample = 20)
l

forecast <-ugarchforecast(l,data=na.omit(growth_rate1),n.ahead=20)
forecast

##################################
######### Rolling Forecast #######
##################################

fit_roll <- ugarchfit(k, data = growth_rate1, out.sample =20)
fore_roll <- ugarchforecast(fit_roll, n.ahead=20, n.roll=20)
fore_roll
f <- fore_roll@forecast
pred <- f$seriesFor[,1]
sum((growth_rate1[(length(growth_rate1)-19):length(growth_rate1)] - pred)^2)/20
par(mfrow=c(1,2))
plot(fore_roll,which=1)
plot(fore_roll,which=2)

#BIC----
#BEST MODEL (3,3,1,1)

k3 <- ugarchspec(mean.model = list(armaOrder=c(3,3)),
                variance.model = list(model="sGARCH",garchOrder=c(1,1)),
                distribution.model="sstd")
l3 <- ugarchfit(data=growth_rate1, spec = k3, out.sample = 20)
coef(l3)

forecast3 <-ugarchforecast(l3,data=na.omit(growth_rate1),n.ahead=20)
forecast3

#################################
######### Rolling Forecast #######
#################################

fit_roll3 <- ugarchfit(k3, data = growth_rate1, out.sample =20)
fore_roll3 <- ugarchforecast(fit_roll3, n.ahead=20, n.roll=20)
fore_roll3
f2 <- fore_roll3@forecast
pred2 <- f2$seriesFor[,1]
sum((growth_rate1[(length(growth_rate1)-19):length(growth_rate1)] - pred2)^2)/20
par(mfrow=c(1,2))
plot(fore_roll3,which=1)
plot(fore_roll3,which=2)
l3






growth_rate1

par(mfrow=c(1,1))



#########   ANALYSIS OF SECOND WAVE   ##########

# Defining growth rate (log-difference of smoothed covid cases)

acf2(wave2)
summary(ur.df(wave2,type="drift"))

log_wave2 <- log(wave2)
growth_rate2 <- diff(log_wave2)
plot(growth_rate2, type='l', xlab = "Days", ylab="Growth Rate",
     main="Growth Rate of Daily Covid Cases in India during Second Wave")
plot(growth_rate2**2, type='l', xlab = "Days", ylab="Growth Rate",
     main="Squared growth Rate of Daily Covid Cases in India during Second Wave")
acf2(growth_rate2)
summary(ur.df(growth_rate2,type="drift"))

# Performing differencing to try to make the data stationary

growth_rate2_diff1 <- diff(growth_rate2,1)
plot(growth_rate2_diff1, type='l',xlab = "Days", ylab="Growth Rate Residuals",
     main="Residuals of Growth Rate of Daily Covid Cases in India during Second Wave")
abline(h = mean(growth_rate2_diff1), col="purple", lwd=2,lty=3)
plot(growth_rate2_diff1**2, type='l',xlab = "Days", ylab="Growth Rate Residual Squared",
     main="Squared Residuals of Growth Rate of Daily Covid Cases in India during Second Wave")

#chart.RollingMean(ts(growth_rate2_diff1),width = 30)
par(mfrow=c(3,1))
chart.RollingMean(ts(growth_rate2_diff1),width = 30)
chart.RollingPerformance(ts(growth_rate2_diff1),width = 30,FUN = 'sd.annualized',main='Rolling 1-month moving average on Variance of Residuals')
chart.RollingPerformance(ts(growth_rate2_diff1**2),width = 30,FUN = 'sd.annualized',main='Rolling 1-month moving average on square Variance of Residuals')

acf2(growth_rate2_diff1)
summary(ur.df(growth_rate2_diff1,type="drift"))

par(mfrow=c(1,1))
chart.Histogram(growth_rate2_diff1,
                methods = c('add.density', 'add.normal'),
                colorset = c("steelblue", 'dark green', 'red'),main='Histogram of growth rate')
legend("topleft",inset=0.1,legend = c("Kernel","Normal Distribution"),col=c("dark green", "red"),lty = 1:1,cex=1.1)
skewness(growth_rate2_diff1)
qqnorm(growth_rate2_diff1)
qqline(growth_rate2_diff1)

shapiro.test(growth_rate2_diff1)
Box.test(growth_rate2_diff1,type =  "Ljung-Box")


# GARCH model

AICC <- function(n,k,AIC){
  return(AIC+ (2*((k)**2+k)/(n-k-1)))
}

NIF_AIC_norm<-c()
NIF_BIC_norm<-c()
NIF_Big_norm <-c()
name_vec<- c()

for (a in seq(0,3)){
  for (b in seq(0,3)){
    for (i in seq(1,4)){
      for (j in seq(1,4)){
        try(k <- ugarchspec(mean.model = list(armaOrder=c(a,b)),variance.model = list(model="sGARCH",garchOrder=c(i,j)),distribution.model="sstd"))
        try(l <-ugarchfit(data=growth_rate2_diff1, spec = k, out.sample = 20))
        try({NIF_AIC_norm <- c(NIF_AIC_norm,infocriteria(l)[1])
        NIF_BIC_norm<- c(NIF_BIC_norm,infocriteria(l)[2])
        name_vec<- c(name_vec,paste(a,b,i,j))})
      }
    }
    
  }
}

WAVE2_results<-data.frame(order=name_vec, AIC=NIF_AIC_norm,BIC=NIF_BIC_norm)
WAVE2_results

for (i in seq(1,length(WAVE2_results$AIC))){
  WAVE2_results$AICC[i] <- AICC(length(growth_rate2_diff1),sum(as.numeric(str_split(WAVE2_results$order[i]," ")[[1]])),WAVE2_results$AIC[i])
}
sum(as.numeric(str_split(WAVE2_results$order[1]," ")[[1]]))
write.csv(WAVE2_results,"/home/srijan/srijan_ds/Sem3/BDA_TSA/Covid_project/wave2_results.csv", row.names = FALSE)

which(WAVE2_results==min(WAVE2_results$AICC),arr.ind=TRUE)
WAVE2_results

#Best : (1,0,1,1)

k2 <- ugarchspec(mean.model = list(armaOrder=c(1,0)),
                variance.model = list(model="sGARCH",garchOrder=c(1,1)),
                distribution.model="sstd")
l2 <- ugarchfit(data=growth_rate2_diff1, spec = k2, out.sample = 20)
l2

forecast2 <-ugarchforecast(l2,data=growth_rate2_diff1,n.ahead=20)
forecast2

##################################
######### Rolling Forecast #######
##################################


fit_roll2 <- ugarchfit(k2, data = growth_rate2_diff1, out.sample =20)
fore_roll2 <- ugarchforecast(fit_roll2, n.ahead=20, n.roll=20)
fore_roll2
f3 <- fore_roll2@forecast
pred3 <- f3$seriesFor[,1]
sum((growth_rate2_diff1[(length(growth_rate2_diff1)-19):length(growth_rate2_diff1)] - pred3)^2)/20
par(mfrow=c(1,2))
plot(fore_roll2,which=1)
plot(fore_roll2,which=2)
