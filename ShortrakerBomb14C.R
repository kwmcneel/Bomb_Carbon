############################################################### 
# Non-linear model of 14C rise 
# Kevin McNeel
###############################################################

library(nlme)
library(lattice)
library(AICcmodavg)
library(nlstools)
library(ggplot2)

### Import data:
rm(X152_bomb_data)
bomb <- read_csv("152_bomb_data.csv")
str(bomb)
# D14C-at-age by year

### Fit model:
Carbon <- function(YMIN, K, R, S, M, Year) {
	YMIN+K*exp(R*(M-Year)+S^2*R^2/2)*pnorm(Year,mean = M+S^2*R,sd=S)
}

### Fit Carbon model across all years:
START <- c(YMIN = -106, K = 185, R = 0, S=2.5, M = 1963)
fit.HB <- nls(D14C ~ Carbon(YMIN, K, R, S, M, Year), data=bomb,subset=Species!="SR", start=START)

summary(fit.HB)
cf1 <- coef(fit.HB)
round(cf1,3);

# Visualize the fitted model:
plot(D14C ~ Year, data=bomb,ylab=expression(paste(Delta^14,"C")),pch=c(3),cex=0)
points(D14C ~ Year, data=bomb, subset=Species==c("HB"),col="grey",pch=2,cex=.5)
points(D14C ~ Year, data=bomb, subset=Species==c("YE"),col="grey",pch=3,cex=.5)
points(D14C ~ Year, data=bomb, subset=Species==c("SR"),col="black",pch=20,cex=1)
#error bars
SRbomb<-subset(bomb,Species=="SR")
Year<-SRbomb$Year
D14C<-SRbomb$D14C
SE14<-SRbomb$SE14
arrows(Year, D14C-SE14, Year, D14C+SE14, length=0.05, angle=90, code=3)
curve(Carbon(cf1[1],cf1[2],cf1[3],cf1[4],cf1[5], x),1940, 1990, col="grey", lwd=2, add=T)
legend("bottomright",c("Halibut","Yelloweye","Shortraker"),pch=c(2,3,20),col=c("grey","grey","black"),bty = "n")

plot(fit.all)

#Predict values for shortraker
predict(object=fit.HB,newdata=data.frame(Year=Year))

##Monte Carlo Condifence Interval 
#install.packages("propagate")
library("propagate")
#pred<-data.frame(predictNLS(fit.HB, newdata=data.frame(Year=seq(1950,1975,by=1))))

  #summary: The mean/error estimates obtained from first-/second-order Taylor expansion and Monte Carlo simulation, together with calculated confidence/prediction intervals based on asymptotic normality. 
  #prop: the complete output from propagate for each value in newdata.
pred<-predictNLS(fit.HB, newdata=data.frame(Year=seq(1950,1975,by=1)))
predsum<-pred$summary
predsum$Year<-seq(1950,1975,by=1)
fitmodbomb<-fitted(fit.HB)
## Plot mean and confidence values from first-/second-order 
## Taylor expansion and Monte Carlo simulation.
plot(bomb$Year, bomb$D14C)
lines(Modbomb$Year, fitted(fit.HB), lwd = 2, col = 1)
points(Year, pred$summary[, 1], col = 2, pch = 16)
lines(Year, pred$summary[, 5], col = 2)
lines(Year, pred$summary[, 6], col = 2)
lines(Year, pred$summary[, 11], col = 4)
lines(Year, pred$summary[, 12], col = 4)

??predictNLS
# Diagnostics:
plot(fit.all)   # residual against fitted values
plot(fit.all, resid(.) ~ fitted(.) | Year)		# Residuals against fitted by year
plot(fit.all, resid(.) ~ Year | Year)		      # Residuals against age by year
plot(fit.all, Year ~ resid(.))						# Boxplots of residuals by year
plot(fit.all, factor(Year) ~ resid(.))           # Boxplots of residuals by age
# Notice that there appears to be some model mis-specification. For example, residuals
# at younger ages are mostly negative before becoming positive on average at ages 7-9
# In other words the Carbon model overestimates length at younger ages and underestimates 
# age at some intermediate ages. To capture the true trend, a more complex growth model
# would be required, for example the 4-parameter Gompertz model 

plot(fit1, factor(Year) ~ resid(.) | Year)    # Boxplots of residuals by age and year


##########mulitplot

#####Set starting parameters
( svGen <- lapply(START,rep,3) )#####replicated to serve are starting values for the General Model with

####General model
bombGen <- function(YMIN, K, R, S, M, Year, Species) {
	YMIN[Species]+K[Species]*exp(R[Species]*(M[Species]-Year)+S[Species]^2*R[Species]^2/2)*pnorm(Year,mean = M[Species]+S[Species]^2*R[Species],sd=S[Species])
}

fitGen <- nls(D14C ~ bombGen(YMIN, K, R, S, M, Year, Species),data=bomb,start=svGen)
summary(fitGen)

#### Testing assumptions
residPlot(fitGen)
hist(residuals(fitGen),main="")

##### Subset models
bomb1M <- function(YMIN, K, R, S, M, Year, Species) {
	YMIN+K*exp(R*(M[Species]-Year)+S^2*R^2/2)*pnorm(Year,mean = M[Species]+S^2*R,sd=S)
}
sv1M <- mapply(rep,START,c(1,1,1,1,3))
bombCom <- function(YMIN, K, R, S, M, Year) {
	YMIN+K*exp(R*(M-Year)+S^2*R^2/2)*pnorm(Year,mean = M+S^2*R,sd=S)
}


##### Fitting subset models
fitM <- nls(D14C ~ bomb1M(YMIN, K, R, S, M, Year, Species),data=bomb,start=sv1M)
fitCom <- nls(D14C ~ bombCom(YMIN, K, R, S, M, Year),data=bomb,start=START)
summary(fitM)
#### Testing the models ####
anova(fitM,fitGen)
anova(fitCom,fitGen)
anova(fitCom,fitM)
AIC(fitGen,fitM,fitCom)
AICc(fitGen,return.K = FALSE, second.ord = TRUE,nobs = NULL)
AICc(fitM,return.K = FALSE, second.ord = TRUE,nobs = NULL)
AICc(fitCom,return.K = FALSE, second.ord = TRUE,nobs = NULL)
cfgen <- coef(fitGen)
cfM <- coef(fitM)
cfCom <- coef(fitCom)

library(ggplot2)
E<-ggplot()+
	geom_point(data=bomb,aes(x=Year,y=D14C),size=1)+
	ylab(expression(paste(Delta^14,"C")))+
	geom_line(data=Eland,aes(x=YEAR,y=(Pollock)),size=1,linetype=2)+
	xlab("Year")+
	theme(panel.background = element_blank())+
	theme(axis.title.y=element_text(vjust=1))+

par(mfrow=c(1,2))
plot(D14C ~ Year,col=bomb$Species, data=bomb,ylab=expression(paste(Delta^14,"C")))
curve(Carbon(cfM[1],cfM[2],cfM[3],cfM[4],cfM[5], x),1940, 1990, col="black", lwd=1, add=T)
curve(Carbon(cfM[1],cfM[2],cfM[3],cfM[4],cfM[6], x),1940, 1990, col=2, lwd=2, add=T)
curve(Carbon(cfM[1],cfM[2],cfM[3],cfM[4],cfM[7], x),1940, 1990, col="green", lwd=1, add=T)
legend("bottomright",c("Halibut","Yelloweye","Shortraker"),lty=1,col=c("black","green","red"))

plot(D14C ~ Year,col=bomb$Species, data=bomb,ylab=expression(paste(Delta^14,"C")),main="All Species Combined")
curve(Carbon(cfCom[1],cfCom[2],cfCom[3],cfCom[4],cfCom[5], x),1940, 1990, col=1, lwd=2, add=T)
legend("bottomright",c("Halibut","Yelloweye","Thornyhead"),pch=1,col=c("black","green","red"))

par(mfrow=c(1,1))
sumGen<-summary(fitGen)
sumM<-summary(fitM)
sumCom<-summary(fitCom)
coefM<-round(sumM$coefficients,3)
CIl<-coefM[,"Estimate"]-(1.96*coefM[,"Std. Error"])
CIH<-coefM[,"Estimate"]+(1.96*coefM[,"Std. Error"])
CIl[13:15];CIH[13:15]
prof <- profile(fitM)
confint(prof, level=0.95)
prof <- profile(fit.HB)
confint(prof, level=0.95)
prof <- profile(fitGen)
confint(prof, level=0.95)
prof <- profile(fit.ST,alpha = 0.5)
confint(prof, level=0.8)
prof <- profile(fitCom)
confint(prof, level=0.5)
?profile
?confint
prof <- profile(fitCom)
confint(prof, level=0.95)
bomb$fitM<-fitted(fitM)
bomb1=bomb
opar <- par(mfrow = c(1,3), oma = c(1.1, 0, 1.1, 0), las = 1)
plot(prof, conf = c(95, 90, 80, 50)/100)
plot(prof, conf = c(95, 90, 80, 50)/100, absVal = FALSE)
mtext("Confidence intervals based on the profile sum of squares",
      side = 3, outer = TRUE)
mtext("BOD data - confidence levels of 50%, 80%, 90% and 95%",
      side = 1, outer = TRUE)
par(opar)
#####Bootstrap########

#mbootTypical <- nlsBoot(fit.ST,niter=1000) # niter should be nearer 1000
#mbootST <- nlsBoot(fit.ST,niter=1000) # niter should be nearer 1000
#ests <- mbootTypical$coefboot
#quantile(ests,c(0.025,0.975))

Bomb.bf <- function(resids, i) {
	bomb1$D14C <- bomb$fitM + resids[i]
	boot.fit <- nls(D14C ~ bomb1M(YMIN, K, R, S, M, Year, Species), data = bomb, 
		start = sv1M)
	coef(boot.fit)
}
BH.boot <- Bomb.bf(resid(boot.fit), 999)





