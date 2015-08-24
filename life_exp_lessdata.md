
```{r}
library(lme4)
library(rjags)

setwd('~/Desktop/Mixed models/bayesianModel')

expec<-read.csv("life_exp_2.csv", header=TRUE)
expec1<-na.omit(expec)
x1<-log10(expec1$pop)
x2<-log10(expec1$gdpPercap)
y<-expec1$lifeExp

ypred<-NULL
pred.data<-NULL
############ LM METHOD
model_lm<-function()
{
m1<-lm(y~x1+x2)

ypred<- coef(m1)[1]+coef(m1)[2]*x1+coef(m1)[3]*x2
resid<-y-ypred

par(mfrow=c(2,3), mar = .1+ c(2,2,2,2))
plot(ypred, y,  xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="lm method")
abline(0,1)
plot(ypred, resid, xlim=c(range(y)), xlab="Predicted", ylab="Residuals", main="lm method")
abline(h=0)

hist(resid)

##predict for 2007
e<-read.csv("life_exp.csv", header=TRUE)
e1<-na.omit(e)


data.subset<-e1[e1$year==1997,]
pop<-log10(data.subset$pop)
gdp<-log10(data.subset$gdpPercap)

pred.data<-NULL

for(i in 1:nrow(data.subset))
{

pred.data[i]<-coef(m1)[1]+ coef(m1)[2]*pop[i]+coef(m1)[3]*gdp[i]
}

y.tmp<-data.subset$lifeExp
plot(pred.data, y.tmp, xlim=c(range(y.tmp)), xlab="Predicted", ylab="Observed", main="lm method 1997" )
abline(0,1)

}




model_lmer<-function()
{
m2<-lmer(y~ x1+x2+(1+x2|country), data=expec1, na.action="na.omit", REML=FALSE)

temp<-data.frame(coef(m2)[1])

ypred<-predict(m2)
resid<-y-ypred

par(mfrow=c(2,3), mar = .1+ c(2,2,2,2))
plot(ypred, y,  xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="lmer method")
abline(0,1)
plot(ypred, resid, xlim=c(range(y)), xlab="Predicted", ylab="Residuals", main="lmer method")
abline(h=0)

hist(resid)


##predict for 1997
e<-read.csv("life_exp.csv", header=TRUE)
e1<-na.omit(e)


data.subset<-e1[e1$year==1997,]
pop<-log10(data.subset$pop)
gdp<-log10(data.subset$gdpPercap)
pred.data<-NULL

for(i in 1:nrow(data.subset))
{
p<-data.subset$country[i]
pred.data[i]<-temp[p, 1]+ temp[p,2]*pop[i]+temp[p,3]*gdp[i]
}

y.tmp<-data.subset$lifeExp
plot(pred.data,y.tmp, xlim=c(range(y.tmp)),xlab="Predicted", ylab="Observed", main="lmer method 1997" )
abline(0,1)

}

expec<-read.csv("life_exp_2.csv", header=TRUE)
expec1<-na.omit(expec)
x1<-log10(expec1$pop)
x2<-log10(expec1$gdpPercap)
y<-expec1$lifeExp










###########################################################################################



bayesian_model_noninformativePrior<-function()
{
cat('model {
for (i in 1:n){
y[i] ~ dnorm (y.hat[i], tau.y)
y.hat[i] <- a[country[i]] + c*x1[i]+b[country[i]]*x2[i]
}
c ~ dnorm (0, .0001)
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif (0, 100)
for (j in 1:J){
a[j] ~ dnorm (mu.a, tau.a)
b[j] ~ dnorm (mu.b, tau.b)
}
mu.a ~ dnorm (0, .0001)
mu.b ~ dnorm (0, .0001)
tau.a <- pow(sigma.a, -2)
sigma.a ~ dunif (0, 100)
tau.b <- pow(sigma.b, -2)
sigma.b ~ dunif (0, 100)
}', file={f <- tempfile()})

  n = 150
 J = 25
#W <- diag (2)
jags_n <-jags.model(f, data = list('y' = y, 'x1' = x1,'x2' = x2, "n" = n, "J"=J,'country'=country),
n.chains = 3, n.adapt = 100)

#as.mcmc(jags_new)

posterior <- coda.samples(jags_n, c("a", "b","c", "mu.a", "mu.b", "sigma.y", "sigma.a", "sigma.b") ,n.iter=10000)
#summary(posterior)
#plot(posterior)

burn.in <- 1000
#summary(posterior, start = burn.in)$quantiles[, c(3,1, 5)]



a.posterior.mean <- summary(posterior, start = burn.in)$statistics[1:25, "Mean"]

b.posterior.mean <- summary(posterior, start = burn.in)$statistics[26:50, "Mean"]
c.posterior.mean <- summary(posterior, start = burn.in)$statistics["c", "Mean"]
sigmay.posterior.mean<-summary(posterior, start = burn.in)$statistics["sigma.y", "Mean"]


ymean<-a.posterior.mean[country]+c.posterior.mean*x1+b.posterior.mean[country]*x2
for(i in 1:n)
ypred[i]<-rnorm(1, ymean[i], sd=sqrt(sigmay.posterior.mean))

yresid=y-ypred

par(mfrow=c(2,3), mar = .1+ c(2,2,2,2))
plot(ypred, y, xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="Bayesian method")
abline(0,1)
plot(ypred, yresid, xlim=c(range(y)), xlab="Predicted", ylab="Residuals", main="Bayesian method")
abline(h=0)

hist(yresid)

#predicion for 2007


e<-read.csv("life_exp.csv", header=TRUE)
e1<-na.omit(e)


data.subset<-e1[e1$year==1997,]

pop<-log10(data.subset$pop)
gdp<-log10(data.subset$gdpPercap)
mean.data<-NULL

for(i in 1:nrow(data.subset))
{
  p<-data.subset$country.r[i]
  mean.data[i]<-a.posterior.mean[p]+c.posterior.mean*pop[i]+b.posterior.mean[p]*gdp[i]
  pred.data[i]<-rnorm(1, mean.data[i],sd=sqrt(sigmay.posterior.mean) )
}


y.tmp<-data.subset$lifeExp
plot(pred.data,y.tmp, xlim=c(range(y.tmp)),xlab="Predicted", ylab="Observed", main="Bayesian method 1997" )
abline(0,1)

}


expec<-read.csv("life_exp_2.csv", header=TRUE)
expec1<-na.omit(expec)
x1<-log10(expec1$pop)
x2<-log10(expec1$gdpPercap)
y<-expec1$lifeExp
country<-factor(expec1$country.r)




####################################################################################################

bayesian_model_informativePrior<-function()
{
cat('model {
for (i in 1:n){
y[i] ~ dnorm(y.hat[i], tau.y)
y.hat[i] <- a[country[i]] + c*x1[i]+ b[country[i]]*x2[i]
}
c ~ dnorm (40, 1)
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif (0, 2)
for (j in 1:J){
a[j] ~ dnorm (a.hat[j], tau.a)
b[j] ~ dnorm (b.hat[j], tau.b)
a.hat[j] <- mu.a
b.hat[j] <- mu.b
}
mu.a ~ dnorm (-217, .00859)
mu.b ~ dnorm (1.82, .155)
tau.a <- pow(sigma.a, -2)
tau.b <- pow(sigma.b, -2)
sigma.a ~ dunif (40, 60)
sigma.b ~ dunif (9, 14)
}', file={f <- tempfile()})

  n = 150
 J = 25
#W <- diag (2)
jags_n <-jags.model(f, data = list('y' = y, 'x1' = x1,'x2' = x2, "n" = n, "J"=J,'country'=country),
n.chains = 3, n.adapt = 100)

#as.mcmc(jags_new)

posterior <- coda.samples(jags_n, c("a", "b","c", "mu.a", "mu.b", "sigma.y", "sigma.a", "sigma.b" ) ,n.iter=10000)
#summary(posterior)
#plot(posterior)

burn.in <- 2000
summary(posterior, start = burn.in)



a.posterior.mean <- summary(posterior, start = burn.in)$statistics[1:25, "Mean"]

b.posterior.mean <- summary(posterior, start = burn.in)$statistics[26:50, "Mean"]
c.posterior.mean <- summary(posterior, start = burn.in)$statistics["c", "Mean"]
sigmay.posterior.mean<-summary(posterior, start = burn.in)$statistics["sigma.y", "Mean"]


ymean<-a.posterior.mean[country]+c.posterior.mean*x1+b.posterior.mean[country]*x2
for(i in 1:n)
ypred[i]<-rnorm(1, ymean[i], sd=sqrt(sigmay.posterior.mean))


yresid=y-ypred

par(mfrow=c(2,3), mar = .1+ c(2,2,2,2))
plot(ypred, y, xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="Bayesian method")
abline(0,1)
plot(ypred, yresid, xlim=c(range(y)), xlab="Predicted", ylab="Residuals", main="Bayesian method")
abline(h=0)

hist(yresid)

#predicion for 2007


e<-read.csv("life_exp.csv", header=TRUE)
e1<-na.omit(e)


data.subset<-e1[e1$year==1997,]

pop<-log10(data.subset$pop)
gdp<-log10(data.subset$gdpPercap)
mean.data<-NULL

for(i in 1:nrow(data.subset))
{
  p<-data.subset$country.r[i]
  mean.data[i]<-a.posterior.mean[p]+c.posterior.mean*pop[i]+b.posterior.mean[p]*gdp[i]
  pred.data[i]<-rnorm(1, mean.data[i],sd=sqrt(sigmay.posterior.mean) )
}


y.tmp<-data.subset$lifeExp
plot(pred.data,y.tmp, xlim=c(range(y.tmp)),xlab="Predicted", ylab="Observed", main="Bayesian method 1997" )
abline(0,1)
}

expec<-read.csv("life_exp_2.csv", header=TRUE)
expec1<-na.omit(expec)
x1<-log10(expec1$pop)
x2<-log10(expec1$gdpPercap)
y<-expec1$lifeExp
country<-factor(expec1$country.r)

model_lm()
model_lmer()
bayesian_model_noninformativePrior()
bayesian_model_informativePrior()












##################################################################
