
```{r}
library(rjags)
library(lme4)

data=read.table("http://bayes.bgsu.edu/multilevel/chapter12/radon.txt", sep="\t",header=TRUE)

y = data$y
x = data$x
county = data$county
u = data$u.full



#classical method
m1<-lm(y~x)

par(mfrow=c(2,2))
plot(fitted(m1), y, xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="lmer method")
abline(0,1)
plot(fitted(m1), resid(m1), xlab="Predicted", ylab="Residuals", main="lmer method")
abline(h=0)
hist(resid(m1))


cat('model {
for (i in 1:n){
y[i] ~ dnorm (y.hat[i], tau.y)
y.hat[i] <- a + b*x[i]
}
b ~ dnorm (0, 1)
a ~ dnorm (0, 1)
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif (0, 10)
}', file={f <- tempfile()})

 n = 919
 J = 85

jags_new <- jags.model(f, data = list('y' = y, 'x' = x, "n" = n),
n.chains = 1, n.adapt = 100)

as.mcmc(jags_new)

posterior <- coda.samples(jags_new, c("a", "b", "sigma.y"), n.iter=5000)
summary(posterior)
plot(posterior)

burn.in <- 1000
summary(posterior, start = burn.in)$quantiles[, c(3,1, 5)]

a.posterior.mean <- summary(posterior, start = burn.in)$statistics["a", 
    "Mean"]
b.posterior.mean <- summary(posterior, start = burn.in)$statistics["b", 
    "Mean"]
sigma.y.posterior.mean <- summary(posterior, start = burn.in)$statistics["sigma.y", 
    "Mean"]

ypred<-a.posterior.mean+b.posterior.mean*x
yresid=y-ypred

par(mfrow=c(2,3))
plot(ypred, y, xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="Bayesian method")
abline(0,1)
plot(ypred, yresid, xlab="Predicted", ylab="Residuals", main="Bayesian method")
abline(h=0)

hist(yresid)

###################################################################################

#mixed effect model

m2<-lmer(y~x+(1+x|county), REML=FALSE)

par(mfrow=c(2,3))
plot(fitted(m2), y, xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="lmer method")
abline(0,1)
plot(fitted(m2), resid(m2), xlab="Predicted", ylab="Residuals", main="lmer method")
abline(h=0)
hist(resid(m2))


#BAYESIAN METHOD

cat('model {
for (i in 1:n){
y[i] ~ dnorm (y.hat[i], tau.y)
y.hat[i] <- a[county[i]] + b*x[i]
}
b ~ dnorm (0, .0001)
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif (0, 100)
for (j in 1:J){
a[j] ~ dnorm (mu.a, tau.a)
}
mu.a ~ dnorm (0, .0001)
tau.a <- pow(sigma.a, -2)
sigma.a ~ dunif (0, 100)
}', file={f <- tempfile()})

 n = 919
 J = 85

jags_new <- jags.model(f, data = list('y' = y, 'x' = x, "n" = n, "J"=J,'county'=county),
n.chains = 1, n.adapt = 100)

#as.mcmc(jags_new)

posterior <- coda.samples(jags_new, c("a", "b", "mu.a", "sigma.y", "sigma.a"), n.iter=5000)
summary(posterior)
plot(posterior)

burn.in <- 1000
#summary(posterior, start = burn.in)$quantiles[, c(3,1, 5)]



aposterior.mean <- summary(posterior, start = burn.in)$statistics[1:85, "Mean"]

b.posterior.mean <- summary(posterior, start = burn.in)$statistics["b", 
    "Mean"]
sigma.y.posterior.mean <- summary(posterior, start = burn.in)$statistics["sigma.y", 
    "Mean"]
}

ypred<-aposterior.mean[county]+b.posterior.mean*x

yresid=y-ypred

par(mfrow=c(2,3))
plot(ypred, y, xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="Bayesian method")
abline(0,1)
plot(ypred, yresid, xlab="Predicted", ylab="Residuals", main="Bayesian method")
abline(h=0)

hist(yresid)

##########################################################################################
#using t-distribution

#BAYESIAN METHOD

cat('model {
for (i in 1:n){
y[i] ~ dnorm (y.hat[i], tau.y)
y.hat[i] <- a[county[i]] + b[county[i]]*x[i]
}
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif (0, 100)
for (j in 1:J){
a[j] <- B[j,1]
b[j] <- B[j,2]
B[j,1:2] ~ dmnorm (B.hat[j,], Tau.B[,])
B.hat[j,1] <- mu.a
B.hat[j,2] <- mu.b
}
mu.a ~ dnorm (0, .0001)
mu.b ~ dnorm (0, .0001)
Tau.B[1:2,1:2] ~ dwish (W[,], df)
df <- 3
Sigma.B[1:2,1:2] <- inverse(Tau.B[,])
Sigma.B[1,1] <- pow(sigma.a, 2)
sigma.a ~ dunif (0, 100)
Sigma.B[2,2] <- pow(sigma.b, 2)
sigma.b ~ dunif (0, 100)
Sigma.B[1,2] <- rho*sigma.a*sigma.b
Sigma.B[2,1] <- Sigma.B[1,2]
rho ~ dunif (-1, 1)
}', file={f <- tempfile()})

 n = 919
 J = 85

jags_new <- jags.model(f, data = list('y' = y, 'x' = x, "n" = n, "J"=J,'county'=county),
n.chains = 1, n.adapt = 100)

#as.mcmc(jags_new)

posterior <- coda.samples(jags_new, c("a", "b", "mu.a", "mu.b", "sigma.y", "sigma.a", "sigma.b","rho" ) ,n.iter=5000)
summary(posterior)
plot(posterior)

burn.in <- 1000
#summary(posterior, start = burn.in)$quantiles[, c(3,1, 5)]



a.posterior.mean <- summary(posterior, start = burn.in)$statistics[1:85, "Mean"]

b.posterior.mean <- summary(posterior, start = burn.in)$statistics[1:85, "Mean"]
sigma.y.posterior.mean <- summary(posterior, start = burn.in)$statistics["sigma.y", 
    "Mean"]
}

ypred<-a.posterior.mean[county]+b.posterior.mean[county]*x

yresid=y-ypred

par(mfrow=c(2,3))
plot(ypred, y, xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="Bayesian method")
abline(0,1)
plot(ypred, yresid, xlab="Predicted", ylab="Residuals", main="Bayesian method")
abline(h=0)

hist(yresid)

#########################################################################################
setwd('~/Desktop/Mixed models/bayesianModel')

expec<-read.csv("life_exp.csv", header=TRUE)
expec1<-na.omit(expec)
x1<-log10(expec1$pop)
x2<-log10(expec1$gdpPercap)
y<-expec1$lifeExp

############ LM METHOD

m1<-lm(y~x1+x2)

ypred<- -15.367-1.154*x1+23.514*x2
resid<-y-ypred

par(mfrow=c(2,3))
plot(ypred, y,  xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="lmer method")
abline(0,1)
plot(ypred, resid, xlab="Predicted", ylab="Residuals", main="lmer method")
abline(h=0)

hist(resid)


########### LMER METHOD #### INTERCEPT ONLY

m2<-lmer(y~ x1+x2+(1|country), data=expec1, na.action="na.omit", REML=FALSE)

ypred<-predict(m2)
resid<-y-ypred

#par(mfrow=c(2,3))
plot(ypred, y,  xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="lmer method")
abline(0,1)
plot(ypred, resid, xlab="Predicted", ylab="Residuals", main="lmer method")
abline(h=0)

hist(resid)

############# BAYESIAN METHOD ### INTERCEPT ONLY
country<-factor(expec1$country.r)


#BAYESIAN METHOD

cat('model {
for (i in 1:n){
y[i] ~ dnorm (y.hat[i], tau.y)
y.hat[i] <- a[country[i]] + b*x1[i]+c*x2[i]
}
b ~ dnorm (40, .25)
c ~ dnorm (0, .25)
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif (0, 100)
for (j in 1:J){
a[j] ~ dnorm (mu.a, tau.a)
}
mu.a ~ dnorm (-220, .016)
tau.a <- pow(sigma.a, -2)
sigma.a ~ dunif (0, 100)
}', file={f <- tempfile()})

 n = 299
 J = 25

jags_new <- jags.model(f, data = list('y' = y, 'x1' = x1,'x2' = x2, "n" = n, "J"=J,'country'=country),
n.chains = 1, n.adapt = 100)

#as.mcmc(jags_new)

posterior <- coda.samples(jags_new, c("a", "b","c", "mu.a", "sigma.y", "sigma.a"), n.iter=5000)
summary(posterior)
plot(posterior)

burn.in <- 1000
#summary(posterior, start = burn.in)$quantiles[, c(3,1, 5)]



a.posterior.mean <- summary(posterior, start = burn.in)$statistics[1:25, "Mean"]

b.posterior.mean <- summary(posterior, start = burn.in)$statistics["b", "Mean"]
c.posterior.mean <- summary(posterior, start = burn.in)$statistics["c", "Mean"]


ypred<-a.posterior.mean[country]+b.posterior.mean*x1+c.posterior.mean*x2

yresid=y-ypred

par(mfrow=c(2,3))
plot(ypred, y, xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="Bayesian method")
abline(0,1)
plot(ypred, yresid, xlab="Predicted", ylab="Residuals", main="Bayesian method")
abline(h=0)

hist(yresid)

############# LMER METHOD ### SLOPE ONLY

m3<-lmer(y~ x1+x2+(1+x2|country), data=expec1, na.action="na.omit", REML=FALSE)

ypred<-predict(m3)
resid<-y-ypred

#par(mfrow=c(2,3))
plot(ypred, y,  xlim=c(range(y)), xlab="Predicted", ylab="Observed", main="lmer method")
abline(0,1)
plot(ypred, resid, xlab="Predicted", ylab="Residuals", main="lmer method")
abline(h=0)

hist(resid)

####### BAYESIAN METHOD ##SLOPE
country<-factor(expec1$country.r)
dwish (W[,], df)
df <- 3

cat('model {
for (i in 1:n){
y[i] ~ dnorm (y.hat[i], tau.y)
y.hat[i] <- a[country[i]] + c*x1[i]+ b[country[i]]*x2[i]
}
c ~ dnorm (0, .0001)
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif (0, 100)
for (j in 1:J){
a[j] ~ dnorm (a.hat[j], tau.a)
b[j] ~ dnorm (b.hat[j], tau.b)
a.hat[j] <- mu.a
b.hat[j] <- mu.b
}
mu.a ~ dnorm (0, .0001)
mu.b ~ dnorm (0, .0001)
tau.a <- pow(sigma.a, -2)
tau.b <- pow(sigma.b, -2)
sigma.a ~ dunif (0, 100)
sigma.b ~ dunif (0, 100)
}', file={f <- tempfile()})

  n = 299
  J = 25
#W <- diag (2)
jags_new <-jags.model(f, data = list('y' = y, 'x1' = x1,'x2' = x2, "n" = n, "J"=J,'country'=country),
n.chains = 1, n.adapt = 100)

#as.mcmc(jags_new)

posterior <- coda.samples(jags_new, c("a", "b","c", "mu.a", "mu.b", "sigma.y", "sigma.a", "sigma.b" ) ,n.iter=10000)

#plot(posterior)

burn.in <- 1000
summary(posterior, start = burn.in)
#summary(posterior, start = burn.in)$quantiles[, c(3,1, 5)]



a.posterior.mean <- summary(posterior, start = burn.in)$statistics[1:25, "Mean"]

b.posterior.mean <- summary(posterior, start = burn.in)$statistics[26:50, "Mean"]
c.posterior.mean <- summary(posterior, start = burn.in)$statistics["c", "Mean"]


ypred<-a.posterior.mean[country]+c.posterior.mean*x1+b.posterior.mean[country]*x2

yresid=y-ypred

par(mfrow=c(2,3))
plot(ypred, y, xlab="Predicted", ylab="Observed", main="Bayesian method")
abline(0,1)
plot(ypred, yresid, xlab="Predicted", ylab="Residuals", main="Bayesian method")
abline(h=0)

hist(yresid)


