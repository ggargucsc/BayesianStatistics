##################################################################
##
## Carlin, Gelfand and Smith (1992): Hierarchical Bayesian Analysis of Changepoint Problems
## Counts of coal mining disasters in Great Britain by year from 1851 to 1962.
##
##################################################################
rm(list=ls(all=TRUE))
set.seed(79861)

y <- c(4,5,4,1,0,4,3,4,0,6,3,3,4,0,2,6,3,3,5,4,5,3,1,4,4,1,5,5,3,4,2,5,2,2,3,4,2,1,3,2,2,1,1,1,1,3,0,0,1,0,1,1,0,0,3,1,0,3,2,2,0,1,1,1,0,1,0,1,0,0,0,2,1,0,0,0,1,1,0,2,3,3,1,1,2,1,1,1,1,2,4,2,0,0,0,1,4,0,0,0,1,0,0,0,0,0,1,0,0,1, 0,1)
n <- length(y)


#update function for lamda
update.lamda<-function(m)
{
  sum.m<-sum(y[1:m])
  
  rgamma(1, sum.m+alpha,m+beta)
  
}

#update function for phi
update.phi<-function(m)
{
  sum.m<-sum(y[1:m])
  sum.n<-sum(y[1:n])
  
  rgamma(1, gam+sum.n-sum.m, n-m+delta)
  
}

#values for m is discrete, using sample function in R to update m 
update.m<-function(ulamda, uphi)
{
  unnormal.den.m<-rep(0,n)
  
  for(j in 1:n)
  {
 sum.j<-sum(y[1:j])
 sum.n<-sum(y[1:n])
 unnormal.den.m[j]<-ulamda^(sum.j)*uphi^(sum.n-sum.j)*exp(-j*ulamda)*exp(-(n-j)*uphi)
  }

#finding probabilities for m 
prob.m<-unnormal.den.m/sum(unnormal.den.m)
sample.m<-sample(1:n, 1, replace=T, prob.m)

return(sample.m)
 
}

#initializing values
alpha<-0.5
beta<-1.5
gam<-0.6
delta<-1.7

#set up MCMC
N.iter<-10000
post<-matrix(0, ncol=3, nrow=N.iter)
colnames(post)<- c("lam", "phi", "m")
post[1,]<-c(1,1,4)

for(i in 2:N.iter)
{
 post[i-1,"m"]<- update.m(post[i-1,"lam"], post[i-1,"phi"]) 
 print(post[i-1,"m"])
 post[i,"lam"]<- update.lamda(post[i-1,"m"])
 post[i,"phi"]<-update.phi(post[i-1,"m"])
 cat("iteration ", i, " of ", N.iter, "    \r")
 flush.console()
}

burnin <- 1:1000

#trace plots for lamda and phi
plot(post[-burnin,"lam"], type="l", main="Trace plot for lamda")
plot(post[-burnin,"phi"], type="l", main="Trace plot for phi" )

#histogram for m
hist(post[-burnin,"m"], xlim=range(20:60), breaks=200, main="Histogram for m")

## conditional posterior density graphs of lamda and phi
plot(density(c(post[-burnin,"lam"], post[-burnin,"phi"])), main="Density for lamda and phi")

## posterior means
apply(post[-burnin,], 2, mean)
#posterior 95% credible intervals
apply(post[-burnin,], 2, quantile, c(0.025, 0.975))

