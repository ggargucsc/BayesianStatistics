# Example from Bayesian Data Analysis, Gelman on hierarchical normal model
# chapter 11 - Basics of Markov Chain Simulation

# Data is based on coagulation time in seconds for blood drawn from 24 animals randomly 
# allocated to four different diets. Different treatments have different number of 
# observations because the randomization was unrestricted.
# Gibbs sampling example

#dataset
y.1<-c(62,60,63,59)
y.2<-c(63,67,71,64,65,66)
y.3<-c(68,66,71,67,68,68)
y.4<- c(56,62,60,61,63,64,63,59) 

J<-4

#creating a matrix to save dataset
y<-matrix(0, nrow=8, ncol=J)
y[1:4,1]<-y.1
y[1:6,2]<-y.2
y[1:6,3]<-y.3
y[,4]<-y.4

n<-c(4,6,6,8)

#update function for theta to compute theta1, theta2, theta3 and theta4
update.thetaj<-function(j, tau2, sig2, mu)
{
  var.thetaj<-1/((1/tau2)+(n[j]/sig2))
  mean.thetaj<-((mu/tau2)+((n[j]*mean(y[1:n[j],j]))/sig2))*var.thetaj
  
  
  rnorm(1, mean=mean.thetaj, sd=sqrt(var.thetaj))
}

#update function for mu
update.mu<-function(theta, tau2)
{
 mean.mu<-mean(theta)
 var.mu<-(tau2)/J
 
 rnorm(1, mean=mean.mu, sd=sqrt(var.mu))
}

#update function for sig2
update.sig2<-function(theta)
{
  sum<-0
  for(j in 1:J)
   {
    sum<-sum+sum((y[1:n[j],j]-theta[j])^2)
   }
  scale.sig2<-sum/2
  
  1/rgamma(1, shape=sum(n)/2, scale=scale.sig2)
}

#update function for tau2
update.tau2<-function(theta, mu)
{
  scale.tau2<-sum((theta - mu)^2)
  
  1/rgamma(1, shape=(J-1)/2, scale=scale.tau2/2)
}

#set up MCMC
N.iter<-5000
#creating a posterior matrix
post<-matrix(0,ncol=7,nrow=N.iter)
colnames(post)<-c("theta1","theta2","theta3", "theta4", "mu","sig2","tau2")
post[1,]<-c(25,26,27,11,32,15,24)


for(i in 2:N.iter)
{
  post[i-1,"theta1"]<-update.thetaj(j=1,tau2=post[i-1,"tau2"],sig2=post[i-1,"sig2"], mu=post[i-1,"mu"])
  post[i-1,"theta2"]<-update.thetaj(j=2,tau2=post[i-1,"tau2"],sig2=post[i-1,"sig2"], mu=post[i-1,"mu"])
  post[i-1,"theta3"]<-update.thetaj(j=3,tau2=post[i-1,"tau2"],sig2=post[i-1,"sig2"], mu=post[i-1,"mu"])
  post[i-1,"theta4"]<-update.thetaj(j=4,tau2=post[i-1,"tau2"],sig2=post[i-1,"sig2"], mu=post[i-1,"mu"])
  theta<-c(post[i-1,"theta1"], post[i-1,"theta2"], post[i-1,"theta3"], post[i-1,"theta4"])
  post[i,"sig2"]<-update.sig2(theta)
  post[i,"tau2"]<-update.tau2(theta, post[i-1,"mu"])
  post[i,"mu"]<-update.mu(theta,post[i,"tau2"])  
}

burnin <- 1:1000

#trace plots
plot(post[1000:3000,"theta1"], type="l")
plot(post[1000:3000,"theta2"], type="l")
plot(post[1000:3000,"theta3"], type="l")
plot(post[1000:3000,"theta4"], type="l")
plot(post[-burnin,"mu"], type="l")
plot(post[-burnin,"sig2"], type="l")
plot(post[-burnin,"tau2"], type="l")

#posterior mean and 90% credible interval
apply(post[,], 2, mean)
apply(post[,], 2, quantile, c(0.05, 0.95))

#density plots
plot(density(post[,"theta1"]), xlim=range(57:64))
plot(density(post[,"mu"]), xlim=range(55:70))
plot(density(post[,"tau2"]), xlim=range(-1:1))
