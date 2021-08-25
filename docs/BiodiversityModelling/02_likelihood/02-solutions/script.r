### Extract species, elevation and coordinates
sutton <- read.csv("/home/local/USHERBROOKE/grad3002/.local/share/Trash/files/bios2.3.github.io/_BiodiversityModelling/02_likelihood/02-solutions/data/sutton.csv", sep=";")
acsa <- sutton$acsa
xy <- sutton[,2:1]
elev<-as.matrix(xy[,1])

### Draw map
par(mar=c(5,0,0,0))
plot(xy,type="n",asp=1,axes=FALSE,xlab="",ylab="",las=1)
points(xy[which(acsa==1),],pch=15,cex=1.25)
points(xy[which(acsa==0),],pch=0,cex=1.25)

###########
# 2.2. Linear regression
###########

### Likelihood function
ll_fn <- function(a,b,sig,E,obs) {

    # Function for the mean
    mu <- a + b*E

    # PDF 
    lik <- dnorm(x=obs, mean = mu, sd = sig)
    loglik <- log(lik)

    # Return loglikelihood
    sum(loglik)
}

### Try different values 
ll_fn(a=25,b=-0.075,sig=10,E=sutton$y,obs=sutton$acsa)
ll_fn(a=25,b=-0.05,sig=10,E=sutton$y,obs=sutton$acsa)
ll_fn(a=30,b=-0.075,sig=10,E=sutton$y,obs=sutton$acsa)

### Run a grid search 
a_vec <- seq(15,45,length.out = 100)
b_vec <- seq(-0.1,-0.01,length.out = 100)
res <- matrix(nr=100,nc=100)
for(i in 1:100) {
    for(j in 1:100) {
        res[i,j] <- ll_fn(a=a_vec[i],b=b_vec[j],sig=10,E=sutton$y,obs=sutton$acsa)
    }
}

### Run a grid search 
a_vec <- seq(15,45,length.out = 100)
b_vec <- seq(-0.1,-0.01,length.out = 100)
res <- matrix(nr=100,nc=100)
for(i in 1:100) {
    for(j in 1:100) {
        res[i,j] <- ll_fn(a=a_vec[i],b=b_vec[j],sig=10,E=sutton$y,obs=sutton$acsa)
    }
}
image(x=a_vec,y=b_vec,z=res)

###########
# 2.2. Logistic regression
###########

### Likelihood function
ll_fn <- function(a,b,sig,E,obs) {

    # Function for the mean
    mu <- a + b*E

    # logit transformation 
    p = exp(mu)/(1+exp(mu))

    # no PDF for logistic regression
    # use the output of the model directly 
    lik <- numeric(length(obs))
    lik[obs==1] = log(p[obs==1])
    lik[obs==0] = log(1-p[obs==0])

    # Return loglikelihood
    sum(log(lik))
}

# Plot the model 
logit_fn <- function(x, a, b) {
    mu <- a + b*x
    exp(mu)/(1+exp(mu))
}

curve(logit_fn(x, a=-0.01, b = -0.01),xlim = c(0, 1000))


###########
# 2.3. Poisson regression
###########

### Likelihood function
ll_fn <- function(a,b,sig,E,obs) {

    # Function for the mean
    mu <- a + b*E

    # PDF 
    lik <- dpois(x=obs, lambda = exp(mu))
    loglik <- log(lik)

    # Return loglikelihood
    sum(loglik)
}



