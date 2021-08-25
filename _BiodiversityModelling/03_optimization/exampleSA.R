# Load the data 
data <- read.table("/home/local/USHERBROOKE/grad3002/Bureau/bios2.github.io/_BiodiversityModelling/02_likelihood/01-likelihood/data/hemlock.txt", header = T)
obs = data[,2]
L = data[,1]

# Likelihood function
h <- function(obs, L, pars) {
    a <- pars[1]
    s <- pars[2]
    sigma <- pars[3]
    mu <- a*L/(a/s + L)
    sum(log(dnorm(obs, mu, sigma)))
}

# Candidate function
c_x <- function(pars_lo, pars_hi) runif(1, pars_lo, pars_hi)

# Set conditions for the simulated annealing sequence
T_fn <- function(T0, alpha, step) T0*exp(alpha*step)

T0 <- 10
alpha = -0.001
nsteps <- 10000

# Prepare an object to store the result
res <- matrix(nr = nsteps, nc = 5)

# Initiate the algorithm
a0 <- 1
s0 <- 1000
sd0 <- 10
pars0 <- c(1, 100, 10)
pars_lo <- c(0, 1, 1)
pars_hi <- c(1000, 1000, 100)

# Main loop
for(step in 1:nsteps) {

    for(j in 1:3) {

        # Draw new value for parameter j
        pars1 <- pars0
        pars1[j] <- c_x(pars_lo[j], pars_hi[j])

        # Evaluate the function
        h1 <- h(obs, L, pars1)
        h0 <- h(obs, L, pars0)

        # Compute the difference
        diff <- h1 - h0

        # Accept if improvement
        if(diff > 0) pars0 <- pars1

        # Accept wrong candidates 
        else {	
            p <- exp(diff/T_fn(T0, alpha, step))
            if(runif(1)<p) pars0 <- pars1
        }
    }
    
    # Record values
    res[step,] <- c(step,pars0,h(obs,L,pars0)) 
}

#Plot the results
plot(c(1:nsteps), res[,5], type = "l", xlab = "Time step", ylab = "h(x)",cex = 2, log = "x")

dev.new()
plot(L, obs, cex = 2, xlab="Light", ylab = "Growth", cex.axis = 2, cex.lab = 2)

add_model<- function(step) {
    Lseq = seq(0,100,0.1)
    a=res[step,2]
    s=res[step,3]
    lines(Lseq,a*Lseq/(a/s+Lseq), lwd = 2)
}
add_model(2)
add_model(10)
add_model(100)
add_model(500)
add_model(1000)
