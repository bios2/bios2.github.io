# Load data
data <- read.table("/home/local/USHERBROOKE/grad3002/Bureau/bios2.github.io/_BiodiversityModelling/01_probabilities/06-problem/data/invasion_NZ.txt", header = T, sep = ",")

# Remove NAs
data <- subset(data, data$BodyMass!="NA" & data$BodyMass!=0)

# Subset the invaders
data_inv <- subset(data, data$Invader_NZ ==1)

# Transform bodymass on a log scale
M <- log(data$BodyMass) 
MI <- log(data_inv$BodyMass)

# Get the parameters of the body size distribution in the regional pool. Looks like the distribution of log body size is itself lognormal 
hist(M)
mean_M <- mean(M)
sd_M <- sd(M)
mean_MI <- mean(MI)
sd_MI <- sd(MI)

# Candidate function for body size selection 
PIgM <- function(obs,mu,sigma) exp(-(mu-obs)^2/2/sigma^2)

# Function for the global body size distribution 
PM <- function(obs,mean_M,sd_M) dnorm(obs,mean_M,sd_M)

# Compute the denominator (integral of the product of candidate function times the global distribution)
PI <- function(mu,sigma,mean_M,sd_M) sigma/(sigma^2+sd_M^2)^0.5*exp(-(mean_M-mu)^2/2/(sigma^2+sd_M^2))	


# Collect everything to get the probability of a certain size
PMgI <- function(obs,mean_M,sd_M,mu,sigma) {
    PIgM(obs,mu,sigma)*
    PM(obs,mean_M,sd_M)/
    PI(mu,sigma,mean_M,sd_M)
} 

# Now evaluate the fit of the model for a certain set of parameters
# Candidate values to propose are the midpoints between the global and the invader distributions
cand_mu <- (mean_M + mean_MI)/2
cand_sigma <- (sd_M + sd_MI)/2

# Now compute the probability of an invader having a mass of 5 
PMgI(obs = 5, mean_M, sd_M, mu = cand_mu, sigma = cand_sigma)

# Do it across the entire dataset and compute its likelihood
p_obs <- PMgI(obs = MI, mean_M, sd_M, mu = cand_mu, sigma = cand_sigma)
sum(log(p_obs))

# Try with a larger optimal size
p_obs <- PMgI(obs = MI, mean_M, sd_M, mu = 6, sigma = cand_sigma)
sum(log(p_obs))


