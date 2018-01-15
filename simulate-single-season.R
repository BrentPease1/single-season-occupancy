#Simulate data for static occupancy model
#Brent Pease @BrentPease1
library(here)
setwd(here())


# Option to set random number generator to get the same data
set.seed(74)


# Set sample sizes
nsites <- 100                                   # Number of sites to be sampled
nsurveys <- 4                                   # Number of surveys (visits)


# Define parameter values
psi <- 0.6                                      # Initial probability of occupancy
p <- 0.3                                        # Probability of detection
alpha <- -0.5                                   # Define intercept for occupancy
beta <- 1.5                                     # Define magnitude of effect of neighbor




# Prepare vector and arrays to hold true and observed values
z <- matrix(NA, nrow=nsites, ncol = 1)          # Expected and perfectly observed occurrence    
y <- matrix(NA, nrow=nsites, ncol=nsurveys)     # Detection / non-detection data
x <- matrix(NA, nrow=nG, ncol=nG)



# generate presence/absence data (the truth)
# Following years
for(i in 1:nG){                             # Loop over sites
    z[i] <- rbinom(n = 1,size = 1,prob = psi)
  }


# Generate detection/non-detection data
for(i in 1:nsites){                                       # Loop over sites                                   # Loop over years
    mu <- z[i] * p                                      # Observation process                            
    for(j in 1:nsurveys){                                 # Loop over surveys
      y[i,j] <- rbinom(1, 1, mu)                        # Populate y array with detection/non-detection data
    }
  }


# `y` is the final simulated detection/non-detection data
print(y[1:6,])                                          # First 6 sites
