#Spatial autologistic single-species, single-season occupancy
#Brent Pease @BrentPease1

#Load libraries
library("R2WinBUGS")
library(here)

setwd(here())

#Read in gridded dataset 
source("create_grid.R")

#Specify bugs directory
BD <- "C:/Users/Brent/Documents/WinBUGS14/"

#Assign objects from grid.data
data<- grid.data
x <- as.vector(data$x)
y<-data$y
J<-apply(!is.na(y),1,sum)
M<-nrow(y)
y<-apply(y,1,sum,na.rm=TRUE)

#Bundle data
win.data <- list ( "y","J","M","x")


setwd(here("bugs_output")) #sink all files into a specific folder
#Specify model
sink("model.txt")
cat("
    model { 
    #Prior distributions
    p ~ dunif(0,1)                
    b0 ~ dnorm(0,.01) 
    b1 ~ dnorm(0,.01)

    #Process and observation model specification
    for(i in 1:M){ 
    z[i] ~ dbin(psi[i],1)      
    logit(psi[i]) <- b0 + b1*x[i]
    mu[i]<-z[i]*p             
    y[i] ~ dbin(mu[i],J[i])  
    } 
    
    #Derived value
    logit(psi0)<-b0
    }
    ",fill=TRUE)
sink()

#Specify initial values
inits <- function()
  list ( z=rbinom(M,1,.4),b0=rnorm(1),p=runif(1),b1=rnorm(1))

#Parameters to monitor
parameters <- c("p","b0","b1","psi0")

#MCMC Settings
nt <- 1; nc <- 2; nb <- 1000; ni <- 11000
#Call BUGS from R
out1<-bugs (win.data, inits, parameters, "model.txt", n.thin=nt,
            n.chains=nc, n.burnin=nb,n.iter=ni,debug=FALSE,
           bugs.directory = BD, working.directory = getwd())
#View output
print(out1$summary, dig = 3)


