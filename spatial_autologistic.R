library(here)
library("R2WinBUGS")
# this function uses the supplied data file "grid.data" and creates
# survey data based on J=5 replicates with p = 0.15, and then fits
# the autologistic model (Panel 9.4) to the observations.
# This is a fairly large spatial data set that can take some time to fit
setwd(here())
source("create_grid.R")
BD <- "C:/Users/Brent/Documents/WinBUGS14/"


p<-0.3
J<-4
data<- grid.data
## supplied data have logit(psi) = -1 + 2.5*x[i]
## where x[i] is the number of neighboring cells occupied
nG<-data$nG
griddim<-data$griddim
z<-data$z
numnn<-data$numnn
N<-data$N


y<-rbinom(length(z),J,p*z)
Ymat<-matrix(NA,nrow=griddim,ncol=griddim)
Ymat[1:nG]<-y

par(cex.axis=1.8,cex.lab=2.0,mar= c(5, 4, 4, 2)*1.2 + 0.1)
image(1:10,1:10,Ymat,col=rev(terrain.colors(10)),xlab="Easting",ylab="Northing")


sink("model.txt")
cat("
    model{
    
    alpha ~ dnorm(0,.01)
    beta  ~ dnorm(0,.01)
    p ~ dunif(0,1)
    
    for(i in 1:nG){
    x[i,1]<-0
    for(j in 1:numnn[i]){
    x[i,j+1]<-x[i,j]+z[N[i,j]]
    }
    logit(psi[i])<- alpha + beta*(x[i,numnn[i]+1]/numnn[i])
    z[i]~dbern(psi[i])
    mu[i]<-z[i]*p
    y[i]~dbin(mu[i],J)
    }
    N.occ <- sum(z[])
    psi.fs <- N.occ/nG
    }
    ",fill=TRUE)
sink()

zst<-ifelse(y>0,1,0)

data <- list ( "numnn","N","nG","y","J")
inits <- function()
  list (z=zst, alpha=rnorm(1),beta=rnorm(1), p=runif(1) )

parameters <- c("alpha","beta","p","N.occ","psi.fs","x")

ni=3000;nb=1000;nt=1;nc=3
out<-bugs (data, inits, parameters, "model.txt", n.thin=nt,n.chains=nc, n.burnin=nb,n.iter=ni,debug=FALSE,
           bugs.directory = BD, working.directory = getwd())

print(out, dig=3)

