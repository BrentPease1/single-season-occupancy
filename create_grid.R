#https://stackoverflow.com/questions/29105175/find-neighbouring-elements-of-a-matrix-in-r

#same dataset while learning
set.seed(74)

#Create matrix (lattice) of discrete, non-overlapping units
mat.row <- 10
mat.col <- 10

mat <- matrix(1:(mat.row*mat.col), mat.row, mat.col)
m2<-cbind(NA,rbind(NA,mat,NA),NA)
addresses <- expand.grid(x = 1:mat.row, y = 1:mat.col)

#Create sets of cells (units) that are neighbors of i, using queen's rule
N<-c()                                           
for(i in 1:-1)
  for(j in 1:-1)
    if(i!=0 || j !=0)
      N<-rbind(N,m2[addresses$x+i+1+nrow(m2)*(addresses$y+j)])
N <- t(N)                                                      #transpose output
N <- t(apply(N, 1, function(x)sort(x, na.last = TRUE)))        #sort output


#Calculate cardinality of each cell
numnn <- apply(N,1,function(x) sum(!(is.na(x))))

#commit the number of grid cells
nG <- mat.row*mat.col

#commit the grid dimensions
griddim <- mat.row

#Define sampling conditions
nsites <- mat.row*mat.col
nsurveys <- 4
nyears <- 1

#define parameter values
alpha <- -1
beta <- 2.5


#define initial occupancy across sites
psi <- 0.4


#define p
p <- 0.3


#Create arrays to hold true occurrence status and sampling realizations 
w <- array(NA, dim = c(nsites,nyears))
z <- array(NA, dim = c(nsites,nyears)) 
y <- matrix(NA, nrow=nsites, ncol = nsurveys)
tmp <- matrix(NA, nrow=nsites, ncol = nsurveys)


#Generate initial presence/absence to build neighborhood structure off of
for(i in unique(N)){
w <- rbinom(mat.row*mat.col,1,psi)}

#Calculate neighborhood autocovariate
x <- array(NA, dim = c(mat.row*mat.col, max(numnn)))
for(i in 1:nG){
  x[i,1]<-0
  for(j in 2:numnn[i]){
    x[i,j]<-x[i,j-1]+w[N[i,j-1]]
  }}

#calculate psi from the initial state layer of presence/absence
for(i in 1:nG){
psi[i] <- plogis(alpha + beta*(x[i,numnn[i]]/numnn[i]))}

#Generate true presence/absence across grid
for(i in 1:nG){
z <- rbinom(nG, 1, psi[i])}

#Generate detection/non-detection data (observed)
for(i in 1:nG){
  mu <- z[i]*p
for(j in 1:nsurveys){
  y[i,j] <- rbinom(1, 1, mu)
}}


#package up objects into a list
grid.data <- list(z=z,y=y,numnn=numnn,N = N, nG = nG, griddim = griddim, x=x)

#Visualize the truth
Zmat<-matrix(NA,nrow=mat.row,ncol=mat.col) #first, create a matrix to hold values of z
z2 <- as.vector(as.integer(z))
Zmat[1:nG]<- z2                             #assign values of z to matrix

#plot
par(cex.axis=1.8,cex.lab=2.0,cex.main=2,mar= c(5, 4, 4, 2)*1.2 + 0.1)
image(1:mat.row,1:mat.col,Zmat,col=rev(terrain.colors(10)),xlab="Easting",ylab="Northing",main = "True presence/absence across grid")

#Visualize detection/non-detection
Ymat<-matrix(NA,nrow=griddim,ncol=griddim)
Ymat[1:nG]<- as.integer(apply(y,1,max))

par(cex.axis=1.8,cex.lab=2.0,cex.main=2,mar= c(5, 4, 4, 2)*1.2 + 0.1)
image(1:mat.row,1:mat.col,Ymat,col=rev(terrain.colors(10)),xlab="Easting",ylab="Northing",main = "Detection/non-detection across grid")

rm(tmp,mu,nsites,nsurveys,nyears,w,s,addresses, m2, mat, N, x, y, Ymat, Zmat, alpha, beta, griddim, i, j, mat.col, mat.row, nG, numnn, p, psi, z, z2)
