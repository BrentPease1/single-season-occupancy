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
N <- t(N)                                        #transpose output

#rm(m2, mat, i, j)                                #clean up environemnt


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
alpha <- -0.05
beta <- .15

#define psi
psi <- array(NA, dim = c(nsites,nyears))
for(i in 1:nsites){
  for(t in 1:nyears){
psi[i,t] <- alpha + beta*(numnn[i]+1/numnn[i])}}
psi <- plogis(psi)

#define p
p <- 0.3


#Create arrays to hold true occurrence status and sampling realizations 
z <- array(NA, dim = c(nsites,nyears)) 
y <- array(NA, dim = c(nsites,nsurveys))


#Generate presence/absence (truth)
for(i in 1:nsites){
  for(t in 1:nyears){
z[i,t] <- rbinom(1,1,psi)
  }
}


#package up objects into a list
grid.data <- list(z=z,numnn=numnn,N = N, nG = nG, griddim = griddim)

#Visualize the truth
Zmat<-matrix(NA,nrow=mat.row,ncol=mat.col) #first, create a matrix to hold values of z
Zmat[1:nG]<-z                              #assign values of z to matrix

#plot
par(cex.axis=1.8,cex.lab=2.0,cex.main=2,mar= c(5, 4, 4, 2)*1.2 + 0.1)
image(1:mat.row,1:mat.col,Zmat,col=rev(terrain.colors(10)),xlab="Easting",ylab="Northing",main = "True presence/absence across grid")

#Visualize detection/non-detection
p = 0.3
y<-rbinom(length(z),1,p*z)
Ymat<-matrix(NA,nrow=griddim,ncol=griddim)
Ymat[1:nG]<-y

par(cex.axis=1.8,cex.lab=2.0,cex.main=2,mar= c(5, 4, 4, 2)*1.2 + 0.1)
image(1:mat.row,1:mat.col,Ymat,col=rev(terrain.colors(10)),xlab="Easting",ylab="Northing",main = "Detection/non-detection across grid")

