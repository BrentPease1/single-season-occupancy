# single-season occupancy using 'unmarked'
# Adapted from Kery & Royle (2016) Applied Hierarchical Modeling in Ecology (chapter 10)
# Intended to be a companion and extended workflow from www.github.com/BrentPease1/capture-history.git, which creates a capture (detection) history from eMammal data
library(unmarked)
library(data.table)
library(AICcmodavg)
library(corrplot) #for plotting correlations


##############################
##### Read in data sets ######
##############################

#capture (detection) history
data <- read.csv("deer_detection_history.csv",header=T) #capture history created in ```Create_Capture_History.R```
str(data)

y <- as.matrix(data[,2:24])    #specify which columns contain capture history

#Occupancy and detection covariates
covariates <- read.csv("detection_covariates.csv",header=TRUE)
str(covariates)

# Unstandardized, original values of occupancy covariates
canopy.orig <- covariates[,"canopy"]
dist.orig <- covariates[,"dist_road"]

#unstandardized, oringial values of detection covariates
date.orig <- as.matrix(covariates[,4:26]) #julian date


###########################################
##### Visualize occupancy covariates ######
###########################################

# Overview of occupancy covariates
covs <- cbind(canopy.orig, dist.orig) #cbind all occupancy covariates
par(mfrow = c(3,3), mar=c(rep(1,4)))                   #define plotting prelims
for(i in 1:2){                        #'1:2' specifies the number of covariates in 'covs'. 
  hist(covs[,i], breaks = 50, col = "grey", main = colnames(covs)[i])
}
pairs(cbind(canopy.orig, dist.orig))  #

#Correlation matrix
cat("Correlation matrix of occupancy covariates:\n")
print(cov.cor <- cor(covs)) # correlation matrix
par(mfrow = c(1,1), mar=rep(0,4)) #reset plotting prelims
source('correlation_pvalue.R')


##################################
##### standardize covariates #####
##################################

#### Center covariates on zero and mean-impute detection covariates
# Compute means and standard deviations
(means <- c(apply(cbind(canopy.orig, dist.orig), 2, mean), date.orig = mean(c(date.orig), na.rm = TRUE)))
(sds <- c(apply(cbind(canopy.orig, dist.orig), 2, sd), date.orig = sd(c(date.orig), na.rm = TRUE)))

# Scale covariates
canopy <- (canopy.orig - means[1]) / sds[1]
dist <- (dist.orig - means[2]) / sds[2]
date <- (date.orig - means[3]) / sds[3]
date[is.na(date)] <- 0


#format data for occuapncy analysis and summarize
umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(canopy = canopy, dist = dist), obsCovs = list(date = date))
summary(umf)

# Fit null model
summary(fm1 <- occu(~1 ~1, data=umf))

#Look at the occupancy and detection estimates on probability scale
#Note that 'backTransform' will only work without intermediate steps if inquiring about the null model
print(occ.null <- backTransform(fm1, type="state"))
print(det.null <- backTransform(fm1, type="det"))


#####################
#### DETECTION ######
#####################

# Fit a series of models for detection first and do model selection
fm2 <- occu(~date ~1, data=umf)
fm3 <- occu(~date+I(date^2) ~1, data=umf)
fm4 <- occu(~date+I(date^2)+I(date^3) ~1, data=umf)


# Put the fitted models in a "fitList"
fms <- fitList("p(.)psi(.)"                     = fm1,
               "p(date)psi(.)"                      = fm2,
               "p(date+date2)psi(.)"                = fm3,
               "p(date+date2+date3)psi(.)"          = fm4)
#Rank detection models by AIC
print(ms <- modSel(fms))

#Look at estimates of top model on probability scale
#Need to first use LinearComb() since there are covariates involved
lc <- linearComb(fm2, c(1, 1), type="det") # Estimate abundance on the log scale when forest=0
backTransform(lc)                           # Abundance on the original scale 


###################
#### OCCUPANCY ####
###################


#Use covariates from top detection model in occupancy model fitting
# Continue with model fitting for occupancy, guided by AIC
fm5 <- occu(~date ~dist, data=umf)
fm6 <- occu(~date ~dist + I(dist^2), data=umf)
fm7 <- occu(~date ~canopy, data=umf)
fm8 <- occu(~date ~canopy + dist, data=umf)
fm9 <- occu(~date ~canopy + dist + I(dist^2), data=umf)

# Put the fitted models in a "fitList"
fms2 <- fitList("p(date)psi(.)"                      = fm2,
               "p(date)psi(dist)"                = fm5,
               "p(date)psi(dist + I(dist^2))"          = fm6,
               "p(date)psi(canopy)"          = fm7,
               "p(date)psi(canopy + dist)"          = fm8,
               "p(date)psi(canopy + dist + I(dist^2))"          = fm9)
#Rank detection models by AIC
(ms2 <- modSel(fms2))

#### 'fm7' is our top model. We will be using this model for assessing GoF and predictions

############################
##### Assess model fit #####
############################
#Conduct a MacKenzie-Bailey Goodness-of-Fit test

#system.time(gof.boot <- mb.gof.test(fm7, nsim = 1000)) #Note: This can take a long time (hours). Best to run it when you are ready, then save the output as .Rdata file to easily load. I suggest once you have the output from this, hash (#) out the mb.gof.test line
#save(gof.boot, file="gof.boot_fm7.Rdata")
#print(gof.boot)

###########################
###### Predictions ########
###########################

# Create new covariates for prediction ('prediction covs')
# Really only need the covariates that were in our top model 
orig.canopy <- seq(0, 100,,87)# 'Create 87 values ranging from 0-100, and R determines the spacing between numbers (,,)'
orig.date <- seq(1, 230,,87)

can.pred <- (orig.canopy - means[1]) / sds[1] # Standardize them like actual covs
date.pred <- (orig.date - means[3]) / sds[3]


# Obtain predictions
newData <- data.frame(canopy=can.pred) 
#Note: if your top model has more than one covariate, do something like this: newData <- data.fream(cov1 = cov1.pred, cov2=0), if you want to predict values of cov1 and vice versa for cov2 
pred.occ.canopy <- predict(fm7, type="state", newdata=newData, appendData=TRUE)

newData <- data.frame(date=date.pred)
pred.det.date <- predict(fm7, type="det", newdata=newData, appendData=TRUE)


# Plot predictions against unstandardized 'prediction covs'
par(mfrow = c(2,2), mar = c(5,5,2,3), cex.lab = 1.2) #set plotting frame

plot(pred.occ.canopy[[1]] ~ orig.canopy, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. occupancy prob.", xlab = "Canopy Cover (%)", frame = F)
matlines(orig.canopy, pred.occ.canopy[,3:4], lty = 1, lwd = 1, col = "grey")

plot(pred.det.date[[1]] ~ orig.date, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Julian Date", frame = F)
matlines(orig.date, pred.det.date[,3:4], lty = 1, lwd = 1, col = "grey")

