# single-season occupancy using 'unmarked'
# Adapted from Kery & Royle (2016) Applied Hierarchical Modeling in Ecology (chapter 10)
# Intended to be a companion and extended workflow from www.github.com/BrentPease1/capture-history.git, which creates a capture (detection) history from eMammal data
library(unmarked)
library(data.table)
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

# Overview of occupancy covariates
covs <- cbind(canopy.orig, dist.orig) #cbind all occupancy covariates
par(mfrow = c(3,3), mar=c(rep(1,4)))                   #define plotting prelims
for(i in 1:2){                        #'1:2' specifies the number of covariates in 'covs'. 
  hist(covs[,i], breaks = 50, col = "grey", main = colnames(covs)[i])
}
pairs(cbind(canopy.orig, dist.orig))  #

#Correlation matrix
cov.cor <- cor(covs) # correlation matrix
par(mfrow = c(1,1), mar=rep(0,4)) #reset plotting prelims
source('correlation_pvalue.R')
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cov.cor, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=TRUE 
)




# Standardize covariates and mean-impute date and duration
# Compute means and standard deviations
(means <- c(apply(cbind(elev.orig, forest.orig), 2, mean), date.orig = mean(c(date.orig), na.rm = TRUE), dur.orig=mean(c(dur.orig), na.rm = TRUE)))
(sds <- c(apply(cbind(elev.orig, forest.orig), 2, sd), date.orig = sd(c(date.orig), na.rm = TRUE), dur.orig=sd(c(dur.orig), na.rm = TRUE)))

# Scale covariates
elev <- (elev.orig - means[1]) / sds[1]
forest <- (forest.orig - means[2]) / sds[2]
date <- (date.orig - means[3]) / sds[3]
date[is.na(date)] <- 0
dur <- (dur.orig - means[4]) / sds[4]
dur[is.na(dur)] <- 0

# Load unmarked, format data and summarize
library(unmarked)
umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(elev = elev, forest = forest), obsCovs = list(time = time, date = date, dur = dur))
summary(umf)

We want to identify a model that is useful for inference, specifically for prediction of squirrel distribution to the whole of Switzerland. We do some step-wise model selection first on the detection part, then on the occupancy part, while keeping the detection part as identified in the first step.

# Fit a series of models for detection first and do model selection
summary(fm1 <- occu(~1 ~1, data=umf))
summary(fm2 <- occu(~date ~1, data=umf))
summary(fm3 <- occu(~date+I(date^2) ~1, data=umf))
summary(fm4 <- occu(~date+I(date^2)+I(date^3) ~1, data=umf))
summary(fm5 <- occu(~dur ~1, data=umf))
summary(fm6 <- occu(~date+dur ~1, data=umf))
summary(fm7 <- occu(~date+I(date^2)+dur ~1, data=umf))
summary(fm8 <- occu(~date+I(date^2)+I(date^3)+dur ~1, data=umf))
summary(fm9 <- occu(~dur+I(dur^2) ~1, data=umf))
summary(fm10 <- occu(~date+dur+I(dur^2) ~1, data=umf))
summary(fm11 <- occu(~date+I(date^2)+dur+I(dur^2) ~1, data=umf))
summary(fm12 <- occu(~date+I(date^2)+I(date^3)+dur+I(dur^2) ~1, data=umf))

# Put the fitted models in a "fitList" and rank them by AIC
fms <- fitList("p(.)psi(.)"                     = fm1,
               "p(date)psi(.)"                      = fm2,
               "p(date+date2)psi(.)"                = fm3,
               "p(date+date2+date3)psi(.)"          = fm4,
               "p(dur)psi(.)"                       = fm5,
               "p(date+dur)psi(.)"                  = fm6,
               "p(date+date2+dur)psi(.)"            = fm7,
               "p(date+date2+date3+dur)psi(.)"      = fm8,
               "p(dur+dur2)psi(.)"                  = fm9,
               "p(date+dur+dur2)psi(.)"             = fm10,
               "p(date+date2+dur+dur2)psi(.)"       = fm11,
               "p(date+date2+date3+dur+dur2)psi(.)" = fm12)

(ms <- modSel(fms))

# Continue with model fitting for occupancy, guided by AIC as we go
# Check effects of elevation
summary(fm13 <- occu(~date+dur+I(dur^2) ~elev, data=umf))
summary(fm14 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2), data=umf))
summary(fm15 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+ I(elev^3), data=umf))
cbind(fm13@AIC, fm14@AIC, fm15@AIC) # model 14 with elev2 best

# Check effects of forest and interactions
summary(fm16 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest, data=umf))
summary(fm17 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest+I(forest^2), data=umf))
summary(fm18 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest+I(forest^2)+elev:forest, data=umf))
summary(fm19 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2), data=umf))
summary(fm20 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2)+I(elev^2):forest, data=umf))
summary(fm21 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2)+I(elev^2):forest+ I(elev^2):I(forest^2), data=umf))
cbind(fm16@AIC, fm17@AIC, fm18@AIC, fm19@AIC, fm20@AIC) # fm20 is best

# Check for some additional effects in detection
summary(fm22 <- occu(~date+dur+I(dur^2)+elev ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2)+I(elev^2):forest, data=umf))
summary(fm23 <- occu(~dur+I(dur^2)+date*(elev+I(elev^2)) ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2)+I(elev^2):forest, data=umf))
summary(fm24 <- occu(~dur+I(dur^2)+date*(elev+I(elev^2))+forest ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2)+I(elev^2):forest, data=umf))
cbind(fm22@AIC, fm23@AIC, fm24@AIC) # None better, hence, stay with model 20


library(AICcmodavg)
system.time(gof.boot <- mb.gof.test(fm20, nsim = 1000))
gof.boot


# Create new covariates for prediction ('prediction covs')
orig.elev <- seq(200, 2500,,100)    # New covs for prediction
orig.forest <- seq(0, 100,,100)
orig.date <- seq(15, 110,,100)
orig.duration <- seq(100, 550,,100)
ep <- (orig.elev - means[1]) / sds[1] # Standardise them like actual covs
fp <- (orig.forest - means[2]) / sds[2]
dp <- (orig.date - means[3]) / sds[3]
durp <- (orig.duration - means[4]) / sds[4]

# Obtain predictions
newData <- data.frame(elev=ep, forest=0)
pred.occ.elev <- predict(fm20, type="state", newdata=newData, appendData=TRUE)
newData <- data.frame(elev=0, forest=fp)
pred.occ.forest <- predict(fm20, type="state", newdata=newData, appendData=TRUE)
newData <- data.frame(date=dp, dur=0)
pred.det.date <- predict(fm20, type="det", newdata=newData, appendData=TRUE)
newData <- data.frame(date=0, dur=durp)
pred.det.dur <- predict(fm20, type="det", newdata=newData, appendData=TRUE)

# Plot predictions against unstandardized 'prediction covs'
par(mfrow = c(2,2), mar = c(5,5,2,3), cex.lab = 1.2)
plot(pred.occ.elev[[1]] ~ orig.elev, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. occupancy prob.", xlab = "Elevation (m)", frame = F)
matlines(orig.elev, pred.occ.elev[,3:4], lty = 1, lwd = 1, col = "grey")
plot(pred.occ.forest[[1]] ~ orig.forest, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. occupancy prob.", xlab = "Forest cover (%)", frame = F)
matlines(orig.forest, pred.occ.forest[,3:4], lty = 1, lwd = 1, col = "grey")
plot(pred.det.date[[1]] ~ orig.date, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Date (1 = 1 April)", frame = F)
matlines(orig.date, pred.det.date[,3:4], lty = 1, lwd = 1, col = "grey")
plot(pred.det.dur[[1]] ~ orig.duration, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Survey duration (min)", frame = F)
matlines(orig.duration, pred.det.dur[,3:4], lty = 1, lwd = 1, col = "grey")


# Predict abundance and detection jointly along two separate covariate gradients
# abundance ~ (forest, elevation) and detection ~ (survey duration, date)
pred.matrix1 <- pred.matrix2 <- array(NA, dim = c(100, 100)) # Define arrays
for(i in 1:100){
  for(j in 1:100){
    newData1 <- data.frame(elev=ep[i], forest=fp[j])       # For abundance
    pred <- predict(fm20, type="state", newdata=newData1)
    pred.matrix1[i, j] <- pred$Predicted
    newData2 <- data.frame(dur=durp[i], date=dp[j])        # For detection
    pred <- predict(fm20, type="det", newdata=newData2)
    pred.matrix2[i, j] <- pred$Predicted
  }
}

par(mfrow = c(1,2), cex.lab = 1.2)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=orig.elev, y=orig.forest, z=pred.matrix1, col = mapPalette(100), axes = FALSE, xlab = "Elevation [m]", ylab = "Forest cover [%]")
contour(x=orig.elev, y=orig.forest, z=pred.matrix1, add = TRUE, lwd = 1.5, col = "blue", labcex = 1.3)
axis(1, at = seq(min(orig.elev), max(orig.elev), by = 250))
axis(2, at = seq(0, 100, by = 10))
box()
title(main = "Expected squirrel occurrence prob.", font.main = 1)
points(data$ele, data$forest, pch="+", cex=1)

image(x=orig.duration, y=orig.date, z=pred.matrix2, col = mapPalette(100), axes = FALSE, xlab = "Survey duration [min]", ylab = "Date (1 = April 1)")
contour(x=orig.duration, y=orig.date, z=pred.matrix2, add = TRUE, lwd = 1.5, col = "blue", labcex = 1.3)
axis(1, at = seq(min(orig.duration), max(orig.duration), by = 50))
axis(2, at = seq(0, 100, by = 10))
box()
title(main = "Expected squirrel detection prob.", font.main = 1)
matpoints(as.matrix(data[, 13:15]), as.matrix(data[, 10:12]), pch="+", cex=1)


# Load the Swiss landscape data from unmarked
data(Switzerland)             # Load Swiss landscape data in unmarked
CH <- Switzerland

# Get predictions of occupancy prob for each 1km2 quadrat of Switzerland
newData <- data.frame(elev = (CH$elevation - means[1])/sds[1], forest = (CH$forest - means[2])/sds[2])
predCH <- predict(fm20, type="state", newdata=newData)

# Prepare Swiss coordinates and produce map
library(raster)
library(rgdal)

# Define new data frame with coordinates and outcome to be plotted
PARAM <- data.frame(x = CH$x, y = CH$y, z = predCH$Predicted)
r1 <- rasterFromXYZ(PARAM)     # convert into raster object

# Mask quadrats with elevation greater than 2250
elev <- rasterFromXYZ(cbind(CH$x, CH$y, CH$elevation))
elev[elev > 2250] <- NA
r1 <- mask(r1, elev)

# Plot species distribution map (Fig. 10-14 left)
par(mfrow = c(1,2), mar = c(1,2,2,5))
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette(100), axes = F, box = F, main = "Red squirrel distribution in 2007")
lakes <- readOGR(".", "lakes")
rivers <- readOGR(".", "rivers")
border <- readOGR(".", "border")
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)

# Plot SE of the species distrbution map (Fig. 10-14 right)
r2 <- rasterFromXYZ(data.frame(x = CH$x, y = CH$y, z = predCH$SE))
r2 <- mask(r2, elev)
plot(r2, col = mapPalette(100), axes = F, box = F, main = "Uncertainty map 2007")
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
points(data$coordx, data$coordy, pch = "+", cex = 0.8)


# Get extent of squirrel occurrence in 2007
sum(predCH$Predicted)                      # All quadrats
sum(predCH$Predicted[CH$elevation < 2250]) # Only at elevations < 2250 m


# Standardise prediction covariate identical to those in analysis
pelev <- (CH$elevation - means[1]) / sds[1]
pforest <- (CH$forest - means[2]) / sds[2]

# Define function that predicts occupancy under model 20
Eocc <- function(fm) {
  betavec <- coef(fm)[1:8]       # Extract coefficients in psi
  DM <- cbind(rep(1,length(pelev)), pelev, pelev^2, pforest, pforest^2, pelev*pforest, pelev*pforest^2, pelev^2*pforest) # design matrix
  pred <- plogis(DM%*%(betavec)) # Prediction = DM * param. vector 
  Eocc <- sum(pred)              # Sum over all Swiss quadrats (no mask)
  Eocc
}

(estimate.of.occurrence <- Eocc(fm20))    # Same as before, without mask
system.time(Eocc.boot <- parboot(fm20, Eocc, nsim=1000, report=10)) # 100 sec
plot(Eocc.boot)         # Plot bootstrap distribution of extent of occurrence
quantile(Eocc.boot@t.star, c(0.025, 0.975))