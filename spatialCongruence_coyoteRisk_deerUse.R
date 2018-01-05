###########################
## R code used to assess the spatial congruency between coyote removal risk
## and deer fawn site selection.
##
## Title:
## Spatial processes decouple management from objectives in a heterogeneous 
## landscape: predator control as a case study
##
## Publication: 
## Ecological Applications, 2018
##
## Code author: 
## Peter Mahoney, 2017
## 
###########################

####################
## Data Management - Rstan
####################

## Coyote Resource Selection
##-------------------------------------------------------
## Availability sampled from within home ranges
## i.e., 3rd Order selection as described in Johnson 1980
##

load('./data/CoyoteSU_20pix_UseAvail.dat')

# Pull Means and SDs, then standardize
meansCoy <- apply(dCoy[, 5:26], 2, mean)
sdsCoy <- apply(dCoy[, 5:26], 2, sd)
dCoy[, 5:26] <- apply(dCoy[, 5:26], 2, function (x) (x - mean(x)) / sd(x))

# Add new covariates
dCoy$ES <- dCoy$E + dCoy$S
dCoy$slope2 <- dCoy$slope^2
dCoy$dist_hw2 <- dCoy$dist_hw^2
dCoy$dist_TerRds2 <- dCoy$dist_TerRds^2
dCoy$dist_grass2 <- dCoy$dist_grass^2

# Data inputs for Stan
Nc <- nrow(dCoy);                 # Number of records
#Yc <- length(unique(dCoy$Year)); # Number of years
Ic <- length(unique(dCoy$ID));    # Number of individuals
#yearsCoy <- as.numeric(as.factor(dCoy$Year)) # Years
indsCoy <- as.numeric(dCoy$ID)    # IDs

# Covariates for coyote resource selection model
# Determined from an independent run of models for coyote RS
covsCoy <- c("dist_con", "dist_mix", "dist_hw", "dist_hw2",
             "dist_shr",
             "dist_grass", "dist_grass2", "dist_roc",
             "vrm", "dist_TerRds", "dist_TerRds2", 'dist_water',
             "E", "S", "W", "NoAsp")
xC <- dCoy[, covsCoy]  # Matrix/data.frame for Stan input
Kc <- ncol(xC);        # Number of parameters/covariates
y_c <- dCoy$Used;      # Response (1/0) for coyote RS model


## Coyote Removal Risk
##-----------------------------------------
## Used <- removal locations
## Avail <- sampled at the study area scale
##

load('./data/AllRemovals_50gpix_UseAvail.dat')

# Capped TRI to allow model convergence
dRem$tri[dRem$tri > 20] <- 20

# Input data used to predict coyote probability of use
# as a covariate in the removal submodel.
# Standardize based on means/sdsCoy so they are on the same scale.
xInt <- dRem
for (n in names(meansCoy)) {
  xInt[, n] <- (xInt[, n] - meansCoy[n]) / sdsCoy[n]
}
xInt$dist_hw2 <- xInt$dist_hw^2
xInt$dist_TerRds2 <- xInt$dist_TerRds^2
xInt$dist_grass2 <- xInt$dist_grass^2
xInt <- xInt[, covsCoy]

# Other removal model covariates
# Pull Means and SDs, then standardize
meansRem <- apply(dRem[, 5:26], 2, mean)
sdsRem <- apply(dRem[, 5:26], 2, sd)
dRem[, 5:26] <- apply(dRem[, 5:26], 2, function (x) (x - mean(x)) / sd(x))

# Modify factors, add new covariates
dRem[which(dRem$StudyArea=='B/S'), 'StudyArea'] <- 'S'  # Not relevant
dRem$StudyArea <- factor(dRem$StudyArea)
dRem$ES <- dRem$E + dRem$S
dRem$slope2 <- dRem$slope^2
dRem$dist_hw2 <- dRem$dist_hw^2
dRem$dist_TerRds2 <- dRem$dist_TerRds^2
dRem$dist_grass2 <- dRem$dist_grass^2
dRem$dist_tc502 <- dRem$dist_tc50^2
dRem$tri2 <- dRem$tri^2

# Data inputs for Stan
Nr <- nrow(dRem);                               # Number of records
Na <- length(unique(dRem$StudyArea));           # Number of Study Areas
areas <- as.numeric(as.factor(dRem$StudyArea)); # Area IDs

# Covariates for coyote removal risk model 
# (+ coyote relative probability of use estimated in overall model)
covsRem <- c('S', 'E', 'dist_tc50', 'dist_tc502', 'tri', 'tri2')
xR <- dRem[, covsRem];  # Matrix/data.frame for Stan input
Kr <- ncol(xR);         # Number of parameters/covariates
y_r <- dRem$Removed;    # Response (1/0) for coyote risk model



## Deer Fawn Site Selection
##---------------------------------
## Availability sampled from within study area
## i.e., 2nd Order selection as described in Johnson 1980
##

load('./data/DeerSU_25pix_UseAvail.dat')

# Pull Means and SDs, then standardize
meansDeer <- apply(d[, 4:26], 2, mean)
sdsDeer <- apply(d[, 4:26], 2, sd)
d[, 4:26] <- apply(d[, 4:26], 2, function (x) (x - mean(x)) / sd(x))

# Add new covariates
d$vrm2 <- d$vrm^2
d$dist_TerRds2 <- d$dist_TerRds^2
d$NDVI2 <- d$NDVI^2
d$dist_shr2 <- d$dist_shr^2
d$dist_grass2 <- d$dist_grass^2
d$dist_asp2 <- d$dist_asp^2
d$dist_pj2 <- d$dist_pj^2
d$dist_con2 <- d$dist_con^2
d$dist_hw2 <- d$dist_hw^2
d$dist_springs2 <- d$dist_springs^2

# Data inputs for Stan
Nd <- nrow(d);                             # Number of records
Yd <- length(unique(d$Year));              # Number of years
yearsDeer <- as.numeric(as.factor(d$Year)) # Years

# Covariates for deer fawn site selection model 
covsDeer <- c('dist_TerRds', 'dist_TerRds2', 'NDVI', 'NDVI2', 'vrm', 'vrm2',
              'dist_shr', 'dist_shr2', 'dist_grass', 'dist_grass2', 'dist_asp', 'dist_asp2',
              'dist_pj', 'dist_pj2', 'dist_con', 'dist_hw', 'dist_hw2',
              'dist_springs', 'dist_springs2', 'E', 'S', 'W', 'NoAsp')
xD <- d[, covsDeer]  # Matrix/data.frame for Stan input
Kd <- ncol(xD);      # Number of parameters/covariates
y_d <- d$Used;       # Response (1/0) for deer model



## Data for prediction and Earth Mover's Distance / Congruence estimaton
##---------------------------------
## Ideally, this will be a data.frame composed of data from every pixel
## within a study area.  Here, for time, I provided data for every 10th pixel
## which has the appearance of coarsening the resolution.
## Peforming this operation is very resource intensive in terms of RAM and storage

load('./data/DataForPrediction_10p.dat')

# Predict Coyote Resource Selection across study area
# Standardize by coyote means and sds defined above
# And create appropriate variables
dC <- d
for (n in names(meansCoy)) {
  dC[, n] <- (dC[, n] - meansCoy[n]) / sdsCoy[n]
}
dC$ES <- dC$E + dC$S
dC$slope2 <- dC$slope^2
dC$dist_hw2 <- dC$dist_hw^2
dC$dist_grass2 <- dC$dist_grass^2
dC$dist_TerRds2 <- dC$dist_TerRds^2

# Predict Coyote Removal Risk across study area
# Standardize by coyote risk means and sds defined above
# And create appropriate variables (+ Coyote Prob. of Use estimated from the model)
dR <- d
for (n in names(meansRem)) {
  dR[, n] <- (dR[, n] - meansRem[n]) / sdsRem[n]
}
dR$dist_tc502 <- dR$dist_tc50^2
dR$tri2 <- dR$tri^2

# Predict Deer Site Selection across study area
# Standardize by deer means and sds defined above
# And create appropriate variables
dD <- d
for (n in names(meansDeer)) {
  dD[, n] <- (dD[, n] - meansDeer[n]) / sdsDeer[n]
}
dD$vrm2 <- dD$vrm^2
dD$dist_TerRds2 <- dD$dist_TerRds^2
dD$NDVI2 <- dD$NDVI^2
dD$dist_shr2 <- dD$dist_shr^2
dD$dist_grass2 <- dD$dist_grass^2
dD$dist_asp2 <- dD$dist_asp^2
dD$dist_pj2 <- dD$dist_pj^2
dD$dist_con2 <- dD$dist_con^2
dD$dist_hw2 <- dD$dist_hw^2
dD$dist_springs2 <- dD$dist_springs^2

# Create Stan inputs
xCoyPred <- dC[, covsCoy]   # Area-wide predictors for coyote model
xRemPred <- dR[, covsRem]   # Area-wide predictors for risk model
xDeerPred <- dD[, covsDeer] # Area-wide predictors for deer model
Npred <- nrow(xRemPred)     # Number of pixels/rows in prediction data



####################
## Hierarchical model - Rstan
####################
library(rstan) # see http://mc-stan.org/users/interfaces/rstan.html
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)  # 14 threads on a 16 thread machine

# Build list of data inputs for Stan
mod.dat <- list(Nc=Nc, Ic=Ic, Kc=Kc, y_c=y_c, indsCoy=indsCoy, xC=xC, #Yc=Yc, 
                Nr=Nr, Na=Na, Kr=Kr, y_r=y_r, areas=areas, xR=xR, xInt=xInt,
                Nd=Nd, Yd=Yd, Kd=Kd, y_d=y_d, yearsDeer=yearsDeer, xD=xD, #)
                xCoyPred=xCoyPred, xRemPred=xRemPred, xDeerPred=xDeerPred, Npred=Npred)

# Fit model
fit <- stan(file='./completeModel.stan', data=mod.dat, 
           iter=450, warmup=300, chains=14, refresh=10, seed=102,  # 1 chain per thread
           control=list(adapt_delta=0.95));

# Save fit
save(fit, file='CompleteModel_wPP.dat')

# Summary output
summ <- summary(fit, fit@sim$pars_oi[c(1:13)], probs = c(0.025, 0.975), digits=3)
write.csv(summ$summary, 'outputFull_Final.csv')



####################
## Posterior Predictive Check
####################
library(reshape2)
yPC <- as.data.frame(extract(fit, 'y_predC')[[1]])
yPD <- as.data.frame(extract(fit, 'y_predD')[[1]])
yPR <- as.data.frame(extract(fit, 'y_predR')[[1]])

yPR$sim <- 1:nrow(yPR)
yp <- yPR
ypd <- melt(yp, id.vars='sim')
ypd <- ypd[order(ypd$sim),]

ggplot(ypd, aes(x=value, group=sim, color=sim)) + 
  geom_density()
  #geom_histogram(bins=10)



####################
## Earth Mover's Distance & Congruence
####################
library(foreach)
library(doParallel)
library(parallel)
library(emdist)
library(raster)

# Pulls predictions from model 'fit'
x <- as.data.frame(extract(fit, 'y_RemPred')[[1]])
y <- as.data.frame(extract(fit, 'y_DeerPred')[[1]])


# Function for calculating Earth Mover's Distance
EMDchains <- function(X, Y, doPar=T, ncores=2) {
  if (!all(dim(X) == dim(Y))) stop('Dimensions of input matrices must be equal!')
  numVals <- ncol(X)
  
  if (doPar) {
    emdClust <- function (r) {
      xn <- density(X[,r], bw='sj', n=499)
      yn <- density(Y[,r], bw='sj', n=499)
      xmat <- cbind(xn$y, xn$x)
      ymat <- cbind(yn$y, yn$x)
      return(emd(xmat, ymat))
    }
    
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, 'emdClust')
    clusterEvalQ(cl, library(emdist))
    clusterEvalQ(cl, 'X')
    clusterEvalQ(cl, 'Y')
    out = unlist(parLapply(cl, 1:numVals, function (x) emdClust(x)))
    stopCluster(cl)
  }
  
  if (!doPar) {
    out <- c()
    for (n in 1:numVals) {
      xn <- density(X[,n], bw='sj')
      yn <- density(Y[,n], bw='sj')
      
      #xmat <- cbind(xn$y, xn$x)
      #ymat <- cbind(yn$y, yn$x)
      
      out <- c(out, emd(cbind(xn$y, xn$x), cbind(yn$y, yn$x)))
    }
  }
  
  return(out)
}

#
# Calculate Earth Mover's Distance
#-----------------------
#

# Run as a single file
#o <- EMDchains(x, y, doPar=T, ncores=4)
#save(o, file='EMD.dat')

# Or break original data into blocks for faster processing
brks <- seq(1, ncol(x), 240)

for (b in 2:length(brks)) {
  subX <- x[, (brks[b-1]+1):brks[b]]
  subY <- y[, (brks[b-1]+1):brks[b]]  
  
  print(paste('Block:', b-1, 'of', length(brks)))
  
  if(b == 2) {
    out <- EMDchains(subX, subY, doPar=T, ncores=15)
  } else 
    out <- c(out, EMDchains(subX, subY, doPar=T, ncores=15))
  
  if (b==100) save(out, file='EMD.dat')
}

save(out, file='EMD.dat')

#
# Estimate congruence
#---------------------

load('./EMD.dat') # loads object called out
emd <- out

# Estimate median deer probability of use
ymd <- apply(y, 2, median)
yl <- 1 / (1 + exp(-ymd))

# Estimate congruence by weighting EMD by median deer prob of use
cong <- (emd^2 * yl^2) / (emd * yl)

# Build raster
library(sp)
library(raster)

xy <- read.csv('./data/DataForPrediction_10p.csv')[,c('Easting', "Northing")]
dd <- as.data.frame(cbind(xy, cong))
east <- seq(min(xy$Easting), max(xy$Easting), 300)    # determines x-resolution, adjust accordingly based on pred data
north <- seq(min(xy$Northing), max(xy$Northing), 300) # determines y-resolution, adjust accordingly based on pred data

# Create SpatialPointsDataFrame from xy-congruence data
spdf <- SpatialPointsDataFrame(coords=dd[, c('Easting','Northing')], data = dd,
                               proj4string = CRS('+proj=utm +zone=12 + datum=NAD83'))

# Create empty raster with east-west resolution
rast <- raster(ncol = length(east), nrow = length(north), 
               ext = extent(spdf))

# Create congruence raster
rastOut <- rasterize(spdf, rast, spdf$cong, fun='last')
plot(rastOut)

# Output raster
writeRaster(rastOut, 'EMD_congruence.tif', 'GTiff')
