###############################################################################
############## General description
# This routine shows the entire pipelines to implement the NNGP-based emulator
# in combination with the active subspace approach for dimension reduction. 
# The details of the methods can be found in the paper titled "Computer Model
# Emulation with High-Dimensional Functional Output in Large-Scale Observing
# System Uncertainty Experiments".
#
############ Instructions for code implementation
#The entire codes have three parts.
# (1) In the first part, R codes are used to implement functional principal
# component analysis (FPCA) so that the functional output from OCO-2 forward 
# model can be represented as a linear combination of basis functions and 
# corresponding FPC scores. 
# (2) In the second part, R codes are used to find the active subspace from 
# the space formed by the state vectors. The key quantity here is the 
# projection matrix that projects a high-dimensional space into a low-dimensional
# space.
# (3) In the third part, MATLAB codes are used to train the NNGP emulator based
# on a set of training data, in which the active variables and viewing geometry
# are treated as input variables, and FPC scores are treated as output. 
###############################################################################
##### set directories and put the data into corresponding directories
# ./data: a directory where data are placed.
# ./code: a directory where all the R codes and MATLAB codes are placed.

rm(list=ls())
setwd("~/Documents/MATLAB/OCO2code")

library(fda)
library(ggplot2)
library(reshape2)
library(ncdf4)

###############################################################################
############################### Part I: FPCA ##############################

setwd("./data")

########################### FPCA for O2 band ##############################
## wavelengths for O2
wvln.o2 = as.numeric(read.csv("o2_wvln.csv",header=F))
## O2-band radiance
o2 = t(read.csv("o2_regular.csv",header=F))

## create spline basis
rangeval.o2 = range(wvln.o2); norder = 4;
breaks.o2 = seq(from=rangeval.o2[1], to=rangeval.o2[2], length.out = 500)
basisfd.o2 = create.bspline.basis(rangeval=rangeval.o2,norder=norder,
                                   breaks=breaks.o2)
## fit spline basis
fd.o2 = Data2fd(argvals=wvln.o2,y=o2,basisobj=basisfd.o2, lambda=1e-25)
func.o2=eval.fd(wvln.o2, fd.o2 ) # functional form

## perform FPCA
fd.pca.o2 = pca.fd(fd.o2,nharm=3)

## mean function in functinal form
m.o2 = mean(fd.o2)
mean.o2 = eval.fd(wvln.o2, m.o2)
## basis functions
basvec.o2=eval.fd(wvln.o2, fd.pca.o2$harmonics)
## reconstruction of radiances
ret.o2=basvec.o2%*%t(fd.pca.o2$scores)

### figures
## eigen value  plot
par(mfrow=c(1,2))
plot(fd.pca.o2$values[1:30], type='o', xlab="index", ylab="Eigenvalues")
grid(lwd=2, nx=10, ny=10)

# proportion of variance explained cumulative plot
plot(cumsum(fd.pca.o2$values[1:30])/sum(fd.pca.o2$values), type='o', xlab="index", ylab="Percentage of variation")
grid(lwd=2, nx=10, ny=10)

## reconstruction using first FPC score
plot(wvln.o2, o2[,1], type='l')
lines(wvln.o2, func.o2[,1], col='red', lwd=2)
lines(wvln.o2, mean.o2+ret.o2[,1], col='green', lwd=2)

legend("bottomright", c("original data", 
                        "functional data", 
                        "retrieved from fPCA"), 
       lty=c(1,1,1), col=c("black", "red", "green"))


####
theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 10
theme_mat$axis.text.x$size = 10
theme_mat$axis.title.y$size = 10
theme_mat$axis.text.y$size = 10
theme_mat$plot.title$size = 12
snd=c(10,20,50,100)
smpidx = rep(wvln.o2,4)
sound = rep(1:4,each=1016)
bdat = data.frame(FrameIdx=1:4064,Soundings=sound,Wavelength=smpidx)
bdat$Soundings = factor(bdat$Soundings,levels=1:4,labels=c("Sounding 10","Sounding 20","Sounding 50", "Sounding 100"))
ysmp=cbind(as.vector(o2[,snd]), as.vector(matrix(rep(as.vector(mean.o2),4),ncol=4)+ret.o2[,snd]))
ymlt = melt(ysmp,varnames=c("FrameIdx","Type"))
ymlt$Type = factor(ymlt$Type, levels=1:2,labels=c("Observed", "Retrieved"))
ymrg = merge(ymlt,bdat)
mplt = ggplot(ymrg,aes(x=Wavelength,y=value, group=Type)) + 
  geom_line(aes(color=Type), lwd=1) + 
  geom_point(aes(color=Type), lwd=0.5) +
  facet_wrap( ~ Soundings,nrow = 2, scale="free") + 
  xlab("Wavelength") + ylab("Radiance") + theme_mat

print(mplt)
###################### End of FPCA for O2-band ############################


########################### FPCA for WCO2 band ############################

wvln.wco2 = as.numeric(read.csv("wco2_wvln.csv",header=F))
wco2 = t(read.csv("wco2_regular.csv",header=F))

rangeval.wco2 = range(wvln.wco2); norder <- 4;
breaks.wco2 = seq(from=rangeval.wco2[1], to=rangeval.wco2[2], length.out = 500)

basisfd.wco2 = create.bspline.basis(rangeval=rangeval.wco2,norder=norder,breaks=breaks.wco2)
fd.wco2 = Data2fd(argvals=wvln.wco2,y=wco2,basisobj=basisfd.wco2, lambda=1e-25)
func.wco2 = eval.fd(wvln.wco2, fd.wco2)
fd.pca.wco2 = pca.fd(fd.wco2,nharm=3)
m.wco2 = mean(fd.wco2)
mean.wco2 = eval.fd(wvln.wco2, m.wco2)
basvec.wco2 = eval.fd(wvln.wco2, fd.pca.wco2$harmonics )
ret.wco2 = basvec.wco2[,1]%*%t(fd.pca.wco2$scores[,1])



## plot eigen values and cumulatives
par(mfrow=c(1,2))
plot(fd.pca.wco2$values[1:30], type='o', xlab="index", ylab="Eigenvalues")
grid(lwd=2, nx=10, ny=10)

plot(cumsum(fd.pca.wco2$values[1:30])/sum(fd.pca.wco2$values), type='o', xlab="index", ylab="Percentage of variation")
grid(lwd=2, nx=10, ny=10)

##ggplot type plots for 4 radiances
theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 10
theme_mat$axis.text.x$size = 10
theme_mat$axis.title.y$size = 10
theme_mat$axis.text.y$size = 10
theme_mat$plot.title$size = 12
snd=c(10,20,50,100)
smpidx = rep(wvln.wco2,4)
sound = rep(1:4,each=1016)
bdat = data.frame(FrameIdx=1:4064,Soundings=sound,Wavelength=smpidx)
bdat$Soundings = factor(bdat$Soundings,levels=1:4,labels=c("Sounding 10","Sounding 20","Sounding 50", "Sounding 100"))
ysmp=cbind(as.vector(wco2[,snd]), as.vector(matrix(rep(as.vector(mean.wco2),4),ncol=4)+ret.wco2[,snd]))
ymlt = melt(ysmp,varnames=c("FrameIdx","Type"))
ymlt$Type = factor(ymlt$Type, levels=1:2,labels=c("Observed", "Retrieved"))

ymrg = merge(ymlt,bdat)
# Plot some radiances
mplt = ggplot(ymrg,aes(x=Wavelength,y=value, group=Type)) + 
  geom_line(aes(color=Type), lwd=1) + 
  geom_point(aes(color=Type), lwd=0.5) +
  facet_wrap( ~ Soundings,nrow = 2, scales='free') + 
  xlab("Wavelength") + ylab("Radiance") + theme_mat

print(mplt)

##################### End of FPCA for WCO2 band ############################


########################## FPCA for SCO2 band ##############################
wvln.sco2 = as.numeric(read.csv("sco2_wvln.csv",header=F))
sco2 = t(read.csv("sco2_regular.csv",header=F))

rangeval.sco2 = range(wvln.sco2); norder <- 4;
breaks.sco2 = seq(from=rangeval.sco2[1], to=rangeval.sco2[2], length.out = 500)

basisfd.sco2 = create.bspline.basis(rangeval=rangeval.sco2,norder=norder,breaks=breaks.sco2)
fd.sco2 = Data2fd(argvals=wvln.sco2,y=sco2,basisobj=basisfd.sco2, lambda=1e-25)
func.sco2 = eval.fd(wvln.sco2, fd.sco2 )
fd.pca.sco2 = pca.fd(fd.sco2,nharm=3)
m.sco2 = mean(fd.sco2)
mean.sco2 = eval.fd(wvln.sco2, m.sco2 )
basvec.sco2 = eval.fd(wvln.sco2, fd.pca.sco2$harmonics )
ret.sco2 = basvec.sco2[,1]%*%t(fd.pca.sco2$scores[,1])

## plot eigne values and cumulative proportion of eigen values
par(mfrow=c(1,2))
plot(fd.pca.sco2$values[1:30], type='o', xlab="index", ylab="Eigenvalues")
grid(lwd=2, nx=10, ny=10)

plot(cumsum(fd.pca.sco2$values[1:30])/sum(fd.pca.sco2$values), type='o', xlab="index", ylab="Percentage of variation")
grid(lwd=2, nx=10, ny=10)

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 10
theme_mat$axis.text.x$size = 10
theme_mat$axis.title.y$size = 10
theme_mat$axis.text.y$size = 10
theme_mat$plot.title$size = 12
snd=c(10,20,50,100)
smpidx = rep(wvln.sco2,4)
sound = rep(1:4,each=1016)
bdat = data.frame(FrameIdx=1:4064,Soundings=sound,Wavelength=smpidx)
bdat$Soundings = factor(bdat$Soundings,levels=1:4,labels=c("Sounding 10","Sounding 20","Sounding 50", "Sounding 100"))
ysmp=cbind(as.vector(sco2[,snd]), as.vector(matrix(rep(as.vector(mean.sco2),4),ncol=4)+ret.sco2[,snd]))
ymlt = melt(ysmp,varnames=c("FrameIdx","Type"))
ymlt$Type = factor(ymlt$Type, levels=1:2,labels=c("Observed", "Retrieved"))

ymrg = merge(ymlt,bdat)
# Plot radiances
mplt = ggplot(ymrg,aes(x=Wavelength,y=value, group=Type)) + 
  geom_line(aes(color=Type), lwd=1) + 
  geom_point(aes(color=Type), lwd=0.5) +
  facet_wrap( ~ Soundings,nrow = 2, scale="free") + 
  xlab("Wavelength") + ylab("Radiance") + theme_mat

print(mplt)

######################## Enf of FPCA for SCO2 band ##########################

############################################################################
###### Save data for Emulation
save(fd.pca.o2, fd.pca.sco2, fd.pca.wco2, 
     file="./data/fPCA_BasisVectors_Means.RData")



#############################################################################
###################### Part II: Active Subspace #############################
# Here the full dataset lnd_nadir_201502_retrieved_state_vector_jacobian_fph.h5 
# needs to be download online at the NASA Goddard Earth Science Data and Information Services Center (GES DISC https://disc.gsfc.nasa.gov/datasets?page=1\&keywords=OCO-2). 
# This dataset from the full-physics forward model is around 9GB. Please use it for caution. 
# Instead, we only provide a subset of the original data, and save them as formats 
# that can be read into R and MATLAB easily.
rm(list=ls())

setwd("~/Documents/MATLAB/OCO2code")

library(xlsx)
library(ggplot2)
source("./code/multiplot.R")

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 16
theme_mat$axis.text.x$size = 16
theme_mat$axis.title.y$size = 16
theme_mat$axis.text.y$size = 16
theme_mat$plot.title$size = 16

### load OCO-2 data
## this dataset is around 9 GB, please download it online
nc1 = nc_open("./data/lnd_nadir_201502_retrieved_state_vector_jacobian_fph.h5") 
rd1 = ncvar_get(nc1,"measured_radiance")
jc1 = ncvar_get(nc1,"jacobian")  ## memory consuming, be cautious
stnms = ncvar_get(nc1,"state_vector_names")
wvln = ncvar_get(nc1,"wavelength")
stsd = ncvar_get(nc1,"state_vector_prior_stddev")
corpr = ncvar_get(nc1,"state_vector_prior_correlation")
stsd_err = ncvar_get(nc1, "measured_radiance_uncert")
qflg = ncvar_get(nc1,"xco2_quality_flag")
# prior mean
mu_x = ncvar_get(nc1, "state_vector_prior_mean")
# X-hat, retrieved state vectors
xhat = ncvar_get(nc1,"retrieved_state_vector")
# Y-hat, fitted radiances: Y-hat = F(X-hat,B)
yhat = ncvar_get(nc1,"modeled_radiance")
# Y, observed radiances: Y
y = ncvar_get(nc1, "measured_radiance")
# solar zenith angle
solzenith = ncvar_get(nc1, "solar_zenith_angle")
# solar azimuth angle
solazimuth = ncvar_get(nc1, "solar_azimuth_angle") 
# instrument zenith angle
instzenith = ncvar_get(nc1, "instrument_zenith_angle")
# instrument azimuth angle
instazimuth = ncvar_get(nc1, "instrument_azimuth_angle")
# latitude and longitude
lat = ncvar_get(nc1, "latitude")
lon = ncvar_get(nc1, "longitude")
nc_close(nc1)

lat = as.numeric(lat)
lon = as.numeric(lon)

###########################################################
source("./code/gradient.f.R")
source("./code/sensitivity.R")
n = 3048
N = dim(yhat)[2]
y = rd1[ ,1:N]
yhat = yhat[ ,1:N]
missingness=rep(0,n,1)
for(i in 1:n){
  missingness[i]=length(which(is.na(y[i,])))
}

## Delete completely missing values
missing.all=which(missingness==N)
Yobs = y[-missing.all, ]
Yhat = yhat[-missing.all, ]
error.vec = stsd_err[ ,1:N]
error.vec = error.vec[-missing.all,]

K = jc1[-missing.all, , ]
rm(jc1)
save("K", file="./jocabian_matrix.RData")
#load("./jocabian_matrix.RData")


m = dim(K)[2]

y.obs = Yobs
y.model = Yhat
out = matrix(NA, nrow=m, ncol=N)

#### calculate active subspace based on detrend and standardized inputs
# notice that the AS is built for xhat.tilde = DInv*(xhat - mu_x)

out = gradient.f(y.obs, K, y.model, stsd, error.vec)

E = sensitivity(out)

#### figures
eigen.E = svd(E)
E.eigenvalue = eigen.E$d
E.eigenvec = eigen.E$u

## plot eigenvalues
df = data.frame(index=seq(1,m), eigenvalue=E.eigenvalue)
h1 = ggplot(df, aes(x=index,y=eigenvalue)) + 
  geom_line(aes(x=index,y=eigenvalue)) + 
  geom_point(aes(x=index,y=eigenvalue)) + 
  xlab("Index") + ylab("Eigenvalue") + theme_mat
plot(h1)

## plot eigenvectors
dfv1 = data.frame(index=seq(1,m), eigenvec=E.eigenvec[ ,1])
h2 = ggplot(dfv1, aes(x=index,y=eigenvec)) + 
  geom_line(aes(x=index,y=eigenvec)) + 
  geom_point(aes(x=index,y=eigenvec)) + 
  xlab("Eigenvector: 1") + ylab("Eigenvector component") + theme_mat
plot(h2)

dfv2 = data.frame(index=seq(1,m), eigenvec=E.eigenvec[ ,2])
h3 = ggplot(dfv2, aes(x=index,y=eigenvec)) + 
  geom_line(aes(x=index,y=eigenvec)) + 
  geom_point(aes(x=index,y=eigenvec)) + 
  xlab("Eigenvector: 2") + ylab("Eigenvector component") + theme_mat
plot(h3)

dfv3 = data.frame(index=seq(1,m), eigenvec=E.eigenvec[ ,3])
h4 = ggplot(dfv3, aes(x=index,y=eigenvec)) + 
  geom_line(aes(x=index,y=eigenvec)) + 
  geom_point(aes(x=index,y=eigenvec)) + 
  xlab("Eigenvector: 3") + ylab("Eigenvector component") + theme_mat
plot(h4)

df1 = data.frame(index=seq(1,m), percentage=cumsum(E.eigenvalue)/sum(E.eigenvalue))
h5 = ggplot(df1, aes(x=index,y=percentage)) + 
  geom_line(aes(x=index,y=percentage)) + 
  geom_point(aes(x=index,y=percentage)) + 
  xlab("Number of eigenvalues") + ylab("Percentage of variation") + theme_mat
plot(h5)

df2 = data.frame(index=seq(1,m-1), gap=E.eigenvalue[-m]-E.eigenvalue[-1])
h6 = ggplot(df2, aes(x=index,y=gap)) + 
  geom_line(aes(x=index,y=gap)) + 
  geom_point(aes(x=index,y=gap)) + 
  xlab("Index") + ylab("Gap") + theme_mat
plot(h6)

#########
#### construct active variable 
# xhat.tilde = DInv*(xhat - mu_x)
xhat.tilde = stsd^(-1)*(xhat - mu_x) # detrended and standardized inputs
# active variable 
s = t(W1)%*%xhat.tilde
save(s, file="./actvars.RData")
write.table(s, file="./actvars.txt", row.names=FALSE, col.names=FALSE)

#############################################################################
#############################################################################

#############################################################################
############################### Part III: NNGP Emulation ####################
##### go to the following MATLAB files
##1. O21_AS_nugget_IndEmulation.m
##2. O22_AS_nugget_IndEmulation.m
##3. O23_AS_nugget_IndEmulation.m
##4. WCO2_AS_nugget_IndEmulation.m
##5. SCO2_AS_nugget_IndEmulation.m
##6. CO2_AS_nugget_SepEmulation.m

#############################################################################
#############################################################################

#############################################################################
############################# Final Results ################################
### go to the R code name Figure_paper.R
