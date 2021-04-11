#######################################################################
#'@title Computation for predictive measures
#'
#'@param z.pred.df a list of principal component scores that are predicted
#' by a Gaussian process
#'
#'@param y.heldout a matrix with each column containing spectrum band at 
#' one sounding 
#' 
#'@param y.fpca a list of mean vector and basis functions obtained from
#' (functional) principal component analysis of the radiances data
#'
#'
#######################################################################

post.summary<-function(z.pred.df, y.heldout, y.fpca){

band.mean = y.fpca$mean
band.basis = y.fpca$basis

n.PC = length(z.pred.df)
n.wvln = dim(y.heldout)[1]

#### 
n.pred = dim(z.pred.df[[1]])[1]
n.samp = dim(z.pred.df[[1]])[2]


post.mean = vector(mode="list", n.PC)
post.sd = vector(mode="list", n.PC)
for(k in 1:n.PC){
  post.mean[[k]] = matrix(apply(z.pred.df[[k]], 1, mean), nrow=1)
  post.sd[[k]] = matrix(apply(z.pred.df[[k]], 1, sd), nrow=1)
}


### concatenate them together
post.mean = do.call(rbind, post.mean)
post.sd = do.call(rbind, post.sd)



if(n.PC==1){
  
  post.pred.perc = t(apply(z.pred.df[[1]], 1, quantile, probs=c(0.025, 0.975)))
  
  #### convert to original scale
  y.pred.mu = t(matrix(rep(c(band.mean), times=n.pred), 
           nrow=n.pred, byrow=TRUE)) +  band.basis[ ,1,drop=FALSE]%*%post.mean
  
  y.pred.CI025 = t(matrix(rep(c(band.mean), times=n.pred), 
          nrow=n.pred, byrow=TRUE)) + band.basis[ ,1,drop=FALSE]%*%post.pred.perc[ ,1]
  
  y.pred.CI975 = t(matrix(rep(c(band.mean), times=n.pred), 
          nrow=n.pred, byrow=TRUE)) + band.basis[ ,1,drop=FALSE]%*%post.pred.perc[ ,2]
  
  # normal distributional assumption
  y.pred.sd = band.basis[ ,1,drop=FALSE]%*%post.sd
  
}else{
  
  y.pred.CI025 = matrix(NA, n.wvln, n.pred)
  y.pred.CI975 = matrix(NA, n.wvln, n.pred)
  #post.samp.tmp = vector(mode="list", n.PC)
  post.samp.mat = matrix(NA, n.PC, n.samp)
  
  ## computing quantiles (most time-consuming part)
  for(j in 1:n.pred){
    
    for(k in 1:n.PC){
      #post.samp.tmp[[k]] = z.pred.df[[k]][j, ]
      post.samp.mat[k, ] = z.pred.df[[k]][j, ]
    }
    #post.samp.mat = do.call(rbind, post.samp.tmp)
    
    
    tmp = band.basis[ ,1:n.PC,drop=FALSE]%*%post.samp.mat
    tmp.org = t(matrix(rep(c(band.mean), times=n.samp), 
                       nrow=n.samp, byrow=TRUE)) + tmp 
    temp = apply(tmp.org, 1, quantile, probs=c(0.025, 0.975))
    y.pred.CI025[ ,j] =  temp[1, ]
    y.pred.CI975[ ,j] = temp[2, ]
  }
  
  #### convert to original scale
  y.pred.mu = t(matrix(rep(c(band.mean), times=n.pred), 
                       nrow=n.pred, byrow=TRUE)) +  band.basis[ ,1:n.PC,drop=FALSE]%*%post.mean
  
  # normal distributional assumption
  y.pred.sd = sqrt(band.basis[ ,1:n.PC,drop=FALSE]^2%*%post.sd^2)
  
}

#### RMSPE
RMSPE = sqrt(mean((y.heldout-y.pred.mu)^2))

#### 95% CI
PCI95 = y.pred.CI025<=y.heldout & y.heldout<=y.pred.CI975
#print(mean(PCI95))

### CRPS
crps = CRPS(y.heldout, y.pred.mu, y.pred.sd)

out = c(RMSPE=signif(RMSPE, 4), signif(mean(PCI95), 4), signif(mean(c(crps)), 4))
names(out) = c("RMSPE", "PCI95", "CRPS")


return(list(measure=out, pred=y.pred.mu, sd=y.pred.sd))

}

######################################################################
######################################################################
