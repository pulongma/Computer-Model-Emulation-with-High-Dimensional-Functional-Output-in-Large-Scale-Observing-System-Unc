gradient.f <- function(y.obs, K, y.model, prior.stsd, error.vec){

  
  N = dim(y.obs)[2] # number of replicates
  for(k in 1:N){
    flag_NA = which(is.na(y.obs[ ,k]))
    out[ ,k] = -t(K[-flag_NA,,k]%*%diag(prior.stsd[ ,k]))%*%(error.vec[-flag_NA,k]^(-1)*(y.obs[-flag_NA,k]-y.model[-flag_NA,k]))
  }
  
  return(out)
}
