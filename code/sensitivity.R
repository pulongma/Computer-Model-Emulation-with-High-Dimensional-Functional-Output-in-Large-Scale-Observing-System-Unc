sensitivity <- function(gradf){
  out <- gradf %*% t(gradf)/dim(gradf)[1]
  
  return(out)
}
