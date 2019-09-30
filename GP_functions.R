#' Power exponential correlation function
Internal_CorrMatPowerExp <- function(x1, x2,theta,return_dCdtheta = FALSE) {
    diffmat =abs(outer(x1,x2,'-'))
    tmax <- 3
    expLS = exp(tmax*theta[1])
    minpower <- 1
    maxpower <- 1.95
    alpha <- minpower + (theta[2]+1)/2 * (maxpower - minpower)
    h = diffmat/expLS
    C = exp(-(h)^alpha)
    if(return_dCdtheta){
      dCdtheta <- cbind(tmax*alpha*C*diffmat^alpha/expLS^alpha,-C*h^alpha*log(h)/2 * (maxpower - minpower))
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
}

Internal_CorrMatCauchySQ <- function(x1, x2,theta, return_dCdtheta = FALSE) {
    diffmat =abs(outer(x1,x2,'-')); 
    
    expLS = exp(3*theta[1])
    expHE = exp(3*theta[2])
    h = diffmat/expLS
    alpha = 2*exp(0+6)/(1+exp(0+6))
    halpha = h^alpha
    pow = -expHE/alpha
    
    C = (1+halpha)^pow
    if(return_dCdtheta){
      dCdtheta = cbind(3*expHE*((1+halpha)^(pow-1))*(halpha),3*C*pow*log(1+halpha))
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }

Internal_neglogpost <- function(theta, ys, Xs) {
  d = dim(Xs)[2]
  n = dim(Xs)[1]
  numpara = 2
  nugget = 10^(-10)
  
  Sigma_t = matrix(1,n,n)
  for (dimlcv in 1:d) { # Loop over dimensions
    V =Internal_CorrMatPowerExp(Xs[,dimlcv], Xs[,dimlcv],theta[(dimlcv-1)*numpara+1:numpara])
    Sigma_t = Sigma_t*V
  }
  Sigma_t = (1-nugget)*Sigma_t+diag(dim(Sigma_t)[1])*nugget
  
  try.chol <- try({Sigma_chol = chol(Sigma_t)}, silent = TRUE)
  if (inherits(try.chol, "try-error")) {
    return(Inf)
  }; rm(try.chol)
  
  output_is_1D <- !is.matrix(ys)
  if(output_is_1D){
    nsamples = length(ys)
  }else{
    nsamples = dim(ys)[1]
  }
  
  Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys,transpose=TRUE))
  sigma2_hat = colSums(as.matrix(ys*Sti_resid))/n
  lDet = 2*sum(log(diag(Sigma_chol)))
  
  neglogpost =  0.1*sum((log(1-theta)-log(theta+1))^2) #start out with prior
  neglogpost =  neglogpost+0.1*sum((log(1-theta)-log(theta+1))^4) #start out with prior
  if(output_is_1D){
    neglogpost = neglogpost+1/2*nsamples*log(sigma2_hat[1])+1/2*lDet#
  }else{
    neglogpost = neglogpost+1/2*(nsamples*mean(log(c(sigma2_hat)))+lDet)
  }
  return(neglogpost)
}

#' Gradient of negative log likelihood posterior
Internal_gneglogpost <- function(theta,ys,Xs) {
    d = dim(Xs)[2]
    n = dim(Xs)[1]
    
    numpara = 2
    nugget = 10^(-10)
    if(is.matrix(ys)){
      ndim = dim(ys)[2]
    }
    
    Sigma_t = matrix(1,n,n)
    dSigma_to = list(matrix(1,n,numpara*n),d) 
    for (dimlcv in 1:d){
      dSigma_to[[dimlcv ]] = matrix(1,n,numpara*n)
    }
    for (dimlcv in 1:d) { # Loop over dimensions
      V = Internal_CorrMatPowerExp(Xs[,dimlcv], Xs[,dimlcv],theta[(dimlcv-1)*numpara+1:numpara],return_dCdtheta=TRUE)
      Sigma_t = Sigma_t*V$C
      TV = rep(V$C, numpara)
      for (dimlcv2 in 1:d){
        if(dimlcv2==dimlcv){
          dSigma_to[[dimlcv2]] =  dSigma_to[[dimlcv2]]*V$dCdtheta
        }else{
          dSigma_to[[dimlcv2]] =  dSigma_to[[dimlcv2]]*TV 
        }
      }
    }
    
   Sigma_t = (1-nugget)*Sigma_t+diag(dim(Sigma_t)[1])*nugget
    for (dimlcv in 1:d) {
      dSigma_to[[dimlcv]] = (1-nugget)*dSigma_to[[dimlcv]]
    }
   
   
    try.chol <- try({Sigma_chol = chol(Sigma_t)}, silent = TRUE)
    if (inherits(try.chol, "try-error")) {
      stop("chol error in gneglogpost #1, this can happen when neglogpost is Inf")
    }; rm(try.chol)
    
    
    tempvec1 = backsolve(Sigma_chol,backsolve(Sigma_chol,ys,transpose = TRUE))
    sigma2_hat_supp = colSums(as.matrix(ys*tempvec1))/dim(Xs)[1]
    if(is.matrix(ys)){
      dsigma2_hat_supp = matrix(0,d*numpara,ncol=dim(ys)[2])
    }else{
      dsigma2_hat_supp = rep(0,d*numpara)
    }
    dlDet_supp=rep(0,d*numpara)
    for (dimlcv in 1:d) {
      for(paralcv in 1:numpara){
        dSigma_supp = as.matrix((dSigma_to[[dimlcv]])[ ,((paralcv-1)*n+1):(paralcv*n)])
        tempvec2= dSigma_supp%*%tempvec1
        if(is.matrix(dsigma2_hat_supp )){
          if(dim(dsigma2_hat_supp)[1]>1.5){
            dsigma2_hat_supp[(dimlcv-1)*numpara+paralcv,] =-colSums(as.matrix(tempvec1*tempvec2))/n
          }else{
            dsigma2_hat_supp[,(dimlcv-1)*numpara+paralcv] =-colSums(as.matrix(tempvec1*tempvec2))/n
          }
        }else{
          dsigma2_hat_supp[(dimlcv-1)*numpara+paralcv] =-colSums(as.matrix(tempvec1*tempvec2))/n
        }
        dlDet_supp[(dimlcv-1)*numpara+paralcv] =
          sum(diag(backsolve(Sigma_chol,backsolve(Sigma_chol,dSigma_supp,transpose = TRUE))))
      }
    }
    lDet_supp = 2*sum(log(diag(Sigma_chol)))
    
    sigma2_hat = sigma2_hat_supp
    dsigma2_hat = dsigma2_hat_supp
    dlDet = dlDet_supp
    lDet = lDet_supp
    
  neglogpost =  0.1*sum((log(1-theta)-log(theta+1))^2) #start out with prior
  gneglogpost = -0.2*(log(1-theta)-log(theta+1))*((1/(1-theta))+1/(1+theta))
  
  neglogpost =  neglogpost+0.1*sum((log(1-theta)-log(theta+1))^4) #start out with prior
  gneglogpost = gneglogpost-0.1*4*((log(1-theta)-log(theta+1))^3)*((1/(1-theta))+1/(1+theta))
  
  output_is_1D <- !is.matrix(ys)
  if(output_is_1D){
    neglogpost =neglogpost +1/2*(n*log(sigma2_hat[1])+lDet)#
   gneglogpost = gneglogpost+1/2*(n*dsigma2_hat / sigma2_hat[1]+dlDet)#n
  }else{
    neglogpost = neglogpost+1/2*(n*mean(log(c(sigma2_hat)))+lDet)
    gneglogpost = gneglogpost+1/2*dlDet
    for(i in 1:ndim){
      gneglogpost = gneglogpost+1/2*1/ndim*n*dsigma2_hat[,i]/sigma2_hat[i]
    }
  }
  
  return(gneglogpost)
}



# Est GP model given data
GPfitting <- function(Xs,Ys,theta0 = rep(0,2*(dim(Xs)[2]))) {
  GP = list()
  nugget = 10^(-10)
  
  n = dim(Xs)[1]
  d = dim(Xs)[2]
  numpara = 2
  
  xoffset = apply(Xs,2,mean)
  Xcenter = sweep(Xs,2,xoffset)
  xscale = apply(abs(Xcenter),2,mean)
  Xrescaled =  sweep(Xcenter,2,xscale,FUN = "/")
  
  if(!is.matrix(Ys)){
    GP$mu = mean(Ys)
    ys = Ys-GP$mu
  }else{
    GP$mu = apply(Ys,2,mean)
    ys <- sweep(Ys, 2, GP$mu)
  }
    opt.out = nlminb(
        theta0,
        objective = Internal_neglogpost,
        gradient = Internal_gneglogpost,
        ys = ys,
        Xs = Xrescaled,
        lower = -rep(1,2*d), upper = rep(1,2*d)
      )
    thetaMAP <- opt.out$par
    
    Sigma_t = matrix(1,n,n)
    for (dimlcv in 1:d) { # Loop over dimensions
      V =Internal_CorrMatPowerExp(Xrescaled[,dimlcv], Xrescaled[,dimlcv],thetaMAP[(dimlcv-1)*numpara+1:numpara])
      Sigma_t = Sigma_t*V
    }
    Sigma_t = (1-nugget)*Sigma_t+diag(dim(Sigma_t)[1])*nugget
    Sigma_chol = chol(Sigma_t)  
      
    GP$pw = backsolve(Sigma_chol,backsolve(Sigma_chol,ys,transpose=TRUE))
    GP$sigma2_hat = colSums(as.matrix(ys*GP$pw))/n
    GP$thetaMAP = thetaMAP
    GP$Sigma_chol = Sigma_chol
    GP$Sigmainv= solve(Sigma_t)
    GP$xoffset = xoffset
    GP$xscale = xscale 
    GP$Xrescaled=Xrescaled
    GP$ismultipleoutputs = is.matrix(ys)
  return(GP)
}

#' Predict with GP object
GPpred <- function(GP, Xp) {
  Xpcenter = sweep(Xp,2,GP$xoffset)
  Xprescaled =  sweep(Xpcenter,2,GP$xscale,FUN = "/")
  numpara = 2
  d = dim(Xp)[2]
  
  n = dim(GP$Sigma_chol)[1]
  np = dim(Xp)[1]
  Sigma_pt = matrix(1,np,n)
  for (dimlcv in 1:d) { # Loop over dimensions
    V =Internal_CorrMatPowerExp(Xprescaled[,dimlcv], GP$Xrescaled[,dimlcv],GP$thetaMAP[(dimlcv-1)*numpara+1:numpara])
    Sigma_pt = Sigma_pt*V
  }
  Xprescaled[,dimlcv]
  Sigma_p = matrix(1,np,np)
  for (dimlcv in 1:d) { # Loop over dimensions
    V =Internal_CorrMatPowerExp(Xprescaled[,dimlcv], Xprescaled[,dimlcv],GP$thetaMAP[(dimlcv-1)*numpara+1:numpara])
    Sigma_p = Sigma_p*V
  }
  
  yhat = sweep(Sigma_pt%*%GP$pw,2,GP$mu, FUN = "+")
  
  if(GP$ismultipleoutputs){
    var_scale = diag(GP$sigma2_hat)
  }else{
    var_scale = GP$sigma2_hat
  }
  
  corrmat_pred = Sigma_p-Sigma_pt%*%GP$Sigmainv%*%t(Sigma_pt)
  
  return(list("pred"=yhat,"corr_mat"=corrmat_pred,"var_scale"=var_scale))
}


