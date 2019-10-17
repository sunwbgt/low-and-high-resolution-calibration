# read simulation data
setwd("E:/Research Projects/Low and High Resolution Calibration/EmuAccord")
source("GP_functions.R")
source("Injury_risk_functions.R")
A = read.csv('Accord_Simulation_Results_LL_Airbag_300_Large_Range.csv', header = TRUE)
# settings
nvar = 3
npara = 7
parameter_default = c(1, 0.45, 1, 1, 0.55, 1, 1)
w_default = rep(50,length(parameter_default))

nres = 4
nsim = 196
XCE = A[,2:(nvar+npara+1)]
YCE = log(as.matrix(A[,(nvar+npara+2):ncol(A)]))

# read field data
B = read.csv('Field_Test_Data_Full_Size_Sedan.csv', header = TRUE)
nfield = nrow(B)

# predict injuries for field data
GP = GPfitting(XCE, YCE)
Xpred_field = cbind(B[, 4 : 6], t(matrix(parameter_default, npara, nfield)))

# read injury data
Amat = cbind(rep(1,length(B[,8])),B[,8])

XF = B[, 4:6]
XFs = B$Age
ZF = 0+(B[,8]>=3)

C <- read.csv('injury_function_data_chest.csv', header = TRUE)
XEs = C$Age
YE = C$ChestD
YE = pmax(YE,10^(-2))
ZE = C$Injury.Indicator

# read experimental data
D <- read.csv('Experimental Data - All Sedans.csv')
exp_mean <- mean(D[, 9])
exp_std <- sqrt(var(D[, 9]))
ZD = D[, 9]
nexp = length(ZD)

# specify data sources: injury test -> 0, field data -> 1, experimental test -> 2
Type1 = rep(0,length(XEs))
Type2 = rep(1,length(XFs))
Type3 = rep(2,length(ZD))

# correlation matrix for GP of delta()
corrmat_x <- function(Xa1,Xa2,Xa3,...,Xb1=NA,Xb2=NA,Xb3=NA){
  if (any(is.na(Xb1))){
    Xb1 = Xa1
    Xb2 = Xa2
    Xb3 = Xa3
  }
  R1 <- exp(-abs(outer(Xa1,Xb1,'-')/12)^1.95)
  R2 <- exp(-abs(outer(Xa2,Xb2,'-')/1)^1.95)
  R3 <- exp(-abs(outer(Xa2,Xb2,'-')/17)^1.95)
  return(R1*R2*R3)
}

# correlation matrix for GP of g()
corrmat <-  function(Xa1,Xa2,Xa3,...,Xb1=NA,Xb2=NA,Xb3=NA){
  if (any(is.na(Xb1))){
    Xb1 = Xa1
    Xb2 = Xa2
    Xb3 = Xa3
  }
  R1 <- exp(-abs(outer(Xa1,Xb1,'-')/40)^1.95)
  # change order from 1.95 to 2 for an easier error-in-variable calculation
  R2 <- exp(-abs(outer(Xa2,Xb2,'-')/0.035)^2)
  Xa3[Xa3 == 2] = 0
  Xb3[Xb3 == 2] = 0  
  R3 <- 0.95+0.05*(1-abs(outer(Xa3,Xb3,'-')))
  return(R1*R2*R3)
}

# correlation matrix for error-in-variable GPs
corrmat_error <-  function(Xa1,Xa2,Xa3,...,Xb1=NA,Xb2=NA,Xb3=NA){
  # specify standard deviation of y in injury tests: 5% of the population standard deviation
  injury_std <- 0.0014
  if (any(is.na(Xb1))){
    Xb1 = Xa1
    Xb2 = Xa2
    Xb3 = Xa3
  }
  R1 <- exp(-abs(outer(Xa1,Xb1,'-')/40)^1.95)
  # change order from 1.95 to 2 for an easier error-in-variable calculation
  R2 <- matrix(1, length(Xa1), length(Xb1))
  for (i in 1 : length(Xa1)) {
	for (j in 1 : length(Xb1)) {
	  # both from injury test
      if ((Xa3[i] == 0) && (Xb3[j] == 0)) {
	    R2[i, j] = exp(-((Xa2[i] - Xb2[j])/0.035/(1+4*injury_std^2/0.035)) ^ 2)/(1+4*injury_std^2/0.035)^(1/2)
      }
      # only one from injury test
      if (((Xa3[i] == 0) && (Xb3[j] > 0)) || ((Xa3[i] > 0) && (Xb3[j] == 0))) {
	    R2[i, j] = exp(-((Xa2[i] - Xb2[j])/0.035/(1+2*injury_std^2/0.035)) ^ 2)/(1+2*injury_std^2/0.035)^(1/2)
      }
      # none from injury test
      if ((Xa3[i] > 0) && (Xb3[j] > 0)) {
        R2[i, j] = exp(-((Xa2[i] - Xb2[j])/0.035) ^ 2)
      }
	}
  }
  # data with labels 0 or 2 come from the same lab environment (share the same type)
  Xa3[Xa3 == 2] = 0
  Xb3[Xb3 == 2] = 0
  R3 <- 0.95+0.05*(1-abs(outer(Xa3,Xb3,'-')))
  return(R1*R2*R3)
}

# derivative of correlation matrix with respect to delta for error-in-variable GPs
dcorrmat_error <-  function(Xa1,Xa2,Xa3,...,Xb1=NA,Xb2=NA,Xb3=NA){
  # specify standard deviation of y in injury tests: 5% of the population standard deviation
  injury_std <- 0.0014
  if (any(is.na(Xb1))){
    Xb1 = Xa1
    Xb2 = Xa2
    Xb3 = Xa3
  }
  R1 <- exp(-abs(outer(Xa1,Xb1,'-')/40)^1.95)
  # change order from 1.95 to 2 for an easier error-in-variable calculation
  R2 <- matrix(1, length(Xa1), length(Xb1))
  for (i in 1 : length(Xa1)) {
	for (j in 1 : length(Xb1)) {
	  # both from injury test
      if ((Xa3[i] == 0) && (Xb3[j] == 0)) {
	    R2[i, j] = exp(-((Xa2[i] - Xb2[j])/0.035/(1+4*injury_std^2/0.035)) ^ 2)/(1+4*injury_std^2/0.035)^(1/2)*2*(Xa2[i]-Xb2[j])/0.035/(1+4*injury_std^2/0.035)
      }
      # only one from injury test
      if (((Xa3[i] == 0) && (Xb3[j] > 0)) || ((Xa3[i] > 0) && (Xb3[j] == 0))) {
	    R2[i, j] = exp(-((Xa2[i] - Xb2[j])/0.035/(1+2*injury_std^2/0.035)) ^ 2)/(1+2*injury_std^2/0.035)^(1/2)*2*(Xa2[i]-Xb2[j])/0.035/(1+2*injury_std^2/0.035)
      }
      # none from injury test
      if ((Xa3[i] > 0) && (Xb3[j] > 0)) {
        R2[i, j] = exp(-((Xa2[i] - Xb2[j])/0.035) ^ 2)*2*(Xa2[i]-Xb2[j])/0.035
      }
	}
  }
  # data with labels 0 or 2 come from the same lab environment (share the same type)
  Xa3[Xa3 == 2] = 0
  Xb3[Xb3 == 2] = 0
  R3 <- 0.95+0.05*(1-abs(outer(Xa3,Xb3,'-')))
  return(R1*R2*R3)
}

# mvec <- function(X,Y){
  # return(-9+0.01*X+15*sqrt(Y))
# }

# use the one provided in the literature - do not use sqrt anywhere else
mvec <- function(X,Y){
  return(-12.597+0.05861*X+1.568*sqrt(Y * 1000))
}

nsamp = 10
normvecs = matrix(rnorm(length(ZF)*nsamp,0,1),nrow=length(ZF),ncol=nsamp) #makes function continous by fixing this

# optimization function for thetaMAP: still use corr_mat because the injury test data should not be included here
PostCalcApprox <- function(x,GPmodel,XF,XFs,ZF,ZD,integrateoveremulator=FALSE) {
  GP_obj = GPpred(GPmodel, cbind(XF, t(matrix(x, npara, nfield))))
  YF = exp(GP_obj$pred)[,1]
  
  #likellihood for field data
  Rv <- corrmat(XFs,YF,Type2)
  if(integrateoveremulator){ #this will be a quick approximation, just to check
    covemuerror = GP_obj$corr_mat*GP_obj$var_scale[1,1] + diag(rep(10^(-8),length(YF)))
    L = chol(covemuerror)
    Rv <- matrix(0,nrow=length(YF),ncol=length(YF))
    
    for(kv in 1:nsamp){ #draw over the emulator error
      YFsamp = exp(GP_obj$pred[,1] + L%*%normvecs[,kv])
      Rv <- Rv+1/nsamp*corrmat(XFs,c(YFsamp),Type2)
    }
  }
  Cv <- 4*Rv+10^(-4)*diag(rep(1,length(YF)))
  
  mv <- mvec(XFs,YF)
  pv <- exp(mv)/(1+exp(mv))
  diagadj = exp(mv)/(1+exp(mv))^2
  Cva = diag(diagadj)%*%Cv%*%diag(diagadj)
  Va = diag(pmax(pv*(1-pv),10^(-4)))
  
  #forgot det term on previous rounds, it is included below
  LikA =0.5*sum(log(eigen(Cva+Va)$values))+ 0.5*t(ZF-pv)%*%solve(Cva+Va,ZF-pv)
  
  #likelihood for exp data
  twoexpvect = matrix((rep(c(25.47,1.75,35,parameter_default),1)),nrow=1,byrow=T)
  pd = rep(exp(GPpred(GPmodel, twoexpvect)$pred[1,1]), length(ZD))
  LikA = LikA + 0.5 * t(ZD - pd)%*%diag(rep(0.005^(-2), length(ZD)))%*%(ZD-pd)
  
  #prior
  Prior = sum(c(w_default)*(x[1:npara]-parameter_default)^2) #pull it closer to true values
  return(LikA+Prior)
}

opt.out1 = optim(
  parameter_default,
  fn = PostCalcApprox,
  XF=XF,XFs=XFs,ZF=ZF,ZD=ZD,GPmodel=GP,integrateoveremulator=TRUE,
  lower = c(apply(XCE[,4:(3+npara)],2,min)), 
  upper = c(apply(XCE[,4:(3+npara)],2,max)),
  control=list('factr' = 10^12, 'maxit'=50), 
  hessian=T,
  method = "L-BFGS-B"
)
opt.out1$par
unname(0+apply(XCE[,4:(3+npara)],2,min))
unname(0+apply(XCE[,4:(3+npara)],2,max))

# MAP estimation for simulation parameters
thetaMAP = opt.out1$par
ApproxCItheta = cbind(thetaMAP-1.96*sqrt(diag(opt.out1$hessian^(-1))),thetaMAP+1.96*sqrt(diag(opt.out1$hessian^(-1))))
print(ApproxCItheta)

# optimization functions
PostCalc <- function(Z,zvec) {
  pv = exp(zvec)/(1+exp(zvec))
  Lik = sum(-(1-Z)*log(1-pv)-Z*log(pv))
  q = zvec-mv
  Prior = 1/2*sum(forwardsolve(t(cholCv),q)^2)
  return(Lik+Prior)
}
gPostCalc <- function(Z,zvec) {
  pv = exp(zvec)/(1+exp(zvec))
  gpv = exp(zvec)/(1+exp(zvec))^2
  Lik = sum(-(1-Z)*log(1-pv)-Z*log(pv))
  gLik = (-(1-Z)*1/(1-pv)*(-gpv)-Z*1/pv*gpv)
  q = zvec-mv
  q1 = forwardsolve(t(cholCv),q)
  Prior = 1/2*sum( q1^2)
  gPrior = backsolve(cholCv, q1)
  return(gLik+gPrior)
}
PostCalcDelta <- function(dvec,Xa1,Xa2,Xa3,etaF,etaD,XFs,ZF,ZD) {
  Cpred <- 4*corrmat_error(Xa1,Xa2,Xa3,Xb1 = XFs,Xb2 = etaF+dvec[1:nfield], Xb3 = rep(1,length(XFs)))
  # have to add abs() to avoid negative y's
  hpred <- mvec(XFs,abs(etaF+dvec[1:nfield]))+c(t(Cpred)%*%wvals)
  pv <- exp(hpred) / (1 + exp(hpred))
  Lik <- sum(-(1-ZF)*log(1-pv)-ZF*log(pv))
  # Liklihood for exp data
  pd <- etaD+rep(dvec[nfield+1], nexp)
  Lik <- Lik + 0.5 * t(ZD - pd)%*%diag(rep(0.005^(-2), length(ZD)))%*%(ZD-pd)
  # Prior for bias
  Prior <- sum(rep(160,nfield+1)*(dvec-rep(0,nfield+1))^2)
  return(Lik+Prior)
}
gPostCalcDelta <- function(dvec,Xa1,Xa2,Xa3,etaF,etaD,XFs,ZF,ZD) {
  Cpred <- 4*corrmat_error(Xa1,Xa2,Xa3,Xb1 = XFs,Xb2 = etaF+dvec[1:nfield], Xb3 = rep(1,length(XFs)))
  # Derivative of the matrix
  Dpred <- 4*dcorrmat_error(Xa1,Xa2,Xa3,Xb1 = XFs,Xb2 = etaF+dvec[1:nfield], Xb3 = rep(1,length(XFs)))
  hpred <- mvec(XFs,abs(etaF+dvec[1:nfield]))+t(Cpred)%*%wvals
  dmvec <- 784/sqrt(abs(etaF+dvec[1:nfield])*1000)
  pv <- exp(hpred) / (1 + exp(hpred))
  gpv <- exp(hpred) / (1 + exp(hpred))
  gLik <- (-(1-ZF)*1/(1-pv)*(-gpv)-ZF*1/pv*gpv)*(dmvec + c(t(Dpred) %*%wvals))
  # Liklihood for exp data - add one more element to the gradient
  gLik <- c(gLik, sum(etaD+rep(dvec[nfield+1], nexp)-ZD)/0.005^2)
  # Prior for bias
  gPrior <- sum(rep(160,nfield+1)*(dvec-rep(0,nfield+1))^2) * (-6) * (dvec - rep(0,nfield+1))
  return(gLik+gPrior)
}

# start iterative procedure to estimate the model discrepancy and injury function
twoexpvect = matrix((rep(c(25.47,1.75,35,thetaMAP),1)),nrow=1,byrow=T)
etaF <- exp((GPpred(GP, cbind(B[, 4 : 6], t(matrix(thetaMAP, npara, nfield))))$pred)[,1])
etaD <- rep(exp(GPpred(GPmodel, twoexpvect)$pred[1,1]), nexp)
bias <- rep(0, nfield + 1)
bias_previous <- rep(999, nfield + 1)
iteration <- 0
while (mean((bias-bias_previous)^2) > 10^(-5)) {
	iteration <- iteration + 1
	bias_previous = bias
	# step 1: estimate g
	YF = exp((GPpred(GP, cbind(B[, 4 : 6], t(matrix(thetaMAP, npara, nfield))))$pred)[,1]) + bias[1 : nfield]
	Cv <- 4*corrmat_error(c(XEs,XFs),c(YE,YF),c(Type1,Type2))+0.01*diag(rep(1,length(c(YE,YF))))
	mv <- mvec(c(XEs,XFs),c(YE,YF))
	cholCv <- chol(Cv)
	opt.out = optim(mv,fn <- PostCalc,gr <- gPostCalc, Z = c(ZE,ZF),method = "L-BFGS-B")
	q1 = forwardsolve(t(cholCv), opt.out$par-mv)
	wvals = backsolve(cholCv, q1)
	# step 2: estimate delta
	opt.out2 = optim(
	  rep(0, nfield + 1),
	  fn <- PostCalcDelta,
	  gr <- gPostCalcDelta,
	  lower = rep(-0.1, nfield + 1), 
      upper = rep(0.5, nfield + 1),  
	  Xa1=c(XEs,XFs),Xa2=c(YE,YF),Xa3=c(Type1,Type2),etaF=etaF,etaD=etaD,XFs=XFs,ZF=ZF,ZD=ZD,
	  method = "L-BFGS-B"
	)
	# export bias for field data
	bias=opt.out2$par
}


Xa1=c(XF[, 1], rep(25.47, nexp))
Xa2=c(XF[, 2], rep(1.75, nexp))
Xa3=c(XF[, 3], rep(35, nexp))
Cv2 <- 4*corrmat_x(Xa1,Xa2,Xa3)+0.01*diag(rep(1,nfield+nexp))
cholCv2 <- chol(Cv2)
q2 = forwardsolve(t(cholCv2), c(bias[1:nfield],rep(bias[nfield+1], nexp)))
wvals2 = backsolve(cholCv2, q2)

#this function will be used to predict
papprox <- function(DataPred){
  #introduce model discrepancy
  Cpred0 <- 4*corrmat_x(Xa1, Xa2, Xa3, Xb1 = DataPred$BMI, Xb2 = DataPred$Stature, DataPred$DeltaV)
  YFpred = exp((GPpred(GP, cbind(DataPred$BMI,DataPred$Stature,DataPred$DeltaV, t(matrix(thetaMAP, npara, nfield))))$pred)[,1]) + t(Cpred0)%*%wvals2

  #calculate injury risks
  Cpred <- 4*corrmat_error(c(XEs,XFs),c(YE,YF),c(Type1,Type2),Xb1 = DataPred$Age,Xb2 = YFpred, Xb3 = rep(1,length(DataPred$Age)))
  zpred = mvec(DataPred$Age,YF)+t(Cpred)%*%wvals
  return(exp(zpred)/(1+exp(zpred)))
}

#show prediction with g() to compare with m()
gapprox <- function(age, ystar){
  #calculate injury risks
  Cpred <- 4*corrmat_error(c(XEs,XFs),c(YE,YF),c(Type1,Type2),Xb1 = rep(age, length(ystar)), Xb2 = ystar, Xb3 = rep(1, length(ystar)))
  zpred = mvec(rep(age, length(ystar)),ystar)+t(Cpred)%*%wvals
  return(exp(zpred)/(1+exp(zpred)))
}

#this is an approximation of the CI for the parameters. 
#No identifiability correction was used
print(ApproxCItheta)

#some sample predictions at our field data
DataPred = data.frame("BMI"=B[,4],"Stature"=B[,5],"DeltaV"=B[,6],"Age"=XFs)

# plot(DataPred$BMI,papprox(DataPred))
# plot(DataPred$Stature,papprox(DataPred))
# plot(DataPred$DeltaV,papprox(DataPred))
# plot(DataPred$Age,papprox(DataPred))

# check prediction errors
InjColor = rep('blue', nrow(DataPred))
InjColor[ZF == 1] = 'red'
plot(papprox(DataPred), col = InjColor, pch = 3, ylim = c(0, 1), ylab = 'Predicted injury risks', xlab = 'Field test index', cex = 1.2, cex.lab = 1.5, cex.axis = 1.5)
abline(h=0.5, lty=2, lwd=2)

# compare m() with g()
ystar <- (1:170) * 0.001
h <- mvec1(rep(age, length(ystar)), ystar)
plot(ystar, exp(h)/(1+exp(h)), type = 'n', lwd = 2, lty = 2, ylim = c(0, 1), xlab = 'Injury values', ylab = 'Predicted injury risks', cex = 1.2, cex.lab = 1.5, cex.axis = 1.5)
i <- 0
for (age in c(30, 50, 70)) {
  i <- i + 1
  h <- mvec1(rep(age, length(ystar)), ystar)
  lines(ystar, gapprox(age, ystar), col = rainbow(3)[i], lty = 1, lwd = 2)
  lines(ystar, exp(h)/(1+exp(h)), col = rainbow(3)[i], lwd = 2, lty = 2)
}
