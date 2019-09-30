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

Type1 = rep(0,length(XEs))
Type2 = rep(1,length(XFs))

# read experimental data
D <- read.csv('Experimental Data - All Sedans.csv')
exp_mean <- mean(D[, 9])
exp_std <- sqrt(var(D[, 9]))
ZD = D[, 9]

corrmat <-  function(Xa1,Xa2,Xa3,...,Xb1=NA,Xb2=NA,Xb3=NA){
  if (any(is.na(Xb1))){
    Xb1 = Xa1
    Xb2 = Xa2
    Xb3 = Xa3
  }
  R1 <- exp(-abs(outer(Xa1,Xb1,'-')/40)^1.95)
  R2 <- exp(-abs(outer(sqrt(Xa2),sqrt(Xb2),'-')/0.12)^1.95)
  R3 <- 0.95+0.05*(1-abs(outer(Xa3,Xb3,'-')))
  return(R1*R2*R3)
}

mvec <- function(X,Y){
  return(-9+0.01*X+15*sqrt(Y))
}

PostCalcApprox <- function(x,GPmodel,XF,XFs,ZF,YE,XEs,ZE,ZD,lambda) {
  GP_obj = GPpred(GPmodel, cbind(XF, t(matrix(x, npara, nfield))))
  YF = exp(GP_obj$pred)[,1]
  
  #likellihood for field data
  Cv <- 4*corrmat(c(XEs,XFs),c(YE,YF),c(Type1,Type2))+0.01*diag(rep(1,length(c(YE,YF))))
  mv <- mvec(c(XEs,XFs),sqrt(c(YE,YF)))
  pv <- exp(mv)/(1+exp(mv))
  diagadj = exp(mv)/(1+exp(mv))^2
  Cva = diag(diagadj)%*%Cv%*%diag(diagadj)
  Va = diag(pmax(pv*(1-pv),10^(-4)))
  LikA = 0.5*t(c(ZE,ZF)-pv)%*%solve(Cva+Va,c(ZE,ZF)-pv)
  
  #likelihood for exp data
  twoexpvect = matrix((rep(c(25.47,1.75,35,parameter_default),2)),nrow=2,byrow=T)
  pd = rep(exp(GPpred(GPmodel, twoexpvect)$pred[1,1]), length(ZD))
  LikA = LikA + 0.5 * t(ZD - pd)%*%diag(rep(0.005^(-2), length(ZD)))%*%(ZD-pd)
  
  #penalization - delta method due to log transformation
  Penal = lambda * sqrt(abs(mean(exp(GP_obj$pred[, 1]) ^ 2 * diag(GP_obj$corr_mat) * GP_obj$var_scale[1, 1])))
  
  #prior
  Prior = sum(c(w_default)*(x[1:npara]-parameter_default)^2) #pull it closer to true values
  return(LikA+Prior-Penal)
}

opt.out1 = optim(
  parameter_default,
  fn = PostCalcApprox,
  XF=XF,XFs=XFs,ZF=ZF,YE=YE,XEs=XEs,ZE=ZE,ZD=ZD,GPmodel=GP,lambda=1000,
  lower = c(apply(XCE[,4:(3+npara)],2,min)), 
  upper = c(apply(XCE[,4:(3+npara)],2,max)),
  control=list('factr' = 10^12, 'maxit'=50), 
  hessian=T,
  method = "L-BFGS-B"
)
opt.out1$par
unname(0+apply(XCE[,4:(3+npara)],2,min))
unname(0+apply(XCE[,4:(3+npara)],2,max))

# 
thetaMAP = opt.out1$par
ApproxCItheta = cbind(thetaMAP-1.96*sqrt(diag(opt.out1$hessian^(-1))),thetaMAP+1.96*sqrt(diag(opt.out1$hessian^(-1))))
print(ApproxCItheta)
YF = exp((GPpred(GP, cbind(B[, 4 : 6], t(matrix(thetaMAP, npara, nfield))))$pred)[,1])
Cv <- 4*corrmat(c(XEs,XFs),c(YE,YF),c(Type1,Type2))+0.01*diag(rep(1,length(c(YE,YF))))
mv <- mvec(c(XEs,XFs),sqrt(c(YE,YF)))

cholCv <- chol(Cv)

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
opt.out = optim(mv,fn <- PostCalc,gr <- gPostCalc, Z = c(ZE,ZF),method = "L-BFGS-B")
q1 = forwardsolve(t(cholCv), opt.out$par-mv)
wvals = backsolve(cholCv, q1)

#this function will be used to predict
papprox <- function(DataPred){
  YFpred = exp((GPpred(GP, cbind(DataPred$BMI,DataPred$Stature,DataPred$DeltaV, t(matrix(thetaMAP, npara, nfield))))$pred)[,1])


  Cpred <- 4*corrmat(c(XEs,XFs),c(YE,YF),c(Type1,Type2),Xb1 = DataPred$Age,Xb2 = YFpred, Xb3 = rep(1,length(DataPred$Age)))

  zpred = mvec(DataPred$Age,sqrt(YF))+t(Cpred)%*%wvals
  return(exp(zpred)/(1+exp(zpred)))
}

#this is an approximation of the CI for the parameters. 
#No identifiability correction was used, but we still have boundary issues
print(ApproxCItheta)

#some sample predictions at our field data
DataPred = data.frame("BMI"=B[,4],"Stature"=B[,5],"DeltaV"=B[,6],"Age"=XFs)

plot(DataPred$BMI,papprox(DataPred))
plot(DataPred$Stature,papprox(DataPred))
plot(DataPred$DeltaV,papprox(DataPred))
plot(DataPred$Age,papprox(DataPred))

# check prediction errors
InjPred = rep(0, nrow(DataPred))
InjPred[papprox(DataPred) > 0.3] = 1
plot(InjPred, col = 'red', pch = 3)
points(ZF, col = 'blue', pch = 2)

InjColor = rep('blue', nrow(DataPred))
InjColor[ZF == 1] = 'red'
plot(papprox(DataPred), col = InjColor, pch = 3, ylim = c(0, 1))