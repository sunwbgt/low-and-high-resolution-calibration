library(MASS)
# read simulation data
setwd("E:/Research Projects/Low and High Resolution Calibration/EmuAccord")
source("GP_functions.R")
source("Injury_risk_functions.R")
A = read.csv('Accord_Simulation_Results_LL_Airbag_300_Large_Range.csv', header = TRUE)
# settings
n_var = 3
n_para = 7
parameter_default = c(1, 0.45, 1, 1, 0.55, 1, 1)
w_default = rep(50,length(parameter_default))
# module 1: emulator
n_res = 4
n_sim = 196
X_sim = A[,2:(n_var+n_para+1)]
Y_sim = log(as.matrix(A[,(n_var+n_para+2):ncol(A)]))
GP_sim = GPfitting(X_sim, Y_sim)

# read field data
B = read.csv('Field_Test_Data_Full_Size_Sedan.csv', header = TRUE)
n_field = nrow(B)
X_field = B[, 4:6]
A_field = B$Age
S_field = cbind(rep(1,length(B[,8])),B[,8]) # indicator for experimental condition (field = 1, experiment = 0)
Z_field = 0+(B[,8]>=3)
# predict injuries for field data
Xpred_field = cbind(B[, 4 : 6], t(matrix(parameter_default, n_para, n_field)))

# read injury data
C <- read.csv('injury_function_data_chest.csv', header = TRUE)
A_injury = C$Age
Y_injury = C$ChestD
Y_injury = pmax(Y_injury,10^(-2))
Z_injury = C$Injury.Indicator
n_injury <- length(A_injury)

# read experimental data
D <- read.csv('Experimental Data - All Sedans.csv')
exp_mean <- mean(D[, 9])
exp_std <- sqrt(var(D[, 9]))
Y_exp = D[, 9] 
X_exp = c(25.47, 1.75, 35)
n_exp = length(Y_exp)

# nugget value for matrices
nugget_val <- 0.0001

# specify data sources: injury test -> 0, field data -> 1, experimental test -> 2
Type1 = rep(0,length(A_injury))
Type2 = rep(1,length(A_field))
Type3 = rep(2,length(Y_exp))

## module 2: Laplace approximation in the injury dataset
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
corrmat_error <-  function(lambda_J, Xa1,Xa2,Xa3,...,Xb1=NA,Xb2=NA,Xb3=NA){
  # specify standard deviation of y in injury tests: 5% of the population standard deviation
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
	    R2[i, j] = exp(-((Xa2[i] - Xb2[j])^2/(0.035+4/lambda_J)))/(1+4/0.035/lambda_J)^(1/2)
      }
      # only one from injury test
      if (((Xa3[i] == 0) && (Xb3[j] > 0)) || ((Xa3[i] > 0) && (Xb3[j] == 0))) {
	    R2[i, j] = exp(-((Xa2[i] - Xb2[j])^2/(0.035+2/lambda_J)))/(1+2/0.035/lambda_J)^(1/2)
      }
      # none from injury test
      if ((Xa3[i] > 0) && (Xb3[j] > 0)) {
        R2[i, j] = exp(-((Xa2[i] - Xb2[j])^2/0.035))
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
dcorrmat_error <-  function(lambda_J, Xa1,Xa2,Xa3,...,Xb1=NA,Xb2=NA,Xb3=NA){
  # specify standard deviation of y in injury tests: 5% of the population standard deviation
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
		R2[i, j] = 0
	    # R2[i, j] = exp(-((Xa2[i] - Xb2[j])^2/(0.035+4/lambda_J)))/(1+4/0.035/lambda_J)^(1/2) * 2*(Xa2[i]-Xb2[j])/(0.035+4/lambda_J)*(-1)
      }
      # only one from injury test
      if (((Xa3[i] == 0) && (Xb3[j] > 0)) || ((Xa3[i] > 0) && (Xb3[j] == 0))) {
	    R2[i, j] = exp(-((Xa2[i] - Xb2[j])^2/(0.035+2/lambda_J)))/(1+2/0.035/lambda_J)^(1/2) * 2*(Xa2[i]-Xb2[j])/(0.035+2/lambda_J)*(-1)
      }
      # none from injury test
      if ((Xa3[i] > 0) && (Xb3[j] > 0)) {
	    R2[i, j] = 0
        # R2[i, j] = exp(-((Xa2[i] - Xb2[j])^2/0.035)) * 2*(Xa2[i]-Xb2[j])/0.035*(-1)
      }
	}
  }
  # data with labels 0 or 2 come from the same lab environment (share the same type)
  Xa3[Xa3 == 2] = 0
  Xb3[Xb3 == 2] = 0
  R3 <- 0.95+0.05*(1-abs(outer(Xa3,Xb3,'-')))
  return(R1*R2*R3)
}

# second derivative of correlation matrix with respect to delta for error-in-variable GPs
ddcorrmat_error <-  function(lambda_J, Xa1,Xa2,Xa3,...,Xb1=NA,Xb2=NA,Xb3=NA){
  # specify standard deviation of y in injury tests: 5% of the population standard deviation
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
		# if (i >= j) {
			# R2[i, j] = -2 / (1 + 4 / 0.035 / lambda_J)^(1/2) / (0.035 + 4 / lambda_J) * exp(-(Xa2[i] - Xb2[j]) ^ 2 / (0.035 + 4 / lambda_J)) * (1 - 2 * (Xa2[i] - Xb2[j]) ^ 2 / (0.035 + 4 / lambda_J))
		# }
		# if (i < j) {
			# R2[i, j] = -2 / (1 + 4 / 0.035 / lambda_J)^(1/2) / (0.035 + 4 / lambda_J) * exp(-(Xa2[i] - Xb2[j]) ^ 2 / (0.035 + 4 / lambda_J)) * (-1 - 2 * (Xa2[i] - Xb2[j]) ^ 2 / (0.035 + 4 / lambda_J))
		# }
		R2[i, j] = 0
	  }
      # only one from injury test # not symmetric for second order
      if ((Xa3[i] > 0) && (Xb3[j] == 0)) {
		R2[i, j] = -2 / (1 + 2 / 0.035 / lambda_J)^(1/2) / (0.035 + 2 / lambda_J) * exp(-(Xa2[i] - Xb2[j]) ^ 2 / (0.035 + 2 / lambda_J)) * (1 - 2 * (Xa2[i] - Xb2[j]) ^ 2 / (0.035 + 2 / lambda_J))
	  }
      if ((Xa3[i] == 0) && (Xb3[j] > 0)) {
		R2[i, j] = -2 / (1 + 2 / 0.035 / lambda_J)^(1/2) / (0.035 + 2 / lambda_J) * exp(-(Xa2[i] - Xb2[j]) ^ 2 / (0.035 + 2 / lambda_J)) * (-1 - 2 * (Xa2[i] - Xb2[j]) ^ 2 / (0.035 + 2 / lambda_J))
	  }
      # none from injury test
      if ((Xa3[i] > 0) && (Xb3[j] > 0)) {
		# if (i >= j) {
			# R2[i, j] = -2 / 0.035 * exp(-(Xa2[i] - Xb2[j]) ^ 2 / 0.035) * (1 - 2 * (Xa2[i] - Xb2[j]) ^ 2 / 0.035)
		# }
		# if (i < j) {
			# R2[i, j] = -2 / 0.035 * exp(-(Xa2[i] - Xb2[j]) ^ 2 / 0.035) * (-1 - 2 * (Xa2[i] - Xb2[j]) ^ 2 / 0.035)
		# }
		R2[i, j] = 0
      }
	}
  }
  # data with labels 0 or 2 come from the same lab environment (share the same type)
  Xa3[Xa3 == 2] = 0
  Xb3[Xb3 == 2] = 0
  R3 <- 0.95+0.05*(1-abs(outer(Xa3,Xb3,'-')))
  return(R1*R2*R3)
}

# mean function - use the one provided in the literature - do not use sqrt anywhere else
mvec <- function(X,Y){
  return(-12.597+0.05861*X+1.568*sqrt(abs(Y * 1000)))
}

dmvec <- function(Y){
  return(0.5*1.568*sqrt(1000)*(abs(Y))^(-1/2))
}

ddmvec <- function(Y) {
  return(0.5*1.568*sqrt(1000)*(-0.5)*(abs(Y))^(-3/2))
}

# psi function
psi_J <- function(h, lambda_J, lambda_h, Z, Y, A) {
	K_J <- 4 * corrmat_error(lambda_J, A, Y, Type1)
	# take abs() of eigen values because there will be some extremely small negative numbers
	return(-sum(Z * log(exp(h) / (1 + exp(h))) + (1 - Z) * log(1 / (1 + exp(h)))) + lambda_h / 2 * t(h - mv_injury) %*% solve(K_J + diag(nugget_val, length(h))) %*% (h - mv_injury))
}

# derivative of psi function
dpsi_J <- function(h, lambda_J, lambda_h, Z, Y, A) {
	K_J <- 4 * corrmat_error(lambda_J, A_injury, Y_injury, Type1)
	return(-Z / (1 + exp(h)) + (1 - Z) * exp(h) / (1 + exp(h)) + lambda_h * t(h - mv_injury) %*% solve(K_J + diag(nugget_val, length(h))))
}

# grid search the best lambda_J and lambda_h
mv_injury <- mvec(A_injury, Y_injury)
lambda_J = 0
lambda_h = 0
max_likelihood = -10^20
for (lambda_J_temp in 5^(1 * {-1:10})) {
	for (lambda_h_temp in 5^(1 * (-1:10))) {
		opt.out = optim(
			rep(0, n_injury),
			fn <- psi_J,
			gr <- dpsi_J,
			lower = rep(-20, n_injury), 
			upper = rep(20, n_injury),  
			lambda_J = lambda_J_temp,
			lambda_h = lambda_h_temp,
			Z = Z_injury,
			Y = Y_injury,
			A = A_injury,
			method = "L-BFGS-B"
		)
		h_hat <- opt.out$par
		W_J <- diag(exp(h_hat) / (1 + exp(h_hat)) ^ 2)
		K_J <- 4 * corrmat_error(lambda_J_temp, A_injury, Y_injury, Type1)
		temp_likelihood <- -lambda_h_temp / 2 * t(h_hat - mv_injury) %*% solve(K_J + diag(nugget_val, length(h_hat))) %*% (h_hat - mv_injury) + sum(log(exp(h_hat * Z_injury) / (1 + exp(h_hat)))) - 0.5 * sum(log(abs(eigen(diag(n_injury) + sqrt(W_J) %*% K_J %*% sqrt(W_J) / lambda_h_temp)$values)))
		if (temp_likelihood >= max_likelihood) {
			max_likelihood <- temp_likelihood
			lambda_J <- lambda_J_temp
			lambda_h <- lambda_h_temp
		}
		# print(c(lambda_h_temp, lambda_J_temp, temp_likelihood, -lambda_h_temp / 2 * t(h_hat - mv_injury) %*% solve(K_J + diag(10 ^ (-4), length(h_hat))) %*% (h_hat - mv_injury), sum(log(exp(h_hat * Z_injury) / (1 + exp(h_hat)))), -0.5 * sum(log(abs(eigen(diag(n_injury) + sqrt(W_J) %*% K_J %*% sqrt(W_J) / lambda_h_temp)$values)))))
	}
}

# optimize h_hat for Laplace approximation
opt.out = optim(
	rep(0, n_injury),
	fn <- psi_J,
	gr <- dpsi_J,
	lower = rep(-20, n_injury), 
	upper = rep(20, n_injury),  
	lambda_J = lambda_J,
	lambda_h = lambda_h,
	Z = Z_injury,
	Y = Y_injury,
	A = A_injury,
	method = "L-BFGS-B"
)
h_hat = opt.out$par
W_J = diag(exp(h_hat) / (1 + exp(h_hat)) ^ 2)

# test: plot the predicted probabilites with the labels in the injury dataset
p_hat = exp(h_hat) / (1 + exp(h_hat))
m_hat = exp(mv_injury) / (1 + exp(mv_injury))
plot(1 : n_injury, Z_injury[order(Z_injury)], col = 'red', cex = 1.5, pch = 3, cex.lab = 1.5, cex.axis = 1.5, xlab = 'sorted samples', ylab = 'probabilities')
points(1 : n_injury, p_hat[order(Z_injury)], col = 'blue', cex = 1.5, pch = 2)
points(1 : n_injury, m_hat[order(Z_injury)], col = 'orange', cex = 1.5, pch = 2)
abline(h = 0.5, lwd = 1.5, lty = 3)

## module 3: Laplace approximation for field data
# covariance matrix for zeta(x)
# correlation matrix for GP of delta()
corrmat_x <- function(Xa1,Xa2,Xa3,...,Xb1=NA,Xb2=NA,Xb3=NA){
  if (any(is.na(Xb1))){
    Xb1 = Xa1
    Xb2 = Xa2
    Xb3 = Xa3
  }
  R1 <- exp(-abs(outer(Xa1,Xb1,'-')/0.4)^1.95)
  R2 <- exp(-abs(outer(Xa2,Xb2,'-')/0.6)^1.95)
  R3 <- exp(-abs(outer(Xa2,Xb2,'-')/0.6)^1.95)
  return(R1*R2*R3)
}

# psi function for zeta
psi_zeta <- function(zeta, lambda_delta, lambda_J, lambda_h, Z_F, X_F, A_F, Type2, Y_J, A_J, Type1, temp_mat_inverse_w, temp_mat_inverse, h_hat) {
	K_zeta <- 4 * corrmat_x(X_F[, 1], X_F[, 2], X_F[, 3])
	mat_cross <- t(4 * corrmat_error(lambda_J, Xa1 = A_F, Xa2 = zeta, Xa3 = Type2, Xb1 = A_J, Xb2 = Y_J, Xb3 = Type1))
	temp_mat <- 4 * corrmat_error(lambda_J, Xa1 = A_F, Xa2 = zeta, Xa3 = Type2)
	mu_hat <- mvec(A_F, zeta) + t(mat_cross) %*% temp_mat_inverse %*% (h_hat - mvec(A_J, Y_J))
	V_hat <- diag(temp_mat - t(mat_cross) %*% temp_mat_inverse_w %*% mat_cross)	
	
	temp_value <- pmax(rep(1, n_field) + pi / 8 / lambda_h * V_hat, 10 ^ (-4)) ^ (-1 / 2) * mu_hat
	lik <- -sum(Z_F * log(exp(temp_value) / (exp(temp_value) + 1)) + (1 - Z_F) * log(1 / (exp(temp_value) + 1)))
	prior <- lambda_delta / 2 * t(zeta - eta_F) %*% solve(K_zeta + diag(nugget_val, n_field)) %*% (zeta - eta_F)
	return(lik + prior)
}

# derivative of psi function for zeta
dpsi_zeta <- function(zeta, lambda_delta, lambda_J, lambda_h, Z_F, X_F, A_F, Type2, Y_J, A_J, Type1, temp_mat_inverse_w, temp_mat_inverse, h_hat) {
	K_zeta <- 4 * corrmat_x(X_F[, 1], X_F[, 2], X_F[, 3])
	mat_cross <- t(4 * corrmat_error(lambda_J, Xa1 = A_F, Xa2 = zeta, Xa3 = Type2, Xb1 = A_J, Xb2 = Y_J, Xb3 = Type1))
	temp_mat <- 4 * corrmat_error(lambda_J, Xa1 = A_F, Xa2 = zeta, Xa3 = Type2)
	V_hat <- temp_mat - t(mat_cross) %*% temp_mat_inverse_w %*% mat_cross
	mu_hat <- mvec(A_F, zeta) + t(mat_cross) %*% temp_mat_inverse %*% (h_hat - mvec(A_J, Y_J))
	dmat <- t(4 * dcorrmat_error(lambda_J, Xa1 = A_F, Xa2 = zeta, Xa3 = Type2))
	dmat_cross <- t(4 * dcorrmat_error(lambda_J, Xa1 = A_F, Xa2 = zeta, Xa3 = Type2, Xb1 = A_J, Xb2 = Y_J, Xb3 = Type1))
	
	# chain rule
	first <- pmax(rep(1, n_field) + pi / 8 / lambda_h * diag(V_hat), nugget_val) ^ (-1 / 2)
	second <- mu_hat
	dfirst_v <- rep(0, n_field)
	dsecond <- rep(0, n_field)
	for (i in 1 : n_field) {
		dfirst_v[i] <- dmat[i, i] - 2 * t(mat_cross[, i]) %*% temp_mat_inverse_w %*% dmat_cross[, i]
		dsecond[i] <- dmvec(zeta[i]) + t(dmat_cross[, i]) %*% temp_mat_inverse %*% (h_hat - mvec(A_J, Y_J))	
	}
	dfirst <- -pi / 16 / lambda_h * pmax(rep(1, n_field) + pi / 8 / lambda_h * diag(V_hat), nugget_val) ^ (-3 / 2) * dfirst_v

	# calculate the derivative
	glik <- -Z_F / (1 + exp(first * second)) * (first * dsecond + second * dfirst) + (1 - Z_F) * exp(first * second) / (1 + exp(first * second)) * (first * dsecond + second * dfirst)
	gprior <- t(lambda_delta * t(zeta - eta_F) %*% solve(K_zeta + diag(nugget_val, n_field)))
	return(glik + gprior)
}

# grid search the best lambda_delta
temp_mat_inverse_w <- solve(4 * corrmat_error(lambda_J, Xa1 = A_injury, Xa2 = Y_injury, Xa3 = Type1) + solve(W_J + diag(nugget_val, n_injury)))
temp_mat_inverse <- solve(4 * corrmat_error(lambda_J, Xa1 = A_injury, Xa2 = Y_injury, Xa3 = Type1) + diag(nugget_val, n_injury))
eta_F <- exp(GPpred(GP_sim, cbind(X_field, t(matrix(parameter_default, n_para, n_field))))$pred)[, 1]
max_likelihood = -10^20
for (lambda_delta_temp in 10^(1 * {0:20})) {
	opt.out = optim(
		eta_F,
		fn <- psi_zeta,
		gr <- dpsi_zeta,
		lower = rep(0.001, n_injury), 
		upper = rep(0.40, n_injury),  
		lambda_delta = lambda_delta_temp,
		lambda_J = lambda_J,
		lambda_h = lambda_h,
		Z_F = Z_field,
		X_F = X_field,
		A_F = A_field,
		Type2 = Type2,
		Y_J = Y_injury,
		A_J = A_injury,
		Type1 = Type1,
		temp_mat_inverse_w = temp_mat_inverse_w,
		temp_mat_inverse = temp_mat_inverse,
		h_hat = h_hat,
		method = "L-BFGS-B"
	)
	zeta_hat = opt.out$par
	# calculate W_zeta for Laplace approximation (second derivative of log p)
	K_zeta <- 4 * corrmat_x(X_field[, 1], X_field[, 2], X_field[, 3])
	mat_cross <- t(4 * corrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2, Xb1 = A_injury, Xb2 = Y_injury, Xb3 = Type1))
	dmat <- t(4 * dcorrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2))
	temp_mat <- 4 * corrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2)
	V_hat <- temp_mat - t(mat_cross) %*% temp_mat_inverse_w %*% mat_cross
	mu_hat <- mvec(A_field, zeta_hat) + t(mat_cross) %*% temp_mat_inverse %*% (h_hat - mvec(A_injury, Y_injury))
	dmat_cross <- t(4 * dcorrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2, Xb1 = A_injury, Xb2 = Y_injury, Xb3 = Type1))
	
	# chain rule
	first <- pmax(rep(1, n_field) + pi / 8 / lambda_h * diag(V_hat), nugget_val) ^ (-1 / 2)
	second <- mu_hat
	dfirst_v <- rep(0, n_field)
	dsecond <- rep(0, n_field)
	for (i in 1 : n_field) {
		dfirst_v[i] <- dmat[i, i] - 2 * t(mat_cross[, i]) %*% temp_mat_inverse_w %*% dmat_cross[, i]
		dsecond[i] <- dmvec(zeta_hat[i]) + t(dmat_cross[, i]) %*% temp_mat_inverse %*% (h_hat - mvec(A_injury, Y_injury))
	}
	dfirst <- -pi / 16 / lambda_h * pmax(rep(1, n_field) + pi / 8 / lambda_h * diag(V_hat), nugget_val) ^ (-3 / 2) * dfirst_v
	
	# second orders
	ddmat_cross <- t(4 * ddcorrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2, Xb1 = A_injury, Xb2 = Y_injury, Xb3 = Type1))
	ddmat <- t(4 * ddcorrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2))
	ddfirst_v <- rep(0, n_field)
	ddsecond <- rep(0, n_field)
	for (i in 1 : n_field) {
		ddfirst_v[i] <- ddmat[i, i] - 2 * t(dmat_cross[, i]) %*% temp_mat_inverse_w %*% dmat_cross[, i] - 2 * t(mat_cross[, i]) %*% temp_mat_inverse_w %*% ddmat_cross[, i]
		ddsecond[i] <- ddmvec(zeta_hat[i]) + t(ddmat_cross[, i]) %*% temp_mat_inverse %*% (h_hat - mvec(A_injury, Y_injury))
	}
	ddfirst <- -pi / 16 / lambda_h * (-pi * 3 / 16 / lambda_h * pmax(rep(1, n_field) + pi / 8 / lambda_h * diag(V_hat), nugget_val) ^ (-5/2) * dfirst_v + (pmax(rep(1, n_field) + pi / 8 / lambda_h * diag(V_hat), nugget_val)) ^ (-3/2) * ddfirst_v)
	W_value <- (first * dsecond + second * dfirst) ^ 2 * exp(first * second) / (1 + exp(first * second)) ^ 2 + (Z_field - (1 - Z_field) * exp(first * second)) / (1 + exp(first * second)) * (ddfirst * second + 2 * dfirst * dsecond + first * ddsecond)
	W_zeta <- diag(as.vector(W_value))
	GP_obj_field <- GPpred(GP_sim, cbind(X_field, t(matrix(parameter_default, n_para, n_field))))
	# calculate marginal likelihood for comparison
	lik <- sum(Z_field * log(exp(first * second) / (exp(first * second) + 1)) + (1 - Z_field) * log(1 / (exp(first * second) + 1)))
	prior1 <- -0.5 * t(zeta_hat - eta_F) %*% solve(K_zeta / lambda_delta_temp + GP_obj_field$corr_mat * GP_obj_field$var_scale[1, 1] + diag(nugget_val, n_field)) %*% (zeta_hat - eta_F)
	prior2 <- -0.5 * sum(log(abs(eigen(K_zeta / lambda_delta_temp + GP_obj_field$corr_mat * GP_obj_field$var_scale[1, 1])$values))) - 0.5 * sum(log(abs(eigen(solve(K_zeta / lambda_delta_temp + GP_obj_field$corr_mat * GP_obj_field$var_scale[1, 1] + diag(nugget_val, n_field)) + W_zeta)$values)))
	temp_likelihood <- lik + prior1 + prior2 
	print(c(lambda_delta_temp, temp_likelihood, lik, prior1, prior2))
	if (temp_likelihood >= max_likelihood) {
		max_likelihood <- temp_likelihood
		lambda_delta <- lambda_delta_temp
	}
}

# optimize zeta_hat for Laplace approximation
opt.out = optim(
	eta_F,
	fn <- psi_zeta,
	gr <- dpsi_zeta,
	lower = rep(0.001, n_injury), 
	upper = rep(0.40, n_injury),  
	lambda_delta = lambda_delta,
	lambda_J = lambda_J,
	lambda_h = lambda_h,
	Z_F = Z_field,
	X_F = X_field,
	A_F = A_field,
	Type2 = Type2,
	Y_J = Y_injury,
	A_J = A_injury,
	Type1 = Type1,
	temp_mat_inverse_w = temp_mat_inverse_w,
	temp_mat_inverse = temp_mat_inverse,
	h_hat = h_hat,
	method = "L-BFGS-B"
)
zeta_hat = opt.out$par

# calculate W_zeta for Laplace approximation (second derivative of log p)
K_zeta <- 4 * corrmat_x(X_field[, 1], X_field[, 2], X_field[, 3])
mat_cross <- t(4 * corrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2, Xb1 = A_injury, Xb2 = Y_injury, Xb3 = Type1))
dmat <- t(4 * dcorrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2))
temp_mat <- 4 * corrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2)
V_hat <- temp_mat - t(mat_cross) %*% temp_mat_inverse_w %*% mat_cross
mu_hat <- mvec(A_field, zeta_hat) + t(mat_cross) %*% temp_mat_inverse %*% (h_hat - mvec(A_injury, Y_injury))
dmat_cross <- t(4 * dcorrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2, Xb1 = A_injury, Xb2 = Y_injury, Xb3 = Type1))
	
# test: plot the predicted probabilites with the labels in the field dataset
p_hat_0 <- exp(mvec(A_field, eta_F)) / (1 + exp(mvec(A_field, eta_F)))
p_hat = exp(mu_hat) / (1 + exp(mu_hat))
plot(1 : n_field, Z_field[order(Z_field)], col = 'red', cex = 1.5, pch = 3, cex.lab = 1.5, cex.axis = 1.5, xlab = 'sorted samples', ylab = 'probabilities')
points(1 : n_field, p_hat[order(Z_field)], col = 'blue', cex = 1.5, pch = 2)
points(1 : n_field, p_hat_0[order(Z_field)], col = 'orange', cex = 1.5, pch = 2)
abline(h = 0.5, lwd = 1.5, lty = 3)	

# chain rule
first <- pmax(rep(1, n_field) + pi / 8 / lambda_h * diag(V_hat), nugget_val) ^ (-1 / 2)
second <- mu_hat
dfirst_v <- rep(0, n_field)
dsecond <- rep(0, n_field)
for (i in 1 : n_field) {
	dfirst_v[i] <- dmat[i, i] - 2 * t(mat_cross[, i]) %*% temp_mat_inverse_w %*% dmat_cross[, i]
	dsecond[i] <- dmvec(zeta_hat[i]) + t(dmat_cross[, i]) %*% temp_mat_inverse %*% (h_hat - mvec(A_injury, Y_injury))
}
dfirst <- -pi / 16 / lambda_h * pmax(rep(1, n_field) + pi / 8 / lambda_h * diag(V_hat), nugget_val) ^ (-3 / 2) * dfirst_v
	
# second orders
ddmat_cross <- t(4 * ddcorrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2, Xb1 = A_injury, Xb2 = Y_injury, Xb3 = Type1))
ddmat <- t(4 * ddcorrmat_error(lambda_J, Xa1 = A_field, Xa2 = zeta_hat, Xa3 = Type2))
ddfirst_v <- rep(0, n_field)
ddsecond <- rep(0, n_field)
for (i in 1 : n_field) {
	ddfirst_v[i] <- ddmat[i, i] - 2 * t(dmat_cross[, i]) %*% temp_mat_inverse_w %*% dmat_cross[, i] - 2 * t(mat_cross[, i]) %*% temp_mat_inverse_w %*% ddmat_cross[, i]
	ddsecond[i] <- ddmvec(zeta_hat[i]) + t(ddmat_cross[, i]) %*% temp_mat_inverse %*% (h_hat - mvec(A_injury, Y_injury))
}
ddfirst <- -pi / 16 / lambda_h * (-pi * 3 / 16 / lambda_h * pmax(rep(1, n_field) + pi / 8 / lambda_h * diag(V_hat), nugget_val) ^ (-5/2) * dfirst_v + (pmax(rep(1, n_field) + pi / 8 / lambda_h * diag(V_hat), nugget_val)) ^ (-3/2) * ddfirst_v)
W_value <- (first * dsecond + second * dfirst) ^ 2 * exp(first * second) / (1 + exp(first * second)) ^ 2 + (Z_field - (1 - Z_field) * exp(first * second)) / (1 + exp(first * second)) * (ddfirst * second + 2 * dfirst * dsecond + first * ddsecond)
W_zeta <- diag(as.vector(W_value))

## Gibbs sampler
n_iter <- 1000
lambda_delta_gibbs <- rep(0, n.iter)
lambda_h_gibbs <- rep(0, n.iter)
lambda_J_gibbs <- rep(0, n.iter)
lambda_E_gibbs <- rep(0, n.iter)
theta_gibbs <- matrix(0, n.iter, n_para)
theta_w <- 50
# initial values
lambda_delta_gibbs[1] <- lambda_delta
lambda_h_gibbs[1] <- lambda_h
lambda_J_gibbs[1] <- lambda_J
lambda_E_gibbs[1] <- 1000000
theta_gibbs[1, ] <- parameter_default
theta_lower <- apply(A[, 5 : 11], 2, min)
theta_upper <- apply(A[, 5 : 11], 2, max)

# calculate some constant matrices
mat_cross <- t(4 * corrmat_x(Xa1 = X_exp[1], Xa2 = X_exp[2], Xa3 = X_exp[3], Xb1 = X_field[, 1], Xb2 = X_field[, 2], Xb3 = X_field[, 3]))
K_zeta <- 4 * corrmat_x(X_field[, 1], X_field[, 2], X_field[, 3])
# start MCMC
for (i in 2 : n_iter) {
	# emulate in X_field
	GP_obj_field <- GPpred(GP_sim, cbind(X_field, t(matrix(theta_gibbs[i - 1, ], n_para, n_field))))
	delta_hat_field <- zeta_hat - exp(GP_obj_field$pred)[, 1]
	# emulate in X_exp
	GP_obj_exp <- GPpred(GP_sim, cbind(matrix(X_exp, n_exp, 3), t(matrix(theta_gibbs[i - 1, ], n_para, n_exp))))
	delta_hat_exp <- t(mat_cross) %*% solve(K_zeta + diag(nugget_val, n_field)) %*% delta_hat_field 
	V_hat_eta_exp <- exp(GP_obj_exp$pred)[1, 1] ^ 2 * GP_obj_exp$corr_mat[1, 1] * GP_obj_exp$var_scale[1,1]
	V_hat_delta_exp <- 4 - t(mat_cross) %*% solve(K_zeta + diag(nugget_val, n_field)) %*% mat_cross
	V_hat_eta_field <- diag(exp(GP_obj_field$pred)[, 1]) %*% (GP_obj_field$corr_mat * GP_obj_field$var_scale[1,1]) %*% diag(exp(GP_obj_field$pred)[, 1])
	# update lambda_delta
	flag <- TRUE
	while (flag) {
		ig <- 1 / rgamma(1, n_exp / 4 + 1 / 2, 1 / 2 * sum((Y_exp - exp(GP_obj_exp$pred)[, 1] - rep(delta_hat_exp, n_exp)) ^ 2))
		lambda_delta_propose <- as.numeric((ig - 1 / lambda_E_gibbs[i - 1] - V_hat_eta_exp) ^ (-1) * V_hat_delta_exp)
		if (lambda_delta_propose > 0) {
			flag <- FALSE
		}
	}	
	det_nom <- 0.5 * sum(log(abs(eigen(solve(K_zeta / lambda_delta_gibbs[i - 1] + V_hat_eta_field + diag(nugget_val, n_field)) + W_zeta)$values))) + 0.5 * sum(log(abs(eigen(K_zeta / lambda_delta_gibbs[i - 1] + V_hat_eta_field)$values)))
	det_denom <- 0.5 * sum(log(abs(eigen(solve(K_zeta / lambda_delta_propose + V_hat_eta_field + diag(nugget_val, n_field)) + W_zeta)$values))) + 0.5 * sum(log(abs(eigen(K_zeta / lambda_delta_propose + V_hat_eta_field)$values)))
	log_nom <- det_nom - 1 / 2 * t(zeta_hat - delta_hat_field) %*% solve(K_zeta / lambda_delta_propose + V_hat_eta_field + diag(nugget_val, n_field)) %*% (zeta_hat - delta_hat_field)
	log_denom <- det_denom - 1 / 2 * t(zeta_hat - delta_hat_field) %*% solve(K_zeta / lambda_delta_gibbs[i - 1] + V_hat_eta_field + diag(nugget_val, n_field)) %*% (zeta_hat - delta_hat_field)
	rho <- min(1, exp(log_nom - log_denom))
	temp <- runif(1, 0, 1)
	if (temp <= rho) {
		lambda_delta_gibbs[i] <- lambda_delta_propose
	}
	if (temp > rho) {
		lambda_delta_gibbs[i] <- lambda_delta_gibbs[i - 1]
	}

	# update lambda_h
	K_J <- 4 * corrmat_error(lambda_J_gibbs[i - 1], A_injury, Y_injury, Type1)
	lambda_h_propose <- rgamma(1, 1, 1 / 2 * t(h_hat - mvec(A_injury, Y_injury)) %*% solve(K_J + diag(nugget_val, n_injury)) %*% (h_hat - mvec(A_injury, Y_injury)))
	first_nom <- pmax(rep(1, n_field) + pi / 8 / lambda_h_propose * diag(V_hat), nugget_val) ^ (-1 / 2)
	first_denom <- pmax(rep(1, n_field) + pi / 8 / lambda_h_propose * diag(V_hat), nugget_val) ^ (-1 / 2)
	log_nom <- 0.5 * sum(log(abs(eigen(diag(n_injury) + sqrt(W_J) %*% K_J %*% sqrt(W_J) / lambda_h_gibbs[i - 1])$values))) + sum(log(exp(first_nom * second) / (1 + exp(first_nom * second)))) + lambda_h_gibbs[i - 1]
	log_denom <- 0.5 * sum(log(abs(eigen(diag(n_injury) + sqrt(W_J) %*% K_J %*% sqrt(W_J) / lambda_h_propose)$values))) + sum(log(exp(first_denom * second) / (1 + exp(first_denom * second)))) + lambda_h_propose
	rho <- min(1, exp(log_nom - log_denom))
	temp <- runif(1, 0, 1)
	if (temp <= rho) {
		lambda_h_gibbs[i] <- lambda_h_propose
	}
	if (temp > rho) {
		lambda_h_gibbs[i] <- lambda_h_gibbs[i - 1]
	}
		
	# update lambda_J
	ig <- runif(1, 10, 15)
	lambda_J_propose <- exp(ig)
	K_J_nom <- 4 * corrmat_error(lambda_J_propose, A_injury, Y_injury, Type1)
	K_J_denom <- 4 * corrmat_error(lambda_J_gibbs[i - 1], A_injury, Y_injury, Type1)
	log_nom <- 0.5 * sum(log(abs(eigen(diag(n_injury) + sqrt(W_J) %*% K_J_denom %*% sqrt(W_J) / lambda_h_gibbs[i])$values))) - lambda_h_gibbs[i] / 2 * t(h_hat - mvec(A_injury, Y_injury)) %*% solve(K_J_nom + diag(nugget_val, n_injury)) %*% (h_hat - mvec(A_injury, Y_injury))
	log_denom <- 0.5 * sum(log(abs(eigen(diag(n_injury) + sqrt(W_J) %*% K_J_nom %*% sqrt(W_J) / lambda_h_gibbs[i])$values))) - lambda_h_gibbs[i] / 2 * t(h_hat - mvec(A_injury, Y_injury)) %*% solve(K_J_denom + diag(nugget_val, n_injury)) %*% (h_hat - mvec(A_injury, Y_injury))
	rho <- min(1, exp(log_nom - log_denom))
	temp <- runif(1, 0, 1)
	if (temp <= rho) {
		lambda_J_gibbs[i] <- lambda_J_propose
	}
	if (temp > rho) {
		lambda_J_gibbs[i] <- lambda_J_gibbs[i - 1]
	}	
	
	# update lambda_E_gibbs
	flag <- TRUE
	while (flag) {
		ig <- 1 / rgamma(1, n_exp / 4 + 1 / 2, 1 / 2 * sum((Y_exp - exp(GP_obj_exp$pred)[, 1] - rep(delta_hat_exp, n_exp)) ^ 2))
		lambda_E_propose <- as.numeric((ig - V_hat_delta_exp / lambda_delta_gibbs[i] - V_hat_eta_exp) ^ (-1))
		if (lambda_E_propose > 0) {
			flag <- FALSE
		}
	}
	rho <- min(1, lambda_E_gibbs[i - 1] / lambda_E_propose)
	temp <- runif(1, 0, 1)
	if (temp <= rho) {
		lambda_E_gibbs[i] <- lambda_E_propose
	}
	if (temp > rho) {
		lambda_E_gibbs[i] <- lambda_E_gibbs[i - 1]
	}
	
	# update theta
	# flag <- TRUE
	# while (flag) {
		# theta_propose <- mvrnorm(1, parameter_default, diag(rep(1 / theta_w, n_para)))
		# if ((min(theta_propose - theta_lower) >= 0) && (min(theta_upper - theta_propose) >= 0)) {
			# flag <- FALSE
		# }
	# }
	theta_propose <- rep(0, 7)
	for (j in 1 : 7) {
		theta_propose[j] <- runif(1, theta_lower[j], theta_upper[j])
	}
	# calculate with theta_propose
	GP_obj_field <- GPpred(GP_sim, cbind(X_field, t(matrix(theta_propose, n_para, n_field))))
	GP_obj_exp <- GPpred(GP_sim, cbind(matrix(X_exp, n_exp, 3), t(matrix(theta_propose, n_para, n_exp))))
	delta_hat_field <- zeta_hat - exp(GP_obj_field$pred)[, 1]	
	V_hat_eta_exp <- exp(GP_obj_exp$pred)[1, 1] ^ 2 * GP_obj_exp$corr_mat[1, 1] * GP_obj_exp$var_scale[1,1]
	V_hat_delta_exp <- 4 - t(mat_cross) %*% solve(K_zeta + diag(nugget_val, n_field)) %*% mat_cross
	V_hat_eta_field <- diag(exp(GP_obj_field$pred)[, 1]) %*% (GP_obj_field$corr_mat * GP_obj_field$var_scale[1,1]) %*% diag(exp(GP_obj_field$pred)[, 1])
	delta_hat_exp <- t(mat_cross) %*% solve(K_zeta + diag(nugget_val, n_field)) %*% delta_hat_field 
	part1_nom <- -1 / 2 * t(zeta_hat - delta_hat_field) %*% solve(K_zeta / lambda_delta_gibbs[i] + V_hat_eta_field + diag(nugget_val, n_field)) %*% (zeta_hat - delta_hat_field)
	part2_nom <- -0.5 * sum(log(abs(eigen(solve(K_zeta / lambda_delta_gibbs[i] + V_hat_eta_field + diag(nugget_val, n_field)) + W_zeta)$values))) - 0.5 * sum(log(abs(eigen(K_zeta / lambda_delta_gibbs[i] + V_hat_eta_field)$values)))
	part3_nom <- -n_exp / 2 * log(1 / lambda_E_gibbs[i] + V_hat_delta_exp / lambda_delta_gibbs[i] + V_hat_eta_exp)
	part4_nom <- -1 / 2 / (1 / lambda_E_gibbs[i] + V_hat_delta_exp / lambda_delta_gibbs[i] + V_hat_eta_exp) * sum((Y_exp - exp(GP_obj_exp$pred)[, 1] - rep(delta_hat_exp, n_exp)) ^ 2)
	part5_nom <- -1 / 2 * t(theta_propose - parameter_default) %*% diag(rep(theta_w, n_para)) %*% (theta_propose - parameter_default)
	# calculate with theta_gibbs
	GP_obj_field <- GPpred(GP_sim, cbind(X_field, t(matrix(theta_gibbs[i - 1, ], n_para, n_field))))
	GP_obj_exp <- GPpred(GP_sim, cbind(matrix(X_exp, n_exp, 3), t(matrix(theta_gibbs[i - 1, ], n_para, n_exp))))
	delta_hat_field <- zeta_hat - exp(GP_obj_field$pred)[, 1]	
	V_hat_eta_exp <- exp(GP_obj_exp$pred)[1, 1] ^ 2 * GP_obj_exp$corr_mat[1, 1] * GP_obj_exp$var_scale[1,1]
	V_hat_delta_exp <- 4 - t(mat_cross) %*% solve(K_zeta + diag(nugget_val, n_field)) %*% mat_cross
	V_hat_eta_field <- diag(exp(GP_obj_field$pred)[, 1]) %*% (GP_obj_field$corr_mat * GP_obj_field$var_scale[1,1]) %*% diag(exp(GP_obj_field$pred)[, 1])	
	delta_hat_exp <- t(mat_cross) %*% solve(K_zeta + diag(nugget_val, n_field)) %*% delta_hat_field 
	part1_denom <- -1 / 2 * t(zeta_hat - delta_hat_field) %*% solve(K_zeta / lambda_delta_gibbs[i] + V_hat_eta_field + diag(nugget_val, n_field)) %*% (zeta_hat - delta_hat_field)
	part2_denom <- -0.5 * sum(log(abs(eigen(solve(K_zeta / lambda_delta_gibbs[i] + V_hat_eta_field + diag(nugget_val, n_field)) + W_zeta)$values))) - 0.5 * sum(log(abs(eigen(K_zeta / lambda_delta_gibbs[i] + V_hat_eta_field)$values)))
	part3_denom <- -n_exp / 2 * log(1 / lambda_E_gibbs[i] + V_hat_delta_exp / lambda_delta_gibbs[i] + V_hat_eta_exp)
	part4_denom <- -1 / 2 / (1 / lambda_E_gibbs[i] + V_hat_delta_exp / lambda_delta_gibbs[i] + V_hat_eta_exp) * sum((Y_exp - exp(GP_obj_exp$pred)[, 1] - rep(delta_hat_exp, n_exp)) ^ 2)
	part5_denom <- -1 / 2 * t(theta_gibbs[i - 1, ] - parameter_default) %*% diag(rep(theta_w, n_para)) %*% (theta_gibbs[i - 1, ] - parameter_default)
	# calculate acceptance ratio
	log_nom <- part1_nom + part2_nom + part3_nom + part4_nom + part5_nom
	log_denom <- part1_denom + part2_denom + part3_denom + part4_denom + part5_denom
	rho <- min(1, exp(log_nom - log_denom))
	temp <- runif(1, 0, 1)
	if (temp <= rho) {
		theta_gibbs[i, ] <- theta_propose
	}
	if (temp > rho) {
		theta_gibbs[i, ] <- theta_gibbs[i - 1, ]
	}		
}

# check traces
burn_in <- 200
par(mfrow = c(5, 1))
plot(lambda_delta_gibbs[(burn_in + 1) : n_iter], type = 'l')
plot(lambda_h_gibbs[(burn_in + 1) : n_iter], type = 'l')
plot(lambda_J_gibbs[(burn_in + 1) : n_iter], type = 'l')
plot(lambda_E_gibbs[(burn_in + 1) : n_iter], type = 'l')
plot(theta_gibbs[(burn_in + 1) : n_iter, 1], type = 'l')

# check posterior density
par(mfrow = c(4, 1))
plot(density(lambda_delta_gibbs[(burn_in + 1) : n_iter]))
plot(density(lambda_h_gibbs[(burn_in + 1) : n_iter]))
plot(density(lambda_J_gibbs[(burn_in + 1) : n_iter]))
plot(density(lambda_E_gibbs[(burn_in + 1) : n_iter]))
