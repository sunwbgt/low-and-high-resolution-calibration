Injury_Risk_ChestD <- function(ChestD, age) {
	return(1 / (1 + exp(12.597 - 0.05861 * age - 1.568 * ((ChestD * 1000) ^ 0.4612))))
}

Injury_Risk_HIC <- function(HIC) {
	return(pnorm(log(HIC), 7.45231, 0.73998))
}

Injury_Risk_FF <- function(FF) {
	return(1 / (1 + exp(5.7949 - 0.5196 * FF / 1)))
}
