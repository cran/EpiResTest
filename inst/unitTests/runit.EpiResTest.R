test.ILR <- function() {

	#require(ADGofTest)

	data(epi)
	
	#set.seed(1)
	
	n <- nrow(epi)
	
	kvalues_rk <- kvalues_wk <- matrix(NA,nrow=n, ncol=n) # kvalues computed from the correct spatial kernel and an arbitrary wrong spatial kernel respectively
	diag(kvalues_rk) <- diag(kvalues_wk) <- 1
	for (i in 1:n){
	for (j in 1:n){
		if (i<j) {
			d_ij <- sqrt((epi$coor_x[i]-epi$coor_x[j])**2 + (epi$coor_y[i]-epi$coor_y[j])**2)
			kvalues_rk[i,j] <- exp(-0.02*(d_ij))
			kvalues_wk[i,j] <- 1.0/(d_ij**1.0)
		}
		if (i>j) {
			kvalues_rk[i,j] <- kvalues_rk[j,i]
			kvalues_wk[i,j] <- kvalues_wk[j,i]
		}
	}
	}
	
	E <- epi$k[which(epi$te!=9e+100 & epi$te!=min(epi$te))]
	I <- epi$k[which(epi$ti!=9e+100)]
	tr <- epi$tr
	ti <- epi$ti
	te <- epi$te
	stb <- epi$stb
	alpha <- 0.002
	beta <- 8
	infected_source <- epi$infected_source

	checkTrue(ad.test(ILR(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb),punif,min=0,max=1)$p.value<=ad.test(ILR(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb))$p.value)
	checkTrue(ad.test(ILR(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb),punif,min=0,max=1)$p.value<=ad.test(ILR(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb))$p.value)
	checkTrue(ad.test(ILR(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb),punif,min=0,max=1)$p.value<=ad.test(ILR(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb))$p.value)
	checkTrue(ad.test(ILR(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb),punif,min=0,max=1)$p.value<=ad.test(ILR(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb))$p.value)
	checkTrue(ad.test(ILR(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb),punif,min=0,max=1)$p.value<=ad.test(ILR(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb))$p.value)
	checkTrue(ad.test(ILR(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb),punif,min=0,max=1)$p.value<=ad.test(ILR(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb))$p.value)

}



test.ILRpath <- function() {

	#require(ADGofTest)

	data(epi)
	
	#set.seed(1)
	
	n <- nrow(epi)
	
	kvalues_rk <- kvalues_wk <- matrix(NA,nrow=n, ncol=n) # kvalues computed from the correct spatial kernel and an arbitrary wrong spatial kernel respectively
	diag(kvalues_rk) <- diag(kvalues_wk) <- 1
	for (i in 1:n){
	for (j in 1:n){
		if (i<j) {
			d_ij <- sqrt((epi$coor_x[i]-epi$coor_x[j])**2 + (epi$coor_y[i]-epi$coor_y[j])**2)
			kvalues_rk[i,j] <- exp(-0.02*(d_ij))
			kvalues_wk[i,j] <- 1.0/(d_ij**1.0)
		}
		if (i>j) {
			kvalues_rk[i,j] <- kvalues_rk[j,i]
			kvalues_wk[i,j] <- kvalues_wk[j,i]
		}
	}
	}
	
	E <- epi$k[which(epi$te!=9e+100 & epi$te!=min(epi$te))]
	I <- epi$k[which(epi$ti!=9e+100)]
	tr <- epi$tr
	ti <- epi$ti
	te <- epi$te
	stb <- epi$stb
	alpha <- 0.002
	beta <- 8
	infected_source <- epi$infected_source

	checkTrue(ad.test(ILR_path(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source),punif,min=0,max=1)$p.value<=ad.test(ILR_path(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source))$p.value)
	checkTrue(ad.test(ILR_path(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source),punif,min=0,max=1)$p.value<=ad.test(ILR_path(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source))$p.value)
	checkTrue(ad.test(ILR_path(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source),punif,min=0,max=1)$p.value<=ad.test(ILR_path(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source))$p.value)
	checkTrue(ad.test(ILR_path(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source),punif,min=0,max=1)$p.value<=ad.test(ILR_path(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source))$p.value)
	checkTrue(ad.test(ILR_path(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source),punif,min=0,max=1)$p.value<=ad.test(ILR_path(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source))$p.value)
	checkTrue(ad.test(ILR_path(kvalues_wk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source),punif,min=0,max=1)$p.value<=ad.test(ILR_path(kvalues_rk, E, I, tr, ti, te, alpha, beta, n, stb, infected_source))$p.value)

}



test.LTR<- function() {

	#require(ADGofTest)

	data(epi)
	
	#set.seed(1)
	
	
	E <- epi$k[which(epi$te!=9e+100 & epi$te!=min(epi$te))]
	I <- epi$k[which(epi$ti!=9e+100)]
	EnI <- E[!(E%in%I)]
	ti <- epi$ti
	te <- epi$te

	tmax <- 60
	shape <- 10 #shape parameter for a Gamma waiting time in class E
	rate <- 2 # rate parameter for a Gamma waiting time in class E
	
	para_gamma <- c(shape, rate)
	para_exp <- shape/rate


	checkTrue(ad.test(LTR(E, I, EnI, ti, te, para_exp, tmax, pexp))$p.value<=ad.test(LTR(E, I, EnI, ti, te, para_gamma, tmax, pgamma))$p.value)
	checkTrue(ad.test(LTR(E, I, EnI, ti, te, para_exp, tmax, pexp))$p.value<=ad.test(LTR(E, I, EnI, ti, te, para_gamma, tmax, pgamma))$p.value)
	checkTrue(ad.test(LTR(E, I, EnI, ti, te, para_exp, tmax, pexp))$p.value<=ad.test(LTR(E, I, EnI, ti, te, para_gamma, tmax, pgamma))$p.value)
	checkTrue(ad.test(LTR(E, I, EnI, ti, te, para_exp, tmax, pexp))$p.value<=ad.test(LTR(E, I, EnI, ti, te, para_gamma, tmax, pgamma))$p.value)
	checkTrue(ad.test(LTR(E, I, EnI, ti, te, para_exp, tmax, pexp))$p.value<=ad.test(LTR(E, I, EnI, ti, te, para_gamma, tmax, pgamma))$p.value)
	checkTrue(ad.test(LTR(E, I, EnI, ti, te, para_exp, tmax, pexp))$p.value<=ad.test(LTR(E, I, EnI, ti, te, para_gamma, tmax, pgamma))$p.value)

}



