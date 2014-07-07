

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>


using namespace Rcpp;
using namespace std;
using namespace arma;


//' @title Infection Link Residual Test (Non-path-explicit)
//' @author Max S.Y. Lau <maxlauhk54@@gmail.com>
//'
//' @description Infection-link Residual Test (Non-path-explicit).
//'
//' @details The residual test specifically designed to test the goodness-of-fit of a spatial kernel.
//'Explicit sources of infections are not required and are imputed in the function.
//' @param kvalues  a nxn symmetric matrix where n is total number of subjects in the population. Entry at position (i,j) is the value of the spatial kernel  \eqn{K(d_{i,j})} where  \eqn{d_{i,j}} is the Euclidean distance between i and j. Diagonal entries are equal to 1.
//' @param E a vector contains the indices of infected cases excludingt the first case. Note that subjects are indexed from 0 (not 1).
//' @param I a vector contains the indices of infectious cases (being infectious or once infectious).
//' @param tr a vector contains the recovery times of all subjects recorded in the order of indexing (entry for an non-recovered should be specified to an arbitrary extreme value such as 9e+10).
//' @param ti a vector contains the times of becoming infectious of all subjects recorded in the order of indexing (entry for an non-recovered should be specified to an arbitrary extreme value such as 9e+10). For a SIR model, {\code{ti}} should set to be {\code{te}}
//' @param te a vector contains the infection times of all subjects recorded in the order of indexing (entry for an non-recovered should be specified to an arbitrary extreme value such as 9e+10).
//' @param alpha the background infection rate
//' @param beta the secondary infection rate
//' @param n the total number of subjects
//' @param stb a vector contains the susceptibility of subjects (entries may set to be 1 if to be ignored).
//' @return a vector contains the set of  imputed residual for infected in the order of indexing.
//' @examples
//'data(epi)
//'
//'set.seed(1)
//'
//'n <- nrow(epi)
//'
//'kvalues_rk <- kvalues_wk <- matrix(NA,nrow=n, ncol=n) 
//'\dontrun{kvalues_wk is to be computed from an arbitrary wrong spatial kernel below}
//'diag(kvalues_rk) <- diag(kvalues_wk) <- 1
//'for (i in 1:n){
//'for (j in 1:n){
//'	if (i<j) {
//'		d_ij <- sqrt((epi$coor_x[i]-epi$coor_x[j])^2 + (epi$coor_y[i]-epi$coor_y[j])^2)
//'		kvalues_rk[i,j] <- exp(-0.02*(d_ij))
//'		kvalues_wk[i,j] <- 1.0/(d_ij^1.0)
//'	}
//'	if (i>j) {
//'		kvalues_rk[i,j] <- kvalues_rk[j,i]
//'		kvalues_wk[i,j] <- kvalues_wk[j,i]
//'	}
//'}
//'}
//'
//'E <- epi$k[which(epi$te!=9e+100 & epi$te!=min(epi$te))]
//'I <- epi$k[which(epi$ti!=9e+100)]
//'tr <- epi$tr
//'ti <- epi$ti
//'te <- epi$te
//'stb <- epi$stb
//'alpha <- 0.002
//'beta <- 8
//'infected_source <- epi$infected_source
//'
//'par(mfrow=c(1,2),mar=c(4,4,4,4))
//'hist(ILR(kvalues_rk,E, I,tr,ti,te,alpha,beta,n,stb), main="Correct Model", xlab="ILR")
//'hist(ILR(kvalues_wk,E,I,epi$tr,ti,te,alpha,beta,n,stb), main="Wrong Model", xlab="ILR")
//' @references Lau, Max SY, George Streftaris, Glenn Marion, Gavin J. Gibson. "New model diagnostics for spatio-temporal systems in epidemiology and ecology." Journal of The Royal Society Interface 11.93 (2014): 20131093.
//' @export
// [[Rcpp::export]]
NumericVector ILR(NumericMatrix kvalues, IntegerVector E, IntegerVector I,  NumericVector tr, NumericVector ti,  NumericVector te, double alpha, double beta,int n,  NumericVector stb){

//ofstream myfile_out; 

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c,iter*1000); // set a seed


std::vector <double> u_imputed(E.size()); // the imputed residuals for all infected, exclude index

for (int i=0; i<=((int)E.size()-1);i++){
//for (int i=0; i<=(7);i++){

	std::vector<double> ic_link; // instantaneous infective challenge of the infection links between infectious/primary infection and the new infected/
	std::vector<double> pr_link; // probabilities of the infection links between infectious/primary infection and the new infected
	
	std::vector<double> ic; // instantaneous infective challenge of the infection links between infectious/primary infection and the new infected/susceptibles
	std::vector<double> pr; // probabilities of the infection links between infectious/primary infection and the new infected/susceptibles
	std::vector<double> cdf; // cumulative probabilities of the infection links between infectious/primary infection and the new infected/susceptibles
	
	double total_ic_link; // the total instantaneous infective challenge to the new infected
	double total_ic; // the total instantaneous infective challenge to the new infected/susceptibles
	
	int link_imputed; // link imputed
	double ic_imputed; // ic of the imputed link
	int rank_imputed; // the rank of p_imputed among all links; 0 means first element

	ic_link.reserve(I.size()+2); // 2 more positions reserved: one for primary infection and one for faciliating computation

	
	ic_link.push_back(0.0); // 1st element in ic
	ic_link.push_back(alpha); // 2nd..

	//-----

	total_ic_link = 0.0; 
	for (int j=0; j<=((int)I.size()-1);j++){
	if ((ti[I[j]]<=te[E[i]]) & (tr[I[j]]>=te[E[i]]) &(I[j]!=E[i])){
	ic_link.push_back(beta*stb[E[i]]*kvalues(E[i],I[j]));
	total_ic_link = total_ic_link + stb[E[i]]*kvalues(E[i],I[j]);
	}
	}
	total_ic_link = alpha + beta*total_ic_link;

	//-----
	pr_link.resize(ic_link.size());

	pr_link.at(0) = 0.0;
	pr_link.at(1) = alpha/total_ic_link; //pr it was infected by primary infection
	
	for (int k=2; k<=(int)(pr_link.size()-1);k++){
	pr_link.at(k) = ic_link.at(k)/total_ic_link;
	}

// 	double *P=&pr_link.at(0);
// 	gsl_ran_discrete_t * g = gsl_ran_discrete_preproc ((int)pr_link.size(),P);
// 	link_imputed= gsl_ran_discrete (r_c, g); // the link imputed
// 	gsl_ran_discrete_free (g);

	IntegerVector x= Range((int) 0, (int)pr_link.size()-1);
	NumericVector prob(pr_link.begin(),pr_link.end());
	link_imputed = Rcpp::RcppArmadillo::sample(x,1, FALSE, prob)[0];
	ic_imputed = ic_link.at(link_imputed);
	//-----

// 	int source = path[E[i]]; // do not need to impute link given transmission path
// 
// 	switch(source==9999){
// 		case 0:{
// 			ic_imputed = beta*stb[E[i]]*kvalues(E[i],source);
// 		break;
// 		}
// 		case 1:{
// 			ic_imputed = alpha;
// 		break;
// 		}
// 	}
	//-----
	
	ic.reserve(ic_link.size()+(I.size()+1)*n);

	total_ic = 0.0; // sum of ic of links from the bunch of infectious to ALL susceptible subject before te.at(E.at(i))

	for (int iu=0; iu<=((int)n-1);iu++){

		double total_ic_iu =0.0; // sum if ic of links from the bunch of infectious to the particular susceptible subject

		switch(te[iu]>te[E[i]]) { // return 1 if the subject is susceptible before te.at(E.at(i))note: this excludes the links connected to xi_E_minus.at(i), where the infection actually happened
		case 1:{
		//ic.push_back(0.0); 
		ic.push_back(alpha); 
	
		for (int j=0; j<=((int)I.size()-1);j++){
		if ((ti[I[j]]<=te[E[i]]) & (tr[I[j]]>=te[E[i]])&(I[j]!=E[i])){ // this ensures we are using same bunch of infectious sources
		ic.push_back(beta*stb[iu]*kvalues(iu,I[j]));
		total_ic_iu = total_ic_iu + stb[iu]*kvalues(iu,I[j]);
		}
		}
		total_ic_iu = alpha + beta*total_ic_iu; 
	
		break;
		}
	
		case 0:{
		break;
		}
		}

		total_ic = total_ic + total_ic_iu;

	}
	
	total_ic = total_ic + total_ic_link; // add the ic from the infectious to infected links

	ic.insert(ic.begin(),ic_link.begin(),ic_link.end()); // combine ic_link and ic

	pr.resize(ic.size());

	for (int k=0; k<=((int)pr.size()-1);k++){
	pr.at(k) = ic.at(k)/total_ic;
	}

	std::sort(pr.begin(),pr.end());

	double p_imputed = ic_imputed / total_ic;
	
	rank_imputed = distance(  pr.begin(), find(pr.begin(),pr.end(),p_imputed) );


	//---------//

//  Rcpp::Rcout << rank_imputed <<"," << pr.size() << "," << p_imputed << endl;

	int num_dup=1; // number of prs same to p_imputed (including p_imputed itself)
 	std::vector<double> dup_test = pr;
 	dup_test.erase(dup_test.begin()+rank_imputed);

	switch(find(dup_test.begin(),dup_test.end(),p_imputed)==dup_test.end()) { // return 1 when no duplicate 
	case 0:{ // with dup
	num_dup = count(pr.begin(),pr.end(),p_imputed);
	break;
	}
	case 1:{
	num_dup=1;
	break;
	}
	}

	//----

	cdf.resize(pr.size());

	cdf.at(0) = 0.0;
	for (int k=1; k<=(int)(cdf.size()-1);k++){
	cdf.at(k) = pr.at(k) + cdf.at(k-1);
	}

	double u_i=0.0; // the u_i to be imputed
	switch(rank_imputed==0){
	case 1:{
	u_i = 0.0;
	break;
	}
	case 0:{
	u_i = R::runif(cdf.at(rank_imputed-1),  cdf.at(rank_imputed-1)+num_dup*p_imputed); // a uniform rv drawn between cdf.at(rank_imputed-1), and cdf.at(rank_imputed-1)+num_dup*p_imputed

	break;
	}
	}

	u_imputed.at(i) =u_i;
	//------------------------


	//pr_link.clear();

	ic_link.clear();

	pr.clear();
	ic.clear();
	cdf.clear();

}

// 	myfile_out.open((string(path4)+string("residual_kernel.txt")).c_str(),ios::app);
// 	myfile_out << endl;
// 	for (int i=0; i<=((int)E.size()-1);i++){
// 	if (i!= (int) E.size()-1) myfile_out << u_imputed.at(i) << ",";
// 	if (i== (int) E.size()-1) myfile_out << u_imputed.at(i);
// 	}
// 	myfile_out.close();

// 	myfile_out.open((string(path4)+string("size_residual.txt")).c_str(),ios::app);
// 	myfile_out <<  u_imputed.size() << endl;
// 	myfile_out.close();

//gsl_rng_free(r_c);


	NumericVector ILR (u_imputed.begin(), u_imputed.end());
	return ILR; 

}




//' @title Infection Link Residual Test (Path Explicit)
//' @author Max S.Y. Lau <maxlauhk54@@gmail.com>
//'
//' @description Infection-link Residual Test (Path Explicit).
//'
//' @details The residual test specifically designed to test the goodness-of-fit of a spatial kernel.
//'Explicit sources of infections are required.
//' @param kvalues  a nxn symmetric matrix where n is total number of subjects in the population. Entry at position (i,j) is the value of the spatial kernel \eqn{K(d_{i,j})} where  \eqn{d_{i,j}} is the Euclidean distance between i and j. Diagonal entries are equal to 1.
//'@param E a vector contains the indices of infected cases excludingt the first case. Note that subjects are indexed from 0 (not 1).
//' @param I a vector contains the indices of infectious cases (being infectious or once infectious).
//' @param tr a vector contains the recovery times of all subjects recorded in the order of indexing (entry for an non-recovered should be specified to an arbitrary extreme value such as 9e+10).
//' @param ti a vector contains the times of becoming infectious of all subjects recorded in the order of indexing (entry for an non-recovered should be specified to an arbitrary extreme value such as 9e+10). For a SIR model, {\code{ti}} set to be {\code{te}}
//' @param te a vector contains the infection times of all subjects recorded in the order of indexing (entry for an non-recovered should be specified to an arbitrary extreme value such as 9e+10).
//' @param alpha the background infection rate
//' @param beta the secondary infection rate
//' @param n the total number of subjects
//' @param stb a vector contains the susceptibility of subjects (entries may set to be 1 if to be ignored).
//' @param path a vector contains the source of infection of the cases (entries for background infections and non-infected should be specified to be 9999 and -99 respectively).
//' @return a vector contains the set of  imputed residual for infected in the order of indexing.
//' @examples
//'data(epi)
//'
//'set.seed(1)
//'
//'n <- nrow(epi)
//'
//'kvalues_rk <- kvalues_wk <- matrix(NA,nrow=n, ncol=n) 
//'\dontrun{kvalues_wk is to be computed from an arbitrary wrong spatial kernel below}
//'diag(kvalues_rk) <- diag(kvalues_wk) <- 1
//'for (i in 1:n){
//'for (j in 1:n){
//'	if (i<j) {
//'		d_ij <- sqrt((epi$coor_x[i]-epi$coor_x[j])^2 + (epi$coor_y[i]-epi$coor_y[j])^2)
//'		kvalues_rk[i,j] <- exp(-0.02*(d_ij))
//'		kvalues_wk[i,j] <- 1.0/(d_ij^1.0)
//'	}
//'	if (i>j) {
//'		kvalues_rk[i,j] <- kvalues_rk[j,i]
//'		kvalues_wk[i,j] <- kvalues_wk[j,i]
//'	}
//'}
//'}
//'
//'E <- epi$k[which(epi$te!=9e+100 & epi$te!=min(epi$te))]
//'I <- epi$k[which(epi$ti!=9e+100)]
//'tr <- epi$tr
//'ti <- epi$ti
//'te <- epi$te
//'stb <- epi$stb
//'alpha <- 0.002
//'beta <- 8
//'infected_source <- epi$infected_source
//'
//'par(mfrow=c(1,2),mar=c(4,4,4,4))
//'hist(ILR_path(kvalues_rk,E,I,tr,ti,te,alpha,beta,n,stb,infected_source),
//'main="Correct Model",xlab="ILR")
//'hist(ILR_path(kvalues_wk,E,I,epi$tr,ti,te,alpha,beta,n,stb,infected_source), 
//'main="Wrong Model",xlab="ILR")
//' @references Lau, Max SY, George Streftaris, Glenn Marion, Gavin J. Gibson. "New model diagnostics for spatio-temporal systems in epidemiology and ecology." Journal of The Royal Society Interface 11.93 (2014): 20131093.
//' @export
// [[Rcpp::export]]
NumericVector ILR_path(NumericMatrix kvalues, IntegerVector E, IntegerVector I,  NumericVector tr, NumericVector ti,  NumericVector te, double alpha, double beta,int n,  NumericVector stb,  IntegerVector path){

std::vector <double> u_imputed(E.size()); // the imputed residuals for all infected, exclude index

for (int i=0; i<=((int)E.size()-1);i++){
//for (int i=0; i<=(7);i++){

	std::vector<double> ic_link; // instantaneous infective challenge of the infection links between infectious/primary infection and the new infected/
	//std::vector<double> pr_link; // probabilities of the infection links between infectious/primary infection and the new infected
	
	std::vector<double> ic; // instantaneous infective challenge of the infection links between infectious/primary infection and the new infected/susceptibles
	std::vector<double> pr; // probabilities of the infection links between infectious/primary infection and the new infected/susceptibles
	std::vector<double> cdf; // cumulative probabilities of the infection links between infectious/primary infection and the new infected/susceptibles
	
	double total_ic_link; // the total instantaneous infective challenge to the new infected
	double total_ic; // the total instantaneous infective challenge to the new infected/susceptibles
	
	//int link_imputed; // link imputed
	double ic_imputed; // ic of the imputed link
	int rank_imputed; // the rank of p_imputed among all links; 0 means first element

	ic_link.reserve(I.size()+2); // 2 more positions reserved: one for primary infection and one for faciliating computation

	
	ic_link.push_back(0.0); // 1st element in ic
	ic_link.push_back(alpha); // 2nd..

	//-----

	total_ic_link = 0.0; 
	for (int j=0; j<=((int)I.size()-1);j++){
	if ((ti[I[j]]<=te[E[i]]) & (tr[I[j]]>=te[E[i]])&(I[j]!=E[i])){
	ic_link.push_back(beta*stb[E[i]]*kvalues(E[i],I[j]));
	total_ic_link = total_ic_link + stb[E[i]]*kvalues(E[i],I[j]);
	}
	}
	total_ic_link = alpha + beta*total_ic_link;

	//-----


	int source = path[E[i]]; // do not need to impute link given transmission path

	switch(source==9999){
		case 0:{
			ic_imputed = beta*stb[E[i]]*kvalues(E[i],source);
		break;
		}
		case 1:{
			ic_imputed = alpha;
		break;
		}
	}
	//-----
	
	ic.reserve(ic_link.size()+(I.size()+1)*n);

	total_ic = 0.0; // sum of ic of links from the bunch of infectious to ALL susceptible subject before te.at(E.at(i))

	for (int iu=0; iu<=((int)n-1);iu++){

		double total_ic_iu =0.0; // sum if ic of links from the bunch of infectious to the particular susceptible subject

		switch(te[iu]>te[E[i]]) { // return 1 if the subject is susceptible before te.at(E.at(i))note: this excludes the links connected to xi_E_minus.at(i), where the infection actually happened
		case 1:{
		//ic.push_back(0.0); 
		ic.push_back(alpha); 
	
		for (int j=0; j<=((int)I.size()-1);j++){
		if ((ti[I[j]]<=te[E[i]]) & (tr[I[j]]>=te[E[i]])&(I[j]!=E[i])){ // this ensures we are using same bunch of infectious sources
		ic.push_back(beta*stb[iu]*kvalues(iu,I[j]));
		total_ic_iu = total_ic_iu + stb[iu]*kvalues(iu,I[j]);
		}
		}
		total_ic_iu = alpha + beta*total_ic_iu; 
	
		break;
		}
	
		case 0:{
		break;
		}
		}

		total_ic = total_ic + total_ic_iu;

	}
	
	total_ic = total_ic + total_ic_link; // add the ic from the infectious to infected links

	ic.insert(ic.begin(),ic_link.begin(),ic_link.end()); // combine ic_link and ic

	pr.resize(ic.size());

	for (int k=0; k<=((int)pr.size()-1);k++){
	pr.at(k) = ic.at(k)/total_ic;
	}

	std::sort(pr.begin(),pr.end());

	double p_imputed = ic_imputed / total_ic;
	
	rank_imputed = distance(  pr.begin(), find(pr.begin(),pr.end(),p_imputed) );
	//---------//

	int num_dup=1; // number of prs same to p_imputed (including p_imputed itself)
 	std::vector<double> dup_test = pr;
 	dup_test.erase(dup_test.begin()+rank_imputed);

	switch(find(dup_test.begin(),dup_test.end(),p_imputed)==dup_test.end()) { // return 1 when no duplicate 
	case 0:{ // with dup
	num_dup = count(pr.begin(),pr.end(),p_imputed);
	break;
	}
	case 1:{
	num_dup=1;
	break;
	}
	}

	cdf.resize(pr.size());

	cdf.at(0) = 0.0;
	for (int k=1; k<=(int)(cdf.size()-1);k++){
	cdf.at(k) = pr.at(k) + cdf.at(k-1);
	}

	double u_i=0.0; // the u_i to be imputed
	switch(rank_imputed==0){
	case 1:{
	u_i = 0.0;
	break;
	}
	case 0:{
	u_i = R::runif(cdf.at(rank_imputed-1),  cdf.at(rank_imputed-1)+num_dup*p_imputed); // a uniform rv drawn between cdf.at(rank_imputed-1), and cdf.at(rank_imputed-1)+num_dup*p_imputed
	break;
	}
	}

	u_imputed.at(i) =u_i;

	ic_link.clear();

	pr.clear();
	ic.clear();
	cdf.clear();

}

	NumericVector ILR (u_imputed.begin(), u_imputed.end());
	return ILR; 

}

//' @title Latent Time Residual
//' @author Max S.Y. Lau <maxlauhk54@@gmail.com>
//'
//' @description  Latent Time Residual.
//'
//' @details The residual test specifically designed to test the goodness-of-fit of the latent period (i.e., waiting time in class E) and the infectious period (i.e., waiting time in class I).
//'@param E a vector contains the indices of infected cases excludingt the first case. Note that subjects are indexed from 0 (not 1).
//' @param I a vector contains the indices of infectious cases (being infectious or once infectious).
//' @param EnI a vector contains the indices of  infected cases which have not yet been infectious (being infectious or once infectious).
//' @param ti a vector contains the times of becoming infectious of all subjects recorded in the order of indexing (entry for an non-recovered should be specified to an arbitrary extreme value such as 9e+10). For testing waiting time in class I, {\code{ti}} should be set to be the time of recovery  {\code{tr}} 
//' @param te a vector contains the infection times of all subjects recorded in the order of indexing (entry for an non-recovered should be specified to an arbitrary extreme value such as 9e+10). . For testing waiting time in class I, {\code{te}} should be set to be the time of recovery  {\code{ti}} 
//' @param para either a one-element numeric vector or two-parameter vector contains the parameter arguments of the waiting time distribution to be specified by the user (see argument {\code{wt}}.
//' @param tmax the maximum duration of the observational period
//' @param wt the cumulative distribution function of the waiting time distribution specified by the user; note that the number of parameter arguments has to match with the length of {\code{para}}.
//' @return a vector contains the set of  imputed residual for infected in the order of indexing.
//' @examples
//'data(epi)
//'
//'set.seed(1)
//'
//'E <- epi$k[which(epi$te!=9e+100 & epi$te!=min(epi$te))]
//'I <- epi$k[which(epi$ti!=9e+100)]
//'EnI <- E[!(E%in%I)]
//'ti <- epi$ti
//'te <- epi$te
//'tmax <- 60
//'shape <- 10 #shape parameter for a Gamma waiting time in class E
//'rate <- 2 # rate parameter for a Gamma waiting time in class E
//'
//'para_gamma <- c(shape, rate)
//'para_exp <- shape/rate # an Exponential distribution to match the mean of the Gamma
//'
//'par(mfrow=c(1,2),mar=c(4,4,4,4))
//'hist(LTR(E, I, EnI, ti, te, para_gamma, tmax, pgamma), main="Correct Model", xlab="LTR")
//'hist(LTR(E, I, EnI, ti, te, para_exp, tmax, pexp), main="Wrong Model", xlab="LTR")
//' @references Lau, Max SY, George Streftaris, Glenn Marion, Gavin J. Gibson. "New model diagnostics for spatio-temporal systems in epidemiology and ecology." Journal of The Royal Society Interface 11.93 (2014): 20131093.
//' @export
// [[Rcpp::export]]
NumericVector LTR ( IntegerVector E, IntegerVector I,  IntegerVector EnI, NumericVector ti, NumericVector te, NumericVector para, double tmax, Function wt){

vector <double> u_imputed; // the imputed residuals for all infected
u_imputed.reserve(E.size());

double cdf=0.0;

for (int i=0; i<=((int)I.size()-1);i++){

	if ((int)para.size()==1) cdf = as<double>(wt( ti[I[i]]-te[I[i]], para[0] ));
	if ((int)para.size()==2) cdf = as<double>(wt( ti[I[i]]-te[I[i]], para[0], para[1] ));

	double u_i = exp(log(1.0-cdf));
	u_imputed.push_back(u_i);
}

for (int i=0; i<=((int)EnI.size()-1);i++){
	if ((int)para.size()==1)  cdf =  as<double>(wt( tmax - te[EnI[i]], para[0]));
	if ((int)para.size()==2)  cdf =  as<double>(wt( tmax - te[EnI[i]], para[0], para[1]));

	double u_i = R::runif(0.0, 1.0-cdf );
	u_imputed.push_back(u_i);
}

	NumericVector LTR (u_imputed.begin(), u_imputed.end());
	return LTR; 

}
