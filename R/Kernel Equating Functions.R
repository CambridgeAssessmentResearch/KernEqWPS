#KERNEL EQUATING FUNCTIONS

#' Function to calculate mu and sigma from a score distribution
#'
#' @param scores Vector of possible scores
#' @param count Number of people achieving each score (weighted numbers, probabilities or proportions also acceptable). Must be same length as scores.
#'
#' @examples
#' x=pmax(0,pmin(10,round(rnorm(200,5,2))))
#' counts=200*(ecdf(x)(0:10)-ecdf(x)(-1:9))#used rather than "table" so all scores included even if zero count
#' DistToMuSig(0:10,counts)
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{mu}{The mean.}
#'   \item{sigma}{The standard deviation.}
#' }
#'
#' @keywords KernEqWPS
#' @export
DistToMuSig=function(scores,count){
	if(length(scores)!=length(count)){print("vector of scores different length from vector of frequencies")
						return(NULL)}
	n=sum(count)
	mu=sum(scores*count)/n

	#check if any non-integers in counts
	anynonint=max(count-round(count))>0.00001
	if(anynonint){var=(1/n)*sum(count*(scores-mu)^2)}#if non-integers in the weights ignore the (n-1) thing as it's now just a distribution
	if(anynonint==FALSE){var=(1/(n-1))*sum(count*(scores-mu)^2)}#if genuine counts then compute sd with an n-1 correction to reduce bias
	sigma=sqrt(var)
	return(list(mu=mu,sigma=sigma))}

#' Continuized density function from a score distribution for a given bandwidth
#'
#' Note: mu and sigma are parameters in here (rather than being calculated from distribution) 
#' to avoid repetitive recalculation when applying this function within the routine to optimise bandwidths
#'
#' @param x Score(s) for which we wish to calculate the density
#' @param fX Number of people achieving each score (weighted numbers, probabilities or proportions also acceptable). Must be same length as scores.
#' @param scores Vector of possible scores.
#' @param mu Mean score in distribution.
#' @param sigma Standard deviation score in distribution.
#' @param h Bandwidth used for continuization.
#'
#' @keywords KernEqWPS
#' @export
fhx=function(x,fX,scores,mu,sigma,h){
	#SUBROUTINE 
	Rjx=function(x,xj,mu,sigma,h){
	ax=sqrt((sigma^2)/(sigma^2+h^2))
	(x-ax*xj-(1-ax)*mu)/(ax*h)}
	
	ax=sqrt((sigma^2)/(sigma^2+h^2))
	f1=function(x,fX,scores,mu,sigma,h){
		sum(fX*dnorm(Rjx(x,scores,mu,sigma,h)))/(ax*h)}
	sapply(x,f1,fX=fX,scores=scores,mu=mu,sigma=sigma,h=h)}

#' Continuized cumulative density function from a score distribution for a given bandwidth
#'
#' Note: mu and sigma are parameters in here (rather than being calculated from distribution) 
#' to avoid repitative recalculation when applying this function within the routine to optimise bandwidths
#'
#' @param x Score(s) for which we wish to calculate the density
#' @param fX Number of people achieving each score (weighted numbers, probabilities or proportions also acceptable). Must be same length as scores.
#' @param scores Vector of possible scores.
#' @param mu Mean score in distribution.
#' @param sigma Standard deviation score in distribution.
#' @param h Bandwidth used for continuization.
#'
#' @keywords KernEqWPS
#' @export
Fhx=function(x,fX,scores,mu,sigma,h){
	#SUBROUTINE 
	Rjx=function(x,xj,mu,sigma,h){
	ax=sqrt((sigma^2)/(sigma^2+h^2))
	(x-ax*xj-(1-ax)*mu)/(ax*h)}
	
	f1=function(x,fX,scores,mu,sigma,h){
		sum(fX*pnorm(Rjx(x,scores,mu,sigma,h)))}
	sapply(x,f1,fX=fX,scores=scores,mu=mu,sigma=sigma,h=h)}

#' Calulate PEN1 criterion for match of continuized and original score distributions
#' 
#' This function calculates the penalty term for a given bandwidth against a given distribution.
#' The definition of PEN1 can be found in Liang and von Davier (2014). The values of PEN1 are not of much interest in their own right.
#' Finding the bandwidth that minimises this penalty (using \emph{optPEN1}) can potentially be useful as identifying a plausible bandwidth to use in subsequent kernel equating.
#'
#' @param h Bandwith being considered.
#' @param scores Vector of possible scores.
#' @param fX Number of people achieving each score (weighted numbers, probabilities or proportions also acceptable). Must be same length as scores.
#' @return The function returns the value of PEN1 for the given bandwidth.
#' @references
#' Liang, T., & von Davier, A. A. (2014). Cross-validation: An alternative bandwidth-selection method in kernel equating. 
#' \emph{Applied Psychological Measurement, 38}(4), 281-295.
#'
#' @examples
#' x1=rowSums(mathsdata)#make some scores
#' tab=tabulate(x1+1)#count how many people got each integer score
#' cbind(0:max(x1),tab)#look at the table
#' PEN1(0.5,0:max(x1),tab)#calculate PEN1 for a given bandwidth
#' optPEN1(0:max(x1),tab)#optimise PEN1
#' @keywords KernEqWPS
#' @export
PEN1=function(h,scores,fX){
	msX=DistToMuSig(scores,fX)
	fhx1=fhx(scores,fX,scores,msX$mu,msX$sigma,h)
	sum((fX-fhx1)^2)}

#' Find bandwidth (between 0.05 and 30) that minimses \emph{PEN1}.
#'
#' This function minimises a defined penalty term for a given bandwidth against a given distribution.
#' The definition of PEN1 can be found in Liang and von Davier (2014). The values of PEN1 are not of much interest in their own right.
#' Finding the bandwidth that minimises this penalty (using \emph{optPEN1}) can potentially be useful as identifying a plausible bandwidth to use in subsequent kernel equating.
#'
#' Note that although this is a typical method of determining the bandwidth it is widely considered to undersmooth. For this reason it is more commonly used when fX has been pre-smoothed with a log-linear model or when the sample size is very large. 
#'
#' @param scores Vector of possible scores.
#' @param fX Number of people achieving each score (weighted numbers, probabilities or proportions also acceptable). Must be same length as scores.
#' @return The function returns a single value (labelled "minimum") representing bandwidth (between 0.05 and 30) that minimises PEN1 and the value on PEN1 (objective) at that value.
#' @references
#' Liang, T., & von Davier, A. A. (2014). Cross-validation: An alternative bandwidth-selection method in kernel equating. 
#' \emph{Applied Psychological Measurement, 38}(4), 281-295.
#'
#' @examples
#' x1=rowSums(mathsdata)#make some scores
#' tab=tabulate(x1+1)#count how many people got each integer score
#' cbind(0:max(x1),tab)#look at the table
#' PEN1(0.5,0:max(x1),tab)#calculate PEN1 for a given bandwidth
#' optPEN1(0:max(x1),tab)#optimise PEN1
#' @keywords KernEqWPS
#' @export
optPEN1=function(scores,fX){optimise(PEN1,c(0.05,30),scores=scores,fX=fX)}

#' One iteration of method to find optimal bandwidth (between 0.05 and 30) that minimses cross-validation error between two half-samples.
#'
#' Based on paper by Liang and von Davier (2014). It works by splitting the data set in half, fitting a continuized distribution with a gaussian kernel 
#' to one half of the data, and then evaluating the log-likelihood of this distribution against the other half of the data.
#'
#' @param scores Vector of possible scores.
#' @param fX Number of people achieving each score (should be whole numbers). Must be same length as scores.
#' @keywords KernEqWPS
#' @examples
#' tot1=rowSums(mathsdata)
#' tab1=table(tot1)
#' FindBestBandwidth1Iter(as.numeric(names(tab1)),as.vector(tab1))
#' @return The function returns a single value representing the optimal bandwidth.
#' @references
#' Liang, T., & von Davier, A. A. (2014). Cross-validation: An alternative bandwidth-selection method in kernel equating. 
#' \emph{Applied Psychological Measurement, 38}(4), 281-295.
#' @export
FindBestBandwidth1Iter=function(scores,fX){
	fXa=rbinom(length(fX),round(fX),0.5)
	fXb=fX-fXa
	msXa=DistToMuSig(scores,fXa)
	#function to find optimal bandwidth based on this split
	BandCrit=function(htry,scores,fXa,fXb,msXa){
		#estimate kernel based on first half-sample
		fhx1=fhx(scores,fX=fXa,scores=scores,mu=msXa$mu,sigma=msXa$sigma,h=htry)
		#renormalise as may not be of exactly equal sample size
		fhx1=fhx1*sum(fXb)/sum(fXa)
		crit=sum(fXb*log(fhx1)-fhx1)
		return(crit)}
	optimise(BandCrit,c(0.05,30),scores=scores,fXa=fXa,fXb=fXb,msXa=msXa,maximum=TRUE)$maximum
	}

#' Method to find optimal bandwidth (between 0.05 and 30) that minimses cross-validation error between two half-samples.
#'
#' Takes the median optimal bandwidth across multiple half-samples. Based on paper by Liang and von Davier (2014).
#'
#' @param scores Vector of possible scores.
#' @param fX Number of people achieving each score (weighted numbers, probabilities or proportions also acceptable). Must be same length as scores.
#' @param niters The number of iterations of applying the algorithm over which to take the median bandwidth (default 100).
#' @param randseed Random seed for splits (default 1234567).
#' @return The function returns a single value representing the optimal bandwidth.
#' @examples
#' x=sample(c(0,0,1,2,2,2,3,4,4,4,4,5,5,6,6,7,8,9,10),200,replace=TRUE)#create x from an unusual true distribution
#' counts=200*(ecdf(x)(0:10)-ecdf(x)(-1:9))#used rather than "table" so all scores included even if zero count
#' FindBestBandwidth(0:10,counts)
#' FindBestBandwidth(0:10,counts,randseed=987654321)
#' 
#' #EXAMPLE OF FINDING BEST BANDWIDTH TO EQUATE FIRST AND SECOND HALVES OF MATHS TEST
#' x1=rowSums(mathsdata[,1:62])
#' y1=rowSums(mathsdata[,63:124])
#' summary(x1)
#' summary(y1)
#' tabx1=table(x1)
#' taby1=table(y1)
#' 
#' BestHx=FindBestBandwidth(as.numeric(names(tabx1)),as.vector(tabx1))
#' BestHx
#' 
#' BestHy=FindBestBandwidth(as.numeric(names(taby1)),as.vector(taby1))
#' BestHy
#' 
#' eq1=KernelEquateFromScoresEG(x1,y1,hX=BestHx,hY=BestHy)
#' plot(sort(unique(x1)),eq1$yx,type='l')
#' #compare to default Andersson and von Davier bandwidth approach
#' eq2=KernelEquateFromScoresEG(x1,y1)
#' lines(sort(unique(x1)),eq2$yx,lty=2,col="purple")
#' #very similar in this case
#' @references
#' Liang, T., & von Davier, A. A. (2014). Cross-validation: An alternative bandwidth-selection method in kernel equating. 
#' \emph{Applied Psychological Measurement, 38}(4), 281-295.
#' @export
#' @keywords KernEqWPS
#' @export
FindBestBandwidth=function(scores,fX,niters=100,randseed=1234567){
	#store the curren random seed
	tmp <- .Random.seed
	#use a set seed so I can replicate results if I need to
	set.seed(randseed)
	BestBandwidth=median(sapply(1:niters,function(x,scores,fX){FindBestBandwidth1Iter(scores,fX)},scores=scores,fX=fX))
	#now put the (global hence "<<-" rather than "<-") random number seed back how it was before so other processes dependent on generating random numbers can continue
	.Random.seed <<- tmp 
	#return result
	BestBandwidth
	}

#' Apply Kernel Equating to Equivalent Groups.
#'
#' Assumes that X and Y are taken by equivalent groups (or by the same group). Requires bandwidths for continuization to be supplied.
#'
#' @param scoresX Vector of possible scores on form X.
#' @param fX Number of people achieving each score on form X (weighted numbers, probabilities or proportions also acceptable). Must be same length as scoresX.
#' @param scoresY Vector of possible scores on form Y.
#' @param fY Number of people achieving each score on form Y (weighted numbers, probabilities or proportions also acceptable). Must be same length as scoresY.
#' @param hX Bandwidth for form X.
#' @param hY Bandwidth for form Y.
#' @keywords KernEqWPS
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{yx}{A vector representing the scores derived from applying yxFunc to scoresX.}
#' }
#'
#' @export
KernelEquateFromDists=function(scoresX,fX,scoresY,fY,hX,hY){

	if(length(scoresX)!=length(fX)){print("vector of X scores different length from vector of frequencies")
						return(NULL)}
	if(length(scoresY)!=length(fY)){print("vector of Y scores different length from vector of frequencies")
						return(NULL)}
	
	msX=DistToMuSig(scoresX,fX)
	msY=DistToMuSig(scoresY,fY)

	#NOW NEED TO FIND THE SCORE POINTS ON Y RELATIVE TO EACH OF THESE CUMULATIVE POINTS
	f2=function(y,fY,scores,muY,sigmaY,hY,cump){
		(Fhx(y,fX=fY,scores=scores,mu=muY,sigma=sigmaY,h=hY)/sum(fY))-cump}
	f4=function(cump,fY,scores,muY,sigmaY,hY){uniroot(f2,c(-1,max(scores)+1),fY=fY,scores=scores,muY=msY$mu,sigmaY=msY$sigma,hY=hY,cump,extendInt="yes")$root}

	yxFunc=function(scores1){sapply(Fhx(scores1,fX=fX,scores=scoresX,mu=msX$mu,sigma=msX$sigma,h=hX)/sum(fX),f4,fY=fY,scores=scoresY,muY=msY$mu,sigmaY=msY$sigma,hY=hY)}
	yx=yxFunc(scoresX)
	return(list(yxFunc=yxFunc,yx=yx))}

#' Apply kernel equating based on two vectors of scores on two forms from equivalent candidates.
#'
#' @param xscores A vector of scores from individual candidates on form X.
#' @param yscores A vector of scores from individual candidates on form Y.
#' @param hX Bandwidth for form X. By default a plug-in estimator from Andersson and von Davier (2014) is used.
#' @param hY Bandwidth for form Y. By default a plug-in estimator from Andersson and von Davier (2014) is used.
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{yx}{A vector representing the scores derived from applying yxFunc to a sorted list of all the unqiue form X scores that actually occur in the data.}
#' }
#'
#' @examples
#' #simple example of equating a roughly normally distributed sets of scores to one with a roughly logistic distribution
#' x=pmax(0,pmin(100,round(rnorm(300,50,15))))
#' y=pmax(0,pmin(100,round(rlogis(300,60,16))))
#' hist(x)
#' hist(y)
#' KE1=KernelEquateFromScoresEG(x,y)
#' plot(0:100,KE1$yxFunc(0:100),type='l')
#'
#' @keywords KernEqWPS
#' @references
#' Andersson, B., & von Davier, A. A. (2014). Improving the bandwidth selection in kernel equating. 
#' \emph{Journal of Educational Measurement, 51}(3), 223-238.
#' @export
KernelEquateFromScoresEG=function(xscores,yscores,hX=NA,hY=NA){

	#could include restriction that xscores and yscores must be the same length (currently not included).
	#if(length(xscores)!=length(yscores)){print("lengths of input scores do not match")
	#				return(NULL)}
	if(any(is.na(xscores))|any(is.na(yscores))){print("Neither xscores nor yscores should contain any missing values")
						return(NULL)}
	#IF NOT SPECIFIED THE USE THE PLUG-IN FORMULA OF Andersson 
	#AND VON DAVIER 2014 ("Improving the Bandwidth Selection in Kernel Equating")
	#TO DETERMINE BANDWIDTH (EQUATION 12)
	if(is.na(hX)){hX=9*sd(xscores)/sqrt(100*(length(xscores)^0.4)-81)}
	if(is.na(hY)){hY=9*sd(yscores)/sqrt(100*(length(yscores)^0.4)-81)}


	tX=table(xscores)
	scoresX=as.numeric(names(tX))
	fX=as.numeric(tX)
	tY=table(yscores)
	scoresY=as.numeric(names(tY))
	fY=as.numeric(tY)
	return(KernelEquateFromDists(scoresX,fX,scoresY,fY,hX,hY))}


#' Apply Levine Observed Score equating to data from a NEAT design.
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param ws Vector of two elements denoting the relative weights of the dx population and the dy population in the synthetic population where equating takes place.
#' @param internal Logical input denoting whether the anchor test in internal or external (default) to the tests being equated.
#'
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{lys}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame combining the sorted unique scores on form X in the data and their equated values on form Y.}
#'   \item{muSx}{Estimated mean on form X within the synthetic population.}
#'   \item{sigmaSx}{Estimated standard deviation on form X within the synthetic population.}
#'   \item{muSy}{Estimated mean on form Y within the synthetic population.}
#'   \item{sigmaSy}{Estimated standard deviation on form Y within the synthetic population.}
#' }
#'
#' @examples
#' #Simulate two data sets with roughly equivalent relationship to underlying "true" scores but a difference in means
#' n1=8000
#' n2=5500
#' t1=rnorm(n1,0.5,1)
#' t2=rnorm(n2,0,1)
#' x=round(pmin(100,pmax(0,50+20*(0.9*t1+rnorm(n1,0,sqrt(1-0.9^2))))))
#' a1=round(pmin(20,pmax(0,10+4*(0.7*t1+rnorm(n1,0,sqrt(1-0.7^2))))))
#' cor(cbind(x,t1,a1))
#' y=round(pmin(100,pmax(0,50+20*(0.9*t2+rnorm(n2,0,sqrt(1-0.9^2))))))
#' a2=round(pmin(20,pmax(0,10+4*(0.7*t2+rnorm(n2,0,sqrt(1-0.7^2))))))
#' cor(cbind(y,t2,a2))
#' LevineObservedEquate(data.frame(x=x,a=a1),data.frame(y=y,a=a2))
#' #equated scores should be close to identity
#'
#' @keywords KernEqWPS
#' @references
#' Andersson, B., & von Davier, A. A. (2014). Improving the bandwidth selection in kernel equating. 
#' \emph{Journal of Educational Measurement, 51}(3), 223-238.
#' @export
LevineObservedEquate=function(dx,dy,ws=NA,internal=FALSE){

	if(!all(sort(names(dx))==c("a","x"))){print("dx should have columns names 'x' and 'a' only")
						return(NULL)}
	if(!all(sort(names(dy))==c("a","y"))){print("dy should have columns names 'y' and 'a' only")
						return(NULL)}

	alldat=rbind(data.frame(x=dx$x,y=rep(NA,dim(dx)[1]),a=dx$a),data.frame(x=rep(NA,dim(dy)[1]),y=dy$y,a=dy$a))
	#unless explicitly defined weights of two population defined by relative frequencies
	if(is.na(ws[1])){ws=as.numeric(prop.table(table(is.na(alldat$x))))}
	w1=ws[1]
	w2=ws[2]
	#from kolen and brennan equations 4.58 and 4.59
	if (internal==FALSE){gam1=(var(dx$x)+cov(dx$x,dx$a))/(var(dx$a)+cov(dx$x,dx$a))
		gam2=(var(dy$y)+cov(dy$y,dy$a))/(var(dy$a)+cov(dy$y,dy$a))}
	#from kolen and brennan equations 4.53 and 4.54
	if (internal==TRUE){gam1=var(dx$x)/cov(dx$x,dx$a)
		gam2=var(dy$y)/cov(dy$y,dy$a)}

	#obvious stuff
	mu1x=mean(dx$x)
	mu2y=mean(dy$y)
	sig1x=sd(dx$x)
	sig2y=sd(dy$y)
	mu1a=mean(dx$a)
	mu2a=mean(dy$a)
	sig1a=sd(dx$a)
	sig2a=sd(dy$a)
	
	#kolen and brennan equations 4.17 to 4.20
	muSx=mu1x-w2*gam1*(mu1a-mu2a)
	muSy=mu2y+w1*gam2*(mu1a-mu2a)
	sigSx=sqrt(
	(sig1x^2)-w2*(gam1^2)*((sig1a^2)-(sig2a^2))+w1*w2*(gam1^2)*((mu1a-mu2a)^2)  
	)
	sigSy=sqrt(
	(sig2y^2)+w1*(gam2^2)*((sig1a^2)-(sig2a^2))+w1*w2*(gam2^2)*((mu1a-mu2a)^2)  
	)

	lys=function(scores){(sigSy/sigSx)*(scores-muSx)+muSy}

	return(list(
		lys=lys
		,EqTable=data.frame(x=sort(unique(dx$x)),yx=lys(sort(unique(dx$x))))
		,muSx=muSx
		,sigSx=sigSx
		,muSy=muSy
		,sigSy=sigSy))}



#' Apply post-stratification Kernel Equating to data from a NEAT design.
#'
#' This method weights the data from the samples taking the two forms being equated
#' (form X and form Y) so that they have an equivalent distribution on the anchor test up to a defined order statistic. 
#' Then, Kernel Equating is applied to the weighted (and thus hopefully equivalent) groups.
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param target A character denoting the synthetic population. Can be "x" to denote weighting to the form X population, "y" to denote weighting to the form Y population or left blank (default) to indicate the data should be weighted to the average across all data.
#' @param order Numeric (integer) input denoting the order up to which the anchor distributions should be matched across populations. For example, 2 would indicate that the means and the means of the squared values (related to the standard deviation) should match. An order of 4 essentially indicates that the means, SDs, skewness and kurtosis of anchor test scores should be matched. The default is 5. Setting this to a lower value (2, 3 or 4) may be useful in the event of error messages (often related to non-convergence). Higher values may be valuable if the distributions of anchor test scores follow an extremely different shape on the two populations.
#' @param hX Bandwidth for form X. By default a plug-in estimator from Andersson and von Davier (2014) is used.
#' @param hY Bandwidth for form Y. By default a plug-in estimator from Andersson and von Davier (2014) is used.
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{KerEquiFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{KerLinFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y based on linear equating in the matched samples.}
#'   \item{EqTable}{A data frame combining the sorted unique scores on form X in the data and their equated values on form Y (both kernel equipercentile and linear equating) on matched data.}
#'   \item{muSx}{Estimated mean on form X within the synthetic population.}
#'   \item{sigmaSx}{Estimated standard deviation on form X within the synthetic population.}
#'   \item{muSy}{Estimated mean on form Y within the synthetic population.}
#'   \item{sigmaSy}{Estimated standard deviation on form Y within the synthetic population.}
#' }
#' @references
#' Andersson, B., & von Davier, A. A. (2014). Improving the bandwidth selection in kernel equating. 
#' \emph{Journal of Educational Measurement, 51}(3), 223-238.
#'
#' @examples
#' #Simulate two data sets with roughly equivalent relationship to underlying "true" scores but a difference in means
#' #This example shows the tendency for the PSE approach to under-adjust for the true difference in ability between groups
#' n1=8000
#' n2=5500
#' t1=rnorm(n1,0.5,1)
#' t2=rnorm(n2,0,1)
#' x=round(pmin(100,pmax(0,50+20*(0.9*t1+rnorm(n1,0,sqrt(1-0.9^2))))))
#' a1=round(pmin(20,pmax(0,10+4*(0.7*t1+rnorm(n1,0,sqrt(1-0.7^2))))))
#' cor(cbind(x,t1,a1))
#' y=round(pmin(100,pmax(0,50+20*(0.9*t2+rnorm(n2,0,sqrt(1-0.9^2))))))
#' a2=round(pmin(20,pmax(0,10+4*(0.7*t2+rnorm(n2,0,sqrt(1-0.7^2))))))
#' cor(cbind(y,t2,a2))
#' PSEObservedEquate(data.frame(x=x,a=a1),data.frame(y=y,a=a2))
#'
#' @keywords KernEqWPS
#' @export
PSEObservedEquate=function(dx,dy,target=NA,order=5,hX=NA,hY=NA){

	if(any(is.na(dx))|any(is.na(dy))){print("Neither dx nor dy should contain any missing values")
						return(NULL)}
	if(!all(sort(names(dx))==c("a","x"))){print("dx should have columns names 'x' and 'a' only")
						return(NULL)}
	if(!all(sort(names(dy))==c("a","y"))){print("dy should have columns names 'y' and 'a' only")
						return(NULL)}

	#IF NOT SPECIFIED THE USE THE PLUG-IN FORMULA OF ANDERSSON 
	#AND VON DAVIER 2014 ("Improving the Bandwidth Selection in Kernel Equating")
	#TO DETERMINE BANDWIDTH (EQUATION 12)
	if(is.na(hX)){hX=9*sd(dx$x)/sqrt(100*(dim(dx)[1]^0.4)-81)}
	if(is.na(hY)){hY=9*sd(dy$y)/sqrt(100*(dim(dy)[1]^0.4)-81)}

	alldat=rbind(data.frame(x=dx$x,y=rep(NA,dim(dx)[1]),a=dx$a),data.frame(x=rep(NA,dim(dy)[1]),y=dy$y,a=dy$a))
	if(is.na(target)){ws=as.numeric(prop.table(table(is.na(alldat$x))))
		targetgroup=NA}
	if(!is.na(target)){
	if(!target%in%c("x","y")){print("target population should be 'x', 'y' or be left blank to averae across")
					return(NULL)}
	if(target=="x"){ws=c(1,0)
		targetgroup=1}
	if(target=="y"){ws=c(0,1)
		targetgroup=2}}
	#FIRST CALCULATE WEIGHTS (WEIGHTING FIRST GROUP TO THE SECOND)
	summary(alldat)
	groups=rep(c(1,2),c(dim(dx)[1],dim(dy)[1]))

	#weight so equivalent anchor test score distribution (first four moments)
	#(standardise anchor test scores before weighting as this seems to help with convergence)
	alldat$a2=(alldat$a-mean(alldat$a))/sd(alldat$a)
	
	form1char="~a2"
	if (order>1){for (ords in 2:order){form1char=paste(form1char,"+I(a2^",ords,")",sep="")}}

	alldat$wt=MDIAwtGroups(formula=as.formula(form1char),data=alldat,groups=groups)

	tX=sum(!is.na(alldat$x))*tapply(alldat$wt[!is.na(alldat$x)],alldat$x[!is.na(alldat$x)],sum)
	scoresX=as.numeric(names(tX))
	fX=as.numeric(tX)
	
	tY=sum(!is.na(alldat$y))*tapply(alldat$wt[!is.na(alldat$y)],alldat$y[!is.na(alldat$y)],sum)
	scoresY=as.numeric(names(tY))
	fY=as.numeric(tY)

	KerEq2=KernelEquateFromDists(scoresX,fX,scoresY,fY,hX,hY)
	
	#Linear equating based on PSE
	msX=DistToMuSig(scoresX,fX)
	msY=DistToMuSig(scoresY,fY)
	linPSEfunc=function(x){msY$mu+(msY$sigma/msX$sigma)*(x-msX$mu)}
	
	return(list(
		KerEquiFunc=KerEq2$yxFunc,
		KerLinFunc=linPSEfunc,
		EqTable=data.frame(x=scoresX,equiyx=KerEq2$yxFunc(scoresX),linyx=linPSEfunc(scoresX)),
		muSx=msX$mu,
		sigSx=msX$sigma,
		muSy=msY$mu,
		sigSy=msY$sigma))}

#' Apply kernel chained equating
#'
#' Chained equating via an anchor test is done by successively applying the \emph{KernelEquateFromScoresEG} function between.
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param hX Bandwidth for form X. By default a plug-in estimator from Andersson and von Davier (2014) is used.
#' @param hY Bandwidth for form Y. By default a plug-in estimator from Andersson and von Davier (2014) is used.
#' @param hA Bandwidth for the anchor form. By default a plug-in estimator from Andersson and von Davier (2014) is used which will differ across populations. If a value is specified then the same bandwidth will be used in both populations.
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{chainedFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame combining the sorted unique scores on form X in the data and their equated values on form Y (chained kernel).}
#' }
#' @references
#' Andersson, B., & von Davier, A. A. (2014). Improving the bandwidth selection in kernel equating. 
#' \emph{Journal of Educational Measurement, 51}(3), 223-238.
#'
#' @examples
#' #Simulate two data sets with roughly equivalent relationship to underlying "true" scores but a difference in means
#' n1=8000
#' n2=5500
#' t1=rnorm(n1,0.5,1)
#' t2=rnorm(n2,0,1)
#' x=round(pmin(100,pmax(0,50+20*(0.9*t1+rnorm(n1,0,sqrt(1-0.9^2))))))
#' a1=round(pmin(20,pmax(0,10+4*(0.7*t1+rnorm(n1,0,sqrt(1-0.7^2))))))
#' cor(cbind(x,t1,a1))
#' y=round(pmin(100,pmax(0,50+20*(0.9*t2+rnorm(n2,0,sqrt(1-0.9^2))))))
#' a2=round(pmin(20,pmax(0,10+4*(0.7*t2+rnorm(n2,0,sqrt(1-0.7^2))))))
#' cor(cbind(y,t2,a2))
#' KernelChainedEquate(data.frame(x=x,a=a1),data.frame(y=y,a=a2))
#'
#' @keywords KernEqWPS
#' @export
KernelChainedEquate=function(dx,dy,hX=NA,hY=NA,hA=NA){
	KExa=KernelEquateFromScoresEG(dx$x,dx$a,hX=hX,hY=hA)
	KEay=KernelEquateFromScoresEG(dy$a,dy$y,hX=hA,hY=hY)
	scoresX=sort(unique(dx$x))
	chainedFunc=function(x){KEay$yxFunc(KExa$yxFunc(x))}
		return(list(chainedFunc=chainedFunc,
		EqTable=data.frame(x=scoresX,equiyx=chainedFunc(scoresX))))}





#' Method to estimate appropriate bandwidth using the Andersson and von Davier plug-in formula.
#' 
#' The formula for suggesting an appropriate bandwidth is very simple and can easily be seen by looking at the function code.
#' 
#' @param scores Vector of possible scores.
#' @param fX Number of people achieving each score. Weighted numbers are also acceptable but must add up to the overall sample size. Must be same length as scores.
#' @return The function returns a single value representing the suggested bandwidth.
#' @examples
#' #EXAMPLE OF FINDING BANDWIDTH TO EQUATE FIRST AND SECOND HALVES OF MATHS TEST
#' x1=rowSums(mathsdata[,1:62])
#' y1=rowSums(mathsdata[,63:124])
#' summary(x1)
#' summary(y1)
#' tabx1=table(x1)
#' taby1=table(y1)
#' 
#' RoughHx=FindAVDBandwidth(as.numeric(names(tabx1)),as.vector(tabx1))
#' RoughHx
#' 
#' RoughHy=FindAVDBandwidth(as.numeric(names(taby1)),as.vector(taby1))
#' RoughHy
#' 
#' eq1=KernelEquateFromScoresEG(x1,y1,hX=RoughHx,hY=RoughHy)
#' plot(sort(unique(x1)),eq1$yx,type='l')
#' 
#' #prove that above code is the same as the default method of equating
#' eq2=KernelEquateFromScoresEG(x1,y1)
#' lines(sort(unique(x1)),eq2$yx,lty=2,col="purple")
#' @references
#' Andersson, B., & von Davier, A. A. (2014). Improving the bandwidth selection in kernel equating. 
#' \emph{Journal of Educational Measurement, 51}(3), 223-238.
#' @export
#' @keywords KernEqWPS
#' @export
FindAVDBandwidth=function(scores,fX){
	mean1=sum(scores*fX)/sum(fX)
	sd1=sqrt(sum(((scores-mean1)^2)*fX)/(sum(fX)-1))
	n=sum(fX)
	hX=9*sd1/sqrt(100*(n^0.4)-81)
	return(hX)
	}



