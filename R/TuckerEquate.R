#' Apply Tucker Observed Score equating to data from a NEAT design.
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param ws Vector of two elements denoting the relative weights of the dx population and the dy population in the synthetic population where equating takes place.
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
#' TuckerEquate(data.frame(x=x,a=a1),data.frame(y=y,a=a2))
#'
#' @export
TuckerEquate=function(dx,dy,ws=NA){

	if(!all(sort(names(dx))==c("a","x"))){print("dx should have columns names 'x' and 'a' only")
						return(NULL)}
	if(!all(sort(names(dy))==c("a","y"))){print("dy should have columns names 'y' and 'a' only")
						return(NULL)}

	alldat=rbind(data.frame(x=dx$x,y=rep(NA,dim(dx)[1]),a=dx$a),data.frame(x=rep(NA,dim(dy)[1]),y=dy$y,a=dy$a))
	#unless explicitly defined weights of two population defined by relative frequencies
	if(is.na(ws[1])){ws=as.numeric(prop.table(table(is.na(alldat$x))))}
	w1=ws[1]
	w2=ws[2]
	#from kolen and brennan equations 4.21 and 4.22
        gam1=cov(dx$x,dx$a)/var(dx$a)
        gam2=cov(dy$y,dy$a)/var(dy$a)

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


