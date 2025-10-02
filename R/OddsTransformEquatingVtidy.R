#tidied code for odds-transform equating
#includes code for the confidence interval

#' Find odds-transform equating functions based on minima, maxima, means and SDs for two tests
#'
#' Rather like circle-arc equating and equating function to intercept with the minima (assumed to 0)
#' , maxima and means of each test. Can optionally also incorporate information on standard deviations.
#'
#' @param meanX Mean score on form X.
#' @param meanY Mean score on form Y.
#' @param maxX Maximum score on form X.
#' @param maxY Maximum score on form Y.
#' @param sdX SD of scores on form X. NA by default in which case estimatek is set to FALSE.
#' @param sdX SD of scores on form X. NA by default in which case estimatek is set to FALSE.
#' @param estimatek Logical value denoting whether the k parameter should be estimated. If FALSE (default) then k is set to 1
#' and the simple approach to Odds-transform equating is applied.
#' @param semeanX Optionally add the standard error of the mean of scores on form X. If both semeanX and semeanY are supplied, and
#' estimatek=FALSE, then an analytic confidence interval will be produced based upon the delta method.
#' @param semeanY Optionally add the standard error of the mean of scores on form X. If both semeanX and semeanY are supplied, and
#' estimatek=FALSE, then an analytic confidence interval will be produced based upon the delta method.
#'
#' Method assumes the minimum available score is zero.
#'
#' @references
#' Benton, T. & Carroll, M. (2025). The biggest equating study in the world ever. Research report. \emph{Cambridge University Press & Assessment}.
#' 
#' @examples
#' ot_eq=OddsTransformFromMeans(20, 30, 50, 50)
#' plot(ot_eq$EqTable$x,ot_eq$EqTable$equiyx,type='l')
#' points(c(0,20,50),c(0,30,50))#show the three points
#'
#' #incorporating SDs
#' ot_eq=OddsTransformFromMeans(20, 45, 50, 50 ,sdX=8,sdY=5,estimatek=TRUE)
#' plot(ot_eq$EqTable$x,ot_eq$EqTable$equiyx,type='l')
#' points(c(0,20,50),c(0,45,50))#show the three points
#' #check reversible
#' ot_eq_rev=OddsTransformFromMeans(45, 20, 50, 50,sdX=5,sdY=8,estimatek=TRUE)
#' plot(ot_eq_rev$EqTable$x,ot_eq_rev$EqTable$equiyx,type='l')
#' points(c(0,45,50),c(0,20,50))#show the three points
#' ot_eq_rev$yxFunc(ot_eq$yxFunc(0:10))
#'
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between 0 and maxX on form X.}
#' }
#'
#' @export
OddsTransformFromMeans=function(meanX,meanY,maxX,maxY,sdX=NA,sdY=NA,estimatek=FALSE,semeanX=NA,semeanY=NA){

	k=1#set default value for k
	#will not overwrite unless both SDs are supplied
	if(is.na(sdX)){estimatek=FALSE}
	if(is.na(sdY)){estimatek=FALSE}

	#find k if required
	if(estimatek==TRUE){
	    k=(sdX*maxX*meanY*(maxY-meanY))/(sdY*maxY*meanX*(maxX-meanX))
	}

	yxFunc=function(x){
		f=(1/k)*log((maxX-meanX)/meanX)+log(meanY/(maxY-meanY))+(1/k)*log(x/(maxX-x))
		y=maxY*exp(f)/(1+exp(f))
		y[x==maxX]=maxY
		return(y)
		}

	scoresX=0:ceiling(maxX)
	EqTable=data.frame(x=scoresX,equiyx=yxFunc(scoresX))

	if(!is.na(semeanX) & !is.na(semeanY) & estimatek==FALSE & k==1){
		x=scoresX
		f=log((maxX-meanX)/meanX)+log(meanY/(maxY-meanY))+log(x/(maxX-x))
		sef=sqrt(
			((maxX*semeanX)/((maxX-meanX)*meanX))^2+
			((maxY*semeanY)/((maxY-meanY)*meanY))^2
			)
		f_lower=f-1.96*sef
		f_upper=f+1.96*sef
		EqTable$equiyx_lower=maxY*exp(f_lower)/(1+exp(f_lower))
		EqTable$equiyx_upper=maxY*exp(f_upper)/(1+exp(f_upper))
		EqTable$equiyx_lower[x==maxX]=maxY
		EqTable$equiyx_upper[x==maxX]=maxY
		}

	return(list(
		yxFunc=yxFunc
		,EqTable=EqTable
		))
	}


#example
#ot_eq=OddsTransformFromMeans(meanX=33.1,meanY=43.2,maxX=80,maxY=120)
#ot_eq$EqTable
#ot_eq$yxFunc(33.1)
#plot(ot_eq$EqTable$x,ot_eq$EqTable$equiyx,type='l')

#with the CI
#ot_eq=OddsTransformFromMeans(meanX=33.1,meanY=43.2,maxX=80,maxY=120,semeanX=2.1,semeanY=2.2)
#ot_eq$EqTable


#' Odds-transform equating based on two vectors of data assuming from equivalent groups
#'
#' Calculates means and sds and applies OddsTransformFromMeans.
#'
#' @param xscores A vector of scores from individual candidates on form X.
#' @param yscores A vector of scores from individual candidates on form Y.
#' @param maxX Maximum score on form X.
#' @param maxY Maximum score on form Y.
#' @param estimatek Logical value denoting whether the k parameter should be estimated. If FALSE (default) then k is set to 1
#' @param CI Logical value denoting whether to attempt to return confidence intervals. This assumes that the data comes from separate independent samples.
#'
#' @references
#' Benton, T. & Carroll, M. (2025). The biggest equating study in the world ever. Research report. \emph{Cambridge University Press & Assessment}.
#' 
#' @examples
#' x=round(runif(500,0,50))
#' y=round(pmin(50,pmax(0,rnorm(500,39,10))))
#' ot_eq=OddsTransformEquate(x,y,50,50,estimatek=TRUE)
#' plot(ot_eq$EqTable$x,ot_eq$EqTable$equiyx,type='l')
#' points(mean(x),mean(y))
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between 0 and maxX on form X.}
#' }
#'
#' @export
OddsTransformEquate=function(xscores,yscores,maxX,maxY,estimatek=FALSE,CI=FALSE){

  xscores=xscores[!is.na(xscores)]
  yscores=yscores[!is.na(yscores)]

  meanX=mean(xscores)
  meanY=mean(yscores)
  sdX=sd(xscores)
  sdY=sd(yscores)
  semeanX=sd(xscores)/sqrt(length(xscores))
  semeanY=sd(yscores)/sqrt(length(yscores))
  if(CI==FALSE){answer=OddsTransformFromMeans(meanX,meanY,maxX,maxY,sdX,sdY,estimatek)}
  if(CI==TRUE){answer=OddsTransformFromMeans(meanX,meanY,maxX,maxY,sdX,sdY,estimatek,semeanX,semeanY)}

	return(answer)
	}

#example code
#x=round(pmin(80,pmax(0,rnorm(500,33.1,8))))
#y=round(pmin(120,pmax(0,rnorm(500,43.2,8))))
#ot_eq=OddsTransformEquate(x,y,80,120,CI=TRUE)
#ot_eq$EqTable
#plot(ot_eq$EqTable$x,ot_eq$EqTable$equiyx,type='l')
#points(mean(x),mean(y))


#' Apply chained equating via the odds-transform method
#'
#' Chained equating via an anchor test is done by successively applying the \emph{OddsTransformEquate} function between anchors and main tests.
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param maxX Maximum score on form X.
#' @param maxY Maximum score on form Y.
#' @param maxA Maximum score on anchor.
#' @param estimatek Logical value denoting whether the k parameter should be estimated. If FALSE (default) then k is set to 1.
#'
#' @references
#' Benton, T. & Carroll, M. (2025). The biggest equating study in the world ever. Research report. \emph{Cambridge University Press & Assessment}.
#' 
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{chainedFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between 0 and maxX on form X.}
#' }
#'
#' @examples
#' #demonstrate method on a 30 item test with an internal 5 item anchor
#' #define 30 rasch item difficulties as equally spread
#' itedifs=rep(seq(-2,2,length=5),6)
#'
#' #simulate population one item scores (and then form scores)
#' n1=300
#' t1=rnorm(n1,0.5,1)
#' ites1=0+(plogis(t1%*%t(rep(1,30))-rep(1,length(t1))%*%t(itedifs))>matrix(runif(n1*30),nrow=n1))
#' scoresX1=rowSums(ites1[,1:30])
#' scoresA1=rowSums(ites1[,26:30])
#' #simulate parallel tests in population two
#' n2=3000
#' t2=rnorm(n2,0,1)
#' ites2=0+(plogis(t2%*%t(rep(1,30))-rep(1,length(t2))%*%t(itedifs))>matrix(runif(n2*30),nrow=n2))
#' scoresY2=rowSums(ites2[,1:30])
#' scoresA2=rowSums(ites2[,26:30])
#'
#' ot_eq=OddsTransformChainedEquate(data.frame(x=scoresX1,a=scoresA1),data.frame(y=scoresY2,a=scoresA2),30,30,5)
#' ot_eq
#'
#' @keywords KernEqWPS
#' @export
OddsTransformChainedEquate=function(dx,dy,maxX,maxY,maxA,estimatek=FALSE){
	OTxa=OddsTransformEquate(dx$x,dx$a,maxX=maxX,maxY=maxA,estimatek=estimatek)
	OTay=OddsTransformEquate(dy$a,dy$y,maxX=maxA,maxY=maxY,estimatek=estimatek)
	scoresX=0:ceiling(maxX)

	chainedFunc=function(x){OTay$yxFunc(OTxa$yxFunc(x))}

	return(list(
		chainedFunc=chainedFunc,
		EqTable=data.frame(x=scoresX,equiyx=chainedFunc(scoresX))
		))
	}

#' Apply odds-transform method based on Tucker equating
#'
#' Odds-transform equating via an anchor test is done by using the Tucker method to estimate the means and SDs on each from in the synthetic population.
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param ws Vector of two elements denoting the relative weights of the dx population and the dy population in the synthetic population where equating takes place.
#' @param maxX Maximum score on form X.
#' @param maxY Maximum score on form Y.
#' @param estimatek Logical value denoting whether the k parameter should be estimated. If FALSE (default) then k is set to 1.
#'
#' @references
#' Benton, T. & Carroll, M. (2025). The biggest equating study in the world ever. Research report. \emph{Cambridge University Press & Assessment}.
#' 
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between 0 and maxX on form X.}
#'   \item{meanX}{Estimated mean form X score in synthetic population.}
#'   \item{meanY}{Estimated mean form Y score in synthetic population.}
#'   \item{sdX}{Estimated SD of form X scores in synthetic population.}
#'   \item{sdY}{Estimated SD form Y scores in synthetic population.}
#' }
#'
#' @examples
#' #demonstrate method on a 30 item test with an internal 5 item anchor
#' #define 30 rasch item difficulties as equally spread
#' itedifs=rep(seq(-2,2,length=5),6)
#'
#' #simulate population one item scores (and then form scores)
#' n1=300
#' t1=rnorm(n1,0.5,1)
#' ites1=0+(plogis(t1%*%t(rep(1,30))-rep(1,length(t1))%*%t(itedifs))>matrix(runif(n1*30),nrow=n1))
#' scoresX1=rowSums(ites1[,1:30])
#' scoresA1=rowSums(ites1[,26:30])
#' #simulate parallel tests in population two
#' n2=3000
#' t2=rnorm(n2,0,1)
#' ites2=0+(plogis(t2%*%t(rep(1,30))-rep(1,length(t2))%*%t(itedifs))>matrix(runif(n2*30),nrow=n2))
#' scoresY2=rowSums(ites2[,1:30])
#' scoresA2=rowSums(ites2[,26:30])
#'
#' ot_eq=OddsTransformTuckerEquate(data.frame(x=scoresX1,a=scoresA1),data.frame(y=scoresY2,a=scoresA2),NA,30,30,5)
#' ot_eq
#'
#'
#' @keywords KernEqWPS
#' @export
OddsTransformTuckerEquate=function(dx,dy,ws=NA,maxX,maxY,estimatek=FALSE){
	te=TuckerEquate(dx,dy,ws)
	meanX=te$muSx
	meanY=te$muSy
	sdX=te$sigSx
	sdY=te$sigSy
	outlist=OddsTransformFromMeans(meanX,meanY,maxX,maxY,sdX,sdY,estimatek)
	outlist$meanX=meanX
	outlist$meanY=meanY
	outlist$sdX=sdX
	outlist$sdY=sdY
	return(outlist)
	}


#' Apply odds-transform method based on Levine equating
#'
#' Odds-transform equating via an anchor test is done by using the Levine method to estimate the means and SDs on each form in the synthetic population.
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param ws Vector of two elements denoting the relative weights of the dx population and the dy population in the synthetic population where equating takes place.
#' @param maxX Maximum score on form X.
#' @param maxY Maximum score on form Y.
#' @param estimatek Logical value denoting whether the k parameter should be estimated. If FALSE (default) then k is set to 1.
#' @param internal Logical input denoting whether the anchor test in internal or external (default) to the tests being equated.
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between 0 and maxX on form X.}
#'   \item{meanX}{Estimated mean form X score in synthetic population.}
#'   \item{meanY}{Estimated mean form Y score in synthetic population.}
#'   \item{sdX}{Estimated SD of form X scores in synthetic population.}
#'   \item{sdY}{Estimated SD form Y scores in synthetic population.}
#' }
#'
#' @references
#' Benton, T. & Carroll, M. (2025). The biggest equating study in the world ever. Research report. \emph{Cambridge University Press & Assessment}.
#' 
#' @examples
#' #demonstrate method on a 30 item test with an internal 5 item anchor
#' #define 30 rasch item difficulties as equally spread
#' itedifs=rep(seq(-2,2,length=5),6)
#'
#' #simulate population one item scores (and then form scores)
#' n1=300
#' t1=rnorm(n1,0.5,1)
#' ites1=0+(plogis(t1%*%t(rep(1,30))-rep(1,length(t1))%*%t(itedifs))>matrix(runif(n1*30),nrow=n1))
#' scoresX1=rowSums(ites1[,1:30])
#' scoresA1=rowSums(ites1[,26:30])
#' #simulate parallel tests in population two
#' n2=3000
#' t2=rnorm(n2,0,1)
#' ites2=0+(plogis(t2%*%t(rep(1,30))-rep(1,length(t2))%*%t(itedifs))>matrix(runif(n2*30),nrow=n2))
#' scoresY2=rowSums(ites2[,1:30])
#' scoresA2=rowSums(ites2[,26:30])
#'
#' ot_eq=OddsTransformLevineEquate(data.frame(x=scoresX1,a=scoresA1),data.frame(y=scoresY2,a=scoresA2),NA,30,30,FALSE,TRUE)
#' ot_eq
#'
#'
#' @keywords KernEqWPS
#' @export
OddsTransformLevineEquate=function(dx,dy,ws=NA,maxX,maxY,estimatek=FALSE,internal=FALSE){
	loe=LevineObservedEquate(dx,dy,ws,internal=internal)
	meanX=loe$muSx
	meanY=loe$muSy
	sdX=loe$sigSx
	sdY=loe$sigSy
	outlist=OddsTransformFromMeans(meanX,meanY,maxX,maxY,sdX,sdY,estimatek)
	outlist$meanX=meanX
	outlist$meanY=meanY
	outlist$sdX=sdX
	outlist$sdY=sdY
	return(outlist)
	}

#' Apply nominal weights equating via the odds-transform method
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param ws Vector of two elements denoting the relative weights of the dx population and the dy population in the synthetic population where equating takes place.
#' @param maxX Maximum score on form X.
#' @param maxY Maximum score on form Y.
#' @param maxA Maximum score on anchor.
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between 0 and maxX on form X.}
#'   \item{meanX}{Estimated mean form X score in synthetic population.}
#'   \item{meanY}{Estimated mean form Y score in synthetic population.}
#' }
#'
#' @references
#' Benton, T. & Carroll, M. (2025). The biggest equating study in the world ever. Research report. \emph{Cambridge University Press & Assessment}.
#' 
#' @examples
#' #demonstrate method on a 30 item test with an internal 5 item anchor
#' #define 30 rasch item difficulties as equally spread
#' itedifs=rep(seq(-2,2,length=5),6)
#'
#' #simulate population one item scores (and then form scores)
#' n1=300
#' t1=rnorm(n1,0.5,1)
#' ites1=0+(plogis(t1%*%t(rep(1,30))-rep(1,length(t1))%*%t(itedifs))>matrix(runif(n1*30),nrow=n1))
#' scoresX1=rowSums(ites1[,1:30])
#' scoresA1=rowSums(ites1[,26:30])
#' #simulate parallel tests in population two
#' n2=3000
#' t2=rnorm(n2,0,1)
#' ites2=0+(plogis(t2%*%t(rep(1,30))-rep(1,length(t2))%*%t(itedifs))>matrix(runif(n2*30),nrow=n2))
#' scoresY2=rowSums(ites2[,1:30])
#' scoresA2=rowSums(ites2[,26:30])
#'
#' ot_eq=OddsTransformNominalEquate(data.frame(x=scoresX1,a=scoresA1),data.frame(y=scoresY2,a=scoresA2),c(0.5,0.5),30,30,5)
#' ot_eq
#'
#' @keywords KernEqWPS
#' @export
OddsTransformNominalEquate=function(dx,dy,ws,maxX,maxY,maxA){

	ws=ws/sum(ws)
	meanX=mean(dx$x)+ws[2]*(maxX/maxA)*(mean(dy$a)-mean(dx$a))
	meanY=mean(dy$y)+ws[1]*(maxY/maxA)*(mean(dx$a)-mean(dy$a))
	outlist=OddsTransformFromMeans(meanX,meanY,maxX,maxY)
	outlist$meanX=meanX
	outlist$meanY=meanY
	return(outlist)
	}

