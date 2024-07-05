#' Find circle-arc equating functions based on minima, maxima and means for two tests
#'
#' Just apply the standard circle-arc equating functions based on the three points on the arc for each form.
#'
#' @param meanX Mean score on form X (or any midpoint for circle arc method).
#' @param meanY Mean score on form Y (or any midpoint for circle arc method).
#' @param minX Minimum score on form X (or lower point for circle arc method).
#' @param minY Minimum score on form Y (or lower point for circle arc method).
#' @param maxX Maximum score on form X (or upper point for circle arc method).
#' @param maxY Maximum score on form Y (or upper point for circle arc method).
#'
#' @examples
#' circ_eq=CircleArcFromMeans(20, 30, 0, 0, 50, 50)
#' plot(circ_eq$EqTable$x,circ_eq$EqTable$equiyx,type='l')
#' points(c(0,20,50),c(0,30,50))#show the three points
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between minX and maxxX on form X.}
#' }
#'
#' @export
CircleArcFromMeans=function(meanX,meanY,minX,minY,maxX,maxY){

	#switch to notation of Livingston and Kim
	#pages 334 to 336 of 
	#Livingston and Kim (2008) from https://www.jstor.org/stable/25651516
	x1=minX
	x2=meanX
	x3=maxX
	y1=minY
	y2=meanY
	y3=maxY

	#if completely linear relationship then centre of circle will be infinitely far away
	#need to deal with this instance first to avoid (rare) errors
	checkcor=cor(c(x1,x2,x3),c(y1,y2,y3))
	if(checkcor==1){
    		yxFunc = function(xscores){y1+(xscores-x1)*(y3-y1)/(x3-x1)}
    		scoresX = floor(x1):ceiling(x2)
    		return(list(yxFunc = yxFunc, EqTable = data.frame(x = scoresX, 
        		equiyx = yxFunc(scoresX))))
	}


	#return to applying main formulae
	Lx2=(y1+((y3-y1)/(x3-x1))*(x2-x1))#L(x2) in the notation of Livingston and Kim
	y2s=y2-Lx2 #y2s is y2* in the notation of Livingston and Kim

	xc=(x3^2-x1^2)/(2*(x3-x1))
	yc=((x1^2)*(x3-x2)-(x2^2+y2s^2)*(x3-x1)+(x3^2)*(x2-x1))/(2*(y2s*(x1-x3)))
	r=sqrt((xc-x1)^2+(yc^2))

	yxFunc=function(xscores){
		if(y2s>=0){yscoresstar=yc+sqrt(r^2-(xscores-xc)^2)}
		if(y2s<0){yscoresstar=yc-sqrt(r^2-(xscores-xc)^2)}
		return(yscoresstar+(y1+((y3-y1)/(x3-x1))*(xscores-x1)))
		}

	scoresX=floor(minX):ceiling(maxX)
	return(list(
		yxFunc=yxFunc
		,EqTable=data.frame(x=scoresX,equiyx=yxFunc(scoresX))
		))
	}

#' Circle-arc based on two vectors of data assuming from equivalent groups
#'
#' Calculates means an applies .
#'
#' @param xscores A vector of scores from individual candidates on form X.
#' @param yscores A vector of scores from individual candidates on form Y.
#' @param maxX Maximum score on form X (or upper point for circle arc method).
#' @param maxY Maximum score on form Y (or upper point for circle arc method).
#' @param minX Minimum score on form X (or lower point for circle arc method). Default of zero.
#' @param minY Minimum score on form Y (or lower point for circle arc method). Default of zero.
#'
#' @examples
#' x=round(runif(500,0,50))
#' y=round(pmin(50,pmax(0,rnorm(500,35,10))))
#' circ_eq=CircleArcEquate(x,y,50,50)
#' plot(circ_eq$EqTable$x,circ_eq$EqTable$equiyx,type='l')
#' circ_eq$EqTable
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between minX and maxxX on form X.}
#' }
#'
#' @export
CircleArcEquate=function(xscores,yscores,maxX,maxY,minX=0,minY=0){

	#use notation of Livingston and Kim
	#pages 334 to 336 of 
	#Livingston and Kim (2008) from https://www.jstor.org/stable/25651516
	x1=minX
	x2=mean(xscores)
	x3=maxX
	y1=minY
	y2=mean(yscores)
	y3=maxY

	Lx2=(y1+((y3-y1)/(x3-x1))*(x2-x1))#L(x2) in the notation of Livingston and Kim
	y2s=y2-Lx2 #y2s is y2* in the notation of Livingston and Kim

	xc=(x3^2-x1^2)/(2*(x3-x1))
	yc=((x1^2)*(x3-x2)-(x2^2+y2s^2)*(x3-x1)+(x3^2)*(x2-x1))/(2*(y2s*(x1-x3)))
	r=sqrt((xc-x1)^2+(yc^2))

	yxFunc=function(xscores){
		if(y2s>=0){yscoresstar=yc+sqrt(r^2-(xscores-xc)^2)}
		if(y2s<0){yscoresstar=yc-sqrt(r^2-(xscores-xc)^2)}
		return(yscoresstar+(y1+((y3-y1)/(x3-x1))*(xscores-x1)))
		}

	scoresX=floor(minX):ceiling(maxX)
	return(list(
		yxFunc=yxFunc
		,EqTable=data.frame(x=scoresX,equiyx=yxFunc(scoresX))
		))
	}

#' Apply chained equating via the circle-arc method
#'
#' Chained equating via an anchor test is done by successively applying the \emph{CircleArcEquate} function between anchors and main tests.
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param maxX Maximum score on form X (or upper point for circle arc method).
#' @param maxY Maximum score on form Y (or upper point for circle arc method).
#' @param maxA Maximum score on anchor (or upper point for circle arc method).
#' @param minX Minimum score on form X (or lower point for circle arc method). Default of zero.
#' @param minY Minimum score on form Y (or lower point for circle arc method). Default of zero.
#' @param minA Minimum score on anchor (or lower point for circle arc method). Default of zero.
#' 
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{chainedFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between minX and maxxX on form X.}
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
#' circ_eq=CircleArcChainedEquate(data.frame(x=scoresX1,a=scoresA1),data.frame(y=scoresY2,a=scoresA2),30,30,5)
#' circ_eq
#' 
#' @keywords KernEqWPS
#' @export
CircleArcChainedEquate=function(dx,dy,maxX,maxY,maxA,minX=0,minY=0,minA=0){
	CAxa=CircleArcEquate(dx$x,dx$a,maxX=maxX,maxY=maxA,minX=minX,minY=minA)
	CAay=CircleArcEquate(dy$a,dy$y,maxX=maxA,maxY=maxY,minX=minA,minY=minY)
	scoresX=floor(minX):ceiling(maxX)

	chainedFunc=function(x){CAay$yxFunc(CAxa$yxFunc(x))}

	return(list(
		chainedFunc=chainedFunc,
		EqTable=data.frame(x=scoresX,equiyx=chainedFunc(scoresX))
		))
	}


#' Apply circle-arc method based on Tucker equating
#'
#' Circle-arc equating via an anchor test is done by using the Tucker method to estimate the means on each from in the synthetic population.
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param ws Vector of two elements denoting the relative weights of the dx population and the dy population in the synthetic population where equating takes place.
#' @param maxX Maximum score on form X (or upper point for circle arc method).
#' @param maxY Maximum score on form Y (or upper point for circle arc method).
#' @param minX Minimum score on form X (or lower point for circle arc method). Default of zero.
#' @param minY Minimum score on form Y (or lower point for circle arc method). Default of zero.
#' 
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between minX and maxxX on form X.}
#'   \item{meanX}{Estimated mean form X score in synthetic population.}
#'   \item{meanY}{Estimated mean form X score in synthetic population.}
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
#' circ_eq=CircleArcTuckerEquate(data.frame(x=scoresX1,a=scoresA1),data.frame(y=scoresY2,a=scoresA2),NA,30,30)
#' circ_eq
#' 
#' @keywords KernEqWPS
#' @export
CircleArcTuckerEquate=function(dx,dy,ws=NA,maxX,maxY,minX=0,minY=0){
	te=TuckerEquate(dx,dy,ws)
	meanX=te$muSx
	meanY=te$muSy
	outlist=CircleArcFromMeans(meanX,meanY,minX,minY,maxX,maxY)
	outlist$meanX=meanX
	outlist$meanY=meanY
	return(outlist)
	}

