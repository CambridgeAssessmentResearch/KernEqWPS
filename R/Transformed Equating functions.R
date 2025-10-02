#' CSEM stabilising transformation followed by linear equating
#'
#' Just linear equating but with a particular atanh transformation before we begin. The transformation is intended to
#' make it so the CSEM will be stable across the score range. It is loosely based on function that is almost linearly related to one suggested in
#' Conditional Standard Errors of Measurement for Scale Scores,
#' Michael J. Kolen, Bradley A. Hanson and Robert L. Brennan,
#' Journal of Educational Measurement , Winter, 1992, Vol. 29, No. 4 (Winter,1992), pp. 285-307.
#' 
#' We assume the minimum available score on each test is zero.
#' 
#' @param xscores A vector of scores from individual candidates on form X.
#' @param yscores A vector of scores from individual candidates on form Y.
#' @param maxX Maximum score on form X.
#' @param maxY Maximum score on form Y.
#'
#' @references
#' Benton, T. & Carroll, M. (2025). The biggest equating study in the world ever. Research report. \emph{Cambridge University Press & Assessment}.
#' 
#' @examples
#' x=round(pmax(0,pmin(100,rnorm(500,50,10))))
#' y=round(pmax(0,pmin(100,rnorm(500,60,10))))
#' tleq=TransformedLinearEquate(x,y,100,100)
#' tleq
#' 
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between 0 and maxxX on form X.}
#' }
#'
#' @export
TransformedLinearEquate=function(xscores,yscores,maxX,maxY){

        #first make the transformaation functions
        transx=function(x){atanh(((x/maxX)-0.5)/0.6)}
        transy=function(y){atanh(((y/maxY)-0.5)/0.6)}
	#...and their inverses
        atransx=function(xt){maxX*(0.5+0.6*tanh(xt))}
        atransy=function(yt){maxY*(0.5+0.6*tanh(yt))}

	xt=transx(xscores)
	yt=transy(yscores)
	lineq=LinearEquate(xt,yt)

	yxFunc=function(xscores){
		atransy(lineq$lys(transx(xscores)))
		}

	scoresX=0:ceiling(maxX)
	return(list(
		yxFunc=yxFunc
		,EqTable=data.frame(x=scoresX,equiyx=yxFunc(scoresX))
		))
	}

#' Chained linear equating after preliminary CSEM stabilising transformation
#'
#' Just chained linear equating but with a particular atanh transformation before we begin. The transformation is intended to
#' make it so the CSEM will be stable across the score range. It is loosely based on function that is almost linearly related to one suggested in
#' Conditional Standard Errors of Measurement for Scale Scores,
#' Michael J. Kolen, Bradley A. Hanson and Robert L. Brennan,
#' Journal of Educational Measurement , Winter, 1992, Vol. 29, No. 4 (Winter,1992), pp. 285-307.
#' 
#' We assume the minimum available score on each test is zero.
#' 
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param maxX Maximum score on form X.
#' @param maxY Maximum score on form Y.
#' @param maxA Maximum score on form A.
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
#' tl_eq=TransformedChainedLinearEquate(data.frame(x=scoresX1,a=scoresA1),data.frame(y=scoresY2,a=scoresA2),30,30,5)
#' tl_eq
#' 
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between 0 and maxxX on form X.}
#' }
#'
#' @export
TransformedChainedLinearEquate=function(dx,dy,maxX,maxY,maxA){

        #first make the transformaation functions
        transx=function(x){atanh(((x/maxX)-0.5)/0.6)}
        transy=function(y){atanh(((y/maxY)-0.5)/0.6)}
        transa=function(a){atanh(((a/maxA)-0.5)/0.6)}
	#...and their inverses
        atransx=function(xt){maxX*(0.5+0.6*tanh(xt))}
        atransy=function(yt){maxY*(0.5+0.6*tanh(yt))}
        atransa=function(at){maxA*(0.5+0.6*tanh(at))}

	dxt=dx
	dxt$x=transx(dxt$x)
	dxt$a=transa(dxt$a)
	dyt=dy
	dyt$y=transy(dyt$y)
	dyt$a=transa(dyt$a)

	lineq=ChainedLinearEquate(dxt,dyt)

	yxFunc=function(xscores){
		atransy(lineq$chainedFunc(transx(xscores)))
		}

	scoresX=0:ceiling(maxX)
	return(list(
		yxFunc=yxFunc
		,EqTable=data.frame(x=scoresX,equiyx=yxFunc(scoresX))
		))
	}

#' Tucker linear equating after preliminary CSEM stabilising transformation
#'
#' Just linear equating using the Tucker method for a NEAT design but with a particular atanh transformation before we begin. The transformation is intended to
#' make it so the CSEM will be stable across the score range. It is loosely based on function that is almost linearly related to one suggested in
#' Conditional Standard Errors of Measurement for Scale Scores,
#' Michael J. Kolen, Bradley A. Hanson and Robert L. Brennan,
#' Journal of Educational Measurement , Winter, 1992, Vol. 29, No. 4 (Winter,1992), pp. 285-307.
#' 
#' We assume the minimum available score on each test is zero.
#' 
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param maxX Maximum score on form X.
#' @param maxY Maximum score on form Y.
#' @param maxA Maximum score on form A.
#' @param ws Vector of two elements denoting the relative weights of the dx population and the dy population in the synthetic population where equating takes place.
#'
#' 
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{yxFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between 0 and maxxX on form X.}
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
#' tl_eq=TransformedTuckerEquate(data.frame(x=scoresX1,a=scoresA1),data.frame(y=scoresY2,a=scoresA2),30,30,5)
#' tl_eq
#' 
#' 
#' @export
TransformedTuckerEquate=function(dx,dy,maxX,maxY,maxA,ws=NA){

	#set up ws if not supplied
	alldat=rbind(data.frame(x=dx$x,y=rep(NA,dim(dx)[1]),a=dx$a),data.frame(x=rep(NA,dim(dy)[1]),y=dy$y,a=dy$a))
	#unless explicitly defined weights of two population defined by relative frequencies
	if(is.na(ws[1])){ws=as.numeric(prop.table(table(is.na(alldat$x))))}

        #first make the transformaation functions
        transx=function(x){atanh(((x/maxX)-0.5)/0.6)}
        transy=function(y){atanh(((y/maxY)-0.5)/0.6)}
        transa=function(a){atanh(((a/maxA)-0.5)/0.6)}
	#...and their inverses
        atransx=function(xt){maxX*(0.5+0.6*tanh(xt))}
        atransy=function(yt){maxY*(0.5+0.6*tanh(yt))}
        atransa=function(at){maxA*(0.5+0.6*tanh(at))}

	dxt=dx
	dxt$x=transx(dxt$x)
	dxt$a=transa(dxt$a)
	dyt=dy
	dyt$y=transy(dyt$y)
	dyt$a=transa(dyt$a)

	lineq=TuckerEquate(dxt,dyt,ws)

	yxFunc=function(xscores){
		atransy(lineq$lys(transx(xscores)))
		}

	scoresX=0:ceiling(maxX)
	return(list(
		yxFunc=yxFunc
		,EqTable=data.frame(x=scoresX,equiyx=yxFunc(scoresX))
		))
	}

