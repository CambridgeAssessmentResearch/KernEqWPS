#' Apply linear equating to data from a EG design.
#'
#' @param xscores A vector of scores from individual candidates on form X.
#' @param yscores A vector of scores from individual candidates on form Y.
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{lys}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{yx}{A vector representing the scores derived from applying yxFunc to a sorted list of all the unqiue form X scores that actually occur in the data.}
#' }
#'
#' @examples
#' #Simulate two vectors with differencein means of 10
#' x=rnorm(500,50,10)
#' y=rnorm(500,60,10)
#' lineq=LinearEquate(x,y)
#' lineq$lys(45:50)
#'
#' @export
LinearEquate=function(xscores,yscores){

	muSx=mean(xscores)
	muSy=mean(yscores)
	sigSx=sd(xscores)
	sigSy=sd(yscores)

	lys=function(scores){(sigSy/sigSx)*(scores-muSx)+muSy}

	return(list(
		lys=lys
		,yx=lys(sort(unique(xscores)))))
	}


#' Apply chained linear equating
#'
#' Chained linear equating via an anchor test is done by successively applying the \emph{LinearEquate} function between anchors and main tests.
#'
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' 
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{chainedFunc}{A function that translates any vector of scores on form X into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame combining the sorted unique scores on form X in the data and their equated values on form Y.}
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
#' lin_eq=ChainedLinearEquate(data.frame(x=scoresX1,a=scoresA1),data.frame(y=scoresY2,a=scoresA2))
#' lin_eq
#' 
#' @keywords KernEqWPS
#' @export
ChainedLinearEquate=function(dx,dy){
	Lxa=LinearEquate(dx$x,dx$a)
	Lay=LinearEquate(dy$a,dy$y)
	scoresX=sort(unique(dx$x))

	chainedFunc=function(x){Lay$lys(Lxa$lys(x))}

	return(list(
		chainedFunc=chainedFunc,
		EqTable=data.frame(x=scoresX,equiyx=chainedFunc(scoresX))
		))
	}
