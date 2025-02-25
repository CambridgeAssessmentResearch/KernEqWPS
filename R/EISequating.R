#' Apply equating based on equivalent expected score on an anchor test
#'
#' At heart, this is the same method as the similar items method suggested by Bramley (2018).
#' However, the name has been changed to reflect the fact it can be used with an actual "anchor" test 
#' rather than simply when we only have "similar" items. Expectes scores on the
#' anchor test for each raw total test score are derived from a loglinear model with a single interaction.
#' Linear interpolation is used to identify expected item scores for non-integer raw scores.
#' Equated scores on each test form are the raw scores that give the same expected anchor test score.
#' For this function to work THE ANCHOR TEST MUST BE INTERNAL TO THE TESTS BEING EQUATED.
#'
#' This function assumes that only integer scores are available on each test form.
#' 
#' The loglinear models behind this approach are fitted to ensure a perfect fit to the marginal distributions of both raw test and anchor test
#' scores. Furthermore, the fitted model will replicate the correlation between raw test and anchor test scores.
#' 
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test.
#' @param maxX Maximum score on form X.
#' @param maxY Maximum score on form Y.
#' @param maxA Maximum score on anchor.
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{EqTable}{A data frame showing the equivalent score on form Y for every integer score between 0 and maxX on form X.}
#' }
#'
#' 
#' @references
#' Bramley, T. (2018, November). \emph{Evaluating the ‘similar items method’ for standard maintaining.}
#' Paper presented at the 19th annual conference of the Association for Educational Assessment in Europe, Arnhem-Nijmegen, The Netherlands. 
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
#' eis_eq=EISEquate(data.frame(x=scoresX1,a=scoresA1),data.frame(y=scoresY2,a=scoresA2),30,30,5)
#' eis_eq
#'
#' @keywords KernEqWPS
#' @export
EISEquate=function(dx, dy, maxX, maxY, maxA){

  #fits a loglinear model for group 1
  count_dat_x <- merge(data.frame(x=0:maxX),data.frame(a=0:maxA))
  count_dat_x$n <- sapply(1:nrow(count_dat_x),function(i) sum((dx$x==count_dat_x$x[i]) & (dx$a==count_dat_x$a[i]), na.rm = T))
  #total score on rest
  count_dat_x$x_minus_a <- count_dat_x$x-count_dat_x$a
  
  #controls if the anchor is internal
  #remove negative cases (will never happen)
  count_dat_x <- count_dat_x[count_dat_x$x_minus_a>=0,]
  #also remove cases where person gets more marks than viable given the anchor
  count_dat_x <- count_dat_x[count_dat_x$x<=(maxX-(maxA-count_dat_x$a)),]
  count_dat_x
  
  #fit a loglinear model (saturated in marginal distribution, linear interaction for relationship)
  loglin_mod_x <- glm(n ~ as.factor(x) + as.factor(a) + I(x*a), data=count_dat_x, family="poisson")
  
  #get fitted values
  count_dat_x$fitted <- loglin_mod_x$fitted 
  count_dat_x$wt.a <- count_dat_x$a*count_dat_x$fitted
  agg_x <- aggregate(count_dat_x[,c("wt.a","fitted")],by=list(x=count_dat_x$x),sum)
  agg_x$Ea_x <- agg_x$wt.a/agg_x$fitted
  
  #same for version y
  #same for version y
  #same for version y
  #same for version y
  #same for version y
  
  #fits a loglinear model for group 1
  count_dat_y <- merge(data.frame(y=0:maxY),data.frame(a=0:maxA))
  count_dat_y$n <- sapply(1:nrow(count_dat_y),function(i) sum((dy$y==count_dat_y$y[i]) & (dy$a==count_dat_y$a[i]), na.rm = T))
  #total score on rest
  count_dat_y$y_minus_a <- count_dat_y$y-count_dat_y$a
  
  #controls if the anchor is internal
  #remove negative cases (will never happen)
  count_dat_y <- count_dat_y[count_dat_y$y_minus_a>=0,]
  #also remove cases where person gets more marks than viable given the anchor
  count_dat_y <- count_dat_y[count_dat_y$y<=(maxY-(maxA-count_dat_y$a)),]
  count_dat_y
  
  #fit a loglinear model (saturated in marginal distribution, linear interaction for relationship)
  loglin_mod_y <- glm(n ~ as.factor(y) + as.factor(a) + I(y*a), data=count_dat_y, family="poisson")
  
  #get fitted values
  count_dat_y$fitted <- loglin_mod_y$fitted 
  count_dat_y$wt.a <- count_dat_y$a*count_dat_y$fitted
  agg_y <- aggregate(count_dat_y[,c("wt.a","fitted")],by=list(y=count_dat_y$y),sum)
  agg_y$Ea_y <- agg_y$wt.a/agg_y$fitted
  
  #make version y into a function (use interpolation)
  eay=approxfun(x=agg_y$y,y=agg_y$Ea_y)
  eay_v2=function(y,targ){eay(y)-targ}

  #get a list of expected scores on anchor associated with each integer score on form X (excluding 0 and the maximum)
  targs=agg_x$Ea_x[2:maxX]

  #numerical search to find scores on form Y that lead to the same expected scores on the anchor  
  eyx=sapply(targs, function(i) uniroot(eay_v2,c(0,maxY),targ=i)$root)
  eyx=c(0,eyx,maxY)
  
  EqTable = data.frame(x = 0:maxX, equiyx = eyx)
  return(EqTable)
  }
