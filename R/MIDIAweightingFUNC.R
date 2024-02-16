#' Function to find weights that can be applied to a data set so that means of columns match assigned targets.
#'
#' Method is from Haberman, 1984.
#'
#' @param datamatrix Numeric matrix (or something that can be quickly converted into one e.g. a data frame).
#' @param targets Vector indicating the desired mean value for each column after weighting.
#' @param tol Numeric input relating to required proximity of weighted means to target means in order for algorithm to be considered converged (default 1e-6).
#' @param maxiter Numeric input denoting the maximum number of iterations that will be applied to find the weights (default 100).
#' @param starttheta Starting parameter vector to begin iterative process (by default a vector of zeros).
#'
#' @return The function returns a list with the following elements:
#' \describe{
#'   \item{wts}{A vector of weights.}
#'   \item{iter}{The number of iterations taken to converge. If this is equal to maxiter then convergenece may not have been achieved.}
#'   \item{FinalWeightedMeans}{The weighted means on each column in the data matrix.}
#'   \item{theta}{Parameters used for weighting.}
#' }
#'
#' @references
#' 
#' Haberman, S. J. (1984). Adjustment by minimum discriminant information. 
#' \emph{The Annals of Statistics, 12}(3), 971-988.
#' 
#' @examples
#' #try weighting one population (with a logistic distribution)
#' #to match E(x), E(x^2),E(x^3),E(x^4), E(x^5), E(x^6) and E(x^7) of another (with a normal distribution)
#' x1=rnorm(50000,100,30)
#' summary(x1)
#' x2=rlogis(50000,120,25)
#' summary(x2)
#' MDIAx2=MDIAweights(cbind(x2,x2^2,x2^3,x2^4,x2^5,x2^6,x2^7),c(mean(x1),mean(x1^2),mean(x1^3),mean(x1^4),mean(x1^5),mean(x1^6),mean(x1^7)))
#' MDIAx2$iter
#' sum(MDIAx2$wts)
#' summary(MDIAx2$wts)
#' summary(MDIAx2$wts/mean(MDIAx2$wts))
#' c(mean(x1),sum(MDIAx2$wts*x2))
#' c(mean(x1^2),sum(MDIAx2$wts*x2^2))
#' c(mean(x1^3),sum(MDIAx2$wts*x2^3))
#' c(mean(x1^4),sum(MDIAx2$wts*x2^4))
#' c(mean(x1^5),sum(MDIAx2$wts*x2^5))
#' c(mean(x1^6),sum(MDIAx2$wts*x2^6))
#' c(mean(x1^7),sum(MDIAx2$wts*x2^7))
#' #see if it worked visually
#' dat1=data.frame(x=c(x1,x2)
#' 	,pop=as.factor(c(rep(1,length(x1)),rep(2,length(x2))))
#' 	,wts=c(rep(1/length(x1),length(x1)),MDIAx2$wts))
#' library(ggplot2)
#' ggplot(data=dat1,aes(x=x,fill=pop))+geom_density(alpha=0.5)#unweighted
#' ggplot(data=dat1,aes(x=x,fill=pop,weight=wts))+geom_density(alpha=0.5)#weighted
#'
#' @keywords KernEqWPS
#' @export
MDIAweights=function(datamatrix,targets,tol=1e-6,maxiter=100,starttheta=NA){

T=as.matrix(datamatrix)
maxs=apply(abs(T),2,max)
t=as.vector(targets)
#divide initial version by fifth of maximum to avoid computational singularity
T=t(t(T)/as.vector(maxs))
t=t/as.vector(maxs)

dimT=dim(T)
nt=length(t)
if(dim(T)[2]!=nt){
  print("Number of targets MUST equal number of variables in input data")
  return(NULL)
}

#initial theta
theta=starttheta*as.vector(maxs)
if(is.na(theta)[1]){theta=rep(0,length(t))}
#calculate current weights given theta coefficient
Ttheta=T%*%theta
Ttheta=Ttheta-max(Ttheta)
w=exp(Ttheta)/sum(exp(Ttheta))
#calculate current weighted mean
m=colSums(T*as.vector(w))

stop=0
if(max(abs(t-m))<tol){stop=1}
iter=0
while(stop==0){
  iter=iter+1
  #calculate current covariance matrix
  #non-looping way
  temp1=t(t(T)-m)
  SIGMA=t(temp1*as.vector(w))%*%(temp1)
  #update theta
  #theta=theta+solve(SIGMA)%*%(t-m)
  theta=theta+MASS::ginv(SIGMA)%*%(t-m)
  #calculate new weights given theta coefficient
  Ttheta=T%*%theta
  Ttheta=Ttheta-max(Ttheta)
  w=exp(Ttheta)/sum(exp(Ttheta))
  #calculate new weighted mean
  m=colSums(T*as.vector(w))
  #print(m)
  #print(cbind(t,m))
  if(max(abs(t-m))<tol|iter>=maxiter){stop=1}
}

#IF METHOD HAS FAILED THEN RESTART AND REFIT WITH A SMALLER STEP SIZE AT EACH ITERATION
#...that is...theta=theta+0.3*MASS::ginv(SIGMA)%*%(t-m)
#...rather than...theta=theta+MASS::ginv(SIGMA)%*%(t-m)

if (iter>=maxiter){
#initial theta
theta=starttheta*as.vector(maxs)
if(is.na(theta)[1]){theta=rep(0,length(t))}
#calculate current weights given theta coefficient
Ttheta=T%*%theta
Ttheta=Ttheta-max(Ttheta)
w=exp(Ttheta)/sum(exp(Ttheta))
#calculate current weighted mean
m=colSums(T*as.vector(w))

stop=0
if(max(abs(t-m))<tol){stop=1}
iter=0
while(stop==0){
  iter=iter+1
  #calculate current covariance matrix
  #non-looping way
  temp1=t(t(T)-m)
  SIGMA=t(temp1*as.vector(w))%*%(temp1)
  #update theta (smaller step size)
  theta=theta+0.3*MASS::ginv(SIGMA)%*%(t-m)
  #calculate new weights given theta coefficient
  Ttheta=T%*%theta
  Ttheta=Ttheta-max(Ttheta)
  w=exp(Ttheta)/sum(exp(Ttheta))
  #calculate new weighted mean
  m=colSums(T*as.vector(w))
  #print(m)
  #print(cbind(t,m))
  if(max(abs(t-m))<tol|iter>=maxiter){stop=1}
}
}

#If convergence still questionable then give an error rather than silently give bad results
if (iter>=maxiter){stop("Error: convergence not achieved. Consider using a lower order of matching.")}

#print(cbind(t,m))

return(list(wts=w,iter=iter,FinalWeightedMeans=maxs*m,theta=theta/as.vector(maxs)))
}

#' Function to find weights to apply to a number of groups in a data set to ensure equivalence on some variables.
#'
#' Method is from Haberman, 1984. 
#'
#' @param formula A formula indicating the variables on which matching should be achieved after weighting. For example, ~x+I(x^2)+y+I(y^2)+I(y^3) would achieve matched values on the mean value of x, the squared mean value of x (essentially equivalent to the variance), and (essentially) the mean, standard deviation and skewness of y.
#' @param data, A data frame.
#' @param groups A vector indicating the group to which each row of the data belongs.
#' @param targetgroup Single value indicating which group defines the population to which all other groups should be weighted. By default all groups are weighted to match the overall mean.
#'
#' @return A vector of weights.
#' @references
#' 
#' Haberman, S. J. (1984). Adjustment by minimum discriminant information. 
#' \emph{The Annals of Statistics, 12}(3), 971-988.
#' 
#' @keywords KernEqWPS
#' @export

MDIAwtGroups=function(formula,data,groups,targetgroup=NA){

ugps=sort(unique(groups))#unique groups
termsDat=sapply(attr(terms(formula),"term.labels"),function(TEMPTERM) with(data,eval(parse(text=TEMPTERM))))
TD2=cbind(rep(1,dim(termsDat)[1]),termsDat)#create a matrix of terms including a column of 1s so can use logistic regression to set starting values

if(!is.na(targetgroup)){targs=colMeans(termsDat[groups==targetgroup,])
		outTD2=TD2[groups==targetgroup,]}
#target group that we are weigthing to is not specified assume it is the 
if(is.na(targetgroup)){targs=colMeans(termsDat)
		outTD2=TD2
		targetgroup=min(groups)-1#(replace target group with unused groups value to help next bit of code)
		}

wts=rep(NA,dim(data)[1])
for(iiiz in 1:length(ugps)){
	if(ugps[iiiz]==targetgroup){wts[groups==ugps[iiiz]]=1/sum(groups==targetgroup)}
	if(ugps[iiiz]!=targetgroup){
		TD2a=TD2[groups==ugps[iiiz],]
		TD3=rbind(TD2a,outTD2)
		outTD3=c(rep(0,dim(TD2a)[1]),rep(1,dim(outTD2)[1]))
		#logistic regression to get starting values for MDIA
		glm1a=glm.fit(outTD3,x=TD3,family=binomial())
		wts[groups==ugps[iiiz]]=MDIAweights(termsDat[groups==ugps[iiiz],],targs,starttheta=glm1a$coefficients[-1])$wts
		}
	}
return(wts)}

