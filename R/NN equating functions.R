#' Create 51x11 crosstab for test by anchor scores
#' 
#' Regardless of the maximum scores available on each test form, it is necessary to reduce the relationship into
#' and matrix of 51 rows and 11 columns in order to feed this into the experimental equating methods using neural networks.
#' This step is completed by this function. If the test is out of 50 and the anchor out of 10 then this function will give very
#' similar results to table(dx$x,dx$a). However, if the total differ from this the data will be squeezed into the correct shape.
#' In particular, if the anchor test has a maximum greater than 10 then the scores will be replaced by categorizations into 11th (0-10)
#' for the purposes of the crosstab. 
#' If the maximum test score differs from 50 then the score range will be broken into 51 equally sized sections
#' and the crosstab will denote the extent of the continuized distribution within each of these ranges.
#' (Continuization is done using the same principles as kernel equating and with a bandwidth of 0.4).
#' 
#' Note that for this function it is assumed that all scores are integers.
#' 
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test. 
#' @param maxx Maximum score available on form X (calculated from the data by default).
#' @param maxa Maximum score available on anchor test (calculated from the data by default).
#' @examples
#' #example simulate a test and an anchor test
#' n1=1000
#' t1=rnorm(n1,0.5,1)
#' dx=data.frame(x=round(pmin(50,pmax(0,25+10*(0.9*t1+rnorm(n1,0,sqrt(1-0.9^2))))))
#' 	,a=round(pmin(10,pmax(0,5+2*(0.7*t1+rnorm(n1,0,sqrt(1-0.7^2)))))))
#' crosstab1=Crosstab51x11(dx)
#' dim(crosstab1)
#' sum(crosstab1)
#' #display matrix as a heatmap
#' image(1:ncol(crosstab1), 1:nrow(crosstab1),t(crosstab1), col = terrain.colors(60), axes = FALSE)
#' axis(1, 1:ncol(crosstab1),0:(ncol(crosstab1)-1))
#' axis(2, 1:nrow(crosstab1),0:(nrow(crosstab1)-1))
#' @keywords KernEqWPS
#' @export
Crosstab51x11=function(dx,maxx=NA,maxa=NA){

if(is.na(maxa)){maxa=max(dx$a)}
if(is.na(maxx)){maxx=max(dx$x)}
if (maxa>10){dx$a=findInterval(dx$a,quantile(dx$a,seq(1/11,10/11,1/11)))}

#if maxx is just 50 then we can simply crosstabulate
if(maxx==50){
	crosstabfunc=function(scorea,dx,maxx){tabulate(match(dx$x[dx$a==scorea],0:maxx),(maxx+1))}
	return(simplify2array(lapply(0:10,crosstabfunc,dx=dx,maxx=maxx)))}

#SUBROUTINE - WORK OUT "NUMBER" WITH EACH SCORE "OUT OF RESCALED 50" USING KERNEL METHODS FOR A PARTICULAR ANCHOR SCORE
#INPUTS ARE dx (A DATA FRAME WITH COLUMNS X AND A), score a (which anchor score)
#ONLY WORKS WITH WHOLE NUMBER SCORES AND ASSUMES MINIMUM AVAILABLE SCORE IS ZERO
WhatProportion1=function(scorea,dx,maxx=NA){
	if(is.na(maxx)){maxx=max(dx$x)}
	if(length(dx$a[dx$a==scorea])==0){return(rep(0,51))}
	fX=tabulate(match(dx$x[dx$a==scorea],0:maxx),(maxx+1))
	msX=DistToMuSig(0:maxx,fX)
	if(is.na(msX$sigma)){msX$sigma=0.1}#if only one observation just set a really small sd
	integrateFunc=function(scorex50){
		lower=((scorex50-0.5)*maxx/50)
		if (scorex50==0){lower=-10}
		upper=((scorex50+0.5)*maxx/50)
		if (scorex50==50){upper=maxx+10}
		integrate(function(tt) fhx(tt,fX=fX,scores=0:maxx,mu=msX$mu,sigma=msX$sigma,h=0.4),lower=lower,upper=upper)$value}
	sapply(0:50,integrateFunc)
	}

simplify2array(lapply(0:10,WhatProportion1,dx=dx,maxx=maxx))}


#' Experimental function for getting percentiles using a standard neural network
#' 
#' This function takes a data frame of test form and anchor scores
#' and estimates the values of the percentiles (1st-99th) for given change in the distribution of anchor scores.
#' 
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test. 
#' @param anchortargettable Table giving distribution of anchor test scores in the target population.
#' @param maxx Maximum score available on form X (calculated from the data by default).
#' @param maxa Maximum score available on anchor test (calculated from the data by default).
#' @param WeightsList A list of neural network parameters used in calculations. Changing this from the default value is not recommended.
#' 
#' @examples
#' #example (compare real and estimated percentiles within a fixed population)
#' n1=1000
#' t1=rnorm(n1,0.5,1)
#' dx=data.frame(x=round(pmin(100,pmax(0,50+20*(0.9*t1+rnorm(n1,0,sqrt(1-0.9^2))))))
#' 	,a=round(pmin(10,pmax(0,5+2*(0.7*t1+rnorm(n1,0,sqrt(1-0.7^2)))))))
#' 
#' percNN=PredPercentileNN(dx,table(dx$a),maxx=100)
#' usualperc=as.vector(quantile(dx$x,seq(0.01,0.99,0.01)))
#' plot(1:99,usualperc,type='l',xlab="Percentile",ylab="Value")
#' lines(1:99,percNN,lty=2)
#' 
#' @keywords KernEqWPS
#' @export
PredPercentileNN=function(dx,anchortargettable=NA,maxx=NA,maxa=NA,WeightsList=EquateNNWeights){
	
	if(any(is.na(dx))){print("input data frame cannot contain missing values")
		return(NULL)}

	if(is.na(anchortargettable[1])){anchortargettable=table(dx$a)}

	if(is.na(maxa)){maxa=max(dx$a)}
	if(is.na(maxx)){maxx=max(dx$x)}
	anctargvec=rep(as.numeric(names(anchortargettable)),anchortargettable)
	nsamp=max(100,min(500,dim(dx)[1]))#training only set up for sample sizes between 100 and 500
	ntarg=max(100,min(500,sum(anchortargettable)))#training only set up for sample sizes between 100 and 500

	if (maxa>10){quanta=quantile(dx$a,seq(1/11,10/11,1/11))
		dx$a=findInterval(dx$a,quanta)
		anctargvec=findInterval(anctargvec,quanta)
		}

	TargsProps11=tabulate(match(anctargvec,0:10),11)
	TargsProps11=TargsProps11/sum(TargsProps11)
	crosstab1=Crosstab51x11(dx)
	crosstab1=crosstab1/sum(crosstab1)
	crosstab1=as.vector((crosstab1))
	X1=c(nsamp,ntarg,TargsProps11,crosstab1)
	X1=matrix(X1,nrow=1)
	layer1=sweep(X1%*%WeightsList$Wfc1,2,WeightsList$bfc1,"+")
	layer1=0.5*layer1*(sign(layer1)+1)
	layer2=sweep(layer1%*%WeightsList$Wfc2,2,WeightsList$bfc2,"+")
	layer2=0.5*layer2*(sign(layer2)+1)
	layer3=sweep(layer2%*%WeightsList$Wfc3,2,WeightsList$bfc3,"+")
	layer3=0.5*layer3*(sign(layer3)+1)
	layer4=sweep(layer3%*%WeightsList$Wfc4,2,WeightsList$bfc4,"+")
	layer4=0.5*layer4*(sign(layer4)+1)
	layer5=sweep(layer4%*%WeightsList$Wfc5,2,WeightsList$bfc5,"+")
	layer5=0.5*layer5*(sign(layer5)+1)
	out=sweep(layer5%*%WeightsList$Wfc6,2,WeightsList$bfc6,"+")
	return((maxx/50)*sort(as.vector(out)))
	}

#' Experimental technique to equate using a pre-trained neural network
#' 
#' This function is designed to work in the NEAT equating design.
#' It assumes all scores are integers and that it is an EXTERNAL anchor test.
#' 
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test. 
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test. 
#' @param anchortargettable Table giving distribution of anchor test scores in the target population.
#' @param maxx Maximum score available on form X (calculated from the data by default).
#' @param maxy Maximum score available on form Y (calculated from the data by default).
#' @param maxa Maximum score available on anchor test (calculated from the data by default).
#' @param WeightsList A list of neural network parameters used in calculations. Changing this from the default value is not recommended.
#'
#' @return The function returns a list with two elements.
#' \describe{
#'   \item{CNNEqFunc}{A function that translates any vector of scores on form X (from 0-maxx) into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame combining scores from 0-maxx in the data and their equated values on form Y.}}
#' 
#' @examples
#' #example where we simulate scores on two parallel tests in two populations using classical test theory
#' n1=400
#' n2=300
#' t1=rnorm(n1,0.5,1)
#' t2=rnorm(n2,0,1)
#' dx1=data.frame(x=round(pmin(100,pmax(0,50+20*(0.9*t1+rnorm(n1,0,sqrt(1-0.9^2)))))),
#' 	a=round(pmin(20,pmax(0,10+4*(0.7*t1+rnorm(n1,0,sqrt(1-0.7^2)))))))
#' dy1=data.frame(y=round(pmin(100,pmax(0,50+20*(0.9*t2+rnorm(n2,0,sqrt(1-0.9^2))))))
#' 	,a=round(pmin(20,pmax(0,10+4*(0.7*t2+rnorm(n2,0,sqrt(1-0.7^2)))))))
#' NNeq1=EquateNN(dx1,dy1,table(dy1$a))
#' chainedeq1=KernelChainedEquate(dx1, dy1)$EqTable
#' pseeq1=PSEObservedEquate(dx1, dy1, target = "y")$EqTable
#' plot(chainedeq1$x,chainedeq1$equiyx,type='l')#chained equating function
#' lines(pseeq1$x,pseeq1$equiyx,lty=3,col="red")#PSE estimated equating function
#' lines(NNeq1$EqTable$x,NNeq1$EqTable$yx,lty=2)#NN estimated equating function
#' lines(0:100,0:100,col="blue",lty=3,lwd=4)#true equating function (identity as set up to be parallel tests)
#' 
#' #example using some real data (but no criterion equate to compare to)
#' dx2=data.frame(x=rowSums(mathsdata[1:250,1:35]),a=rowSums(mathsdata[1:250,41:50]))
#' summary(dx2)
#' dy2=data.frame(y=rowSums(mathsdata[251:500,51:90]),a=rowSums(mathsdata[251:500,41:50]))
#' summary(dy2)
#' NNeq2=EquateNN(dx2,dy2,table(dy1$a))
#' chainedeq=KernelChainedEquate(dx2, dy2)$EqTable
#' plot(chainedeq$x,chainedeq$equiyx,type='l')#chained equating function
#' lines(NNeq2$EqTable$x,NNeq2$EqTable$yx,lty=2)#NN estimated equating function
#' 
#' @export
EquateNN=function(dx,dy,anchortargettable,maxx=NA,maxy=NA,maxa=NA,WeightsList=EquateNNWeights){

if(is.na(maxa)){maxa=max(c(dx$a,dy$a))}
anctargvec=rep(as.numeric(names(anchortargettable)),anchortargettable)

if (maxa>10){quanta=quantile(c(dx$a,dy$a),seq(1/11,10/11,1/11))
		dx$a=findInterval(dx$a,quanta)
		dy$a=findInterval(dy$a,quanta)
		anctargvec=findInterval(anctargvec,quanta)
		}
if(is.na(maxx)){maxx=max(dx$x)}
if(is.na(maxy)){maxy=max(dy$y)}

percX=PredPercentileNN(dx,anchortargettable=table(anctargvec),maxx=maxx,maxa=10,WeightsList=WeightsList)
percY=PredPercentileNN(dx=data.frame(x=dy$y,a=dy$a),anchortargettable=table(anctargvec),maxx=maxy,maxa=10,WeightsList=WeightsList)

x1=as.numeric(percX)
x1=c(min(0,min(x1)),x1,max(maxx,max(x1)))
y1=as.numeric(percY)
y1=c(min(0,min(y1)),y1,max(maxy,max(y1)))
interp=approxfun(x1,y1)
EqRes=interp(0:maxx)
return(list(
	NNEqFunc=interp
	,EqTable=data.frame(x=0:maxx,yx=EqRes))
	)
}


#' Experimental function for getting percentiles using an approximation to a trained convolutional neural network
#' 
#' This function takes a data frame of test form and anchor scores
#' and estimates the values of the percentiles (1st-99th) for given change in the distribution of anchor scores.
#' See Benton (2017) for more details.
#' 
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test. 
#' @param anchortargettable Table giving distribution of anchor test scores in the target population.
#' @param maxx Maximum score available on form X (calculated from the data by default).
#' @param maxa Maximum score available on anchor test (calculated from the data by default).
#' @param WeightsList A list of neural network parameters used in calculations. Changing this from the default value is not recommended.
#' 
#' @references
#' 
#' Benton, T. (2017). Can AI learn to equate?, 
#' \emph{presented at the International Meeting of the Psychometric Society, Zurich, 2017}. Cambridge, UK: Cambridge Assessment.
#' 
#' @examples
#' #example (compare real and estimated percentiles within a fixed population)
#' n1=1000
#' t1=rnorm(n1,0.5,1)
#' dx=data.frame(x=round(pmin(100,pmax(0,50+20*(0.9*t1+rnorm(n1,0,sqrt(1-0.9^2))))))
#' 	,a=round(pmin(10,pmax(0,5+2*(0.7*t1+rnorm(n1,0,sqrt(1-0.7^2)))))))
#' 
#' percCNN=PredPercentileCNN(dx,table(dx$a),maxx=100)
#' usualperc=as.vector(quantile(dx$x,seq(0.01,0.99,0.01)))
#' plot(1:99,usualperc,type='l',xlab="Percentile",ylab="Value")
#' lines(1:99,percCNN,lty=2)
#' 
#' @keywords KernEqWPS
#' @export
PredPercentileCNN=function(dx,anchortargettable=NA,maxx=NA,maxa=NA,WeightsList=ApproxCNNWeights){
	
	if(any(is.na(dx))){print("input data frame cannot contain missing values")
		return(NULL)}

	if(is.na(anchortargettable[1])){anchortargettable=table(dx$a)}

	if(is.na(maxa)){maxa=max(dx$a)}
	if(is.na(maxx)){maxx=max(dx$x)}
	anctargvec=rep(as.numeric(names(anchortargettable)),anchortargettable)
	nsamp=max(100,min(500,dim(dx)[1]))#training only set up for sample sizes between 100 and 500
	ntarg=max(100,min(500,sum(anchortargettable)))#training only set up for sample sizes between 100 and 500

	if (maxa>10){quanta=quantile(dx$a,seq(1/11,10/11,1/11))
		dx$a=findInterval(dx$a,quanta)
		anctargvec=findInterval(anctargvec,quanta)
		}

	TargsProps11=tabulate(match(anctargvec,0:10),11)
	TargsProps11=TargsProps11/sum(TargsProps11)
	crosstab1=Crosstab51x11(dx)
	crosstab1=crosstab1/sum(crosstab1)
	crosstab1=as.vector((crosstab1))

	X0=matrix(crosstab1,nrow=1)
	layer0=X0%*%WeightsList$Wfc0

	X1=c(nsamp,ntarg,TargsProps11,layer0)
	X1=matrix(X1,nrow=1)
	layer1=sweep(X1%*%WeightsList$Wfc1,2,WeightsList$bfc1,"+")
	layer1=0.5*layer1*(sign(layer1)+1)
	out=sweep(layer1%*%WeightsList$Wfc2,2,WeightsList$bfc2,"+")
	return((maxx/50)*sort(as.vector(out)))
	}

#' Experimental technique to equate using an approximation to a pre-trained convolutional neural network
#' 
#' This function is designed to work in the NEAT equating design.
#' It assumes all scores are integers and that it is an EXTERNAL anchor test.
#' 
#' @param dx Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test. 
#' @param dy Data frame with variables "y" and "a" representing scores for individual candidates on form Y and on the anchor test. 
#' @param anchortargettable Table giving distribution of anchor test scores in the target population.
#' @param maxx Maximum score available on form X (calculated from the data by default).
#' @param maxy Maximum score available on form Y (calculated from the data by default).
#' @param maxa Maximum score available on anchor test (calculated from the data by default).
#' @param WeightsList A list of neural network parameters used in calculations. Changing this from the default value is not recommended.
#'
#' @return The function returns a list with two elements.
#' \describe{
#'   \item{CNNEqFunc}{A function that translates any vector of scores on form X (from 0-maxx) into equivalent scores on form Y.}
#'   \item{EqTable}{A data frame combining scores from 0-maxx in the data and their equated values on form Y.}}
#' 
#' @examples
#' #example where we simulate scores on two parallel tests in two populations using classical test theory
#' n1=400
#' n2=300
#' t1=rnorm(n1,0.5,1)
#' t2=rnorm(n2,0,1)
#' dx1=data.frame(x=round(pmin(100,pmax(0,50+20*(0.9*t1+rnorm(n1,0,sqrt(1-0.9^2)))))),
#' 	a=round(pmin(20,pmax(0,10+4*(0.7*t1+rnorm(n1,0,sqrt(1-0.7^2)))))))
#' dy1=data.frame(y=round(pmin(100,pmax(0,50+20*(0.9*t2+rnorm(n2,0,sqrt(1-0.9^2))))))
#' 	,a=round(pmin(20,pmax(0,10+4*(0.7*t2+rnorm(n2,0,sqrt(1-0.7^2)))))))
#' CNNeq1=EquateCNN(dx1,dy1,table(dy1$a))
#' chainedeq1=KernelChainedEquate(dx1, dy1)$EqTable
#' pseeq1=PSEObservedEquate(dx1, dy1, target = "y")$EqTable
#' plot(chainedeq1$x,chainedeq1$equiyx,type='l')#chained equating function
#' lines(pseeq1$x,pseeq1$equiyx,lty=3,col="red")#PSE estimated equating function
#' lines(CNNeq1$EqTable$x,CNNeq1$EqTable$yx,lty=2)#CNN estimated equating function
#' lines(0:100,0:100,col="blue",lty=3,lwd=4)#true equating function (identity as set up to be parallel tests)
#' 
#' #example using some real data (but no criterion equate to compare to)
#' dx2=data.frame(x=rowSums(mathsdata[1:250,1:35]),a=rowSums(mathsdata[1:250,41:50]))
#' summary(dx2)
#' dy2=data.frame(y=rowSums(mathsdata[251:500,51:90]),a=rowSums(mathsdata[251:500,41:50]))
#' summary(dy2)
#' CNNeq2=EquateCNN(dx2,dy2,table(dy1$a))
#' chainedeq=KernelChainedEquate(dx2, dy2)$EqTable
#' plot(chainedeq$x,chainedeq$equiyx,type='l')#chained equating function
#' lines(CNNeq2$EqTable$x,CNNeq2$EqTable$yx,lty=2)#CNN estimated equating function
#' 
#' @export
EquateCNN=function(dx,dy,anchortargettable,maxx=NA,maxy=NA,maxa=NA,WeightsList=ApproxCNNWeights){

if(is.na(maxa)){maxa=max(c(dx$a,dy$a))}
anctargvec=rep(as.numeric(names(anchortargettable)),anchortargettable)

if (maxa>10){quanta=quantile(c(dx$a,dy$a),seq(1/11,10/11,1/11))
		dx$a=findInterval(dx$a,quanta)
		dy$a=findInterval(dy$a,quanta)
		anctargvec=findInterval(anctargvec,quanta)
		}
if(is.na(maxx)){maxx=max(dx$x)}
if(is.na(maxy)){maxy=max(dy$y)}

percX=PredPercentileCNN(dx,anchortargettable=table(anctargvec),maxx=maxx,maxa=10,WeightsList=WeightsList)
percY=PredPercentileCNN(dx=data.frame(x=dy$y,a=dy$a),anchortargettable=table(anctargvec),maxx=maxy,maxa=10,WeightsList=WeightsList)

x1=as.numeric(percX)
x1=c(min(0,min(x1)),x1,max(maxx,max(x1)))
y1=as.numeric(percY)
y1=c(min(0,min(y1)),y1,max(maxy,max(y1)))
interp=approxfun(x1,y1)
EqRes=interp(0:maxx)
return(list(
	CNNEqFunc=interp
	,EqTable=data.frame(x=0:maxx,yx=EqRes))
	)
}


