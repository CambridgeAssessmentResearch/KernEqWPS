###EXAMPLE OF WHY THE LOG-LINEAR SMOOTHING STEP MAY POSSIBLY BE UNECESSARY
library(KernEqWPS)#without pre-smoothing
library(kequate)#with pre-smoothing
#use IRT to simulate a large data set (20,000 cases)
library(mirt)
mod1=mirt(mathsdata,1)
thetas=as.matrix(rnorm(20000,0,1))
sim1=simdata(model=mod1,Theta=thetas)

###WILL COMPARE MULTIPLE EQUATING METHODS
###COULD ITERATE THE NEXT BIT OF THE PROCESS HERE
#for (iterations in 1:100){

#break larger simulated maths data into 2 pseudo-tests
dim(sim1)
itesX=SampleItemsToHitTarget(sim1,50)#randomly select form X to have maximum of 50
itesY=SampleItemsToHitTarget(sim1[,-itesX],50)#randomly select form Y to have maximum of 50
#make form scores for equivalent groups
XP=rowSums(sim1[,itesX])
YP=rowSums(sim1[,-itesX][,itesY])
length(YP)
#Equate form X to form Y with low bandwidth to get a criterion equate
eqCRIT=KernelEquateFromScoresEG(XP,YP,hX=0.3,hY=0.3)
plot(0:50,eqCRIT$yxFunc(0:50)-(0:50),type='l',xlab="x",ylab="y(x)-x",ylim=c(-15,15))

#TAKE A SMALLER SAMPLE (300 CASES)
samp1=sample(1:dim(sim1)[1],300)
Xsamp=XP[samp1]
Ysamp=YP[samp1]
#Equating without pre-smoothing and with plug-in bandwidth
eqWPS=KernelEquateFromScoresEG(Xsamp,Ysamp)
lines(0:50,eqWPS$yxFunc(0:50)-(0:50),lty=2,col="red")

#Equating without pre-smoothing and bandwidth found via cross-validation
hX1=FindBestBandwidth(sort(unique(Xsamp)),table(Xsamp))
hY1=FindBestBandwidth(sort(unique(Ysamp)),table(Ysamp))
eqWPS2=KernelEquateFromScoresEG(Xsamp,Ysamp,hX=hX1,hY=hY1)
lines(0:50,eqWPS2$yxFunc(0:50)-(0:50),lty=2,col="purple")


#Equating with pre-smoothing
freqX <- kefreq(Xsamp,0:50)
freqY <- kefreq(Ysamp,0:50)
# Fit the log-linear models with increasing order until adding more terms leads to two successive increase in AIC
form1="frequency ~ I(X)"
glmX <-glm(as.formula(form1), family='poisson', data=freqX, x=TRUE)
stop=0
order=2
while(stop<1){
	form1=paste(form1,"+I(X^",order,")")
	newglm<-glm(as.formula(form1), family='poisson', data=freqX, x=TRUE)
	if(summary(newglm)$aic>=summary(glmX)$aic){stop=stop+0.5}
	if(summary(newglm)$aic<summary(glmX)$aic){
		glmX=newglm
		stop=0
		}
	order=order+1}

form1="frequency ~ I(X)"
glmY <-glm(as.formula(form1), family='poisson', data=freqY, x=TRUE)
stop=0
order=2
while(stop<1){
	form1=paste(form1,"+I(X^",order,")")
	newglm<-glm(as.formula(form1), family='poisson', data=freqY, x=TRUE)
	if(summary(newglm)$aic>=summary(glmY)$aic){stop=stop+0.5}
	if(summary(newglm)$aic<summary(glmY)$aic){
		glmY=newglm
		stop=0
		}
	order=order+1}

eqKEQUATE <- kequate("EG", 0:50, 0:50, glmX, glmY)
lines(0:50,eqKEQUATE@equating[,1]-(0:50),lty=3,col="blue")

errorWPS=eqWPS$yxFunc(0:50)-eqCRIT$yxFunc(0:50)
errorWPS2=eqWPS2$yxFunc(0:50)-eqCRIT$yxFunc(0:50)
errorKEQUATE=eqKEQUATE@equating[,1]-eqCRIT$yxFunc(0:50)

XPcount=kefreq(XP,0:50)[,2]

#print the weighted average absolute errors of equating of each method
print(c(sum(abs(errorWPS)*XPcount)/sum(XPcount),
	sum(abs(errorWPS2)*XPcount)/sum(XPcount),
	sum(abs(errorKEQUATE)*XPcount)/sum(XPcount)))
#running this full simulation many times usually shows very little difference in average error
#}

