###THIS CODE GIVES AN EXAMPLE OF HOW THE ANALYSES PRESENTED IN BENTON (2017) WERE COMPLETED
#(SEE http://cambridgeassessment.org.uk/Images/424229-can-ai-learn-to-equate-.pdf)
library(KernEqWPS)

#break maths data into 3 pseudo-tests
dim(mathsdata)
itesX=SampleItemsToHitTarget(mathsdata,50)#randomly select form X to have 50 marks
itesY=SampleItemsToHitTarget(mathsdata[,-itesX],50)#randomly select form Y to have 50 marks
itesA=SampleItemsToHitTarget(mathsdata[,-itesX][,-itesY],10)#randomly select form A to have 10 marks
#make form scores on the same items in two populations
#population P (typical ability)
XP=rowSums(mathsdata[,itesX])
YP=rowSums(mathsdata[,-itesX][,itesY])
AP=rowSums(mathsdata[,-itesX][,-itesY][,itesA])
summary(XP)
summary(YP)
summary(AP)
cor(cbind(XP,YP,AP))
#population Q (high ability)
XQ=rowSums(mathsdata2[,itesX])
YQ=rowSums(mathsdata2[,-itesX][,itesY])
AQ=rowSums(mathsdata2[,-itesX][,-itesY][,itesA])
summary(XQ)
summary(YQ)
summary(AQ)
cor(cbind(XQ,YQ,AQ))

#SINGLE GROUP DESIGN
#Equate form X to form Y in population P
eqEG=KernelEquateFromScoresEG(XP,YP)
plot(min(XP):max(XP),eqEG$yxFunc(min(XP):max(XP)),type='l',xlab="x",ylab="y(x)")#make a plot of equating line
plot(sort(unique(XP)),eqEG$yx,type='l',xlab="x",ylab="y(x)")#alternative way of getting plot

#NEAT DESIGN
#(imagine that different populations took different forms but both took an external anchor)
#Equate form X in population P to form Y in population Q

#CHAINED EQUATINNG
eqCHAIN=KernelChainedEquate(data.frame(x=XP,a=AP),data.frame(y=YQ,a=AQ))
plot(min(XP):max(XP),eqCHAIN$chainedFunc(min(XP):max(XP)),type='l',xlab="x",ylab="y(x)")#make a plot of equating line
#compare to criterion equate from single group design
plot(min(XP):max(XP),eqCHAIN$chainedFunc(min(XP):max(XP))-eqEG$yxFunc(min(XP):max(XP)),type='l',xlab="x",ylab="Difference from criterion",ylim=c(-3,3))

#POST-STRATIFICATION (PSE) EQUATINNG
eqPSE=PSEObservedEquate(data.frame(x=XP,a=AP),data.frame(y=YQ,a=AQ),target="x")
plot(min(XP):max(XP),eqPSE$KerEquiFunc(min(XP):max(XP)),type='l',xlab="x",ylab="y(x)")#make a plot of equating line
#compare to criterion equate from single group design
plot(min(XP):max(XP),eqPSE$KerEquiFunc(min(XP):max(XP))-eqEG$yxFunc(min(XP):max(XP)),type='l',xlab="x",ylab="Difference from criterion",ylim=c(-3,3))

#LEVINE LINEAR EQUATING
eqLEV=LevineObservedEquate(data.frame(x=XP,a=AP),data.frame(y=YQ,a=AQ),ws=c(1,0))
plot(min(XP):max(XP),eqLEV$lys(min(XP):max(XP)),type='l',xlab="x",ylab="y(x)")#make a plot of equating line
#compare to criterion equate from single group design
plot(min(XP):max(XP),eqLEV$lys(min(XP):max(XP))-eqEG$yxFunc(min(XP):max(XP)),type='l',xlab="x",ylab="Difference from criterion",ylim=c(-3,3))

#HYBRID PSE-LEVINE EQUATING
eqHYB=HybridEquate(data.frame(x=XP,a=AP),data.frame(y=YQ,a=AQ),target="x")
plot(min(XP):max(XP),eqHYB$hybridFunc(min(XP):max(XP)),type='l',xlab="x",ylab="y(x)")#make a plot of equating line
#compare to criterion equate from single group design
plot(min(XP):max(XP),eqHYB$hybridFunc(min(XP):max(XP))-eqEG$yxFunc(min(XP):max(XP)),type='l',xlab="x",ylab="Difference from criterion",ylim=c(-3,3))

#EXPERIMENTAL EQUATING USING A TRAINED NEURAL NETWORK
eqNN=EquateNN(data.frame(x=XP,a=AP),data.frame(y=YQ,a=AQ),anchortargettable=table(AP))
plot(min(XP):max(XP),eqNN$NNEqFunc(min(XP):max(XP)),type='l',xlab="x",ylab="y(x)")#make a plot of equating line
#compare to criterion equate from single group design
plot(min(XP):max(XP),eqNN$NNEqFunc(min(XP):max(XP))-eqEG$yxFunc(min(XP):max(XP)),type='l',xlab="x",ylab="Difference from criterion",ylim=c(-3,3))
#a little askew at lowest of the score range (good everywhere else)

#EXPERIMENTAL EQUATING USING An APPROXIMATION TO A TRAINED CONVOLUTIONAL NEURAL NETWORK
eqCNN=EquateCNN(data.frame(x=XP,a=AP),data.frame(y=YQ,a=AQ),anchortargettable=table(AP))
plot(min(XP):max(XP),eqCNN$CNNEqFunc(min(XP):max(XP)),type='l',xlab="x",ylab="y(x)")#make a plot of equating line
#compare to criterion equate from single group design
plot(min(XP):max(XP),eqCNN$CNNEqFunc(min(XP):max(XP))-eqEG$yxFunc(min(XP):max(XP)),type='l',xlab="x",ylab="Difference from criterion",ylim=c(-3,3))
#a little askew at lowest of the score range (good everywhere else)

#weighted mean absolute errors of equating
distX=tabulate(XP+1)#count of XP from 1:max(XP)
sum(distX*abs(eqCHAIN$chainedFunc(0:max(XP))-eqEG$yxFunc(0:max(XP))))/sum(distX)
sum(distX*abs(eqPSE$KerEquiFunc(0:max(XP))-eqEG$yxFunc(0:max(XP))))/sum(distX)
sum(distX*abs(eqLEV$lys(0:max(XP))-eqEG$yxFunc(0:max(XP))))/sum(distX)
sum(distX*abs(eqHYB$hybridFunc(0:max(XP))-eqEG$yxFunc(0:max(XP))))/sum(distX)
sum(distX*abs(eqNN$NNEqFunc(0:max(XP))-eqEG$yxFunc(0:max(XP))))/sum(distX)
sum(distX*abs(eqCNN$CNNEqFunc(0:max(XP))-eqEG$yxFunc(0:max(XP))))/sum(distX)




