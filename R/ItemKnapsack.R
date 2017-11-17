#' Function to find a set of items within a test with given available maximum score
#'
#' This function finds a random selection of indices within a matrix of items scores so that the sum of the maximum scores available on these items meets a predefined total.
#' This is done using brute computing power to solve an extremely simple version of the famous Knapsack problem where all the values we need to consider (the maximum item scores) are (fairly) small integers.
#' The function will print the number of selected items with each number of marks in the chosen solution.
#' 
#' @param itemat Matrix of item scores.
#' @param targtot The target sum of maximum scores on selected items.
#' @param chooseclosestproportional Logical value indicating whether the proportion of items with each number of available marks be chosen so as to match the proportions for the overall test as closely as possible (default=TRUE).
#'
#' @return The function returns a vector denoting which items should be included.
#'
#' @examples
#' library(KernEqWPS)
#' summary(mathsdata)
#' maxes=apply(mathsdata,2,max,na.rm=TRUE)
#' sum(maxes)
#' #make form X with 50 marks
#' itesX=SampleItemsToHitTarget(as.matrix(mathsdata),50)
#' formX=mathsdata[,itesX]
#' #note difference in selected item totals if we set chooseclosestproportional=FALSE
#' SampleItemsToHitTarget(as.matrix(mathsdata),50,chooseclosestproportional=FALSE)
#' #make form Y with 50 marks from the remainder
#' remain=mathsdata[,-itesX]
#' itesY=SampleItemsToHitTarget(as.matrix(remain),50)
#' formY=remain[,itesY]
#' #make form A with 10 marks from the remainder
#' remain=remain[,-itesY]
#' itesA=SampleItemsToHitTarget(as.matrix(remain),10)
#' formA=remain[,itesA]
#' 
#' length(itesX)
#' length(itesY)
#' length(itesA)
#' 
#' #check
#' x=rowSums(formX)
#' y=rowSums(formY)
#' a=rowSums(formA)
#' summary(x)
#' summary(y)
#' summary(a)
#' cor(cbind(x,y,a))
#' 
#' @keywords KernEqWPS knapsack
#' @export
SampleItemsToHitTarget=function(itemat,targtot,chooseclosestproportional=TRUE){

maxes=as.vector(apply(itemat,2,max,na.rm=TRUE))
nites=length(maxes)
tabmax=table(maxes)
#tabmax
tabmaxes=as.numeric(names(tabmax))
#remove zeros
tabmax=tabmax[tabmaxes>0]
tabmaxes=as.numeric(names(tabmax))
#LOOP THROUGH ALL COMBINATIONS OF THESE POSSIBILITIES 
#(CAN'T TAKE MORE ITEMS THAN ARE AVAILABLE)
#(ALSO DON'T NEED MORE ITEMS THAN TARGET TOTAL DIVIDED BY MINIMUM AVAILABLE ITEM SCORE)
comb1=data.frame(0:min(tabmax[1],floor(targtot/tabmaxes[1])))
names(comb1)[1]=paste("N",tabmaxes[1],sep="")
if(length(tabmax)>1){
for(jjz in 2:length(tabmax)){
comb1=merge(comb1,data.frame(0:min(tabmax[jjz],floor(targtot/tabmaxes[jjz]))))
names(comb1)[jjz]=paste("N",tabmaxes[jjz],sep="")}}

comb1$nchose=rowSums(comb1)
comb1$tot=0
for(jjz in 1:length(tabmax)){comb1$tot=comb1$tot+tabmaxes[jjz]*comb1[,jjz]}
targcombs=comb1[comb1$tot==targtot,]
dim(targcombs)
if(dim(targcombs)[1]==0){
	print("No available combination! Try a different target total")
	return(NULL)
	}

#work out how many ways there are of doing each combination
targcombs$ncomb=1
for(jjz in 1:length(tabmax)){targcombs$ncomb=targcombs$ncomb*choose(tabmax[jjz],targcombs[,jjz])}

#keep as close as possible to expected proportions of item max distributions
proptab=prop.table(tabmax)
targcombs$diffromprop=0
for(jjz in 1:length(tabmax)){targcombs$diffromprop=targcombs$diffromprop+abs(proptab[jjz]-(targcombs[,jjz]/targcombs$nchose))}


if(chooseclosestproportional==TRUE){targcombs=targcombs[targcombs$diffromprop==min(targcombs$diffromprop),]}


targcombs$wt=targcombs$ncomb/sum(targcombs$ncomb)
targcombs$cumwt=cumsum(targcombs$wt)

#now pick just one combination
rand1=runif(1)
selection=max(1,findInterval(rand1,targcombs$cumwt))
neach=targcombs[selection,1:length(tabmax)]
print(neach)
chosenites=NULL
for(jjz in 1:length(tabmax)){
	possibles=(1:nites)[maxes==tabmaxes[jjz]]
	if(length(possibles)==1){possibles=c(possibles,possibles)}#deal with a bug in how "sample" works if sampling from one element
	newchoice=sample(possibles,as.numeric(neach[jjz]))
	chosenites=c(chosenites,newchoice)
	}
sum(maxes[chosenites])
return(chosenites)}
