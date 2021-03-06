% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NN equating functions.R
\name{PredPercentileNN}
\alias{PredPercentileNN}
\title{Experimental function for getting percentiles using a standard neural network}
\usage{
PredPercentileNN(
  dx,
  anchortargettable = NA,
  maxx = NA,
  maxa = NA,
  WeightsList = EquateNNWeights
)
}
\arguments{
\item{dx}{Data frame with variables "x" and "a" representing scores for individual candidates on form X and on the anchor test.}

\item{anchortargettable}{Table giving distribution of anchor test scores in the target population.}

\item{maxx}{Maximum score available on form X (calculated from the data by default).}

\item{maxa}{Maximum score available on anchor test (calculated from the data by default).}

\item{WeightsList}{A list of neural network parameters used in calculations. Changing this from the default value is not recommended.}
}
\description{
This function takes a data frame of test form and anchor scores
and estimates the values of the percentiles (1st-99th) for given change in the distribution of anchor scores.
}
\examples{
#example (compare real and estimated percentiles within a fixed population)
n1=1000
t1=rnorm(n1,0.5,1)
dx=data.frame(x=round(pmin(100,pmax(0,50+20*(0.9*t1+rnorm(n1,0,sqrt(1-0.9^2))))))
	,a=round(pmin(10,pmax(0,5+2*(0.7*t1+rnorm(n1,0,sqrt(1-0.7^2)))))))

percNN=PredPercentileNN(dx,table(dx$a),maxx=100)
usualperc=as.vector(quantile(dx$x,seq(0.01,0.99,0.01)))
plot(1:99,usualperc,type='l',xlab="Percentile",ylab="Value")
lines(1:99,percNN,lty=2)

}
\keyword{KernEqWPS}
