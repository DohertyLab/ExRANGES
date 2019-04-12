\name{ExRANGES}
\alias{calc.slope}
\alias{sample.pval.calc}
\alias{ExRANGES}
\alias{RANGES}
\title{ExRANGES Functions}
\usage{
#calc.slopes(time.series, cycle=F, last.time.step)

#sample.pval.calc(slopes, sample.size=10000)
}
\description{
#The following two functions are used to calcuate the RANGES values to be used to calculate ExRANGES.
#calc.slopes finds the change in expression between consecutive timepoints for each gene. 
#sample.pval.calc finds the significance of the change in expression between timepoints for each gene.
}
\arguments{
#calc.slopes()

#time.series - A matrix with rows as genes and columns as sample names. Samples must be labeled as XX_samplename, XX=Numeric number for time.

#cycle - If data is cyclical then cycle equals T. 

#last.time.step (numeric) - If cycle equals T, then set last.time.step equal to differnce in time between the start and end of the time series.

#########################################################################

#sample.pval.calc()

#slopes - Output of calc.slopes(). Should be transposed pvalues of slopes between time points.

#sample.size - How many times should the slopes be sampled for each gene to calculate a pvalue.
  
#use_density_function - TRUE: pval estimation by cumulative distribution function followed by kernal density estimate.
#                       FALSE: pval estimation by empirical cumulative distribution function.
}
\examples{

#time.series<-matrix(sample(1:1000,100), nrow = 10, ncol = 10)
#colnames(time.series)<-paste(1:10,"time",sep = "_")
#slopes=calc.slopes(time.series = time.series)
#significance_of_slope=sample.pval.calc(slopes = slope, sample.size = 1000, use_density_function = T)
  

}
