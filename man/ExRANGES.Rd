\name{ExRANGES}
\alias{calc.slope}
\alias{sample.pval.calc}
\alias{ExRANGES}
\alias{RANGES}
\title{ExRANGES Functions}
\usage{
calc.slopes(time.series, cycle=F, last.time.step)

sample.pval.calc(slopes, sample.size=10000)
}
\description{
The following two functions are used to calcuate the RANGES values to be used to calculate ExRANGES.

}
\arguments{
calc.slopes()

time.series - A matrix with rows as genes and columns as sample names. Samples must be labeled as XX_samplename, XX=Numeric number for time.

cycle - If data is cyclical then define the last time step form the last column to the first.

last.time.step - Time step from last column to the first column

#########################################################################

sample.pval.calc()

slopes - Output of calc.slopes(). Should be transposed pvalues of slopes between time points.

sample.size - How many time should the slopes be sampled for each gene to calculate a pvalue.

}
\examples{

#time.series<-matrix(sample(1:1000,100), nrow = 10, ncol = 10)
#colnames(time.series)<-paste(1:10,"time",sep = "_")
#slopes=calc.slopes(time.series = time.series)
#signifance_of_slope=sample.pval.calc(slopes = slopes, sample.size = 1000)

}