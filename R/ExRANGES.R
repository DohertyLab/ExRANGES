#' A ExRANGES Functions
#'
#' The following two functions are used to calcuates the RANGES values to be used to calculate ExRANGES values.
#' time.series a matrix with rows as genes and column as sample names. Samples must be labeled as XX_samplename, XX=Numeric number for time.
#' cycle If data is cyclical then define the last time step form the last column to the first.
#' last.time.step time step from last column to the fist column
#' Time series, slope, ExRANGES
#' calc.slopes(time.series.matrix, cycle=F)
#' calc.slopes(time.series.matrix, cycle=T, last.time.step=3.5)

require(pbapply)
calc.slopes<-function(time.series, cycle=F, last.time.step){
  print("Calculating Slopes")
  slopes<-pbapply(time.series,1,function(x) getrowdiff(rw = x, cycle = cycle, last.time.step = last.time.step))
  if(cycle==F){
    rownames(slopes)<-paste(paste("Previous Timepoint", rownames(slopes), sep = "->"))
  }
  else{
    rownames(slopes)<-paste(ifelse(test = as.numeric(sapply(strsplit(rownames(slopes), split = "_"), FUN = function(x) x[1]))==min(as.numeric(sapply(strsplit(rownames(slopes), split = "_"), FUN = function(x) x[1]))),
                                   yes = "Last Timepoint",
                                   no =  "Previous Timepoint"), rownames(slopes), sep = "->")
  }
  return(slopes)
}

getrowdiff<-function(rw, cycle, last.time.step){
  splitnames<-strsplit(names(rw), split="_", fixed=TRUE)
  samples<-sapply(splitnames, function(x) as.numeric(x[1]))
  if(cycle==F){
    diffnums<-lapply(1:length(rw), function(x) GRD2(x, samples, rw))
  }
  else{
    diffnums<-lapply(1:length(rw), function(x) GRD2.C(x, samples, rw, last.time.step))
  }
  diffnums<-unlist(diffnums)
  return(diffnums)
}
GRD2<-function(x, samples, rw){
  if(length((samples[which(samples[x] < samples)])) > 0) {
    ret<-(rw[which(samples == min(samples[which(samples[x]<samples)]))]-rw[x])/abs((min(samples[which(samples[x]<samples)])-samples[x]))
  }
  else{ret<-NULL}
  return(ret)
}
GRD2.C<-function(x, samples, rw, last.time.step){
  if ( samples[x] == max(samples)){
    ret<-(rw[which(samples == min(samples))]-rw[x])/last.time.step
  }
  else{
    ret<-(rw[which(samples == min(samples[which(samples[x]<samples)]))]-rw[x])/abs((min(samples[which(samples[x]<samples)])-samples[x]))
  }
  return(ret)
}

#' The following functions is used to calcuates the RANGES values to be used to calculate ExRANGES values.
#' slopes from calc.slopes, should be slopes between time points.
#' sample.size How many time should the slopes be sampled for each gene to calculate a pvalue.
#' Time series, slope, ExRANGES
#' sample.pval.calc(slopes=matrix.of.slopes, sample.size=10000)
                     
sample.pval.calc<-function(slopes, sample.size=10000){
  print("Calculating Slope Distributions")
  distributions<-pbapply(slopes,2,function(x) sample(x,sample.size,replace=T))
  print("Calculating per Gene p.value")
  list.of.ecdf<-pblapply(1:length(distributions[1,]),function(x) ecdf(distributions[,x]))
  gene.change<-(slopes) # change when colnames on an object with less than two dimensions
  colnames(gene.change)<-colnames(distributions)
  pvals<-lapply(1:length(gene.change[1,]), function(x) list.of.ecdf[[x]](gene.change[,x]))
  print("Formating Matrix...")
  pvals<-sapply(pvals, unlist)
  colnames(pvals)<-colnames(distributions)
  pvals<-t(pvals)
  pvals<-pvals.transform(pvals, sample.size)
  colnames(pvals)<-rownames(slopes)
 # if(length(colnames(pvals)) > length(colnames(unique(pvals)))){
  slopenames=unique(colnames(pvals))
  pvals_ret=pbsapply(unique(colnames(pvals)),function(x) apply(pvals[,which(colnames(pvals)==x)],1,mean))
 # } else{
 #   pvals_ret=pvals
 # }
  return(pvals_ret)
}

###Log transform pvals
pvals.transform<-function(pvals, sample.size){
  test<-pvals
  test.up<-test
  test.down<-1-test
  test.up[test.up==0]<-(1/sample.size)
  test.down[test.down==0]<-(1/sample.size)
  test.up<--log(test.up, 10)
  test.down<--log(test.down, 10)
  test<-ifelse(test.up<test.down, -(test.down), test.up)
  pvals<-test
  pvals<--(pvals)
  return(pvals)
}





