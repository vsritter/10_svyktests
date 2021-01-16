Rcpp::cppFunction("
NumericMatrix Rcpp_matrix_List_sum(List x) {
    int n = x.size();
    NumericMatrix result = as<NumericMatrix>(x[0]);
    for ( int i = 1; i < n; ++i ) {
        result += as<NumericMatrix>(x[i]);
    }
    return result;
}")


onestage_fast<-function(x, strata, clusters, nPSU, fpc, lonely.psu=getOption("survey.lonely.psu"),stage=0, cal){
  if (NROW(x)==0)
    return(matrix(0,NCOL(x),NCOL(x)))
  stratvars<-tapply(1:NROW(x), list(factor(strata)), function(index){
    myonestrat(x[index,,drop=FALSE], clusters[index],
             nPSU[index][1], fpc[index], ##changed from fpc[index][1], to allow pps(brewer)
             lonely.psu=lonely.psu,stratum=strata[index][1], stage=stage,cal=cal)
  })
  p<-NCOL(x)
  nstrat<-length(unique(strata))
  nokstrat<-sum(sapply(stratvars,function(m) !any(is.na(m))))
  # apply(array(unlist(stratvars),c(p,p,length(stratvars))),1:2,sum,na.rm=TRUE)*nstrat/nokstrat
  Rcpp_matrix_List_sum(stratvars)*nstrat/nokstrat
}
environment(onestage_fast) <- asNamespace('survey')



