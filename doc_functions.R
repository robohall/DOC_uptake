

## some funtions necessary for tracer DOC analysis, such as converting del to fractions, integration routines etc


nutint<-function(time,nut, bkg){
  
  ### remove all cases with no data
  
  nutd<-na.omit(data.frame(TIME=time,NUT=nut))
  
  #####subtract background
  nutbkg<-nutd$NUT-bkg
  
  ##correct negative values to zero
  nutcorr<-numeric(length(nutbkg))
  for (i in 1:length(nutbkg)){
    nutcorr[i]<-if (nutbkg[i]>0) (nutbkg[i]) else (0)}
  
  ##below routine integrates
  ydiff<- nutcorr[-1]+ nutcorr[-length(nutcorr)]
  sum(diff(nutd$TIME)*ydiff/2)
  
}



deltoAF<- function (del) {
  Rst<-0.0112372
  R<-(del/1000+1)*Rst
  AF<-R/(1+R)
  AF
  
}


AFtodel <- function (AF) {
  Rst<-0.0112372
  R<- AF/(1-AF)
  del<- (R/Rst-1)*1000
  
  del
}
