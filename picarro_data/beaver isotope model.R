
library(dplyr)
library(readr)
library(lubridate)

dociso<-read_csv("/Users/bobhall/Dropbox/DOC_expts/picarro_data/2018_07_23_Beaver13C_DOC/2018_07_23_beaver.csv")

beaver23jul_field<- read.csv("/Users/bobhall/Dropbox/DOC_expts/13CAdditionDatasheets.BeaverCreek.20180723.csv")

plot(dociso$samplenum,dociso$delCO2)



beaver23jul_field$dtime<- lubridate::parse_date_time(paste(beaver23jul_field$Date, beaver23jul_field$Time, sep=" "), "mdY HMS", tz='America/Denver' )
beaver23jul_field$dtime <- as.POSIXct(as.character(beaver23jul_field$dtime ), format="%Y-%m-%d %H:%M:%S")

beaver23jul<- left_join(beaver23jul_field, dociso, by="Syringe")




beaver23jul_down<-beaver23jul[beaver23jul$Site=="Downstream",]

plot(beaver23jul_down$dtime,beaver23jul_down$delCO2)
plot(beaver23jul_down$dtime,beaver23jul_down$CH4)
plot(beaver23jul_down$dtime,beaver23jul_down$delCO2)


beaver23jul_up<-beaver23jul[beaver23jul$Site=="Upstream",]

plot(beaver23jul_up$dtime,beaver23jul_up$delCO2)

head(beaver23jul_down)
plot(beaver23jul_down$dtime, beaver23jul_down$EC.us.cm)

##Tracer data

beaver23jul_down$rtime<- (as.numeric(beaver23jul_down$dtime)-as.numeric(as.POSIXct("2018-07-23 10:00:00"))) /60


beaver23jul_down<-beaver23jul_down[beaver23jul_down$rtime>70,]
beaver23jul_down<- beaver23jul_down[order(beaver23jul_down$rtime),] 


plot(beaver23jul_down$rtime,beaver23jul_down$delCO2)
plot(beaver23jul_down$rtime,beaver23jul_down$EC.us.cm)
plot(beaver23jul_down$rtime,beaver23jul_down$delCH4)

###integrate

#integrate function
nutint<-function(time,nut, bkg){
  

  #####subtract background
  nutbkg<-nut-bkg
  
  ##correct negative values to zero
  nutcorr<- ifelse(nutbkg>0, nutbkg,0)
  
  ##below routine integrates
  ydiff<- nutcorr[-1]+ nutcorr[-length(nutcorr)]
  sum(diff(time)*ydiff/2)
  
}

bkg=mean(beaver23jul_down$EC.us.cm[beaver23jul_down$rtime<80])

nutint(time=beaver23jul_down$rtime, nut=beaver23jul_down$EC.us.cm, bkg=bkg)

q<- 2000*2100/478.6

8775/60

tracermle<-function(P, nut, t, dist, M, Q, bkg){
  #p1 is velocity
  #P2 is dispersion
  
  A<-Q/P[1]
  
  d1<- M/(2*A*(pi*P[2]*t)^0.5)
  d2<-   -1*((dist-P[1]*t)^2)
  d3<- 4*P[2]*t
  
  Cxt<- bkg+(d1*exp(d2/d3))
  
  
  sqdiff<-(nut-Cxt)^2
  length(t)*(log(((sum(sqdiff, na.rm=T)/length(t))^0.5)) +0.5*log(6.28))   + ((2*sum(sqdiff, na.rm=T)/length(t))^-1)*sum(sqdiff, na.rm=T)
  
}


nut.mle<-nlm(tracermle, p=c(5, 50), nut=beaver23jul_down$EC.us.cm, t=beaver23jul_down$rtime, dist=483, M=2000*2100, Q=8775, bkg=bkg)

tracerplot<-function(U, D, nut, t,  dist, M, Q, bkg, lamb){
  
  
  A<-Q/U
  
  d1<- M/(2*A*(pi*D*t)^0.5)
  d2<-   -1*((dist-U*t)^2)
  d3<- 4*D*t
  
  
  
  
  Cxt<- bkg+(d1*exp((d2/d3)-lamb*t))
  
  plot(t,nut, xlim=c(0,t[length(t)]), ylab="Conductivity", xlab="Time from start of addition (min)", pch=16, col='blue')
  lines(t,Cxt)
}



tracerplot(U=nut.mle$estimate[1], D=nut.mle$estimate[2], nut=beaver23jul_down$EC.us.cm, t=beaver23jul_down$rtime, dist=483, M=2000*2100, Q=8775, bkg=bkg, lamb=0)

tracerplot(U=nut.mle$estimate[1], D=nut.mle$estimate[2], nut=beaver23jul_down$EC.us.cm, t=beaver23jul_down$rtime, dist=483, M=2000*2100, Q=8775, bkg=bkg, lamb=0)


Rst<-0.0112

deltoAF<- function (del) {
  Rst<-Rst
  R<-(del/1000+1)*Rst
  AF<-R/(1+R)
  AF
  
}

AFtodel <- function (AF) {
  Rst<-(Rst)
  R<- AF/(1-AF)
  del<- (R/Rst-1)*1000
  
  del
}


beaver23jul_down$isoconc<- deltoAF(beaver23jul_down$delCO2)*2.5*12
plot(beaver23jul_down$rtime,beaver23jul_down$isoconc)


#beaver background

beaver.del.bkg1<-beaver23jul$delCO2[beaver23jul$Site=="Downstream" & beaver23jul$dtime <as.numeric(as.POSIXct("2018-07-23 11:20:00")) ]  #tracer hit after 11:20

isobkg1<-deltoAF(beaver.del.bkg1)*2.5*12
mean(isobkg1)
mean(beaver.del.bkg1)
##we added 6 g glucose, 186 g/mol for 13C-6
(6*13)*1000*6/186

beaver.del.bkg2<-beaver23jul_down$delCO2[beaver23jul_down$rtime >140 ]

mean(deltoAF(beaver.del.bkg2)*2.5*12)

lm(c(0.32622,0.32657)~c(80,175))

slantbkg<-beaver23jul_down$rtime*3.684e-6 +3.259e-1

isotracerplot<-function(U, D, nut, t,  dist, M, Q, bkg, lamb){
  
  
  A<-Q/U
  
  d1<- M/(2*A*(pi*D*t)^0.5)
  d2<-   -1*((dist-U*t)^2)
  d3<- 4*D*t
  
  
  
  Cxt.0<- bkg+(d1*exp((d2/d3)))
  Cxt<- bkg+(d1*exp((d2/d3)-lamb*t))
  increase<- -(Cxt.0-Cxt)+bkg
  
  plot(t,nut, xlim=c(0,t[length(t)]), ylab="13C-DIC (mg/L)", xlab="Time from start of addition (min)", pch=16, col='blue')
  lines(t,increase)
}



isotracerplot(U=nut.mle$estimate[1], D=nut.mle$estimate[2], nut=beaver23jul_down$isoconc, t=beaver23jul_down$rtime, dist=483, M=2516, Q=8775, bkg=slantbkg,
              lamb=-0.0015)



0.15*60/(4.8*7)



##Calc for well addition of glucose.



well_iso_bkg<-deltoAF(mean(beaver.del.bkg1))*2.5*12

#add 10 per mil

0.01*well_iso_bkg




##calc isotope tracer for wally side channel.  Guessing at D to get a 2 h wide curve

tracerplot(U=600/120, D=nut.mle$estimate[2], nut=NA, t=beaver23jul_down$rtime, dist=483, M=2000*2100, Q=8775, bkg=bkg, lamb=0)


Q<-3.6*60
D<-2
M<-1000
bkg<-0
lamb<-0
dist<-41
U<-41/60
A<-Q/U


t<-seq(1:360)

d1<- M/(2*A*(pi*D*t)^0.5)
d2<-   -1*((dist-U*t)^2)
d3<- 4*D*t

Cxt<- bkg+(d1*exp((d2/d3)-lamb*t))
plot(t,Cxt)


plot(t,nut, xlim=c(0,t[length(t)]), ylab="Conductivity", xlab="Time from start of addition (min)", pch=16, col='blue')
lines(t,Cxt)



isotracertestplot<-function(U, D, nut, t,  dist, M, Q, bkg, lamb){
  
  
  A<-Q/U
  
  d1<- M/(2*A*(pi*D*t)^0.5)
  d2<-   -1*((dist-U*t)^2)
  d3<- 4*D*t
  
  
  
  Cxt.0<- bkg+(d1*exp((d2/d3)))
  Cxt<- bkg+(d1*exp((d2/d3)-lamb*t))
  increase<- -(Cxt.0-Cxt)+bkg
  
  delincrease<-((increase/bkg)-1)*1000
  
  plot(t,delincrease, xlim=c(0,t[length(t)]), ylab="del increase", xlab="Time from start of addition (min)", pch=16, col='blue')
 
}





isotracertestplot(U=600/120, D=70, t=seq(1:360), dist=600, M=2516, Q=360*60, bkg=0.3262,
              lamb=-0.0025)

#what could lambda be?  glucose is 1-5 mm/min.  Assume z=0.4

1*0.001/0.1  ##0.0025
5*0.001/0.4 #0.0125



##Blaine


isotracertestplot(U=41/60, D=2, t=seq(1:360), dist=10, M=1000, Q=60*6, bkg=0.3262,
                  lamb=-0.0025)


isotracertestpeak<-function(U, D, nut, t,  dist, M, Q, bkg, lamb){
  
  
  A<-Q/U
  
  d1<- M/(2*A*(pi*D*t)^0.5)
  d2<-   -1*((dist-U*t)^2)
  d3<- 4*D*t
  
  
  
  Cxt.0<- bkg+(d1*exp((d2/d3)))
  Cxt<- bkg+(d1*exp((d2/d3)-lamb*t))
  increase<- -(Cxt.0-Cxt)+bkg
  
  delincrease<-((increase/bkg)-1)*1000
  
  max(delincrease)
  #plot(t,delincrease, xlim=c(0,t[length(t)]), ylab="del increase", xlab="Time from start of addition (min)", pch=16, col='blue')
  
}



peak<-NA

for (i in 1:200){
  
  peak[i]<- isotracertestpeak(U=41/60, D=2, t=seq(1:360), dist=i, M=1000, Q=60*6, bkg=0.3262,
                              lamb=-0.001)

}


plot(seq(1:200), peak)


###Blaine iso model

blainecondup<-read_csv("/Users/bobhall/Dropbox/DOC_expts/cond_data/20190809_Creston.Up_Leachate_Cond_QC.csv")
blainefield<-read_csv("/Users/bobhall/Dropbox/DOC_expts/13CAdditionData_CrestonLeachate_20190809.csv")
blaineiso<-read_csv("/Users/bobhall/Dropbox/DOC_expts/picarro_data/2019_08_09_Blaine/2019_08_09_Blaine.csv")

blainefield$time<- mdy_hms(paste(blainefield$SampleDate, "", blainefield$SampleTime) )
blainefield$time<-force_tz( blainefield$time,tzone="US/Mountain")


slugstart<-force_tz(ymd_hms("2019-08-09 10:16:00"),tzone="US/Mountain") 
blainefield$slugtime<-blainefield$time - slugstart


blainecondup$time<- dmy_hms(paste(blainecondup$SampleDate, "", blainecondup$SampleTime) )
blainecondup$time<-force_tz( blainecondup$time,tzone="US/Mountain")


blaine<-left_join(blainefield, blaineiso, by="Syringe")
blaineup<-blaine[blaine$SiteID=="Up",]
blaineup<-blaineup[blaineup$SampleDate=="8/9/19",]

blainecondup$slugtime<-as.numeric(blainecondup$time - as.numeric((ymd_hms("2019-08-09 10:16:00 MDT"))))/60

plot(blainecondup$slugtime, blainecondup$SpC.us.cm)

plot(blaine$slugtime,blaine$delCO2)

mean(blainecondup$SpC.us.cm[1:8])

plot(blainecondup$slugtime, (blainecondup$SpC.us.cm-700)/max((blainecondup$SpC.us.cm-700)) )
points(blaine$slugtime, (blaine$delCO2+22)/ max(blaine$delCO2+22, na.rm=T), col="blue"   )

###calcQ
nutint(time=blainecondup$slugtime, nut=blainecondup$SpC.us.cm, bkg=700) #7045

Q<- 400*2100/(7045*60)  

##calc transport

nut.mle<-nlm(tracermle, p=c(0.5, 1), 
      nut=blainecondup$SpC.us.cm, t=blainecondup$slugtime, dist=40, M=400*2100, Q=119, bkg=700)

tracerplot(U=nut.mle$estimate[1], D=nut.mle$estimate[2], nut=blainecondup$SpC.us.cm, t=blainecondup$slugtime, dist=40, M=400*2100, Q=119, bkg=700, lamb=0)


tracerplot(U=0.5, D=2, nut=blainecondup$SpC.us.cm, t=blainecondup$slugtime, dist=40, M=400*2100, Q=119, bkg=700, lamb=0)




####Analyze glucose addition at Blaine

blainegcond<-read_csv("/Users/bobhall/Dropbox/DOC_expts/cond_data/20190808_Creston.ALL_13CGlucose_Cond.csv")
blainegfield<-read_csv("/Users/bobhall/Dropbox/DOC_expts/13CAdditionData_CrestonGlucose_20190808.csv")
blainegiso<-read_csv("/Users/bobhall/Dropbox/DOC_expts/picarro_data/2019_08_08_Blaine/2019_08_08_Blaine.csv")


blainegfield$time<- mdy_hms(paste(blainegfield$Date, "", blainegfield$SampleTime) )
blainegfield$time<-force_tz( blainegfield$time,tzone="US/Mountain")


slugstart<-force_tz(ymd_hms("2019-08-08 10:11:00"),tzone="US/Mountain") 
blainegfield$slugtime<-blainegfield$time - slugstart


blainegcond$time<- dmy_hms(paste(blainegcond$SampleDate, "", blainegcond$SampleTime) )
blainegcond$time<-force_tz( blainegcond$time,tzone="US/Mountain")
blainegcond$slugtime<-blainegcond$time - slugstart

blainegcondup<-blainegcond[blainegcond$SiteID=="Up", ]

blaineg<-left_join(blainegfield, blainegiso, by="Syringe")
blainegup<-blaineg[blaineg$SiteID=="Up",]
#blainegup<-blainegup[blainegup$SampleDate=="8/8/19",]



#plot(blainecondup$slugtime, blainecondup$SpC.us.cm)

plot(blainegup$slugtime,blainegup$delCO2)

mean(blainegcondup$SpC.us.cm[1:8])

plot(blainegcondup$slugtime, (blainegcondup$SpC.us.cm-689)/max((blainegcondup$SpC.us.cm-689)) )
points(blainegup$slugtime, (blainegup$delCO2+22)/ max(blainegup$delCO2+22, na.rm=T), col="blue", pch=16  )

###calcQ
nutint(time=as.numeric(blainegcondup$slugtime), nut=blainegcondup$SpC.us.cm, bkg=689) #8235

Q<- 400*2100/(8235*60)  

##calc transport

nut.mle<-nlm(tracermle, p=c(0.5, 1), 
             nut=blainegcondup$SpC.us.cm, t=as.numeric(blainegcondup$slugtime), dist=40, M=400*2100, Q=102, bkg=689)

tracerplot(U=nut.mle$estimate[1], D=nut.mle$estimate[2], nut=blainegcondup$SpC.us.cm, t=as.numeric(blainegcondup$slugtime), dist=40, M=400*2100, Q=102, bkg=689, lamb=0)


#now isotope
3/186

##we added 3 g glucose, 186 g/mol for 13C-6
(3*13)*1000*6/186

blainegup$isoconc <- deltoAF(blainegup$delCO2)*6*12

blainegup.bkg<-0.7805


isotracerplot(U=nut.mle$estimate[1], D=nut.mle$estimate[2], nut=blainegup$isoconc, t=as.numeric(blainegup$slugtime), dist=41, M=1258, Q=102, bkg=blainegup.bkg,
              lamb=-0.0012)

blainegup$isoconc[4]<-blainegup$isoconc[3]

nutint(time=as.numeric(blainegup$slugtime), nut=blainegup$isoconc, bkg=0.7805)

deltoAF(-18)*1.5*12

deltoAF(-20)*1.5*12

