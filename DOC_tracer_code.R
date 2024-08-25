
#Generic DOC tracer code based for initial exploration of "blaine0808_code"

library(dplyr)
library(lubridate)
library(zoo)

source("./doc_functions.R")



####First, Blaine 8 August

blaine0808_solute<-read.csv("./DOC_runs/13CAdditionData_CrestonGlucose_20190808_analyzed.csv")
blaine0808_co2<-read.csv("./picarro_data/2019_08_08_Blaine/2019_08_08_Blaine.csv")


blaine0808<-left_join(blaine0808_solute,blaine0808_co2, by="Syringe")


blaine0808$dtime <- dmy_hms(paste(blaine0808$Date, blaine0808$SampleTime), tz="America/Denver")
blaine0808$MinFrom0 <- (as.numeric(blaine0808$dtime) - as.numeric(dmy_hms("08/08/19 10:11:00", tz="America/Denver")))/60

##here below call it data
#use only down for now
blaine0808_down<- blaine0808[blaine0808$SiteID=="Down", ]


data<-blaine0808_down[!is.na(blaine0808_down$CO2), ]  #delete rows with NA and call it data


##release info
sugar<-3*(72/180)*1000  #mg added C
condmass<- 400*2100
add_ratio<-sugar/condmass  #0.001428571
DIC<-  5 #meq/L
reach_length<-61
travel_time<-134
velocity<- reach_length/travel_time

##SpC and hydrology
#blaine0808_up$SpC.us.cm_1[1]<-690
bkg<-mean( data$SpC.us.cm_1[1])
plot(data$MinFrom0, (data$SpC.us.cm_1-bkg))


spcint<- nutint(data$MinFrom0, data$SpC.us.cm_1, bkg=bkg )
data_Q<- condmass/ spcint  #171 L/min

delbkg<-mean(data$Delta.13C[1])
plot(data$MinFrom0, (data$Delta.13C))


##########
##Now CO2

plot(data$MinFrom0, (data$delCO2))
bkg_delCO2<-mean((data$delCO2[1:2]))
data_co2AF_corr<-deltoAF(data$delCO2)-deltoAF(bkg_delCO2)
data_13DIC<-5*13*data_co2AF_corr #mole to mass here
plot(data$MinFrom0, data_13DIC, pch=16, col="red", ylab= "13C mass of extra DIC", 
     xlab="Time from solute addition (min)")

int13co2<- nutint(data$MinFrom0, data_13DIC, bkg=0 )
int13co2*data_Q


######



AF_corr<-deltoAF(data$Delta.13C)-deltoAF(delbkg)
data_13C<-data$DOC_conc*AF_corr

mean(data$DOC_conc)/12 #DOC in mmol/L


plot(data$MinFrom0, (data_13C), ylim=c(0,0.2), pch=16, xlab="Time from solute addition (min)", 
     ylab= "13C mass of added DOC", col="blue")
points(data$MinFrom0, (0.00143*(data$SpC.us.cm_1-bkg)),
       pch=16)

##plot that which was taken up
removed_c<- add_ratio*(data$SpC.us.cm_1-bkg)-data_13C
removed_c<-na.approx(removed_c)

int_removed_C<- nutint(data$MinFrom0, removed_c, bkg=0)
plot(data$MinFrom0,  removed_c, ylim=c(0,0.065), pch=16, xlab="Time from solute addition (min)", 
     ylab= "13C taken up (black) or mineralized (red)")
points(data$MinFrom0, data_13DIC[1:24], pch=16, col="red")

int13DIC<- nutint(data$MinFrom0, data_13DIC[1:24], bkg=0 ) 

int13DIC_flux<-int13DIC*data_Q

int13c<- nutint(data$MinFrom0, data_13C, bkg=0 ) 

int13c_flux<-int13c*data_Q  #that which passed by the downstream station
int_removed_C_flux<-int_removed_C*data_Q #that whic got taken up

#uptake length
K<-  -(log(int13c/spcint) - log(sugar/condmass))/reach_length  #0.036
1/K
#another way
-(log(int13c_flux) -log(sugar))/reach_length





#############
####model
#dB/dx = -k1*B

#integrate that which was taken up
#U= (B_0/-K) * exp(-K*x)
#(B_0/-K) * exp(-K*40)
#rearrange to solve for B_0
U<-int_removed_C_flux

B_0<- K*U/(1-exp(-K*40))  ##correct

x<-seq(1:40)
plot(x, B_0*exp(x*-K))

#C is modeled flux of 13CDIC
Kc<-0.0018  ##fit by eye here
C<- (B_0*Kc/-K)*(exp(-K*x)-1)
plot(x,C)


##immediate release
f1<-0.20  #fit by eye
removed_c_flux<- removed_c*data_Q
pred_imm<- f1* removed_c_flux
DIC_flux<- data_13DIC*data_Q
plot(data$MinFrom0,  pred_imm, ylim=c(0,3), pch=16, xlab="Time from solute addition (min)", 
     ylab= "pred 13Creleased (black) or data (red)")
points(data$MinFrom0, DIC_flux, pch=16, col="red")



#cumulative uptake
#removed_c_flux_diff<- removed_c_flux[-1]+ removed_c_flux[-length(removed_c_flux)]
#U_g<- cumsum(diff(data$MinFrom0)*removed_c_flux_diff/2)
#U_g<-c(0,U_g) # adding back the deleted term from differencing, check this
#B_0<- K*U_g/(1-exp(-K*40))
#B_0<- K*U_g/(1-exp(-K*40))- ( (K*U_g/(1-exp(-K*40))) * exp(-Kc*blaine0808_up$MinFrom0)) 

U_g<-NA
U_g[1]<-0
for (i in 2: length (data$MinFrom0) ){
  int_u<- (data$MinFrom0[i]-data$MinFrom0[i-1])*(removed_c_flux[i])
  U_g[i]<- U_g[i-1]+int_u
}
U_g<-U_g*(1-f1)
U_g


B_0<- K*U_g/(1-exp(-K*reach_length))  * exp(-Kc*velocity*data$MinFrom0)
pred_delayed<- (B_0*Kc/-K)*(exp(-K*reach_length)-1)


pd<-(U_g*Kc* exp(-Kc*velocity*data$MinFrom0)) #simplification of above


pred_DIC<- pred_imm  +  pred_delayed

plot(data$MinFrom0,   DIC_flux, ylim=c(0,3), pch=16, xlab="Time from solute addition (min)", 
     ylab= "pred 13Creleased (black) or data (red)")
lines(data$MinFrom0,pred_delayed , pch=16, col="red")
lines(data$MinFrom0, pred_imm, pch=16, col="blue")
lines(data$MinFrom0, pred_DIC, pch=16, col="purple")

lines(data$MinFrom0, pred_imm/data_Q, pch=16, col="blue")


intpred_DIC<- nutint(data$MinFrom0, pred_DIC, bkg=0 )  #125
intpred_DIC/sugar






########
###now for leachate
########




blaine0809_solute<-read.csv("./DOC_runs/13CAdditionData_CrestonLeachate_20190809_analyzed.csv")

blaine0809_co2<-read.csv("./picarro_data/2019_08_09_Blaine/2019_08_09_Blaine.csv")


blaine0809<-left_join(blaine0809_solute,blaine0809_co2, by="Syringe")


blaine0809$dtime <- mdy_hms(paste(blaine0809$Date, blaine0809$SampleTime), tz="America/Denver")
blaine0809$MinFrom0 <- (as.numeric(blaine0809$dtime) - as.numeric(mdy_hms("08/09/19 10:15:00", tz="America/Denver")))/60

##here below call it data
#use only down for now
blaine0809_up<- blaine0809[blaine0809$SiteID=="Up", ]


data<-blaine0809_up  #delete rows with NA and call it data


##release info
sugar<-1160  #mg added C
condmass<- 400*2100
add_ratio<-sugar/condmass  #0.001428571
DIC<-  5 #meq/L
reach_length<-41
velocity = 0.433

bkg<-mean( data$SpC.us.cm_1[1:5])
plot(data$MinFrom0[1:41], (data$SpC.us.cm_1[1:41]-bkg))


spcint<- nutint(data$MinFrom0[1:39], data$SpC.us.cm_1[1:39], bkg=bkg )
data_Q<- condmass/ spcint  #126 L/min

delbkg<-mean(data$Delta.13C[1:5])
plot(data$MinFrom0, (data$Delta.13C))

AF_corr<-deltoAF(data$Delta.13C)-deltoAF(delbkg)
data_13C<-data$DOC_conc*AF_corr

plot(data$MinFrom0[1:41], (data_13C[1:41]), ylim=c(0,0.15), pch=16, xlab="Time from solute addition (min)", 
     ylab= "13C mass of added DOC",, col="blue")
points(data$MinFrom0, (add_ratio*(data$SpC.us.cm_1-bkg)),
       pch=16)

nutint(data$MinFrom[1:39],data_13C[1:39],bkg=0)




##########
##Now CO2

plot(data$MinFrom0, (data$delCO2))
bkg_delCO2<-mean((data$delCO2[1:2]))
data_co2AF_corr<-deltoAF(data$delCO2)-deltoAF(bkg_delCO2)
data_13DIC<-5*13*data_co2AF_corr #mole to mass here
plot(data$MinFrom0, data_13DIC, pch=16, col="red", ylab= "13C mass of extra DIC", 
     xlab="Time from solute addition (min)")

int13co2<- nutint(data$MinFrom0, data_13DIC, bkg=0 )
int13co2*data_Q
data$DIC_flux<-data_13DIC*data_Q


##plot that which was taken up
removed_c<- add_ratio*(data$SpC.us.cm_1-bkg)-data_13C
removed_c<-na.approx(removed_c)

int_removed_C<- nutint(data$MinFrom0[1:39], removed_c[1:39], bkg=0)

plot(data$MinFrom0,  (add_ratio*(data$SpC.us.cm_1-bkg)-data_13C), ylim=c(0,0.1), pch=16, xlab="Time from solute addition (min)", 
     ylab= "13C taken up (black) or mineralized (red)")
#points(data$MinFrom0, data_13DIC, pch=16, col="red")

int13c<- nutint(data$MinFrom0[1:39], data_13C[1:39], bkg=0 ) 

int13DIC<- nutint(data$MinFrom0[1:39], data_13C[1:39], bkg=0 ) 

int13c_flux<-int13c*data_Q  #that which passed by the downstream station
int_removed_C_flux<-int_removed_C*data_Q #that whic got taken up

#uptake length
K<-  -(log(int13c/spcint) - log(sugar/condmass))/reach_length  #0.036
1/K
#another way
-(log(int13c_flux) -log(sugar))/reach_length




#############
####model
#dB/dx = -k1*B

#integrate that which was taken up
#U= (B_0/-K) * exp(-K*x)
#(B_0/-K) * exp(-K*40)
#rearrange to solve for B_0
U<-int_removed_C_flux

B_0<- K*U/(1-exp(-K*reach_length))  ##correct

x<-seq(1:40)
plot(x, B_0*exp(x*-K))

#C is modeled flux of 13CDIC
Kc<-0.002  ##fit by eye here
C<- (B_0*Kc/-K)*(exp(-K*x)-1)
plot(x,C)


##immediate release
f1<-0.30  #fit by eye
removed_c_flux<- removed_c[1:39]*data_Q
pred_imm<- f1* removed_c_flux
DIC_flux<- data_13DIC[1:39]*data_Q
plot(data$MinFrom0[1:39],  pred_imm[1:39], ylim=c(0,3), pch=16, xlab="Time from solute addition (min)", 
     ylab= "pred 13Creleased (black) or data (red)")
points(data$MinFrom0[1:32], DIC_flux[1:32], pch=16, col="red")



#cumulative uptake
#removed_c_flux_diff<- removed_c_flux[-1]+ removed_c_flux[-length(removed_c_flux)]
#U_g<- cumsum(diff(data$MinFrom0)*removed_c_flux_diff/2)
#U_g<-c(0,U_g) # adding back the deleted term from differencing, check this
#B_0<- K*U_g/(1-exp(-K*40))
#B_0<- K*U_g/(1-exp(-K*40))- ( (K*U_g/(1-exp(-K*40))) * exp(-Kc*blaine0808_up$MinFrom0)) 


U_g<-NA
U_g[1]<-0
for (i in 2: length (data$MinFrom0) ){
  int_u<- (data$MinFrom0[i]-data$MinFrom0[i-1])*(removed_c_flux[i])
  U_g[i]<- U_g[i-1]+int_u
}
U_g[40:45]<-rep(U_g[39],6)
U_g<-(1-f1)*U_g


B_0<- K*U_g/(1-exp(-K*40))  * exp(-Kc*velocity*data$MinFrom0)
pred_delayed<- (B_0*Kc/-K)*(exp(-K*40)-1)




pd<-(U_g*Kc* exp(-Kc*data$MinFrom0)) #simplification of above


pred_DIC<- pred_imm  +  pred_delayed[!is.na(pred_delayed)]

plot(data$MinFrom0,   data$DIC_flux, ylim=c(0,3), pch=16, xlab="Time from solute addition (min)", 
     ylab= "pred 13Creleased (black) or data (red)")
lines(data$MinFrom0,pred_delayed , pch=16, col="red")
lines(data$MinFrom0[1:39], pred_imm, pch=16, col="blue")
lines(data$MinFrom0, pred_DIC, pch=16, col="purple")


intpred_DIC<- nutint(data$MinFrom0, pred_DIC, bkg=0 )  #125
intpred_DIC/sugar



########
#######
#Leachate Down

blaine0809_down<- blaine0809[blaine0809$SiteID=="Down", ]


data<-blaine0809_down  #delete rows with NA and call it data


##release info
sugar<-1460  #mg added C
condmass<- 400*2100
add_ratio<-sugar/condmass 
DIC<-  5.4 #meq/L
reach_length<-65

ntt<-150
velocity<-reach_length/ntt

##2 points dorked up from not shaking probe

data$SpC.us.cm_1[21:23]<-c(686, 683,678)
bkg<-mean( data$SpC.us.cm_1[1:3])
plot(data$MinFrom0[1:28], (data$SpC.us.cm_1[1:28]-bkg))

lm(c(0,-2.6)~c(70,345))

condcorr<- data$SpC.us.cm_1[1:28] - (0.662-0.009455*data$MinFrom0[1:28])

spcint<- nutint(data$MinFrom0[1:28], condcorr, bkg=bkg )
data_Q<- condmass/ spcint  #211 L/min


ntt<-150
velocity<-65/150
depth<- (data_Q*0.001)/(velocity*1.4)

delbkg<-mean(data$Delta.13C[1:2])
plot(data$MinFrom0, (data$Delta.13C))


AF_corr<-deltoAF(data$Delta.13C)-deltoAF(delbkg)
data_13C<-data$DOC_conc*AF_corr

plot(data$MinFrom0[1:28], (data_13C[1:28]), ylim=c(0,0.15), pch=16, xlab="Time from solute addition (min)", 
     ylab= "13C mass of added DOC", col="blue")
points(data$MinFrom0[1:28], (add_ratio*(condcorr-bkg)),
       pch=16)

nutint(data$MinFrom0[1:28],data_13C[1:28],bkg=0)*data_Q

##plot that which was taken up
removed_c<- add_ratio*(condcorr-bkg)-data_13C[1:28]
int_removed_C<- nutint(data$MinFrom0[1:28], removed_c[1:28], bkg=0)
plot(data$MinFrom0[1:28],  (add_ratio*(condcorr-bkg)-data_13C[1:28]), ylim=c(-0.02,0.1), pch=16, xlab="Time from solute addition (min)", 
     ylab= "13C taken up (black) or mineralized (red)")
#points(data$MinFrom0, data_13DIC, pch=16, col="red")

int13c<- nutint(data$MinFrom0[1:28], data_13C[1:28], bkg=0 ) 

int13c_flux<-int13c*data_Q  #that which passed by the downstream station
int_removed_C_flux<-int_removed_C*data_Q #that whic got taken up

#uptake length
K<-  -(log(int13c/spcint) - log(sugar/condmass))/reach_length  
1/K
#another way
(log(int13c_flux) -log(sugar))/reach_length
 vf<-1440*K*data_Q*0.001/1.4

# WHat would DOC have to ve to support ER of -6 g O2
 ER.c<- 6*0.4
 ER.c/vf

##########
##Now CO2

plot(data$MinFrom0, (data$delCO2))
bkg_delCO2<-mean((data$delCO2[1:2]))
data_co2AF_corr<-deltoAF(data$delCO2)-deltoAF(bkg_delCO2)
data_13DIC<-5*13*data_co2AF_corr #mole to mass here
plot(data$MinFrom0, data_13DIC, pch=16, col="red", ylab= "13C mass of extra DIC", 
     xlab="Time from solute addition (min)")

int13co2<- nutint(data$MinFrom0, data_13DIC, bkg=0 )
int13co2*data_Q
#just the release
nutint(data$MinFrom0[1:28], data_13DIC[1:28], bkg=0 )*data_Q
nutint(data$MinFrom0[1:28], data_13DIC[1:28], bkg=0 )*177


#############
####model
#dB/dx = -k1*B

#integrate that which was taken up
#U= (B_0/-K) * exp(-K*x)
#(B_0/-K) * exp(-K*40)
#rearrange to solve for B_0
U<-int_removed_C_flux

B_0<- K*U/(1-exp(-K*40))  ##correct

x<-seq(1:40)
plot(x, B_0*exp(x*-K))

#C is modeled flux of 13CDIC
Kc<-0.005  ##fit by eye here
C<- (B_0*Kc/-K)*(exp(-K*x)-1)
plot(x,C)


##immediate release
f1<-0.4  #fit by eye
removed_c_flux<- removed_c*data_Q
pred_imm<- f1* removed_c_flux
DIC_flux<- data_13DIC*data_Q
plot(data$MinFrom0[1:28],  pred_imm, ylim=c(0,3), pch=16, xlab="Time from solute addition (min)", 
     ylab= "pred 13Creleased (black) or data (red)")
points(data$MinFrom0, DIC_flux, pch=16, col="red")



#cumulative uptake
#removed_c_flux_diff<- removed_c_flux[-1]+ removed_c_flux[-length(removed_c_flux)]
#U_g<- cumsum(diff(data$MinFrom0[1:28])*removed_c_flux_diff/2)
#U_g<-c(0,U_g) # adding back the deleted term from differencing, check this
#B_0<- K*U_g/(1-exp(-K*40))
#B_0<- K*U_g/(1-exp(-K*40))- ( (K*U_g/(1-exp(-K*40))) * exp(-Kc*blaine0808_up$MinFrom0)) 
U_g<-NA
U_g[1]<-0
for (i in 2: length (data$MinFrom0) ){
  int_u<- (data$MinFrom0[i]-data$MinFrom0[i-1])*(removed_c_flux[i])
  U_g[i]<- U_g[i-1]+int_u
}
U_g[29:32]<-rep(U_g[27],4)
U_g<-(1-f1)*U_g


B_0<- K*U_g/(1-exp(-K*reach_length))  * exp(-Kc*velocity*data$MinFrom0)
pred_delayed<- (B_0*Kc/-K)*(exp(-K*reach_length)-1)

pd<-(U_g*Kc* exp(-Kc*data$MinFrom0)) #simplification of above


pred_DIC<- c(pred_imm,rep(0,4))  +  pred_delayed

plot(data$MinFrom0,   DIC_flux, ylim=c(0,3), pch=16, xlab="Time from solute addition (min)", 
     ylab= "pred 13Creleased (black) or data (red)")
lines(data$MinFrom0,pred_delayed , pch=16, col="red")
lines(data$MinFrom0[1:28], pred_imm, pch=16, col="blue")
lines(data$MinFrom0, pred_DIC, pch=16, col="purple")


intpred_DIC<- nutint(data$MinFrom0, pred_DIC, bkg=0 )  #125
intpred_DIC/sugar




######
#####
####second Creston glucose experiment

#no joining CO2 and DOC here.  Instead we joined by hand because we resused some CO2 syringes, this disabing the ability to use left_join based on syringes

blaine0815<-read.csv("./DOC_runs/13CAdditionData_CrestonGlucose_20190815_analyzed.csv")



blaine0815$dtime <- mdy_hms(paste(blaine0815$SampleDate, blaine0815$SampleTime), tz="America/Denver")
blaine0815$MinFrom0 <- (as.numeric(blaine0815$dtime) - as.numeric(mdy_hms("08/15/19 09:00:00", tz="America/Denver")))/60



blaine0815_up<- blaine0815[blaine0815$SiteID=="Up", ]

plot(blaine0815_up$MinFrom0,blaine0815_up$SpC.us.cm, xlim=c(0,90))

data<-blaine0815_up[!is.na(blaine0815_up$CO2), ]  #delete rows with NA and call it data


##release info
sugar<-3*(72/180)*1000  #mg added C
condmass<- 400*2100
add_ratio<-sugar/condmass  #0.001428571
DIC<-  mean(blaine0815$Alkalinity..meq.L., na.rm=T) #meq/L
reach_length<-45
travel_time<-74
velocity<- reach_length/travel_time


##SpC and hydrology

bkg<-mean( data$SpC.us.cm[1:4])
plot(data$MinFrom0, (data$SpC.us.cm-bkg))


spcint<- nutint(data$MinFrom0[1:28], data$SpC.us.cm[1:28], bkg=bkg )
data_Q<- condmass/ spcint  #113 L/min

delbkg<-mean(data$Delta.13C[1:4])
plot(data$MinFrom0, (data$Delta.13C))

## missing  [23].  Interpolate using zoo
data$Delta.13C[1:33]<- na.approx(data$Delta.13C[1:33])


####CO2

plot(data$MinFrom0, (data$delCO2))
bkg_delCO2<-mean((data$delCO2[1:4]))
data_co2AF_corr<-deltoAF(data$delCO2)-deltoAF(bkg_delCO2)
data_13DIC<-5*13*data_co2AF_corr #mole to mass here
plot(data$MinFrom0, data_13DIC, pch=16, col="red", ylab= "13C mass of extra DIC", 
     xlab="Time from solute addition (min)")

int13co2<- nutint(data$MinFrom0, data_13DIC, bkg=0 )
int13co2*data_Q


####DOC

AF_corr<-deltoAF(data$Delta.13C)-deltoAF(delbkg)
data_13C<-data$DOC_conc*AF_corr

mean(data$DOC_conc, na.rm=T)/12 #DOC in mmol/L


plot(data$MinFrom0, (data_13C), ylim=c(0,0.2), pch=16, xlab="Time from solute addition (min)", 
     ylab= "13C mass of added DOC", col="blue")
points(data$MinFrom0, (add_ratio*(data$SpC.us.cm-bkg)),
       pch=16)

##plot that which was taken up
removed_c<- add_ratio*(data$SpC.us.cm-bkg)-data_13C
removed_c<-na.approx(removed_c)

int_removed_C<- nutint(data$MinFrom0[1:28], removed_c[1:28], bkg=0)  ##only integrate addition
plot(data$MinFrom0[1:length( removed_c)],  removed_c, ylim=c(-0.005,0.065), pch=16, xlab="Time from solute addition (min)", 
     ylab= "13C taken up (black) or mineralized (red)")
points(data$MinFrom0, data_13DIC, pch=16, col="red")

int13DIC<- nutint(data$MinFrom0, data_13DIC, bkg=0 ) 

int13DIC_flux<-int13DIC*data_Q

int13c<- nutint(data$MinFrom0, data_13C, bkg=0 ) 

int13c_flux<-int13c*data_Q  #that which passed by the downstream station
int_removed_C_flux<-int_removed_C*data_Q #that whic got taken up

#uptake length
K<-  -(log(int13c/spcint) - log(sugar/condmass))/reach_length  #0.036
1/K
#another way
-(log(int13c_flux) -log(sugar))/reach_length



#############
####model
#dB/dx = -k1*B

#integrate that which was taken up
#U= (B_0/-K) * exp(-K*x)
#(B_0/-K) * exp(-K*40)
#rearrange to solve for B_0
U<-int_removed_C_flux

B_0<- K*U/(1-exp(-K*reach_length))  ##correct

x<-seq(1:reach_length)
plot(x, B_0*exp(x*-K))

#C is modeled flux of 13CDIC
Kc<-0.0018  ##fit by eye here
C<- (B_0*Kc/-K)*(exp(-K*x)-1)
plot(x,C)


##immediate release
f1<-0.20  #fit by eye
removed_c_flux<- removed_c[1:28]*data_Q
pred_imm<- f1* removed_c_flux
DIC_flux<- data_13DIC*data_Q
plot(data$MinFrom0[1:28],  pred_imm, ylim=c(-1,3),xlim=c(0,1440), pch=16, xlab="Time from solute addition (min)", 
     ylab= "pred 13Creleased (black) or data (red)")
points(data$MinFrom0, DIC_flux, pch=16, col="red")



#cumulative uptake
#removed_c_flux_diff<- removed_c_flux[-1]+ removed_c_flux[-length(removed_c_flux)]
#U_g<- cumsum(diff(data$MinFrom0)*removed_c_flux_diff/2)
#U_g<-c(0,U_g) # adding back the deleted term from differencing, check this
#B_0<- K*U_g/(1-exp(-K*40))
#B_0<- K*U_g/(1-exp(-K*40))- ( (K*U_g/(1-exp(-K*40))) * exp(-Kc*blaine0808_up$MinFrom0)) 

U_g<-NA
U_g[1]<-0
for (i in 2: length (data$MinFrom0) ){
  int_u<- (data$MinFrom0[i]-data$MinFrom0[i-1])*(removed_c_flux[i])
  U_g[i]<- U_g[i-1]+int_u
}
U_g<-U_g*(1-f1)
U_g


B_0<- K*U_g/(1-exp(-K*reach_length))  * exp(-Kc*velocity*data$MinFrom0)
pred_delayed<- (B_0*Kc/-K)*(exp(-K*reach_length)-1)


pd<-(U_g*Kc* exp(-Kc*data$MinFrom0)) #simplification of above


pred_DIC<- pred_imm  +  pred_delayed[1:28]

plot(data$MinFrom0,   DIC_flux, ylim=c(0,3), pch=16, xlab="Time from solute addition (min)", 
     ylab= "pred 13Creleased (black) or data (red)")
lines(data$MinFrom0,pred_delayed , pch=16, col="red")
lines(data$MinFrom0, pred_imm, pch=16, col="blue")
lines(data$MinFrom0, pred_DIC, pch=16, col="purple")

lines(data$MinFrom0, pred_imm/data_Q, pch=16, col="blue")


intpred_DIC<- nutint(data$MinFrom0, pred_DIC, bkg=0 )  #125
intpred_DIC/sugar






#####
#####
##CO2 conc
plot(data$MinFrom0, (data$CO2))

kH_CO2<-29* exp( 2400 * ((1/(22+273.15)) - 1/(298.15)) ) 
CO2_mol_water<- data$CO2  / kH_CO2
CO2_mol_air<- data$CO2/ ((22+273.15)*0.08026/0.9)
CO2_conc<-(CO2_mol_water*70 + CO2_mol_air*70)/100 

CO2sat<-2400/kH_CO2

CO2calc2<-function(x) {
  x$kH_CO2<-29* exp( 2400 * ((1/(x$equiltemp+273.15)) - 1/(298.15)) ) #temp correction for Henry's law
  x$CO2_mol_water<- x$CO2_vol_air  / x$kH_CO2
  x$CO2_mol_air<- x$CO2_vol_air/ ((x$equiltemp+273.15)*0.08026/x$pressure)
  
  x$CO2_conc<-(x$CO2_mol_water*x$volwater + x$CO2_mol_air*x$volair)/100 # µmol/L or mmol/m3
  return(x)
  
}

khstream<-29* exp( 2400 * ((1/(15+273.15)) - 1/(298.15)) ) 
temp<-seq(1:25)
29* exp( 2400 * ((1/(temp+273.15)) - 1/(298.15)) ) 

3.4e-2* exp( 2400 * ((1/(temp+273.15)) - 1/(298.15)) ) 
co2sat<-400e-6*0.9 * 3.4e-2* exp( 2400 * ((1/(temp+273.15)) - 1/(298.15)) ) 

## Need to set certain columns to desired values if all the same
## If each sample is different, simple set the column to a column specified in the metadata file
names(sum_file)

csat<-16
co2flux<-(CO2_conc-csat)*13*0.3  # flux in mmol m-2 d-1
(CO2_conc-csat)*13*0.3*12*0.001 #flux in g m-2 d-1
#ER is 9 g O2 0r
9*1000/32 #mmol O or C

co2time<-  as_datetime(  data$MinFrom0*60+ as.numeric(ymd_hms("2019-08-10 10:15:00")) )

sum_file$CO2_vol_air <- sum_file$CO2
sum_file$equiltemp <- sum_file$eq_temp
sum_file$volair <- sum_file$air_mL
sum_file$volwater <- sum_file$H2O_mL
sum_file$pressure <- 0.9

