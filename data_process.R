library(dplyr)
library(lubridate)
library(zoo)

source("./doc_functions.R")

site_info<-read.csv("./site_data.csv")
site_info$condmass<-site_info$salt_mass*2100
site_info$add_ratio<- site_info$doc_mass_13C_mg/site_info$condmass


#more site data transform
site_info$velocity<- site_info$length/site_info$travel_time

#some means
mean(site_info$alkalinity[1:6])




#########Possible to skip running this code because the file it derives is saved.
#### Blaine 8 August down

b_g1_down_info<-site_info[site_info$code=="b_g1_down",]

blaine0808_solute<-read.csv("./DOC_runs/13CAdditionData_CrestonGlucose_20190808_analyzed.csv")
blaine0808_co2<-read.csv("./picarro_data/2019_08_08_Blaine/2019_08_08_Blaine.csv")
blaine0808<-left_join(blaine0808_solute,blaine0808_co2, by="Syringe")
blaine0808$dtime <- dmy_hms(paste(blaine0808$Date, blaine0808$SampleTime), tz="America/Denver")
blaine0808$MinFrom0 <- (as.numeric(blaine0808$dtime) - as.numeric(dmy_hms("08/08/19 10:11:00", tz="America/Denver")))/60


blaine0808_down<- blaine0808[blaine0808$SiteID=="Down", ]
data<-blaine0808_down[!is.na(blaine0808_down$CO2), ]  #delete rows with NA and call it data


##SpC and hydrology
bkg<-mean( data$SpC.us.cm_1[1])
plot(data$MinFrom0, (data$SpC.us.cm_1-bkg))


spcint<- nutint(data$MinFrom0, data$SpC.us.cm_1, bkg=bkg )
data_Q_b_g1_down<- b_g1_down_info$condmass / spcint  #171 L/min




##CO2
plot(data$MinFrom0, (data$delCO2))
bkg_delCO2<-mean((data$delCO2[1:2]))
data_co2AF_corr<-deltoAF(data$delCO2)-deltoAF(bkg_delCO2)
data_13DIC<-b_g1_down_info$alkalinity*13*data_co2AF_corr #mole to mass here
plot(data$MinFrom0, data_13DIC, pch=16, col="red", ylab= "13C mass of extra DIC", 
     xlab="Time from solute addition (min)")

int13co2<- nutint(data$MinFrom0, data_13DIC, bkg=0 )
int13co2*data_Q

##DOC
delbkg<-mean(data$Delta.13C[1])
plot(data$MinFrom0, (data$Delta.13C))
AF_corr<-deltoAF(data$Delta.13C)-deltoAF(delbkg)
data_13C<-data$DOC_conc*AF_corr

removed_c<- b_g1_down_info$add_ratio*(data$SpC.us.cm_1-bkg)-data_13C


b_g1_down<-data.frame(code=rep("b_g1_down", length(data$MinFrom0)), MinFrom0 = data$MinFrom0, SpC_corr=data$SpC.us.cm_1-bkg, DOC_conc=data$DOC_conc, Delta.13C= data$Delta.13C,
                      AF_corr=AF_corr,doc_13C=data_13C, delCO2=data$delCO2, co2AF_corr=data_co2AF_corr, DIC13=data_13DIC, removed_c=removed_c)




###Blaine 8 aug up



b_g1_up_info<-site_info[site_info$code=="b_g1_up",]

blaine0808_solute<-read.csv("./DOC_runs/13CAdditionData_CrestonGlucose_20190808_analyzed.csv")
blaine0808_co2<-read.csv("./picarro_data/2019_08_08_Blaine/2019_08_08_Blaine.csv")
blaine0808<-left_join(blaine0808_solute,blaine0808_co2, by="Syringe")
blaine0808$dtime <- dmy_hms(paste(blaine0808$Date, blaine0808$SampleTime), tz="America/Denver")
blaine0808$MinFrom0 <- (as.numeric(blaine0808$dtime) - as.numeric(dmy_hms("08/08/19 10:11:00", tz="America/Denver")))/60


blaine0808_up<- blaine0808[blaine0808$SiteID=="Up", ]
data<-blaine0808_up[!is.na(blaine0808_up$CO2), ]  #delete rows with NA and call it data


##SpC and hydrology
data$SpC.us.cm_1[1]<-data$SpC.us.cm_1[2]  #mising first spc, sub in second one
bkg<-mean( data$SpC.us.cm_1[1:4])
plot(data$MinFrom0, (data$SpC.us.cm_1-bkg))


spcint<- nutint(data$MinFrom0, data$SpC.us.cm_1, bkg=bkg )
data_Q_b_g1_down<- b_g1_down_info$condmass / spcint 




##CO2
plot(data$MinFrom0, (data$delCO2))
bkg_delCO2<-mean((data$delCO2[1:4]))
data_co2AF_corr<-deltoAF(data$delCO2)-deltoAF(bkg_delCO2)
data_13DIC<-b_g1_up_info$alkalinity*13*data_co2AF_corr #mole to mass here
plot(data$MinFrom0, data_13DIC, pch=16, col="red", ylab= "13C mass of extra DIC", 
     xlab="Time from solute addition (min)")

int13co2<- nutint(data$MinFrom0, data_13DIC, bkg=0 )
int13co2*data_Q

##DOC
delbkg<-mean(data$Delta.13C[1:3])
plot(data$MinFrom0, (data$Delta.13C))
AF_corr<-deltoAF(data$Delta.13C)-deltoAF(delbkg)
data_13C<-data$DOC_conc*AF_corr

removed_c<- b_g1_up_info$add_ratio*(data$SpC.us.cm_1-bkg)-data_13C



b_g1_up<-data.frame(code=rep("b_g1_up", length(data$MinFrom0)), MinFrom0 = data$MinFrom0, SpC_corr=data$SpC.us.cm_1-bkg, DOC_conc=data$DOC_conc, Delta.13C= data$Delta.13C,
                      AF_corr=AF_corr,doc_13C=data_13C, delCO2=data$delCO2, co2AF_corr=data_co2AF_corr, DIC13=data_13DIC, removed_c=removed_c)



#####
#####
###Blaine 15 aug up  'glucose 2'



b_g2_up_info<-site_info[site_info$code=="b_g2_up",]

blaine0815<-read.csv("./DOC_runs/13CAdditionData_CrestonGlucose_20190815_analyzed.csv")

blaine0815$dtime <- mdy_hms(paste(blaine0815$SampleDate, blaine0815$SampleTime), tz="America/Denver")
blaine0815$MinFrom0 <- (as.numeric(blaine0815$dtime) - as.numeric(mdy_hms("08/15/19 09:00:00", tz="America/Denver")))/60

blaine0815_up<- blaine0815[blaine0815$SiteID=="Up", ]

plot(blaine0815_up$MinFrom0,blaine0815_up$SpC.us.cm, xlim=c(0,90))

data<-blaine0815_up[!is.na(blaine0815_up$CO2), ] 




##SpC and hydrology

bkg<-mean( data$SpC.us.cm[1:4])
plot(data$MinFrom0, (data$SpC.us.cm-bkg), xlim=c(0,600))
data$SpC.us.cm[30]<-678 ##this low conductivity had to be error.  Interpolated
##slant backgrouns, need to correct

slantslope<-(data$SpC.us.cm[29]-bkg) / data$MinFrom0[29]

condcorr<-(data$SpC.us.cm[1:28]-bkg)- slantslope*data$MinFrom0[1:28]
condcorr2<- (data$SpC.us.cm[29:37]-bkg)- slantslope*data$MinFrom0[28]
condcorr<-c(condcorr,condcorr2)

plot(data$MinFrom0, condcorr, xlim=c(0,600))

spcint<- nutint(data$MinFrom0, condcorr, bkg=0 )
data_Q_b_g2_up<- b_g2_up_info$condmass / spcint 




##CO2
plot(data$MinFrom0, (data$delCO2))
bkg_delCO2<-mean((data$delCO2[1:4]))
data_co2AF_corr<-deltoAF(data$delCO2)-deltoAF(bkg_delCO2)
data_13DIC<-b_g2_up_info$alkalinity*13*data_co2AF_corr #mole to mass here
plot(data$MinFrom0, data_13DIC, pch=16, col="red", ylab= "13C mass of extra DIC", 
     xlab="Time from solute addition (min)")

int13co2<- nutint(data$MinFrom0, data_13DIC, bkg=0 )


##DOC
delbkg<-mean(data$Delta.13C[1:4])
plot(data$MinFrom0, (data$Delta.13C))
AF_corr<-deltoAF(data$Delta.13C)-deltoAF(delbkg)

data_13C<-data$DOC_conc*AF_corr

removed_c<- b_g2_up_info$add_ratio*condcorr-data_13C
removed_c[34:37]<-c(0,0,0,0)
removed_c<- na.approx(removed_c)

plot(data$MinFrom0, removed_c)

b_g2_up<-data.frame(code=rep("b_g2_up", length(data$MinFrom0)),MinFrom0 = data$MinFrom0, SpC_corr=condcorr, DOC_conc=data$DOC_conc, Delta.13C= data$Delta.13C,
                      AF_corr=AF_corr,doc_13C=data_13C, delCO2=data$delCO2, co2AF_corr=data_co2AF_corr, DIC13=data_13DIC, removed_c=removed_c)



####
###Blaine G2 down
####

b_g2_down_info<-site_info[site_info$code=="b_g2_down",]

blaine0815<-read.csv("./DOC_runs/13CAdditionData_CrestonGlucose_20190815_analyzed.csv")

blaine0815$dtime <- mdy_hms(paste(blaine0815$SampleDate, blaine0815$SampleTime), tz="America/Denver")
blaine0815$MinFrom0 <- (as.numeric(blaine0815$dtime) - as.numeric(mdy_hms("08/15/19 09:00:00", tz="America/Denver")))/60


blaine0815_down<- blaine0815[blaine0815$SiteID=="Down", ]

plot(blaine0815_down$MinFrom0,blaine0815_down$SpC.us.cm, xlim=c(0,300))

data<-blaine0815_down[!is.na(blaine0815_down$CO2), ] 




##SpC and hydrology

bkg<-mean( data$SpC.us.cm[1:6])
plot(data$MinFrom0, (data$SpC.us.cm-bkg), xlim=c(0,600))

##slant backgrouns, need to correct

slantslope<-(data$SpC.us.cm[39]-bkg) / data$MinFrom0[39]

condcorr<-(data$SpC.us.cm-bkg)- slantslope*data$MinFrom0



plot(data$MinFrom0, condcorr, xlim=c(0,600))

spcint<- nutint(data$MinFrom0, condcorr, bkg=0 )
data_Q_b_g2_down<- b_g2_down_info$condmass / spcint 




##CO2
plot(data$MinFrom0, (data$delCO2))
bkg_delCO2<-mean((data$delCO2[1:5]))
data_co2AF_corr<-deltoAF(data$delCO2)-deltoAF(bkg_delCO2)
data_13DIC<-b_g2_down_info$alkalinity*13*data_co2AF_corr #mole to mass here
plot(data$MinFrom0, data_13DIC, pch=16, col="red", ylab= "13C mass of extra DIC", 
     xlab="Time from solute addition (min)")

int13co2<- nutint(data$MinFrom0, data_13DIC, bkg=0 )


##DOC
delbkg<-mean(data$Delta.13C[1:2])
plot(data$MinFrom0, (data$Delta.13C))
AF_corr<-deltoAF(data$Delta.13C)-deltoAF(delbkg)

data_13C<-data$DOC_conc*AF_corr
data_13C[38] <- -1.5e-3  #interpolate here


removed_c<- b_g2_down_info$add_ratio*condcorr-data_13C
removed_c[41:45]<-c(0,0,0,0,0)


plot(data$MinFrom0, removed_c)

b_g2_down<-data.frame(code=rep("b_g2_down", length(data$MinFrom0)), MinFrom0 = data$MinFrom0, SpC_corr=condcorr, DOC_conc=data$DOC_conc, Delta.13C= data$Delta.13C,
                    AF_corr=AF_corr,doc_13C=data_13C, delCO2=data$delCO2, co2AF_corr=data_co2AF_corr, DIC13=data_13DIC, removed_c=removed_c)



###Blaine 9 Aug leachate up

b_l_up_info<-site_info[site_info$code=="b_l_up",]

blaine0809_solute<-read.csv("./DOC_runs/13CAdditionData_CrestonLeachate_20190809_analyzed.csv")
blaine0809_co2<-read.csv("./picarro_data/2019_08_09_Blaine/2019_08_09_Blaine.csv")

blaine0809<-left_join(blaine0809_solute,blaine0809_co2, by="Syringe")

blaine0809$dtime <- mdy_hms(paste(blaine0809$Date, blaine0809$SampleTime), tz="America/Denver")
blaine0809$MinFrom0 <- (as.numeric(blaine0809$dtime) - as.numeric(mdy_hms("08/09/19 10:15:00", tz="America/Denver")))/60

##here below call it data
#use only down for now
blaine0809_up<- blaine0809[blaine0809$SiteID=="Up", ]


data<-blaine0809_up  #delete rows with NA and call it data



##SpC and hydrology
data$SpC.us.cm_1[45]<-data$SpC.us.cm_1[44]  #mising last spc, sub in penultimate one
bkg<-mean( data$SpC.us.cm_1[1:8])
plot(data$MinFrom0, (data$SpC.us.cm_1-bkg))


slantslope<-(data$SpC.us.cm[42]-bkg) / data$MinFrom0[42]

condcorr1<-(data$SpC.us.cm[1:42]-bkg)- slantslope*data$MinFrom0[1:42]
condcorr2<- data$SpC.us.cm[43:45]-bkg
condcorr<- c(condcorr1,condcorr2)

spcint<- nutint(data$MinFrom0, condcorr, bkg=0 )
data_Q_b_l_up<- b_l_up_info$condmass / spcint 




##CO2
plot(data$MinFrom0, (data$delCO2))
bkg_delCO2<-mean((data$delCO2[1:5]))
data_co2AF_corr<-deltoAF(data$delCO2)-deltoAF(bkg_delCO2)
data_13DIC<-b_l_up_info$alkalinity*13*data_co2AF_corr #mole to mass here
plot(data$MinFrom0, data_13DIC, pch=16, col="red", ylab= "13C mass of extra DIC", 
     xlab="Time from solute addition (min)")

int13co2<- nutint(data$MinFrom0, data_13DIC, bkg=0 )


##DOC
data$Delta.13C<- na.approx(data$Delta.13C)  # missing 3 DOC:  interpolated 
delbkg<-mean(data$Delta.13C[1:5])
plot(data$MinFrom0, (data$Delta.13C))
AF_corr<-deltoAF(data$Delta.13C)-deltoAF(delbkg)
data_13C<-data$DOC_conc*AF_corr
data_13C<-na.approx(data_13C) ##missing 3, interpolate

removed_c<- b_l_up_info$add_ratio*(condcorr)-data_13C
plot(data$MinFrom0, condcorr)


b_l_up<-data.frame(code=rep("b_l_up", length(data$MinFrom0)), MinFrom0 = data$MinFrom0, SpC_corr=condcorr, DOC_conc=data$DOC_conc, Delta.13C= data$Delta.13C,
                    AF_corr=AF_corr,doc_13C=data_13C, delCO2=data$delCO2, co2AF_corr=data_co2AF_corr, DIC13=data_13DIC, removed_c=removed_c)


######
###Blaine 9 Aug leachate down

b_l_down_info<-site_info[site_info$code=="b_l_down",]

blaine0809_solute<-read.csv("./DOC_runs/13CAdditionData_CrestonLeachate_20190809_analyzed.csv")
blaine0809_co2<-read.csv("./picarro_data/2019_08_09_Blaine/2019_08_09_Blaine.csv")

blaine0809<-left_join(blaine0809_solute,blaine0809_co2, by="Syringe")

blaine0809$dtime <- mdy_hms(paste(blaine0809$Date, blaine0809$SampleTime), tz="America/Denver")
blaine0809$MinFrom0 <- (as.numeric(blaine0809$dtime) - as.numeric(mdy_hms("08/09/19 10:15:00", tz="America/Denver")))/60

##here below call it data

blaine0809_down<- blaine0809[blaine0809$SiteID=="Down", ]


data<-blaine0809_down  #delete rows with NA and call it data



##SpC and hydrology


##
plot(data$MinFrom0, data$SpC.us.cm_1)
data$SpC.us.cm_1[21:23]<-c(686, 683,678)
bkg<-mean( data$SpC.us.cm_1[1:2])
plot(data$MinFrom0[1:28], (data$SpC.us.cm_1[1:28]-bkg))

lm(c(0,-2.6)~c(70,345))

condcorr<- (data$SpC.us.cm_1[1:28] - (0.662-0.009455*data$MinFrom0[1:28]) ) -bkg
plot(data$MinFrom0[1:28], condcorr)

condcorr<-c(condcorr, c(0,0,0,0)) # make late time cond 0.  Won't affect integration

spcint<- nutint(data$MinFrom0, condcorr, bkg=0 )
data_Q_b_l_down<- b_l_down_info$condmass / spcint  #184 L/min

##


##CO2
plot(data$MinFrom0, (data$delCO2))
bkg_delCO2<-mean((data$delCO2[1:4]))
data_co2AF_corr<-deltoAF(data$delCO2)-deltoAF(bkg_delCO2)
data_13DIC<-b_l_down_info$alkalinity*13*data_co2AF_corr #mole to mass here
plot(data$MinFrom0, data_13DIC, pch=16, col="red", ylab= "13C mass of extra DIC", 
     xlab="Time from solute addition (min)")

int13co2<- nutint(data$MinFrom0, data_13DIC, bkg=0 )
int13co2*170


##DOC

delbkg<-mean(data$Delta.13C[1:2])
plot(data$MinFrom0, (data$Delta.13C))
AF_corr<-deltoAF(data$Delta.13C)-deltoAF(delbkg)
data_13C<-data$DOC_conc*AF_corr
plot(data$MinFrom0, data_13C)
plot(data$MinFrom0, b_l_down_info$add_ratio*(condcorr))

removed_c<- b_l_down_info$add_ratio*(condcorr)-data_13C
plot(data$MinFrom0, removed_c)

170*nutint(time=data$MinFrom0, data_13C, bkg=0)
170*nutint(time=data$MinFrom0, removed_c, bkg=0)

b_l_down<-data.frame(code=rep("b_l_down", length(data$MinFrom0)), MinFrom0 = data$MinFrom0, SpC_corr=condcorr, DOC_conc=data$DOC_conc, Delta.13C= data$Delta.13C,
                   AF_corr=AF_corr,doc_13C=data_13C, delCO2=data$delCO2, co2AF_corr=data_co2AF_corr, DIC13=data_13DIC, removed_c=removed_c)




#stick the 6 releases together

data_doc<- rbind(b_g1_up, b_g1_down, b_g2_up, b_g2_down, b_l_up, b_l_down)

write.csv(data_doc,"./data_doc.csv")



###
####



