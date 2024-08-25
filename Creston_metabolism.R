

library(streamMetabolizer)
library(lubridate)



# the below function estimates bp (in mm Hg) based on altitude and local sealevel corrected barmetric pressure.  
#In my case I assummed the same bp for each time step.  If you have actual bp dat you can merge later
bpcalc<- function(bpst, alt) {
  bpst*25.4*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+15)))
}
#end of function

#this function claculate oxygen saturation based on Garcia and Gordon 1995.  Alternatively you can use a function in streamMetabolizer.
#BP is units of absolute pressure in mm Hg
osat<- function(temp, bp) {
  
  tstd<-log((298.15-temp) / (273.15 + temp))
  
  a0<-2.00907
  a1<-3.22014
  a2<-4.0501
  a3<-4.94457
  a4<- -0.256847
  a5<- 3.88767
  
  u<-10^(8.10765-(1750.286/(235+temp)))
  
  sato<-(exp(a0 + a1*tstd + a2*tstd^2 + a3*tstd^3 + a4*tstd^4+a5*tstd^5))*((bp-u)/(760-u))*1.42905
  sato
}

bpcalc(29.92,2900/3.28)

crest_press<- bpcalc(29.92,2900/3.28)/760 #(atm)

osat<- function(temp, bp) {
  
  tstd<-log((298.15-temp) / (273.15 + temp))
  
  a0<-2.00907
  a1<-3.22014
  a2<-4.0501
  a3<-4.94457
  a4<- -0.256847
  a5<- 3.88767
  
  u<-10^(8.10765-(1750.286/(235+temp)))
  
  sato<-(exp(a0 + a1*tstd + a2*tstd^2 + a3*tstd^3 + a4*tstd^4+a5*tstd^5))*((bp-u)/(760-u))*1.42905
  sato
}

blaine_65<-read.csv("./blaine_65_oxy.csv") #this is a chopped down PME file from sensor "George" located at 61 m
blaine_65$oxysat<- osat(blaine_65$temp, bpcalc(29.92,2900/3.28))
blaine_65$dtime<- as_datetime(blaine_65$unixtime)
blaine_65$local_time<-with_tz(blaine_65$dtime, tz="America/Denver")

plot(blaine_65$dtime, blaine_65$oxy/blaine_65$oxysat)


blaine_65$solar.time<-convert_UTC_to_solartime(blaine_65$dtime, longitude= -114.133, time.type="mean solar")

blaine_65$light<- calc_light(blaine_65$solar.time, latitude=48.186, longitude=-114.133, max.PAR =2326, attach.units = F)


#MOdle with so called `normal' pooling becuase discharge did not change`
blaine_name <- mm_name(type='bayes', pool_K600='normal', err_obs_iid = T, err_proc_iid =T)
blaine_specs <- specs(blaine_name, K600_daily_meanlog_meanlog=2, K600_daily_meanlog_sdlog=0.7, K600_daily_sdlog_sigma=0.05, burnin_steps=1000, 
                  saved_steps=1000)


blaine_65_sm<-data.frame(DO.obs=blaine_65$oxy, DO.sat=blaine_65$oxysat, 
                             temp.water=blaine_65$temp, depth=rep(0.28,length(blaine_65$temp)), 
                             light=blaine_65$light, solar.time=blaine_65$solar.time)


#blaine_65_fit <- metab(blaine_specs, data=blaine_65_sm, info=c(site='Blaine65', source='Bob Hall'))

#save(blaine_65_fit, file="./blaine_65_fit.RData")

load(file="./blaine_65_fit.RData")

plot_DO_preds(predict_DO(blaine_65_fit))

plot_metab_preds(predict_metab(blaine_65_fit))

blaine_params<-get_params(blaine_65_fit , uncertainty='ci')
write.csv(blaine_params, "blaine_params.csv")



########

blaine_params<-read.csv( "blaine_params.csv")

plot(blaine_params$K600.daily, blaine_params$ER.daily, ylim=c(-15,0), xlim=c(0,16), pch=16)

plot(blaine_params$GPP.daily, blaine_params$ER.daily, ylim=c(-15,0), xlim=c(0,16), pch=16)


mean(blaine_params$ER.daily, na.rm=T)
mean(blaine_params$GPP.daily, na.rm=T)
mean(blaine_params$K600.daily, na.rm=T)

sd(blaine_params$ER.daily, na.rm=T)
sd(blaine_params$GPP.daily, na.rm=T)
sd(blaine_params$K600.daily, na.rm=T)
k<-seq(1:30)
plot(k, dlnorm(k,2,0.7)) #checking out my prior here on gas exchnage.  You will want to set your own.

Blaine_hr<- mean(blaine_params$ER.daily, na.rm=T) + 0.45*mean(blaine_params$GPP.daily, na.rm=T)

plot(blaine_65$solar.time[100:244], 100*blaine_65$oxy[100:244]/blaine_65$oxysat[100:244],
     pch=16, col="blue", xlab="Solar time", ylab="% oxygen saturation")


###
##Estimate time specific NEP based on oxygen for the 6 days

lightinaday<-blaine_65$light[100:990]/(sum(blaine_65$light[100:990])/6  )

NEP<- mean(blaine_params$GPP.daily, na.rm=T)*lightinaday +  mean(blaine_params$ER.daily, na.rm=T)*(10/1440)




NEP_mmol<- (1440/10)*NEP*1000/32

NEP_gC<- (1440/10)*NEP*12/32
plot(NEP_gC)

plot(blaine_65$solar.time[100:990], -NEP_gC, type="l", ylim=c(-300,600), xlab="Solar time", ylab= "-NEP mmol m-2 d-1" )
points (co2time,co2flux, pch=16, col="dark green")


plot(blaine_65$local_time[100:990], -NEP_gC, type="l", ylim=c(-4,8), xlab="Solar time", ylab= "-NEP mmol m-2 d-1" )



###Big switch to CO2 here
##CO2 emission

pco2<-read.csv("./pco2_blaine.csv")
pco2$dtime<-mdy_hms(pco2$time, tz="America/Denver")

plot(pco2$dtime, pco2$white_co2)
points (pco2$dtime, pco2$bela_co2, col="red")


plot(pco2$dtime, pco2$white_temp, pch=16, cex=0.8)
points(pco2$dtime, pco2$bela_temp,pch=16, cex=0.8, col="red")

#Henryys law mol/ (l *atm)
KH.CO2 <- function(tempC){
  tempK <- tempC + 273.15
  KH.CO2 <- exp( -58.0931 + (90.5069*(100/tempK)) +22.2940*log((tempK/100), base=exp(1)) )
  KH.CO2
}


KCO2fromK600<- function (temp,K600) {
  K600/(600/(1742-(temp*91.24)+(2.208*temp^2)-(0.0219*temp^3)))^-0.5
}

pco2$bela_co2_conc<- 0.001*12*pco2$bela_co2*(1/crest_press)*KH.CO2(pco2$bela_temp) #mg C/L
#units for the untrusting
# mg/ug * 12 ug/umol * uatm/L * umol / (L * utam)
pco2$satconc<-crest_press*0.001*12*400*KH.CO2(pco2$bela_temp)



k_blaine<- 13.5* 0.195 #K from sM, 0.195 is z Units m/d

pco2$bela_co2_flux<- KCO2fromK600(pco2$bela_temp,k_blaine) *(pco2$bela_co2_conc-pco2$satconc)

plot(pco2$dtime,pco2$bela_co2_flux, ylim=c(0,5))
lines(blaine_65$local_time[100:990], -NEP_gC)
pco2$neplight<- calc_light(force_tz(pco2$dtime, "UTC"), latitude=48.186, longitude=-114.133, max.PAR =2326, attach.units = F)


#NEP<- mean(blaine_params$GPP.daily, na.rm=T)*neplight/sum +  mean(blaine_params$ER.daily, na.rm=T)*(15/1440)


(sum(pco2$bela_co2_flux)*10/1440)/7  # 3.1 g C/d CO2 emission


par(mai=c(0.7,0.7,0.1,0.1), mgp=c(2,1,0))
plot(pco2$dtime,pco2$bela_co2_flux, ylim=c(-7,5), ylab="Flux gC m-2 d-1", xlab="", pch=16, cex=0.8)
lines(blaine_65$local_time[100:990], -NEP_gC, col="blue")
lines(pco2$dtime,rep(0, length(pco2$dtime)))
lines(pco2m$dtime,pco2m$C.mod.flux, col="red", lwd=1.8)



                                                     

                                                     