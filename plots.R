

library(ggplot2)
library(scales)

library(wesanderson)

data_doc<- read.csv("./data_doc.csv")

#dic_pred<-read.csv("./dic_pred.csv")
#dic_params<-read.csv("./dic_params.csv")

###DOC uptake plot

data_doc <- data_doc %>% 
  left_join(select(KQ, code, q), by="code")

data_doc$DIC_flux = data_doc$DIC13 * data_doc$q


######
####### Integrate DOC taken up and CO2 released for a table
####
sum_c<- data_doc %>% group_by(code) %>% summarize(dicsum=(nutint(MinFrom0,DIC13,0)), docsum=(nutint(MinFrom0,removed_c,0)))
sum_c<- left_join(sum_c, KQ, by="code")
sum_c$doc_flux<- sum_c$docsum*sum_c$q
sum_c$dic_flux<- sum_c$dicsum*sum_c$q
sum_c$ratio<-sum_c$dic_flux/sum_c$doc_flux




##removed and produced plot

##datafram with labels


##way better plot

d<-pred_steps[pred_steps$group==1,]

ploterr<-function(d, flux) {
  alpha<-0.05
  data<- d[d$step==1,]
  plot(data$MinFrom0,data$pred_tot, col=alpha("purple",alpha), type="l",  xlab="", ylab="", xlim=c(5,1500), ylim=c(0,3), log="x")
  lines(data$MinFrom0,data$pred_delayed, col=alpha("red",alpha))
  lines(data$MinFrom0,data$pred_imm, col=alpha("blue",alpha))
  
  for (i in 2:200) { 
    data<- d[d$step==i,]
 
  lines(data$MinFrom0,data$pred_tot, col=alpha("purple",alpha))
  lines(data$MinFrom0,data$pred_delayed, col=alpha("red",alpha))
  lines(data$MinFrom0,data$pred_imm, col=alpha("blue",alpha))
  }
  points(flux$MinFrom0,flux$DIC_flux, pch=16)
}

#ploterr(d, flux=dic_pred[dic_pred$code=="b_g1_up",])

##the big model fit figure.  It needs some color help!
par(mfrow=(c(3,2)), mai=c(0.3,0.35,0.03,0.03), mgp=c(1.8,1,0), omi=c(0.4,0.4,0.1,0.1))
ploterr(d=pred_steps[pred_steps$group==1,],flux=data_doc[data_doc$code=="b_g1_up",])
text("g_1_up", x=10, y=2.5)
ploterr(d=pred_steps[pred_steps$group==2,],flux=data_doc[data_doc$code=="b_g1_down",])
text("g_1_down", x=10, y=2.5)
legend(5, 1.5, legend=c("Immediate", "Delayed", "Total"),
       col=c( "blue", "red","purple"), lty=1, cex=1, box.lty=0)
ploterr(d=pred_steps[pred_steps$group==3,],flux=data_doc[data_doc$code=="b_g2_up",])
text("g_2_up", x=10, y=2.5)
ploterr(d=pred_steps[pred_steps$group==4,],flux=data_doc[data_doc$code=="b_g2_down",])
text("g_2_down", x=10, y=2.5)
ploterr(d=pred_steps[pred_steps$group==5,],flux=data_doc[data_doc$code=="b_l_up",])
text("l_up", x=10, y=2.5)
ploterr(d=pred_steps[pred_steps$group==6,],flux=data_doc[data_doc$code=="b_l_down",])
text("l_down", x=10, y=2.5)

mtext("Time from start of addition (min)", side = 1, outer = T)
mtext("13C DIC flux (mg/min)", side = 2, outer = T)


###Joint distributions of parameters

ggplot( doc_steps) +
  geom_point( aes(x=1/kc, y=fi, color = as.factor(i)),alpha=0.1, size =0.1) + ylim(0,0.8) +xlim(0,1)+
  theme_classic()+
  scale_colour_manual(values=c("#EBCC2A","#F2AD00","#D69C4E", "#F98400", "#ABDDDE","#5BBCD6"),name="", breaks=c(1,2,3,4,5,6),
                      labels=c("g_1_up", "g_1_down", "g_2_up", "g_2_down", "l_up", "l_down"))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  labs(x = expression(paste(italic(1/k[c])," (d)")), y = expression(italic(f[i])))+
  geom_point(data=med_doc_steps, aes(x=1/kc, y=fi))


###CO2 flux plot (Run CO2 code in metabolism plus CO2 model to ge this to work)

pdf(file="figures/co2_model.pdf", height=3, width=6)

layout(matrix(c(1,1,1,2,2,1,1,1,2,2), nrow=2, byrow=TRUE) )

par(mai=c(0.5,0.5,0.1,0.1), mgp=c(2,1,0))
plot(pco2$dtime,pco2$bela_co2_flux, ylim=c(-7,5), ylab="Flux gC m-2 d-1", xlab="", pch=16, cex=0.8)
lines(blaine_65$local_time[100:990], -NEP_gC, col="blue")
lines(pco2$dtime,rep(0, length(pco2$dtime)))
lines(pco2m$dtime,pco2m$C.mod.flux, col="red", lwd=1.8)


plot(DIC_pred-DIC_sat, (blaine65_trim$oxy-blaine65_trim$oxysat)/32, xlim=c(-0.05,0.35), ylim=c(-0.2,0.2), pch=16, col="red" ,
     xlab="DIC departure (mmol/L)", ylab="O2 departure (mmol/L)")
lines(c(-0.05,0.35), c(0,0))
lines( c(0,0), c(-0.2,0.2))
lines( c(0.17,0.35), c(0.08,-0.1), lwd=2)
points( (dat$bela_co2_conc/12)-csat(dat$bela_temp, bp), (blaine65_trim$oxy-blaine65_trim$oxysat)/32, 
        col="lightblue" )

legend(0.05, -0.15, legend=c("CO2", "DIC"),
       col=c( "lightblue", "red"),  cex=1, pch=16,box.lty=0)

dev.off()







#####
###Old plotting code.  Will probbaly delete