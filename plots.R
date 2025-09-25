
#----Load Packages-----

library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)
library(wesanderson)

data_doc<- read.csv("./data_doc.csv")


##################################################################
# Figure 1: DOC Uptake and DIC Production Fluxes By Sampling Site
# Dependencies: Run entire "modeling_rates_stan.R" script
##################################################################

#---Join necessary dataframes

data_doc <- data_doc %>% 
  left_join(select(KQ, code, q), by="code")

#---Calculate DOC uptake and DIC production fluxes at each time point

data_doc$DIC_flux = data_doc$DIC13 * data_doc$q
data_doc$DOC_flux = data_doc$removed_c * data_doc$q

#---Integrate DOC taken up and CO2 released at each site

sum_c<- data_doc %>% group_by(code) %>% summarize(dicsum=(nutint(MinFrom0,DIC13,0)), docsum=(nutint(MinFrom0,removed_c,0)))
sum_c<- left_join(sum_c, KQ, by="code")
sum_c$doc_flux<- sum_c$docsum*sum_c$q
sum_c$dic_flux<- sum_c$dicsum*sum_c$q
sum_c$ratio<-sum_c$dic_flux/sum_c$doc_flux
sum_c$ratio<-round(sum_c$ratio, 2)

#---Dictate order by sampling site for facet_wrap (this is for the breakthrough curves)

data_doc$order = factor(data_doc$code, 
                        levels=c('b_g1_up','b_g1_down','b_g2_up','b_g2_down','b_l_up','b_l_down'),
                        labels=c(expression("Glucose"["1,up"]),expression("Glucose"["1,down"]), 
                                 expression("Glucose"["2,up"]),expression("Glucose"["2,down"]), 
                                 expression("Leachate"["up"]),expression("Leachate"["down"])))

data_doc <- data_doc %>%
  arrange(order)

#---Dictate order by sampling site for facet_wrap (this is for the summary labels)

sum_c$order = factor(sum_c$code, 
                     levels=c('b_g1_up','b_g1_down','b_g2_up','b_g2_down','b_l_up','b_l_down'),
                     labels=c(expression("Glucose"["1,up"]),expression("Glucose"["1,down"]), 
                              expression("Glucose"["2,up"]),expression("Glucose"["2,down"]), 
                              expression("Leachate"["up"]),expression("Leachate"["down"])))

sum_c <- sum_c %>%
  arrange(order)

#---Generate summary labels for each sampling site

sum_c$doc_flux_label <- c("DOC = 713 mg", "DOC = 1069 mg", 
                          "DOC = 885 mg","DOC = 1115 mg", 
                          "DOC = 614 mg", "DOC = 478 mg")

sum_c$dic_flux_label <- c("DIC = 196 mg", "DIC = 416 mg", 
                          "DIC = 437 mg","DIC = 493 mg", 
                          "DIC = 513 mg", "DIC = 983 mg")

sum_c$ratio_label <- c("Ratio = 0.28", "Ratio = 0.39", 
                       "Ratio = 0.49","Ratio = 0.44", 
                       "Ratio = 0.84", "Ratio = 2.05")

#---Generate Figure 1 (DOC uptake and resulting DIC production from each tracer experiment)

fig1plot <- ggplot(data = data_doc)+ theme_bw()+
  geom_point(aes(x=MinFrom0, y=DOC_flux), fill = "#F98400", size =3, shape = 21, alpha=0.8)+
  geom_point(aes(x=MinFrom0, y=DIC_flux), fill = "#005AB5", size =3, shape = 21, alpha=0.8)+
  scale_x_log10(limits = c(2, 2000))+
  scale_y_continuous(limits = c(-2.02, 14.16))+
  facet_wrap(~order,  ncol=2, labeller = label_parsed)+
  labs(x = expression("Time from Addition Start (min)"), y = expression(""^13*"C Flux (mg/min)"))+
  geom_text(data = sum_c, aes(x=10, y=12, label = doc_flux_label, size = 10), color = "#F98400")+
  geom_text(data = sum_c, aes(x=10, y=10, label = dic_flux_label,  size = 10), color = "#005AB5")+
  geom_text(data = sum_c, aes(x=10, y=8, label = ratio_label,  size = 10))+
  theme(plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 10), 
        strip.text = element_text(color = "black", size=12))
fig1plot  

ggsave(fig1plot, filename = "doc_dic.pdf", width = 6.5, height = 6.5, units = "in", dpi = 300)


####################################################################
# Figure 2: Immediate, Delayed, and Total Respiration Model Results
# Dependencies: Run entire "modeling_rates_stan.R" script
####################################################################

#---Build function to generate plot 

d<-pred_steps[pred_steps$group==1,]

ploterr<-function(d, flux) {
  alpha<-0.05
  data<- d[d$step==1,]
  plot(data$MinFrom0,data$pred_tot, col=alpha("#005AB5",alpha), type="l",  xlab="", ylab="", xlim=c(5,1500), ylim=c(0,3), log="x")
  lines(data$MinFrom0,data$pred_delayed, col=alpha("#DC267F",alpha))
  lines(data$MinFrom0,data$pred_imm, col=alpha("#F98400",alpha))
  
  for (i in 2:200) { 
    data<- d[d$step==i,]
    
    lines(data$MinFrom0,data$pred_tot, col=alpha("#005AB5",alpha))
    lines(data$MinFrom0,data$pred_delayed, col=alpha("#DC267F",alpha))
    lines(data$MinFrom0,data$pred_imm, col=alpha("#F98400",alpha))
  }
  points(flux$MinFrom0,flux$DIC_flux, pch=21)
}

#---Generate Figure 2 (Immediate, Delayed, and Total Respiration Modeling Results) 

par(mfrow=(c(3,2)), mai=c(0.3,0.35,0.03,0.03), mgp=c(1.8,1,0), omi=c(0.4,0.4,0.1,0.1))
ploterr(d=pred_steps[pred_steps$group==1,],flux=data_doc[data_doc$code=="b_g1_up",])
text(expression("Glucose"["1,up"]), x=12, y=2.5, cex = 1.1)
ploterr(d=pred_steps[pred_steps$group==2,],flux=data_doc[data_doc$code=="b_g1_down",])
text(expression("Glucose"["1,down"]), x=13.5, y=2.5, cex = 1.1)
legend(5, 1.5, legend=c("Immediate", "Delayed", "Total"),
       col=c( "#F98400", "#DC267F","#005AB5"), lty=1, cex=1, box.lty=0)
ploterr(d=pred_steps[pred_steps$group==3,],flux=data_doc[data_doc$code=="b_g2_up",])
text(expression("Glucose"["2,up"]), x=12, y=2.5, cex = 1.1)
ploterr(d=pred_steps[pred_steps$group==4,],flux=data_doc[data_doc$code=="b_g2_down",])
text(expression("Glucose"["2,down"]), x=13.5, y=2.5, cex = 1.1)
ploterr(d=pred_steps[pred_steps$group==5,],flux=data_doc[data_doc$code=="b_l_up",])
text(expression("Leachate"["up"]), x=12, y=2.5, cex = 1.1)
ploterr(d=pred_steps[pred_steps$group==6,],flux=data_doc[data_doc$code=="b_l_down",])
text(expression("Leachate"["down"]), x=13.5, y=2.5, cex = 1.1)

mtext("Time from start of addition (min)", side = 1, outer = T)
mtext(expression(paste(''^13*"C DIC flux (mg/min)",sep="")), side = 2, outer = T)


########################################################################################################################
# Figure 3: Joint Parameter Distribution of Immediate Fraction of DOC Respired and Residence Time of Delayed Respiration
# Dependencies: Run entire "modeling_rates_stan.R" script
########################################################################################################################

#---Generate Figure 3 (Joint Parameter Distribution of Immediate Fraction of DOC Respired and Residence Time of Delayed Respiration) 

fig3plot <- ggplot(doc_steps) +
  geom_point( aes(x=1/kc, y=fi, color = as.factor(i)),alpha=0.4, size =0.1) + ylim(0,0.8) +xlim(0,1)+
  theme_classic()+
  scale_colour_manual(values=c("#EBCC2A","#F2AD00","#D69C4E", "#F98400", "#ABDDDE","#5BBCD6"), name=NULL, breaks=NULL, labels = NULL)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  labs(x = expression(paste(italic(1/k[c])," (d)")), y = expression(italic(f[i])))+
  geom_point(data=med_doc_steps, aes(x=1/kc, y=fi, fill = as.factor(i)), size =2.5, shape = 21, stroke = 1)+
  scale_fill_manual(values=c("#EBCC2A","#F2AD00","#D69C4E", "#F98400", "#ABDDDE","#5BBCD6"),name="",
                    labels=c(expression("Glucose"["1,up"]),expression("Glucose"["1,down"]), 
                             expression("Glucose"["2,up"]),expression("Glucose"["2,down"]), 
                             expression("Leachate"["up"]),expression("Leachate"["down"])))+ 
  guides(fill=guide_legend(title=expression(''^13~C~Addition)))+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 10))
fig3plot

ggsave(fig3plot, filename = "param.pdf", width = 5, height = 4, units = "in", dpi = 300)


########################################################################################################################
# Figure 4: CO2 Emission Fluxes and C vs. O Departure Plots
# Dependencies: Run entire "modeling_rates_stan.R", "data_process.R", "Creston_metabolism.R", "CO2 model.R" scripts
########################################################################################################################

#---OPTIONAL: Build Figure 4 dataframe after you have run all dependent scripts 

#fig4 <- data.frame(dtime = pco2m$dtime,
#bela_co2_flux = dat$bela_co2_flux,
#C.mod.flux = pco2m$C.mod.flux, 
#NEP_gC = NEP_gC, 
#DIC_pred = DIC_pred, 
#DIC_sat = DIC_sat, 
#oxy = blaine65_trim$oxy, 
#oxysat = blaine65_trim$oxysat, 
#bela_co2_conc = dat$bela_co2_conc,
#bela_co2_sat = csat(dat$bela_temp, bp))

#write.csv(fig4, file = "fig4.csv")

#---Load in Figure 4 data

fig4 <- read.csv("fig4.csv")
fig4$dtime <- as.POSIXct(fig4$dtime, format = "%Y-%m-%d %H:%M:%S")


#---Generate Figure 4 (CO2 Emission Fluxes and C vs. O Departure Plots)

fig4a<-ggplot()+
  theme_bw()+
  coord_cartesian(xlim = c(as.POSIXct("2019-08-10 04:00:00"), as.POSIXct("2019-08-17 0:00:00")), ylim = c(-7,5))+
  geom_point(data = fig4, aes(x = dtime, y = bela_co2_flux), size = 2, color = "#ae5c00", fill = "#fb9b31", shape = 21, alpha = 0.7)+
  geom_line(data = fig4, aes(x = dtime, y = -NEP_gC), color = "#006435", size = 1)+
  geom_line(data = fig4, aes(x = dtime, y = C.mod.flux), , color = "#643500", size = 1)+
  geom_line(data = fig4, aes(x = dtime, y = -NEP_gC_EQ1.7), color = "#00A357", size = 1, linetype = "dashed")+
  geom_line(data = fig4, aes(x = dtime, y = C.mod.flux_EQ1.7), , color = "#A35700", size = 1, linetype = "dashed")+
  geom_hline(yintercept = 0, color = "black")+
  annotate("text", x = as.POSIXct("2019-08-16 16:00:00"), y = 4.8, label = (expression(bold(CO['2,obs']))), color = "#F98400", size = 4)+
  annotate("text", x = as.POSIXct("2019-08-16 13:00:00"), y = -3, label = (expression(bold('-NEP'))), color = "#006435", size = 4)+
  annotate("text", x = as.POSIXct("2019-08-16 16:00:00"), y = 2.2, label = (expression(bold(CO['2,pred']))), color = "#643500", size = 4)+
  labs(y = expression(CO['2']~~Emission~Flux~~'(g'~C~m^-2~d^-1*')'))+ 
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 10))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, 
                                    linetype="solid"))
fig4a

fig4b<-ggplot()+
  theme_bw()+
  coord_cartesian(xlim=c(-0.05,0.35), ylim=c(-0.2,0.2))+
  geom_point(data = fig4, aes(x = (DIC_pred-DIC_sat), y = (oxy-oxysat)/32), size = 2, shape = 21, color = "#005AB5", fill = "#669cd3", alpha = 0.7)+
  geom_point(data = fig4, aes(x = (bela_co2_conc/12)-bela_co2_sat, y = (oxy-oxysat)/32), size = 2, shape = 21, color = "#ae5c00", fill = "#fb9b31", , alpha = 0.7)+
  geom_segment(data = fig4, aes(x = 0.17, y = 0.08, xend = 0.35, yend = -0.1), size = 1, color = "black")+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  annotate("text", x = 0.05, y = 0.188, label = (expression(bold(CO['2']))), color = "#F98400", size = 4)+
  annotate("text", x = 0.22, y = 0.19, label = (expression(bold(DIC))), color = "#005AB5", size = 4)+
  labs(x = expression(C~departure~'(mmol/L)'), y = expression(O['2']~departure~'(mmol/L)'))+ 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 10), 
        legend.position = c(0.05, -0.15))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, 
                                    linetype="solid"))
fig4b

fig4plot<-fig4a + fig4b + plot_layout(widths = c(2, 1))
fig4plot

fig4plot <- fig4plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 14, face = "bold", vjust = -1))
fig4plot

ggsave(fig4plot, filename = "co2_model.pdf", width = 9, height = 4, units = "in", dpi = 300)




#---Generate Supplemental Figure with EQ = 1.7 (EQ as in Diamond et al., 2025, L&O)

fig4$NEP_gC_EQ1.7 <- fig4$NEP_gC * 0.58
fig4$C.mod.flux_EQ1.7 <- fig4$C.mod.flux * 0.58

suppfig<-ggplot()+
  theme_bw()+
  coord_cartesian(xlim = c(as.POSIXct("2019-08-10 04:00:00"), as.POSIXct("2019-08-17 0:00:00")), ylim = c(-7,5))+
  geom_point(data = fig4, aes(x = dtime, y = bela_co2_flux), size = 2, color = "#ae5c00", fill = "#fb9b31", shape = 21, alpha = 0.7)+
  geom_line(data = fig4, aes(x = dtime, y = -NEP_gC, linetype = "EQ = 1.0"), color = "#006435", size = 0.75)+
  geom_line(data = fig4, aes(x = dtime, y = C.mod.flux, linetype = "EQ = 1.0"), color = "#643500", size = 0.75)+
  geom_line(data = fig4, aes(x = dtime, y = -NEP_gC_EQ1.7, linetype = "EQ = 1.7"), color = "#00A357", size = 0.75)+
  geom_line(data = fig4, aes(x = dtime, y = C.mod.flux_EQ1.7, linetype = "EQ = 1.7"), color = "#A35700", size = 0.75)+
  geom_hline(yintercept = 0, color = "black")+
  annotate("text", x = as.POSIXct("2019-08-16 16:00:00"), y = 4.8, label = (expression(bold(CO['2,obs']))), color = "#F98400", size = 4)+
  annotate("text", x = as.POSIXct("2019-08-16 13:00:00"), y = -3, label = (expression(bold('-NEP'))), color = "#006435", size = 4)+
  annotate("text", x = as.POSIXct("2019-08-16 16:00:00"), y = 2.2, label = (expression(bold(CO['2,pred']))), color = "#643500", size = 4)+
  scale_linetype_manual(values = c("EQ = 1.0" = "solid", "EQ = 1.7" = "dashed")) +
  guides(linetype = guide_legend(override.aes = list(color = "black", size = 0.75), keywidth = unit(2, "lines"))) +
  labs(y = expression(CO['2']~~Emission~Flux~~'(g'~C~m^-2~d^-1*')'))+ 
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 10),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 10),
        legend.position = c(0.9, 0),
        legend.justification = c(0.5, 0),
        legend.background = element_blank(), 
        legend.spacing.x = unit(1, "cm"))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, 
                                    linetype="solid"))
suppfig

ggsave(suppfig, filename = "supp_EQcompare.pdf", width = 6.5, height = 4, units = "in", dpi = 300)
