mcwdTypes <- c("CWD_Pr","WDfixed_T","WDaet_T","WDaet2_T",
              "WDpet_T", "WDpet2_T", "WDpet_T*", "WDpet2_T*",
              "CWD_Cr", "WDfixed_C", "WDaet_C", "WDaet2_C",	
              "WDpet_C", "WDpet2_C", "WDpet_C*", "WDpet2_C*")

#absDeviations FROM MEAN (in mm)----
absDeviations <- anomalies %>%
                    select(year, lon, lat,
                           "CWD_Pr"	= aet_pet_r,			
                           "WDfixed_T" = preTrmm_et100_r,	
                           "WDaet_T"		= preTrmm_aet_r,	
                           "WDaet2_T"	= preTrmm_aetAvg_r,
                           "WDpet_T"		= preTrmm_pet_r,		
                           "WDpet2_T"	= preTrmm_petAvg_r,		
                           "WDpet_T*"	= preTrmm_pet_nr,	
                           "WDpet2_T*"	= preTrmm_petAvg_nr,
                           "CWD_Cr"		= aet_petCru_r,		
                           "WDfixed_C"	= preCru_et100_r,	
                           "WDaet_C"		= preCru_aet_r,		
                           "WDaet2_C"	= preCru_aetAvg_r,	
                           "WDpet_C"		= preCru_petCru_r,		
                           "WDpet2_C"	= preCru_petCruAvg_r,	
                           "WDpet_C*"	= preCru_petCru_nr,	
                           "WDpet2_C*"	= preCru_petCruAvg_nr)
absDeviations <- gather(absDeviations, key = "mcwd", value = "deviation", 4:ncol(absDeviations))

absDeviations$mcwd <- factor(absDeviations$mcwd, mcwdTypes)

deviationLevels <- c(-200, -150, -100, -50, -25, 0)
deviationLabels <- c("< -200", "< -150", "< -100", "< -50", "< -25", "< 0", "> 0")
deviationColors <- rev(c("greenyellow", "yellow", "orange", "indianred2", "red", "firebrick4"))

absDeviations <- absDeviations %>%
  ungroup() %>%
  mutate(deviationFromMean = ifelse(deviation < deviationLevels[1], deviationLabels[1],
                                    ifelse(deviation < deviationLevels[2], deviationLabels[2],
                                           ifelse(deviation < deviationLevels[3], deviationLabels[3],
                                                  ifelse(deviation < deviationLevels[4], deviationLabels[4],
                                                         ifelse(deviation < deviationLevels[5], deviationLabels[5],
                                                                ifelse(deviation < deviationLevels[6], deviationLabels[6], "> 0"))))))) %>%
  filter(deviationFromMean != "> 0")

absDeviations$deviationFromMean <- factor(absDeviations$deviationFromMean,  deviationLabels)

absDeviations %>%
  regmap(region = "ama") +
  geom_tile(aes(fill = deviationFromMean)) +
  scale_fill_manual(values = deviationColors, name = "Annual MCWD anomaly in mm") +
  facet_grid(mcwd ~ year) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "left",
        legend.direction = "vertical",
        legend.key.height = unit(2, "cm"),
        legend.text = element_text(size = 14, angle = 90),
        legend.title = element_text(size = 14, angle = 90),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0))

## ANOMALIES (speiCru, speiTrmm, spiTrmm) ----
anomaliesSdBackUp <- anomaliesSd %>%
  select(year, lon, lat,
         "CWD_Pr"	= aet_pet_r,			
         "WDfixed_T" = preTrmm_et100_r,	
         "WDaet_T"		= preTrmm_aet_r,	
         "WDaet2_T"	= preTrmm_aetAvg_r,
         "WDpet_T"		= preTrmm_pet_r,		
         "WDpet2_T"	= preTrmm_petAvg_r,		
         "WDpet_T*"	= preTrmm_pet_nr,	
         "WDpet2_T*"	= preTrmm_petAvg_nr,
         "CWD_Cr"		= aet_petCru_r,		
         "WDfixed_C"	= preCru_et100_r,	
         "WDaet_C"		= preCru_aet_r,		
         "WDaet2_C"	= preCru_aetAvg_r,	
         "WDpet_C"		= preCru_petCru_r,		
         "WDpet2_C"	= preCru_petCruAvg_r,	
         "WDpet_C*"	= preCru_petCru_nr,	
         "WDpet2_C*"	= preCru_petCruAvg_nr)

#anomalies correlation plot---
corrplot(cor.fk(
  as.data.frame(
    anomaliesSdBackUp %>% select(-year, -lon, -lat))),
  method = 'color',
  type = 'lower',
  addCoef.col = "black",
  diag = FALSE,
  tl.col="black",
  tl.srt=45,
  tl.cex = 1,
  cl.cex = 1,
  cl.lim = c(0, 1),
  col=colorRampPalette(c("firebrick4", "red", "pink",
                         "lightsteelblue1","lightsteelblue1","lightsteelblue", "lightsteelblue3"))(200),
  number.cex = 1) #Text label color and rotation

#AnomaliesSd maps
anomaliesSdBackUp <- gather(anomaliesSdBackUp, key = "mcwd", value = "anomalySd", 4:ncol(anomaliesSdBackUp))

# #Anomalies distribution----
# distributionsPalette <- c("chartreuse", "mediumorchid3", "red1", "seagreen1",
#                           "cyan", "darkgoldenrod1", "navy", "dodgerblue3")
# anomaliesDistribution <- ggdensity(anomaliesSd, x = "anomalySd", color = "mcwd", alpha = 0,
#                                    linetype = "solid", palette = distributionsPalette, size = 1) +
#   labs(x = "Anomalies", y = "Density")
# anomaliesQqplot <- ggqqplot(anomaliesSd, x = "anomalySd", color = "mcwd", palette = distributionsPalette, size = 0.2) + ylim(-4,4)
# 
# ggarrange(anomaliesDistribution, anomaliesQqplot, nrow = 1, ncol = 2, common.legend = TRUE)

anomaliesSdBackUp$mcwd <- factor(anomaliesSdBackUp$mcwd, mcwdTypes)

anomalyLevels <- c(-2, -1.5, -1, 1, 1.5, 2)
anomalyLables <- c("< -2", "< -1.5", "< -1", "-1 to 1", "< 1", "< 1.5", "< 2")

anomaliesSdBackUp <- anomaliesSdBackUp %>%
  ungroup() %>%
  mutate(sdLevel = ifelse(anomalySd < anomalyLevels[1], anomalyLables[1],
                          ifelse(anomalySd < anomalyLevels[2], anomalyLables[2],
                                 ifelse(anomalySd < anomalyLevels[3], anomalyLables[3],
                                        ifelse(anomalySd <  anomalyLevels[4], anomalyLables[4],
                                               ifelse(anomalySd < anomalyLevels[5], anomalyLables[5],
                                                      ifelse(anomalyLevels < anomalyLevels[6], anomalyLables[6], anomalyLables[7])))))))
anomaliesSdBackUp$sdLevel <- factor(anomaliesSdBackUp$sdLevel, rev(anomalyLables))

anomaliesSdBackUp %>%
  regmap(region = "ama") +
  geom_tile(aes(fill = anomalySd)) +
  scale_fill_gradient2(low = "#FF3300", midpoint = 0, mid = "ivory", high = "#0000FF", name = "SD",
                       breaks = c(min(anomaliesSdBackUp$anomalySd), seq(-2, -1, 0.5), seq(1, 2, 0.5)),
                       labels = c(min(anomaliesSdBackUp$anomalySd), seq(-2, -1, 0.5), seq(1, 2, 0.5))) +
  facet_grid(mcwd ~ year) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "left",
        legend.direction = "vertical",
        legend.key.height = unit(3, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0))





