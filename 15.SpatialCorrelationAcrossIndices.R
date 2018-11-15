#Spatial correlation function. Accomadates both monthly and annual anomalies as the @dataset is defining-----
spatialCorrelation <- function(dataset, var1 = "", var2 = "", title = "", period =  1998:2005, ylabTitle = "")
{
  dataset %>%
    filter(year %in% period) %>%
    select(lon, lat, !!as.name(var1), !!as.name(var2))  %>%
    filter(complete.cases(.)) %>%
    group_by(lon,lat) %>%
    summarise(spatialCor = cor(!!as.name(var1), !!as.name(var2), method = "kendall")) %>%
    regmap(region = "ama") +
    geom_tile(aes(fill = spatialCor)) +
    scale_fill_gradientn(name = "Kendall tau", 
                         colours = c("tomato3", "lightgoldenrod1", "royalblue1"),
                         values = scales::rescale(c(-1, -0.6, 0, 0.6, 1))) +
    labs(title = title, y = ylabTitle) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.key.height = unit(1, "cm"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          strip.text = element_text(size = 18))
}

#Plotting of monthly anomalies of MCWD computed with TRMM against SPEI and SPI------
SCPlots <- list() #spatial correlation plots

SCPairs <- list(c("cwd_preTrmm_et100_r_monthlySd", "cwd_preTrmm_petAvg_r_monthlySd", "WDfixed_T vs WDpet_T"),
                c("cwd_preTrmm_et100_r_monthlySd", "cwd_preTrmm_aetAvg_r_monthlySd", "WDfixed_T vs WDaet_T"),
                c("cwd_preTrmm_et100_r_monthlySd", "cwd_aet_pet_r_monthlySd", "WDfixed_T vs CWD_Pr"),
                c("cwd_preTrmm_et100_r_monthlySd", "spiTrmm03", "WDfixed_T vs SPI03_T"),
                c("cwd_preTrmm_et100_r_monthlySd", "speiTrmm03", "WDfixed_T vs SPEI03_T"),
                c("spiTrmm03", "speiTrmm03", "SPI03_T vs SPEI03_T"))

for (pairIndex in 1:length(SCPairs)) {
  var1 = SCPairs[[pairIndex]][1]
  var2 = SCPairs[[pairIndex]][2]
  title = SCPairs[[pairIndex]][3]
  
  plot = spatialCorrelation(amazonWD, var1, var2, title)
  SCPlots[[pairIndex]] <- eval(substitute(plot))
}

do.call(ggarrange, c(plotlist = SCPlots, nrow = 2, ncol = 3, common.legend = TRUE, legend = "right"))

#Plotting of monthly anomalies of MCWD/SPEI/SPI computed with TRMM against those computed with CRU-------

pl <- list()
SCPairsTrCru <- list(c("speiTrmm03", "speiCru03", "SPEI", period = 1999:2004),
                     c("cwd_preCru_petCru_r_monthlySd", "cwd_preTrmm_pet_r_monthlySd", "WDpet", period = 1999:2004),
                     c("cwd_preCru_aet_r_monthlySd", "cwd_preTrmm_aet_r_monthlySd", "WDaet", period = 1999:2004),
                     c("cwd_preCru_et100_r_monthlySd", "cwd_preTrmm_et100_r_monthlySd", "WDfixed", period = 1999:2004),
                     c("speiTrmm03", "speiCru03", "SPEI", period = c(1998, 2005)),
                     c("cwd_preCru_petCru_r_monthlySd", "cwd_preTrmm_pet_r_monthlySd", "WDpet", period = c(1998, 2005)),
                     c("cwd_preCru_aet_r_monthlySd", "cwd_preTrmm_aet_r_monthlySd", "WDaet", period = c(1998, 2005)),
                     c("cwd_preCru_et100_r_monthlySd", "cwd_preTrmm_et100_r_monthlySd", "WDfixed", period = c(1998, 2005)))

for (pairIndex in 1:length(SCPairsTrCru)) {
  var1 = SCPairsTrCru[[pairIndex]][1]
  var2 = SCPairsTrCru[[pairIndex]][2]
  title = SCPairsTrCru[[pairIndex]][3]
  period = SCPairsTrCru[[pairIndex]][4]
  
  ylabTitle <- ""
  if (pairIndex == 1) {
    ylabTitle = "1999:2004 correlation"
  } else if (pairIndex == 5) {
    ylabTitle = "1998 and 2005 correlation"
  }
  
  plot = spatialCorrelation(amazonWD, var1, var2, title, period, ylabTitle)
  pl[[pairIndex]] <- eval(substitute(plot))
}

do.call(ggarrange, c(plotlist = pl, nrow = 2, ncol = 4, common.legend = TRUE, legend = "right"))


#Annual MCWD anomalies correlated between each other------
#Plotting of annual anomalies of MCWD computed with TRMM against each other------
p <- list() #spatial correlation plots

pairsTrmm <- list(c("preTrmm_et100_r", "preTrmm_pet_r", "WDfixed_T vs WDpet_T"),
              c("preTrmm_et100_r", "preTrmm_aet_r", "WDfixed_T vs WDaet_T"),
              c("aet_pet_r", "preTrmm_et100_r", "WDfixed_T vs CWD_Pr"),
              c("preTrmm_et100_r", "preTrmm_petAvg_r", "WDfixed_T vs WDpet2_T"),
              c("preTrmm_et100_r", "preTrmm_aetAvg_r", "WDfixed_T vs WDaet2_T"),
              c("preTrmm_et100_r", "preTrmm_pet_nr", "WDfixed_T vs WDpet_T*"))

for (pairIndex in 1:length(pairsTrmm)) {
  var1 = pairsTrmm[[pairIndex]][1]
  var2 = pairsTrmm[[pairIndex]][2]
  title = pairsTrmm[[pairIndex]][3]
  
  plot = spatialCorrelation(anomaliesSd, var1, var2, title)
  p[[pairIndex]] <- eval(substitute(plot))
}

do.call(ggarrange, c(plotlist = p, nrow = 2, ncol = 3, common.legend = TRUE, legend = "right"))

##Plotting of annual anomalies of MCWD computed with TRMM and CRU against each other------
pTrCr <- list()
pairsTrmmVsCru <- list(c("aet_pet_r", "aet_petCru_r", "CWD"),
                       c("preTrmm_aet_r", "preCru_aet_r", "WDaet"),
                       c("preTrmm_pet_r", "preCru_petCru_r", "WDpet"),
                       c("preTrmm_et100_r", "preCru_et100_r", "WDfixed"))

for (pairIndex in 1:length(pairsTrmmVsCru)) {
  var1 = pairsTrmmVsCru[[pairIndex]][1]
  var2 = pairsTrmmVsCru[[pairIndex]][2]
  title = pairsTrmmVsCru[[pairIndex]][3]
  
  plot = spatialCorrelation(anomaliesSd, var1, var2, title)
  pTrCr[[pairIndex]] <- eval(substitute(plot))
}

do.call(ggarrange, c(plotlist = pTrCr, nrow = 1, ncol = 4, common.legend = TRUE, legend = "right"))

