#Time series of drought categories with all indices
droughtCategories <- c("Extremely dry/ < -1.6", "Severely dry/ < -1.3", "Moderately dry/ < -0.8", "Abnormally dry/ < -0.5")
droughtLevels <- c(-1.6, -1.3, -0.8, -0.5)

colNames = c("speiTrmm03", "cwd_preTrmm_et100_r_monthlySd", "cwd_preTrmm_aet_r_monthlySd", "cwd_preTrmm_pet_r_monthlySd", "cwd_preTrmm_pet_nr_monthlySd",
                 "spiTrmm03", "cwd_aet_pet_r_monthlySd","cwd_preTrmm_aetAvg_r_monthlySd", "cwd_preTrmm_petAvg_r_monthlySd", "cwd_preTrmm_petAvg_nr_monthlySd",
                 "speiCru03", "cwd_preCru_et100_r_monthlySd", "cwd_preCru_aet_r_monthlySd",  "cwd_preCru_petCru_r_monthlySd", "cwd_preCru_petCru_nr_monthlySd", 
                 "spiCru03", "cwd_aet_petCru_r_monthlySd",  "cwd_preCru_aetAvg_r_monthlySd",  "cwd_preCru_petCruAvg_r_monthlySd", "cwd_preCru_petCruAvg_nr_monthlySd")

titles = c("SPEI03_T", "WDfixed_T", "WDaet_T",  "WDpet_T",  "WDpet_T*",
           "SPI03_T",  "CWD_Pr",  "WDaet2_T", "WDpet2_T", "WDpet2_T*",
           "SPEI03_C", "WDfixed_C", "WDaet_C", "WDpet_C", "WDpet_C*",       
           "SPI03_C", "CWD_Cr", "WDaet2_C", "WDpet2_C", "WDpet2_C*")

droughtLevelsFun <- function(colName = "", title = ""){
  definePalette <- "Reds"
  topHline <- list("firebrick", 2, "dashed", 50)
  lowHline <- list("chocolate1", 2, "dashed", 20)
  
  severityLevels <- amazonWD[complete.cases(amazonWD),] %>%
    dplyr::select(year, month, lon, lat, indexName = !!as.name(colName)) %>%
    mutate(date = make_date(year, month, 01),
           severityLevel =          ifelse(indexName < droughtLevels[1], droughtCategories[1],
                                           ifelse(indexName < droughtLevels[2], droughtCategories[2],
                                                  ifelse(indexName < droughtLevels[3], droughtCategories[3],
                                                         ifelse(indexName < droughtLevels[4], droughtCategories[4], "Near normal or wet"))))) %>%
    group_by(date, severityLevel) %>%
    summarise(percentage = round((n()/1945)*100), 0) %>%
    filter(severityLevel != "Near normal or wet")
  
  severityLevels$severityLevel <- factor(severityLevels$severityLevel, droughtCategories)
  
  sev <- expand.grid(date = unique(severityLevels$date), severityLevel = unique(severityLevels$severityLevel))
  
  severityLevels <- full_join(severityLevels, sev, by = c("date" = "date", "severityLevel" = "severityLevel"))
  
  severityLevels$percentage[which(is.na(severityLevels$percentage))] <- 0
  # severityLevels <- severityLevels %>%
  #   mutate(percentage = ifelse(is.na(percentage), 0, percentage))
  
  severityLevelsPlot <- ggplot(severityLevels, aes(x=date, y=percentage, fill=severityLevel, color)) +
    geom_area(colour="black", size=.2, alpha=.4) +
    scale_fill_brewer(palette=definePalette, direction = -1, breaks=rev(levels(severityLevels$severityLevel))) +
    ggtitle(title) +
    geom_hline(yintercept=topHline[[4]], linetype=topHline[[3]], size=topHline[[2]], color = topHline[[1]]) +
    labs(y = "Drought extent (in %)") +
    #geom_hline(yintercept=lowHline[[4]], linetype=lowHline[[3]], size=lowHline[[2]], color = lowHline[[1]]) +
    theme(axis.title.y.left = element_text(size = 12),
          axis.title.x.bottom = element_blank(),
          axis.text.x.bottom = element_text(size = 18),
          axis.text.y.left = element_text(size = 18),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.width = unit(2, "cm"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 16),
          title = element_text(size = 24),
          plot.title = element_text(size=14),
          strip.text = element_text(size = 18))
  
  return(severityLevelsPlot)
}
  
  plots <- list()
  
  for (varIndex in 1:length(colNames)) {
    var = colNames[varIndex]
    title = titles[varIndex]
    plot = droughtLevelsFun(var, title)
    plots[[varIndex]] <- eval(substitute(plot))
  }
  
do.call(ggarrange, c(plotlist = plots, nrow = 4, ncol = 5, common.legend = TRUE, legend = "bottom"))

#Drought extent (in %) during the 2005 dry season (Jul-Sep)-------
guideLine60 <- list(60, "dotted", 2, "orange")
guideLine50 <- list(50, "dotted", 2, "yellow")
guideLine40 <- list(40, "dotted", 2, "white")

extent2005 <- amazonWD[complete.cases(amazonWD),] %>%
  filter(year == 2005 & month %in% 7:9) %>%
  select(starts_with("spei"), starts_with("spi"), ends_with("monthlySd"), -ends_with("nr_monthlySd")) %>%
  gather() %>%
  mutate(source = ifelse(grepl("cwd_aet", key, fixed = TRUE), "AET-PET",
                        ifelse(grepl("Cru", key, fixed = TRUE), "CRU","TRMM")),
         severityLevel = ifelse(value < droughtLevels[1], droughtCategories[1],
                               ifelse(value < droughtLevels[2], droughtCategories[2],
                                      ifelse(value < droughtLevels[3], droughtCategories[3],
                                             ifelse(value < droughtLevels[4], droughtCategories[4], "Near normal or wet"))))) %>%
  filter(severityLevel != "Near normal or wet") %>% 
  group_by(key, severityLevel, source) %>%
  summarise(percentage = round((n()/1945)*100 /3,0)) %>%
  group_by(key) %>%
  mutate(sumPercentages = sum(percentage)) %>%
  ungroup() %>%
  mutate(labelPercentage = ifelse(percentage > 2, paste(percentage, "%"), ""),
         severityLevel = factor(severityLevel, droughtCategories))

  extent2005$key[which(extent2005$key == "cwd_aet_petCru_r_monthlySd")] <- "CWD_Cr"
  extent2005$key[which(extent2005$key == "cwd_aet_pet_r_monthlySd")] <- "CWD_Pr"
  extent2005$key[which(extent2005$key == "speiCru03" | extent2005$key == "speiTrmm03")] <- "SPEI03"
  extent2005$key[which(extent2005$key == "spiCru03" | extent2005$key == "spiTrmm03")] <- "SPI03"
  extent2005$key[which(extent2005$key == "cwd_preCru_aet_r_monthlySd" | extent2005$key == "cwd_preTrmm_aet_r_monthlySd")] <- "WDaet"
  extent2005$key[which(extent2005$key == "cwd_preCru_aetAvg_r_monthlySd" | extent2005$key == "cwd_preTrmm_aetAvg_r_monthlySd")] <- "WDaet2"
  extent2005$key[which(extent2005$key == "cwd_preCru_petCru_r_monthlySd" | extent2005$key == "cwd_preTrmm_pet_r_monthlySd")] <- "WDpet"
  extent2005$key[which(extent2005$key == "cwd_preCru_petCruAvg_r_monthlySd" | extent2005$key == "cwd_preTrmm_petAvg_r_monthlySd")] <- "WDpet2"
  extent2005$key[which(extent2005$key == "cwd_preCru_et100_r_monthlySd" | extent2005$key == "cwd_preTrmm_et100_r_monthlySd")]  <- "WDfixed"
  
 extent2005 %>%
  filter(source != "AET-PET") %>%
  ggplot(aes(x=key, y=percentage, fill=severityLevel)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label = labelPercentage), position = position_stack(vjust = 0.5), size = 7) +
  xlab("Indices (SPEI/SPI/MCWD)\n") +
  ylab("Drought extent for the 2005 dry season - JAS (in %)") +
  scale_fill_brewer(palette = "Reds", direction = -1) +
  theme(axis.text.x.bottom = element_text(size = 16, hjust = 1),
        axis.text.y.left = element_text(size = 16),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y.left = element_text(size = 18),
        title = element_text(size = 28),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(2.5, "cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 0, hjust = 0),
        strip.text.x = element_text(size = 28)) +
  # geom_hline(yintercept=guideLine60[[1]], linetype=guideLine60[[2]], size=guideLine60[[3]], color = guideLine60[[4]]) +
  # geom_hline(yintercept=guideLine50[[1]], linetype=guideLine50[[2]], size=guideLine50[[3]], color = guideLine50[[4]]) +
  # geom_hline(yintercept=guideLine40[[1]], linetype=guideLine40[[2]], size=guideLine40[[3]], color = guideLine40[[4]]) +
  coord_flip() +
  facet_wrap(source~.)
 
 extent2005 %>%
   filter(source == "AET-PET") %>%
   ggplot(aes(x=key, y=percentage, fill=severityLevel)) + 
   geom_bar(stat="identity") +
   geom_text(aes(label = labelPercentage), position = position_stack(vjust = 0.5), size = 9) +
   xlab("Indices (MCWD)") +
   ylab("Drought extent for the 2005 dry season - JAS (in %)") +
   guides(fill=guide_legend(nrow = 2, byrow = TRUE)) +
   scale_fill_brewer(palette = "Reds", direction = -1) +
   theme(axis.text.x.bottom = element_text(size = 12, hjust = 0.5),
         axis.text.y.left = element_text(size = 12),
         axis.title.x.bottom = element_text(size = 14),
         axis.title.y.left = element_text(size = 14),
         legend.position = "bottom",
         legend.direction = "horizontal",
         legend.key.width = unit(1, "cm"),
         legend.text = element_text(size = 12),
         legend.title = element_text(size = 0, hjust = 0),
         title = element_text(size = 28)) #+
   #geom_hline(yintercept=guideLine60[[1]], linetype=guideLine60[[2]], size=guideLine60[[3]], color = guideLine60[[4]]) +
   #geom_hline(yintercept=guideLine40[[1]], linetype=guideLine40[[2]], size=guideLine40[[3]], color = guideLine40[[4]])

#absDeviations (in mm)
 severityExtent <- absDeviations %>%
   filter(year == 2005 & deviationFromMean != "> 0" & !grepl("*",mcwd, fixed = TRUE)) %>%
   mutate(source = ifelse(grepl("CWD_", mcwd, fixed = TRUE), "CWD",
                          ifelse(grepl("_C", mcwd, fixed = TRUE), "CRU","TRMM"))) %>%
   group_by(mcwd, deviationFromMean, source) %>%
   summarise(percentage = round((n()/1945)*100, 0)) %>%
   filter(percentage != 0)

 severityExtent$mcwd <- as.character(severityExtent$mcwd)

severityExtent$mcwd[which(severityExtent$mcwd == "aet_pet_r")] <- "CWD_Pr"
severityExtent$mcwd[which(severityExtent$mcwd == "aet_petCru_r")] <- "CWD_Cr"
severityExtent$mcwd[which(severityExtent$mcwd == "WDfixed_T" | severityExtent$mcwd == "WDfixed_C")] <- "WDfixed"
severityExtent$mcwd[which(severityExtent$mcwd == "WDpet_T" | severityExtent$mcwd == "WDpet_C")] <- "WDpet"
severityExtent$mcwd[which(severityExtent$mcwd == "WDpet2_T" | severityExtent$mcwd == "WDpet2_C")] <- "WDpet2"
severityExtent$mcwd[which(severityExtent$mcwd == "WDaet_T" | severityExtent$mcwd == "WDaet_C")] <- "WDaet"
severityExtent$mcwd[which(severityExtent$mcwd == "WDaet2_T" | severityExtent$mcwd == "WDaet2_C")] <- "WDaet2"

severityExtent %>%
  filter(source != "CWD") %>%
  mutate(labelPercentage = ifelse(percentage > 3, paste(percentage, "%"), "")) %>%
  ggplot(aes(x=mcwd, y=percentage, fill=deviationFromMean)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label = labelPercentage), position = position_stack(vjust = 0.5), size = 5) +
  xlab("Annual MCWD values for 2005\n") +
  ylab("Drought extent (in %)") +
  scale_fill_manual(values = deviationColors[1:length(deviationColors)]) +
  theme(
        axis.text.x.bottom = element_text(size = 18),
        axis.text.y.left = element_text(size = 18),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(3, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 24)) +
  # geom_hline(yintercept=guideLine60[[1]], linetype=guideLine60[[2]], size=guideLine60[[3]], color = guideLine60[[4]]) +
  # geom_hline(yintercept=guideLine50[[1]], linetype=guideLine50[[2]], size=guideLine50[[3]], color = guideLine50[[4]]) +
  # geom_hline(yintercept=guideLine40[[1]], linetype=guideLine40[[2]], size=guideLine40[[3]], color = guideLine40[[4]]) +
  coord_flip() +
  facet_wrap(source~.)

severityExtent %>%
  filter(source == "CWD") %>%
  mutate(labelPercentage = ifelse(percentage > 3, paste(percentage, "%"), "")) %>%
  ggplot(aes(x=mcwd, y=percentage, fill=deviationFromMean)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label = labelPercentage), position = position_stack(vjust = 0.5), size = 7) +
  xlab("\nAnnual MCWD values for 2005") +
  ylab("Drought extent (in %)") +
  scale_fill_manual(values = deviationColors[2:length(deviationColors)]) +
  theme(axis.text.x.bottom = element_text(size = 14),
    axis.text.y.left = element_text(size = 14),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 0, hjust = 0),
    strip.text = element_text(size = 24)) # +
  # geom_hline(yintercept=60, linetype=guideLine40[[2]], size=guideLine40[[3]], color = guideLine40[[4]])

#Temporal development of absDeviations ----

  absDeviationLevelsFun <- function(wdName = "", title = ""){
    definePalette <- deviationColors
    hline <- list("firebrick", 2, "dashed", 40)
    colorScheme <- c(deviationColors[1:6])
    
    #the next condition takes into account that the time series
    #in quastion have not registered any pixels of the highest categories
    if (wdName == "CWD_Pr") {
      colorScheme <- c(deviationColors[2:6])
    } else if (wdName == "CWD_Cr"){
      colorScheme <- c(deviationColors[3:6])
    }
    
    dep <-  absDeviations %>%
      filter(mcwd == wdName) %>%
      filter(deviationFromMean != "> 0") %>%
      group_by(mcwd, deviationFromMean, year) %>%
      summarise(percentage = round((n()/1945)*100, 0))
    
    d <- expand.grid(year = unique(dep$year), deviationFromMean = unique(dep$deviationFromMean))
    dep <- full_join(d, dep, by = c("year" = "year", "deviationFromMean" = "deviationFromMean")) %>%
      mutate(percentage = ifelse(is.na(percentage), 0, percentage),
             mcwd = wdName) %>%
      arrange(year, deviationFromMean)
    
    
    dep <- dep %>%
      filter(mcwd == wdName) %>%
      ggplot(aes(x=year, y=percentage, fill=deviationFromMean)) +
      geom_area(stat = "identity", colour="black", size=.2, alpha=.4) +
      scale_fill_manual(values=colorScheme) +
      ggtitle(title) +
      labs(y = "Drought extent (in %)") +
      guides(fill=guide_legend(title="Deviation from\n mean (in mm)")) + 
      #geom_hline(yintercept=hline[[4]], linetype=hline[[3]], size=hline[[2]], color = hline[[1]]) +
      theme(axis.title.y.left = element_text(size = 12),
            axis.title.x.bottom = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.width = unit(3, "cm"),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 18))
    
    return(dep)
  }

  titlesAnnual <- c("WDfixed_T", "WDaet_T", "WDpet_T", "WDpet_T*",   
                    "CWD_Pr", "WDaet2_T", "WDpet2_T",  "WDpet2_T*",
                    "WDfixed_C", "WDaet_C", "WDpet_C", "WDpet_C*",
                    "CWD_Cr", "WDaet2_C", "WDpet2_C", "WDpet2_C*")
  
  
  absDeviationsPlots <- list()
  for (wdIndex in 1:length(titlesAnnual)) {
    var = as.character(titlesAnnual[wdIndex])
    plot = absDeviationLevelsFun(var, var)
    absDeviationsPlots[[wdIndex]] <- eval(substitute(plot))
  }
  
  do.call(ggarrange, c(plotlist = absDeviationsPlots, nrow = 4, ncol = 4, common.legend = TRUE, legend = "bottom"))
  