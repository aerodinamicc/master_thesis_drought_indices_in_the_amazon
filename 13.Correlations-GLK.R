#Correlations -----
indicesTimeSeries <- amazonWD[complete.cases(amazonWD),] %>%
  select(starts_with("spei"), starts_with("spi"), ends_with("monthlySd"))

indicesTimeSeries <- indicesTimeSeries %>%
  select(CWD_Pr = cwd_aet_pet_r_monthlySd,
         SPEI03_T = speiTrmm03,
         SPI03_T = spiTrmm03,
         "WDfixed_T" = cwd_preTrmm_et100_r_monthlySd,
         "WDaet_T" = cwd_preTrmm_aet_r_monthlySd,
         "WDaet2_T" = cwd_preTrmm_aetAvg_r_monthlySd,
         "WDpet_T" = cwd_preTrmm_pet_r_monthlySd,
         "WDpet2_T" = cwd_preTrmm_petAvg_r_monthlySd,
         "WDpet_T*" = cwd_preTrmm_pet_nr_monthlySd,
         "WDpet2_T*" = cwd_preTrmm_petAvg_nr_monthlySd,
         CWD_Cr = cwd_aet_petCru_r_monthlySd,
         SPEI03_C = speiCru03,
         SPI03_C = spiCru03,
         "WDfixed_C" = cwd_preCru_et100_r_monthlySd,
         "WDaet_C" = cwd_preCru_aet_r_monthlySd,
         "WDaet2_C" = cwd_preCru_aetAvg_r_monthlySd,
         "WDpet_C" = cwd_preCru_petCru_r_monthlySd,
         "WDpet2_C" = cwd_preCru_petCruAvg_r_monthlySd,
         "WDpet_C*" = cwd_preCru_petCru_nr_monthlySd,
         "WDpet2_C*" = cwd_preCru_petCruAvg_nr_monthlySd)

# #Normality check ----
# monthlyAnomaliesDistribution <- gather(indicesTimeSeries) %>%
#   select("Anomaly" = key, "SD" = value)
# 
# den <- ggdensity(monthlyAnomaliesDistribution, x = "SD", color = "Anomaly", alpha = 0,
#                  linetype = "solid", palette = c(distributionsPalette,"orange", "yellow", "blue"), size = 1) +
#   labs(x = "Anomalies", y = "Density")
# qq <- ggqqplot(monthlyAnomaliesDistribution, x = "SD", color = "Anomaly", palette = c(distributionsPalette,"orange", "yellow", "blue"), size = 0.2) + ylim(-4,4)
# ggarrange(den, qq, nrow = 1, ncol = 2, common.legend = TRUE)

# Nice visualization of correlations
corrplot.mixed(cor.fk(indicesTimeSeries), order = "AOE", addrect = 3,
               tl.cex = 1, tl.pos = "lt", cl.cex = 1, number.cex = 0.8)

#gleichlaufigkeit -----
source('glk.R')

gleichlaufigkeit <- glk(indicesTimeSeries)

corrplot(gleichlaufigkeit,
         is.corr = FALSE,
         method = 'color',
         cl.lim = c(min(gleichlaufigkeit),max(gleichlaufigkeit)),
         col=colorRampPalette(c("firebrick4", "red", "pink", "darkseagreen1","darkseagreen2","darkseagreen3", "darkseagreen"))(50),
         type = 'lower',
         addCoef.col = "black",
         diag = FALSE,
         tl.col="black", tl.srt=45, tl.cex = 1, cl.cex = 1, number.cex = 0.8) #Text label color and rotation
