#dry months count
drySeasonLengthPlot <- monthly_trmm %>%
  filter(year %in% 1998:2005) %>%
  group_by(lon, lat, year) %>%
  summarise(seasonLength = sum(monthly_prec < 100)) %>%
  regmap(region = "ama") +
  geom_tile(aes(x = lon, y = lat, fill = seasonLength)) +
  scale_fill_gradient2(low="darkgreen", midpoint = 5, mid = "yellow", high = "red", name="Dry months", 
                      labels = c(12:0),
                      breaks = c(12:0)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(1.6, "cm"))

annualPrecipitation <- monthly_trmm %>%
  filter(year %in% 1998:2005) %>%
  group_by(lon, lat, year) %>%
  summarise(annualPrec = sum(monthly_prec))

annualPrecipitation$labelPrec <- cut(
  annualPrecipitation$annualPrec,
  breaks = c(seq(0, 2800, by = 400), max(annualPrecipitation$annualPrec)),
  labels = c("0 - 400", "400-800", "800 - 1200", "1200 - 1600", "1600 - 2000", "2000 - 2400", "2400 - 2800", "2800 <"),
  right  = FALSE,
  include.lowest = TRUE)

annualPrecPlot <- annualPrecipitation  %>%
  regmap(region = "ama") +
  geom_tile(aes(x = lon, y = lat, fill = labelPrec)) +
  scale_fill_brewer(palette = "Blues", name = "Rainfall in mm/year")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(0.6, "cm"))
rm("annualPrecipitation")

cycleLine1 <- list("black", 1, "dotted", 50)
cycleLine2 <- list("black", 1, "dotted", 100)
cycleLine3 <- list("black", 1, "dotted", 150)
divisions <- c(-64.5, -7)
regions <- c("NW", "NE", "SW", "SE")

pPetRibbon <- amazonWD %>%
  mutate(region = ifelse(lon < divisions[1] & lat > divisions[2], regions[1],
                         ifelse(lon < divisions[1] & lat <= divisions[2], regions[3],
                                ifelse(lon >= divisions[1] & lat > divisions[2], regions[2],
                                       ifelse(lon >= divisions[1] & lat <= divisions[2], regions[4], NA))))) %>%
  group_by(region, month) %>%
  summarise(pet = mean(pet), pre = mean(preTrmm)) %>%
  ungroup() %>%
  mutate(date = as.Date(paste("2005", as.character(month), "01"), format = "%Y%m%d"),
         region = factor(region, regions)) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = pet), size = 2, color = "red") +
  geom_line(aes(y = pre), size = 2, color = "blue") +
  geom_ribbon(aes(ymin=pre,ymax=pet, fill= ifelse(pre < pet, TRUE, NA), alpha=0.6)) +
  scale_fill_manual("",values=c("red"))+
  scale_x_date(date_labels = "%b") +
  theme(axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),
        legend.position = "bottom") +
  geom_hline(yintercept=cycleLine1[[4]], linetype=cycleLine1[[3]], size=cycleLine1[[2]], color = cycleLine1[[1]]) +
  geom_hline(yintercept=cycleLine2[[4]], linetype=cycleLine2[[3]], size=cycleLine2[[2]], color = cycleLine2[[1]]) +
  geom_hline(yintercept=cycleLine3[[4]], linetype=cycleLine3[[3]], size=cycleLine3[[2]], color = cycleLine3[[1]]) +
  facet_wrap(region~.)

cwdRibbon <- amazonWD %>%
  mutate(region = ifelse(lon < divisions[1] & lat > divisions[2], regions[1],
                         ifelse(lon < divisions[1] & lat <= divisions[2], regions[3],
                                ifelse(lon >= divisions[1] & lat > divisions[2], regions[2],
                                       ifelse(lon >= divisions[1] & lat <= divisions[2], regions[4], NA))))) %>%
  group_by(region, month) %>%
  summarise(pet = mean(pet), aet = mean(aet)) %>%
  ungroup() %>%
  mutate(date = as.Date(paste("2005", as.character(month), "01"), format = "%Y%m%d"),
         region = factor(region, regions)) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = pet), size = 2, color = "red") +
  geom_line(aes(y = aet), size = 2, color = "springgreen3") +
  geom_ribbon(aes(ymin=aet,ymax=pet, fill= ifelse(pet > aet, TRUE, NA), alpha=0.6)) +
  scale_fill_manual("",values=c("red"))+
  scale_x_date(date_labels = "%b") +
  theme(axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),
        legend.position = "bottom") +
  geom_hline(yintercept=cycleLine2[[4]], linetype=cycleLine2[[3]], size=cycleLine2[[2]], color = cycleLine2[[1]]) +
  geom_hline(yintercept=cycleLine3[[4]], linetype=cycleLine3[[3]], size=cycleLine3[[2]], color = cycleLine3[[1]]) +
  facet_wrap(region~.)

grid.arrange(pPetRibbon, cwdRibbon, ncol=2)

grid.arrange(annualPrecPlot, drySeasonLengthPlot, ncol=2)
# #miGHT BE DELETED-----
# amazonWD %>%
#   filter(year == 2005) %>%
#   filter(complete.cases(.)) %>%
#   mutate(region = ifelse(lon < divisions[1] & lat > divisions[2], regions[1],
#                          ifelse(lon < divisions[1] & lat <= divisions[2], regions[3],
#                                 ifelse(lon >= divisions[1] & lat > divisions[2], regions[2],
#                                        ifelse(lon >= divisions[1] & lat <= divisions[2], regions[4], NA))))) %>%
#   group_by(region, month) %>%
#   summarise(speiTrmm = mean(speiTrmm03), speiCru = mean(speiCru03), spiTrmm = mean(spiTrmm03), spiCru = mean(spiCru03)) %>%
#   ungroup() %>%
#   mutate(date = as.Date(paste("2005", as.character(month), "01"), format = "%Y%m%d"),
#          region = factor(region, regions)) %>%
#   ggplot(aes(x = date)) +
#   geom_line(aes(y = speiTrmm), size = 2, color = "red") +
#   geom_line(aes(y = speiCru), size = 2, color = "blue") +
#   geom_line(aes(y = spiTrmm), size = 2, color = "green") +
#   geom_line(aes(y = spiCru), size = 2, color = "gray") +
#   facet_wrap(region~.) +
#   geom_hline(yintercept=0, linetype="solid", size=1, color = "black") 
  
