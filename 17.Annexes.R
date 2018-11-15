#Agreement map between WDfixed_T and WDfixed_C in regards to pixels below -50 mm annual deviation of MCWD
anomalies %>%
  filter(year == 2005) %>%
  select(lon, lat, year, preTrmm_et100_r, preCru_et100_r) %>%
  mutate(agreement = ifelse(preTrmm_et100_r < -50 & preCru_et100_r < -50, "Confirmed by both",
                            ifelse(preTrmm_et100_r < -50 | preCru_et100_r < -50, "Not confirmed",0))) %>%
  filter(agreement != 0) %>%
  mutate(agreement = factor(agreement, sort(unique(agreement), decreasing = FALSE))) %>%
  select(lon, lat, agreement) %>%
  regmap(region = "ama") +
  geom_tile(aes(fill = agreement)) +
  scale_fill_manual(values = c("firebrick1", "orange"), name = "Agreement") +
  ggtitle("Agreement map between WDfixed_T and WDfixed_C\nin regards to pixels below -50 mm annual MCWD deviation")+
  theme(title = element_text(size = 18),
        axis.title.x.bottom = element_text(size = 14),
        axis.title.y.left = element_text(size = 14))

#"Monthly values (rainfall and P) averaged over the whole Amazon basin (in mm)"
pre <- amazonWD %>%
  select(year, month, preTrmm, preCru) %>%
  mutate(date = dmy(paste("01/", month, "/", year, sep = ""))) %>%
  select(-year, -month) %>%
  group_by(date) %>%
  summarise("TRMM" = sum(preTrmm)/n(),
            "CRU" = sum(preCru)/n()) %>%
  gather(Precipitation, value, 2:3) %>%
  ggplot() +
  geom_line(aes(x = date, y = value, color = Precipitation), size = 2) +
  ggtitle("Monthly rainfall values averaged over the whole Amazon basin (in mm)")+
  theme(title = element_text(size = 18),
        axis.title.x.bottom = element_text(size = 14),
        axis.title.y.left = element_text(size = 14))
  
pet <- amazonWD %>%
  select(year, month, pet, petCru) %>%
  mutate(date = dmy(paste("01/", month, "/", year, sep = ""))) %>%
  select(-year, -month) %>%
  group_by(date) %>%
  summarise("Princeton University" = sum(pet)/n(),
            "CRU" = sum(petCru)/n()) %>%
  gather(PET, value, 2:3) %>%
  ggplot() +
  geom_line(aes(x = date, y = value, color = PET), size = 2) +
  ggtitle("Monthly PET values averaged over the whole Amazon basin (in mm)")+
  theme(title = element_text(size = 18),
        axis.title.x.bottom = element_text(size = 14),
        axis.title.y.left = element_text(size = 14))

ggarrange(pre, pet, nrow = 2)

