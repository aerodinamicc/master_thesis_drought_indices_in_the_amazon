amazonWD %>%
  select(lon, lat, year, month, SPEI03_T = speiTrmm03, SPI03_T = spiTrmm03, SPEI03_C = speiCru03, SPI03_C = spiCru03) %>%
  filter(year == 2005 & month %in% c(6, 7, 8, 9, 10)) %>%
  mutate(month = factor(ifelse(month == 6, "June",
                               ifelse(month == 7, "July",
                                      ifelse(month == 8, "August",
                                             ifelse(month == 9, "September", "October")))), c("June", "July", "August", "September", "October"))) %>%
  gather(source, value, c("SPEI03_T", "SPI03_T", "SPEI03_C", "SPI03_C")) %>%
  mutate(source = factor(source, c("SPEI03_T", "SPI03_T", "SPEI03_C", "SPI03_C"))) %>%
  regmap(region = "ama") +
  geom_tile(aes(x = lon, y = lat, fill = value)) +
  theme(strip.text.y = element_text(size = 16),
        strip.text.x = element_text(size = 16),
        axis.text.x.bottom = element_text(size = 12),
        axis.text.y.left = element_text(size = 12),
        axis.title.y.left = element_text(size = 12),
        axis.title.x.bottom = element_text(size = 12)) +
  scale_fill_gradient2() +
  facet_grid(source~month)
