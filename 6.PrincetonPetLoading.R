ncfile <- nc_open("pe_penman_rnet_empirical_monthly_1948-2008.nc")

pe <- readFromNcFile(ncfile, "pe", "longitude", "latitude", moreDims = TRUE)

petForSpei <- pe %>%
  filter(year > 1997) %>%
  mutate(daysInMonth = days_in_month(date)) %>%
  mutate(pet = value*daysInMonth) %>%
  dplyr::select(-daysInMonth, -value, -date)

pe <- pe %>%
  filter(year %in% 1998:2005) %>%
  mutate(daysInMonth = days_in_month(date)) %>%
  mutate(pet = value*daysInMonth) %>%
  dplyr::select(-daysInMonth, -value, -date)

pe <- refineResolution(pe, "pe")
petForSpei <- refineResolution(petForSpei, "pe")

pe <- cropTibbleToAmazonExtent(pe)
petForSpei <- cropTibbleToAmazonExtent(petForSpei) %>%
  arrange(lon, lat, year, month)

speiComputed <- add_column(speiComputed, pet = petForSpei$pet)
rm("petForSpei", 'ncfile')
