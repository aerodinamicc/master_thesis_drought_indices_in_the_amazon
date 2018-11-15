et_raw <- nc_open("LandFluxEVAL.merged.89-05.monthly.all.nc")
amazon_et <- readFromNcFile(et_raw, "ET_mean")
rm("et_raw")

#Going from daily to monthly ET and mean values computation-----------------------------
amazon_et <- amazon_et %>%
  filter(year %in% 1998:2005) %>%
  mutate(daysInMonth = days_in_month(date), aet = daysInMonth*value) %>%
  dplyr::select(-value, -date, -daysInMonth)

amazon_et <- refineResolution(amazon_et, "totalEt")
amazon_et <- cropTibbleToAmazonExtent(amazon_et)
amazon_et$aet<-round(as.numeric(amazon_et$aet), 3)