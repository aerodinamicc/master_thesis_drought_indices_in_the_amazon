amazonPREC <- monthly_trmm %>%
  filter(year > 1997 & year < 2006) %>%
  select(-station_id, -record_number) %>%
  mutate(stays = (!(lon == -76 & lat == -12.0)),
         preTrmm = monthly_prec) %>%
  filter(stays == TRUE) %>%
  arrange(lon, lat, year, month)%>%
  select(-stays, -monthly_prec)

amazon98_05meanET <- amazon_et %>%
  group_by(lon,lat, month) %>%
  mutate(meanET = mean(aet))

amazonWD <- add_column(amazonPREC, aet = amazon_et$aet)
amazonWD <- add_column(amazonWD, preCru = speiCruComputed$pre)
amazonWD <- add_column(amazonWD, aetAvg = amazon98_05meanET$meanET)
amazonWD <- add_column(amazonWD, pet = pe$pet)
amazonWD <- add_column(amazonWD, petCru = speiCruComputed$pet)
rm("amazonPREC", "amazon98_05meanET", "amazon_et")

amazonWD <- amazonWD %>%
  group_by(lon, lat, month) %>%
  mutate(et100 = 100, petAvg = mean(pet), petCruAvg = mean(petCru)) %>%
  arrange(lon, lat, year, month)

# Computing CWD --------------------------
amazonWD <- amazonWD %>%
  mutate(cwd_preTrmm_et100_r = as.double(NA),
         cwd_preTrmm_aet_r = as.double(NA),
         cwd_preTrmm_aetAvg_r = as.double(NA),
         cwd_preTrmm_pet_r = as.double(NA),
         cwd_preTrmm_petAvg_r = as.double(NA),
         cwd_preTrmm_pet_nr = as.double(NA),
         cwd_preTrmm_petAvg_nr = as.double(NA),
         cwd_aet_pet_r = as.double(NA),
         cwd_preCru_et100_r = as.double(NA),
         cwd_preCru_aet_r = as.double(NA),
         cwd_preCru_aetAvg_r = as.double(NA),
         cwd_preCru_petCru_r = as.double(NA),
         cwd_preCru_petCruAvg_r = as.double(NA),
         cwd_preCru_petCru_nr = as.double(NA),
         cwd_preCru_petCruAvg_nr = as.double(NA),
         cwd_aet_petCru_r = as.double(NA))

cwdVector <- colnames(amazonWD[14:ncol(amazonWD)])

amazonWD <- computeCwd(amazonWD, cwdVector)

amazonWD <- cbind(amazonWD, speiTrmm03 = speiComputed$spei3, spiTrmm03 = speiComputed$spi3)
amazonWD <- cbind(amazonWD, speiCru03 = speiCruComputed$spei3, spiCru03 = speiCruComputed$spi3)

amazonWD <- amazonWD %>%
  mutate(season = ifelse(month %in% c(4:9), "dry", "wet"))
