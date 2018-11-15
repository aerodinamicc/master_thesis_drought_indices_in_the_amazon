dir(pattern = ".nc")
cruPreNc <- nc_open("cru_ts4.01.1901.2016.pre.dat.nc")
cruPetNc <- nc_open("cru_ts4.01.1901.2016.pet.dat.nc")

cruPre <- readFromNcFile(cruPreNc, "pre")
cruPre <- cruPre %>%
  filter(year %in% 1998:2005)
cruPre <- cropTibbleToAmazonExtent(cruPre)

cruPet <- readFromNcFile(cruPetNc, "pet")
cruPet <- cruPet %>%
  filter(year %in% 1998:2005)
cruPet <- cropTibbleToAmazonExtent(cruPet)

speiCruComputed <- cruPre %>%
  arrange(lat, lon, year, month) %>%
  select(-date, pre = value)

speiCruComputed <- cbind(speiCruComputed, cruPet %>%
  arrange(lat, lon, year, month) %>%
  mutate(pet = days_in_month(date) * value) %>%
  select(pet))

rm('cruPreNc', 'cruPetNc', 'cruPre', 'cruPet')
  