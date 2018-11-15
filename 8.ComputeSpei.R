source('tidyspei.R')

speiComputed <- speiComputed %>%
  group_by(lat, lon) %>%
  arrange(year, month) %>%
  nest() %>%
  mutate(
    spei3 = purrr::map(data, tidyspei, scale = 3),
    spi3 = purrr::map(data, tidyspei, scale = 3, mode = "spi")
  ) %>%
  unnest(data, spei3, spi3)

speiComputed <- speiComputed %>%
  arrange(lon, lat, year, month) %>%
  filter(year > 1997 & year < 2006)

#sapply(speiComputed, function(x){sum(is.infinite(x))})
#there are 7 Inf values of SPEI
speiComputed$spei3[speiComputed$spei3 %in% c(Inf, -Inf)] <- NA

#sapply(speiComputed, function(x){sum(is.na(x))})
#SPEI 3993, SPI 4049

#CRU
speiCruComputed <- speiCruComputed %>%
  group_by(lat, lon) %>%
  arrange(year, month) %>%
  nest() %>%
  mutate(
    spei3 = purrr::map(data, tidyspei, scale = 3),
    spi3 = purrr::map(data, tidyspei, scale = 3, mode = "spi")
  ) %>%
  unnest(data, spei3, spi3)

speiCruComputed <- speiCruComputed %>%
  arrange(lon, lat, year, month) %>%
  filter(year %in% 1998:2005)

#sapply(speiCruComputed, function(x){sum(is.infinite(x))})
#there are 97 Inf values of SPEI
speiCruComputed$spei3[speiCruComputed$spei3 %in% c(Inf, -Inf)] <- NA

#sapply(speiCruComputed, function(x){sum(is.na(x))})
#unreasonably high number of NA values epsecially with SPI - 22378, SPEI - 9987
#input data is checked - no NA value present there
