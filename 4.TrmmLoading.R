# TRMM Data loading ----
trmm <- nc_open("TRMM_Amazon_1998_April_2017.nc")

trmm_tidy <- ncdf_to_tidy(trmm, .var = "prec", .start = "1998-01-01")

monthly_trmm <- trmm_tidy %>%
  group_by(year, month, lon, lat) %>%
  summarise(monthly_prec = sum(.var))

rm("trmm_tidy")

# introducing record_number for the purpose of ordering
monthly_trmm <- monthly_trmm %>%
  group_by(lon, lat) %>%
  mutate(record_number = row_number(lon))

# introducing station_id for the purpose of ordering
monthly_trmm <- monthly_trmm %>%
  group_by(year, month) %>%
  mutate(station_id = row_number(year))

# ordering
monthly_trmm <- monthly_trmm %>%
  arrange(station_id, record_number)

speiComputed <- monthly_trmm %>%
  mutate(stays = (!(lon == -76 & lat == -12.0))) %>%
  filter(stays == TRUE) %>%
  dplyr::select(-stays, -record_number, -station_id) %>%
  filter(year > 1997 & year < 2009) %>%
  arrange(lon, lat, year, month)

colnames(speiComputed)[ncol(speiComputed)] <- "pre" 

#Cumulative water deficit
trmm_cwd <- monthly_trmm %>% 
  mutate(wd = NA) %>% 
  group_by(station_id) %>% 
  nest() %>% 
  mutate(cwd_df = purrr::map(data, cwd)) %>% 
  unnest(cwd_df)
