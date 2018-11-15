anomalies <- tibble()
anomaliesSd <- tibble()

for (cwd in 1:length(cwdVector)) {
  result <- computeMcwdAndAnomalies(amazonWD, anomalies, anomaliesSd, cwdVector[cwd])
  amazonWD <- result [[1]]
  anomalies <- result[[2]]
  anomaliesSd <- result[[3]]
}

anomalies <- anomalies %>%
  arrange(lon, lat, year)

anomaliesSd <- anomaliesSd %>%
  arrange(lon, lat, year)

amazonWD <- amazonWD %>%
  arrange(lon, lat, year, month)