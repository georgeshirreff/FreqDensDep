library(tidyverse) 

# stem = "ward_nonlinear_bayesIntensity_grouped"
# n_files = 15

stem = "ward_nonlinear_bayesIntensity_multiTarget"
n_files = 240
# stem = "ward_nonlinear_bayesIntensity_catHosp"
# n_files = 1248
# stem = "ward_nonlinear_bayesIntensity_sensAnal"
# n_files = 240


tab <- read_csv(paste0("analysis/pieces/", stem, 1, ".csv"))

for(i in 2:n_files){
  print(i)
  tab <- bind_rows(tab, read_csv(paste0("analysis/pieces/",stem, i, ".csv"), col_types = cols()))
}

tab %>% write_csv(paste0("analysis/", stem, ".csv"))
