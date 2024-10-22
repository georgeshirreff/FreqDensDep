# prefix = "ward_nonlinear_bayesIntensity_multiTarget"
prefix = "ward_nonlinear_bayesIntensity_catHosp"
prefix = "ward_nonlinear_bayesIntensity_sensAnal"


files = system(paste0("ls analysis/pieces/", prefix, "*"), intern = T)

res = read_csv(files[1])
for(f in files[-1]){
  res = rbind(res, read_csv(f))
}

res %>% write_csv(paste0("analysis/", prefix, ".csv"))
