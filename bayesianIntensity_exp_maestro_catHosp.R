library(tidyverse)
library(FME)
library(lubridate)
library(tictoc)
library(modi)

tic()

ward_rates <- read_csv("FreqDens_ward_rates_catHosp.csv"
                       , locale = locale(tz = "CET")) %>% 
  filter(ward_id != 15)

analysis_index = commandArgs(trailingOnly = T)[1] %>% as.numeric #1248
# analysis_index = 

#### analysis table ####
cutoff_vec = 1
ward_vec = c(0, ward_rates$ward_id %>% unique)
cH_vec = ward_rates$this_catHosp %>% unique
per_vec = c("All", "Day", "Night")
analysis_grid <- rbind(expand.grid(list(ward_id = ward_vec
                                  , cutoff = cutoff_vec
                                  , this_catHosp = cH_vec
                                  , that_catHosp = "ALL"
                                  , this_per = per_vec))
                       , expand.grid(list(ward_id = ward_vec
                                  , cutoff = cutoff_vec
                                  , this_catHosp = cH_vec
                                  , that_catHosp = "patient"
                                  , this_per = per_vec))) %>% 
  unique
analysis_grid <- analysis_grid[analysis_grid$ward_id %>% order, ]
rownames(analysis_grid) = NULL
analysis_grid %>% nrow


### PARAMETERS ###

this_cutoff = analysis_grid$cutoff[analysis_index]
this_ward_id = analysis_grid$ward_id[analysis_index]
this_cH = analysis_grid$this_catHosp[analysis_index]
that_cH = analysis_grid$that_catHosp[analysis_index]
this_per = analysis_grid$this_per[analysis_index]

print(c(this_ward_id, this_cutoff, as.character(this_cH), as.character(that_cH), as.character(this_per)))

### SUBSET ###
this_ward_rates = ward_rates %>% 
  filter(Nthat > 0) %>% 
  filter(form_minutes > 0) %>% 
  # mutate(contactIntensity = pmax(contactIntensity, 1/10000)) %>% 
  mutate(Period = ifelse(hour(time) %in% 6:17, "Day", "Night")) %>% 
  filter(Period == this_per | this_per == "All") %>% 
  filter(this_catHosp == this_cH, that_catHosp == that_cH) %>% 
  filter(ward_id == this_ward_id|this_ward_id == 0) %>% 
  group_by(ward_id) %>% 
  filter(as.numeric(difftime(time, min(time), units = "hours")) >= this_cutoff) %>% 
  transmute(ward_id, Nthat, form_minutes, break_minutes, contactIntensity) %>% 
  ungroup

this_tib = this_ward_rates

# this_ward_rates %>% 
#   mutate(z = break_minutes/form_minutes) %>% 
#   pull(z)
#   summarise(sum(z))

total_breakmin = sum(this_ward_rates$break_minutes)
total_formmin = sum(this_ward_rates$form_minutes)
mean_intensity = total_breakmin/total_formmin
total_lines = nrow(this_ward_rates)


### MCMC ###

# source("freqdens_bayesianIntensity_MCMC.R", local = T)
# source("freqdens_bayesianIntensity_MCMC_norm.R", local = T)
source("freqdens_bayesianIntensity_MCMC_exp.R", local = T)
# source("freqdens_bayesianIntensity_MCMC_gamma.R", local = T)

### OUTPUT ###

ward_nonlinear <- tibble(ward_id = this_ward_id
                         , cutoff = this_cutoff
                         , this_catHosp = this_cH
                         , that_catHosp = that_cH
                         , Period = this_per

                         , eta_bayes = eta_bayes
                         , eta_lo = eta_lo
                         , eta_hi = eta_hi

                         , lines = total_lines
                         , formmin = total_formmin
                         , breakmin = total_breakmin
                         , mean_intensity = mean_intensity
)

ward_nonlinear %>%  write_csv(paste0("output/ward_nonlinear_bayesIntensity_exp_catHosp", analysis_index, ".csv"))
# ward_nonlinear %>%  write_csv(paste0("output/ward_nonlinear_bayesIntensity_gamma", analysis_index, ".csv"))

toc()
