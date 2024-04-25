library(tidyverse)
library(FME)
library(lubridate)
library(tictoc)
library(modi)

tic()

ward_rates <- read_csv("FreqDens_ward_rates_catHosp.csv"
                       , locale = locale(tz = "CET")) %>%
  filter(ward_id != 15)

analysis_index = commandArgs(trailingOnly = T)[1] %>% as.numeric #1-1248
# analysis_index = 82
print(analysis_index)

#### analysis table ####
popcut_vec = 10
ward_vec = c(0, ward_rates$ward_id %>% unique)
cH_vec = ward_rates$this_catHosp %>% unique
# stat_vec = c("ALL", "PA", "V", "PE")
# stat_vec = "ALL"
per_vec = c("All", "Day", "Night")
boundary_vec = c(8)
analysis_grid <- rbind(expand.grid(list(ward_id = ward_vec
                                  , popcut = popcut_vec
                                  , this_catHosp = cH_vec
                                  , that_catHosp = c("ALL", "patient")
                                  , this_boundary = boundary_vec
                                  , this_per = c("Day", "Night")))
                       , expand.grid(list(ward_id = ward_vec
                                          , popcut = popcut_vec
                                          , this_catHosp = cH_vec
                                          , that_catHosp = c("ALL", "patient")
                                          , this_boundary = -1
                                          , this_per = "All"))) %>% 
  unique 
analysis_grid <- analysis_grid[analysis_grid$ward_id %>% order, ]
rownames(analysis_grid) = NULL
analysis_grid %>% nrow
# View(analysis_grid)

### PARAMETERS ###

this_popcut = analysis_grid$popcut[analysis_index]
this_ward_id = analysis_grid$ward_id[analysis_index]
this_cH = analysis_grid$this_catHosp[analysis_index]
that_cH = analysis_grid$that_catHosp[analysis_index]
this_boundary = analysis_grid$this_boundary[analysis_index]
this_per = analysis_grid$this_per[analysis_index]

print(c(this_ward_id, as.character(this_popcut), as.character(this_cH), as.character(that_cH), as.character(this_boundary), as.character(this_per)))

### PREPARE DATASET ###
ward_rates_edit = ward_rates %>% 
  group_by(ward_id) %>%  #cut off the start and end of the data according to total population size
  mutate(cum_max = map_dbl(.x = seq_along(Nall), 
                           ~ max(Nall[1:.x])), 
         cum_max_reverse = map_dbl(.x = seq_along(Nall), 
                                   ~ max(Nall[.x:n()]))) %>% 
  filter(cum_max >= this_popcut & cum_max_reverse >= this_popcut) %>% 
  ungroup %>% 
  mutate(Period = ifelse(hour(time) %in% this_boundary:(this_boundary+11), "Day", "Night")) %>%  # identify the time periods
  mutate(Period = case_when(this_boundary == -1 ~ "All",
                            T ~ Period)) %>% # and set to all if this_boundary == -1
  filter(Nthat > 0, form_minutes > 0) # do the standard filtering for negligible time periods

### SUBSET ###
this_ward_rates <- ward_rates_edit %>% 
  # now do the major filtering for analysis (or not, in the case of post-processing)
  filter(this_catHosp == this_cH, that_catHosp == that_cH) %>% #crucial difference if this comes after!
  filter(ward_id == this_ward_id|this_ward_id == 0) %>% 
  filter(Period == this_per) %>% 
  transmute(ward_id, Nthat, form_minutes, break_minutes, contactIntensity) %>% 
  ungroup

this_tib = this_ward_rates



total_breakmin = sum(this_ward_rates$break_minutes)
total_formmin = sum(this_ward_rates$form_minutes)
mean_intensity = total_breakmin/total_formmin
total_lines = nrow(this_ward_rates)


### MCMC ###

# source("freqdens_bayesianIntensity_MCMC.R", local = T)
# source("freqdens_bayesianIntensity_MCMC_norm.R", local = T)
# source("freqdens_bayesianIntensity_MCMC_exp.R", local = T)
# source("freqdens_bayesianIntensity_MCMC_gamma.R", local = T)
source("freqdens_bayesianIntensity_MCMC_multiTarget.R", local = T)

### OUTPUT ###

ward_nonlinear <- tibble(ward_id = this_ward_id
                         , popcut = this_popcut
                         , this_catHosp = this_cH
                         , that_catHosp = that_cH
                         , this_boundary = this_boundary
                         , Period = this_per

                         
                         , lines = total_lines
                         , formmin = total_formmin
                         , breakmin = total_breakmin
                         , mean_intensity = mean_intensity
                         
                         , phi_median 
                         , phi_best 
                         
                         , phi_lo
                         , phi_hi
                         
) %>% 
  cbind(tibble(s = 1:length(phi_vec100), phi = phi_vec100) %>% 
          pivot_wider(names_prefix = "phi_sample", names_from = "s", values_from = "phi"))

ward_nonlinear %>%  write_csv(paste0("output/ward_nonlinear_bayesIntensity_catHosp", analysis_index, ".csv"))
# ward_nonlinear %>%  write_csv(paste0("output/ward_nonlinear_bayesIntensity_gamma", analysis_index, ".csv"))

toc()
