library(tidyverse)


ward_intervals <- read_csv("data/FreqDens_ward_intervals.csv"
                           , locale = locale(tz = "CET"))

ward_intervals_catHosp <- read_csv("data/FreqDens_ward_intervals_catHosp.csv"
                               , locale = locale(tz = "CET"))


# typeID <- read_csv("data/typeID.csv")
typeID <- read_csv("data/typeID_hospAbb.csv")

newID_ref <- typeID %>% 
  arrange(newID) %>%
  mutate(newID_hospAbb = paste0(hospAbb, " - ", newID)) %>% 
  select(ward_id, newID, newID_hospAbb) %>% 
  rbind(tibble(ward_id = 0, newID = "All wards", newID_hospAbb = "All wards")) %>% 
  mutate(newID = factor(newID, levels = newID))


# newID_ref <- typeID %>% 
#   select(ward_id, newID) %>% 
#   arrange(newID) %>% 
#   rbind(tibble(ward_id = 0, newID = "All wards")) %>% 
#   mutate(newID = factor(newID, levels = newID))


# read results of Bayesian analysis
ward_nonlinear_main <- read_csv(paste0("analysis/ward_nonlinear_bayesIntensity_multiTarget.csv"))

ward_nonlinear_catHosp <- read_csv(paste0("analysis/ward_nonlinear_bayesIntensity_catHosp.csv"))

ward_nonlinear_sensAnal <- read_csv(paste0("analysis/ward_nonlinear_bayesIntensity_sensAnal.csv")) %>% 
  rename(boundary = this_boundary)


ward_nonlinear_plottable = ward_nonlinear_main %>% 
  left_join(newID_ref) %>% 
  filter(!is.na(newID)) %>% 
  select(ward_id, newID, this_status, that_status, Period, phi_median, phi_lo, phi_hi)

ward_nonlinear_catHosp_plottable = ward_nonlinear_catHosp %>% 
  left_join(newID_ref) %>% 
  filter(!is.na(newID)) %>% 
  select(ward_id, newID, this_catHosp, that_catHosp, Period, phi_median, phi_lo, phi_hi)



### produce usable dataset ###

this_popcut = ward_nonlinear_main$popcut %>% unique
this_boundary = ward_nonlinear_main$this_boundary %>% unique %>% {.[. > 0]}

ward_intervals_edit = ward_intervals %>% 
  group_by(ward_id) %>%  #cut off the start and end of the data according to total population size
  mutate(cum_max = map_dbl(.x = seq_along(Nall), 
                           ~ max(Nall[1:.x])), 
         cum_max_reverse = map_dbl(.x = seq_along(Nall), 
                                   ~ max(Nall[.x:n()]))) %>% 
  filter(cum_max >= this_popcut & cum_max_reverse >= this_popcut) %>% 
  ungroup %>% 
  mutate(Period = ifelse(hour(time) %in% this_boundary:(this_boundary+11), "Day", "Night")) %>%  # identify the time periods
  mutate(Period = case_when(this_boundary == -1 ~ "All",
                            T ~ Period)) %>% # this should never apply in the post-processing
  filter(Nthat > 0, form_minutes > 0) # do the standard filtering for negligible time periods


ward_intervals_loglik <- ward_intervals_edit %>% 
  {rbind(., mutate(., ward_id = 0))} %>% 
  left_join(newID_ref) %>% 
  filter(!is.na(newID)) %>% 
  left_join(ward_nonlinear_main %>% filter(Period == "All") %>% transmute(ward_id, this_status, that_status, phi_nondiurnal = phi_median, phi_nondiurnal_lo = phi_lo, phi_nondiurnal_hi = phi_hi)) %>%
  left_join(ward_nonlinear_main %>% filter(Period != "All") %>% transmute(ward_id, Period, this_status, that_status, phi_diurnal = phi_median, phi_diurnal_lo = phi_lo, phi_diurnal_hi = phi_hi)) %>%
  mutate(phi_diurnal = ifelse(is.na(phi_diurnal), phi_nondiurnal, phi_diurnal), 
         phi_diurnal_lo = ifelse(is.na(phi_diurnal_lo), phi_nondiurnal_lo, phi_diurnal_lo), 
         phi_diurnal_hi = ifelse(is.na(phi_diurnal_hi), phi_nondiurnal_hi, phi_diurnal_hi)) %>% 
  # group_by(ward_id, this_status, that_status) %>% 
  # filter(form_minutes > 0) %>% 
  # mutate(var = var(contactIntensity)) %>% 
  ## calculate all the scaling parameters and corresponding rates for each model
  group_by(ward_id, this_status, that_status) %>%
  mutate(w_constant = sum(break_minutes)/sum(form_minutes)
         , w_linear = sum(break_minutes)/sum(form_minutes*Nthat)
         , w_nonlin = sum(break_minutes)/sum(form_minutes*Nthat^phi_nondiurnal), 
         , w_nonlin_lo = sum(break_minutes)/sum(form_minutes*Nthat^phi_nondiurnal_lo), 
         , w_nonlin_hi = sum(break_minutes)/sum(form_minutes*Nthat^phi_nondiurnal_hi), 
         , sd_constant = sd(break_minutes/form_minutes)
         , sd_linear = sd(break_minutes/form_minutes*Nthat)
         , sd_nonlin = sd(break_minutes/form_minutes*Nthat^phi_nondiurnal)
         , sd_nonlin_lo = sd(break_minutes/form_minutes*Nthat^phi_nondiurnal_lo)
         , sd_nonlin_hi = sd(break_minutes/form_minutes*Nthat^phi_nondiurnal_hi)
         ) %>% 
  group_by(ward_id, this_status, that_status, Period) %>% #group by day/night to calculate the day and night parameters
  mutate(w_cdiurnal = sum(break_minutes)/sum(form_minutes)
         , w_ldiurnal = sum(break_minutes)/sum(form_minutes*Nthat)
         , w_ndiurnal = sum(break_minutes)/sum(form_minutes*Nthat^phi_diurnal), 
         , w_ndiurnal_lo = sum(break_minutes)/sum(form_minutes*Nthat^phi_diurnal_lo), 
         , w_ndiurnal_hi = sum(break_minutes)/sum(form_minutes*Nthat^phi_diurnal_hi), 
         , sd_cdiurnal = sd(break_minutes/form_minutes)
         , sd_ldiurnal = sd(break_minutes/form_minutes*Nthat)
         , sd_ndiurnal = sd(break_minutes/form_minutes*Nthat^phi_diurnal)
         , sd_ndiurnal_lo = sd(break_minutes/form_minutes*Nthat^phi_diurnal_lo)
         , sd_ndiurnal_hi = sd(break_minutes/form_minutes*Nthat^phi_diurnal_hi)
         ) %>% 
  mutate(rate_constant = w_constant #calculate the rates for both non-diurnal and diurnal models
         , rate_linear = w_linear*Nthat
         , rate_nonlin = w_nonlin*Nthat^phi_nondiurnal
         , rate_cdiurnal = w_cdiurnal
         , rate_ldiurnal = w_ldiurnal*Nthat
         , rate_ndiurnal = w_ndiurnal*Nthat^phi_diurnal
         ) %>%
  
  ## calculate corresponding likelihoods for each time window
  mutate(across(starts_with("rate")
                , .fns = list(loglik = function(r) dexp(x = contactIntensity, rate = 1/r, log = T))
                , .names = "{.fn}_{.col}")) %>% 
  rename_all(function(x) gsub("loglik_rate_", "loglik_", x))

ward_intervals_loglik %>% names


### and for catHosp

ward_intervals_catHosp_edit = ward_intervals_catHosp %>% 
  group_by(ward_id) %>%  #cut off the start and end of the data according to total population size
  mutate(cum_max = map_dbl(.x = seq_along(Nall), 
                           ~ max(Nall[1:.x])), 
         cum_max_reverse = map_dbl(.x = seq_along(Nall), 
                                   ~ max(Nall[.x:n()]))) %>% 
  filter(cum_max >= this_popcut & cum_max_reverse >= this_popcut) %>% 
  ungroup %>% 
  mutate(Period = ifelse(hour(time) %in% this_boundary:(this_boundary+11), "Day", "Night")) %>%  # identify the time periods
  mutate(Period = case_when(this_boundary == -1 ~ "All",
                            T ~ Period)) %>% # this should never apply in the post-processing
  filter(Nthat > 0, form_minutes > 0) # do the standard filtering for negligible time periods


ward_intervals_catHosp_loglik <- ward_intervals_catHosp_edit %>% 
  {rbind(., mutate(., ward_id = 0))} %>% 
  left_join(newID_ref) %>% 
  filter(!is.na(newID)) %>% 
  left_join(ward_nonlinear_catHosp %>% filter(Period == "All") %>% transmute(ward_id, this_catHosp, that_catHosp, phi_nondiurnal = phi_median)) %>%
  left_join(ward_nonlinear_catHosp %>% filter(Period != "All") %>% transmute(ward_id, Period, this_catHosp, that_catHosp, phi_diurnal = phi_median)) %>%
  mutate(phi_diurnal = ifelse(is.na(phi_diurnal), phi_nondiurnal, phi_diurnal)) %>% 
  # group_by(ward_id, this_status, that_status) %>% 
  # filter(form_minutes > 0) %>% 
  # mutate(var = var(contactIntensity)) %>% 
  ## calculate all the scaling parameters and corresponding rates for each model
  group_by(ward_id, this_catHosp, that_catHosp) %>%
  mutate(w_constant = sum(break_minutes)/sum(form_minutes)
         , w_linear = sum(break_minutes)/sum(form_minutes*Nthat)
         , w_nonlin = sum(break_minutes)/sum(form_minutes*Nthat^phi_nondiurnal)) %>% 
  group_by(ward_id, this_catHosp, that_catHosp, Period) %>% #group by day/night to calculate the day and night parameters
  mutate(w_cdiurnal = sum(break_minutes)/sum(form_minutes)
         , w_ldiurnal = sum(break_minutes)/sum(form_minutes*Nthat)
         , w_ndiurnal = sum(break_minutes)/sum(form_minutes*Nthat^phi_diurnal)) %>% 
  mutate(rate_constant = w_constant #calculate the rates for both non-diurnal and diurnal models
         , rate_linear = w_linear*Nthat
         , rate_nonlin = w_nonlin*Nthat^phi_nondiurnal
         , rate_cdiurnal = w_cdiurnal
         , rate_ldiurnal = w_ldiurnal*Nthat
         , rate_ndiurnal = w_ndiurnal*Nthat^phi_diurnal
  ) %>%
  
  ## calculate corresponding likelihoods for each time window
  mutate(across(starts_with("rate")
                , .fns = list(loglik = function(r) dexp(x = contactIntensity, rate = 1/r, log = T))
                , .names = "{.fn}_{.col}")) %>% 
  rename_all(function(x) gsub("loglik_rate_", "loglik_", x))

### and for sensitivity analysis 

popcut_vec = ward_nonlinear_sensAnal$popcut %>% unique
boundary_vec = ward_nonlinear_sensAnal$boundary %>% unique %>% {.[. > 0]} %>% sort
