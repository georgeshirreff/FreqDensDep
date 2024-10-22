source("2.1_load_analysis.R")

ward_aic = ward_intervals_loglik %>% 
  ## sum likelihoods
  group_by(ward_id, newID, this_status, that_status) %>% 
  summarise(across(starts_with("loglik"), sum)
            # , eta = eta_median[1]
  ) %>% 
  pivot_longer(-c("ward_id", "newID", "this_status", "that_status"), names_sep = "_", names_to = c(".value", "Model")) %>% 
  mutate(K = case_when(Model == "constant" ~ 1
                       , Model == "linear" ~ 1
                       , Model == "nonlin" ~ 2
                       , Model == "cdiurnal" ~ 2
                       , Model == "ldiurnal" ~ 2
                       , Model == "ndiurnal" ~ 4
  )) %>% 
  mutate(AIC = 2*K - 2*loglik) %>% 
  group_by(ward_id, newID, this_status, that_status) %>% 
  mutate(dAIC = AIC - min(AIC)) %>% 
  pivot_wider(id_cols = c("ward_id", "newID", "this_status", "that_status"), names_from = "Model", values_from = "dAIC") %>% 
  ungroup


tab <- ward_aic %>% 
  filter(ward_id != 15) %>% 
  filter(this_status == "ALL", that_status == "ALL") %>% 
  arrange(newID) %>% 
  transmute(Ward = newID, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal)

# tab_this <- ward_aic %>% 
#   filter(ward_id != 15) %>% 
#   filter(this_status != "ALL", that_status == "ALL") %>% 
#   arrange(newID) %>% 
#   transmute(Ward = newID, this_status, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
#   mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))
# 
# tab_that <- ward_aic %>% 
#   filter(ward_id != 15) %>% 
#   filter(this_status == "ALL", that_status != "ALL") %>% 
#   arrange(newID) %>% 
#   transmute(Ward = newID, that_status, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
#   mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))

tab_thisthat <- ward_aic %>% 
  filter(ward_id != 15) %>% 
  filter(this_status != "ALL", that_status != "ALL") %>% 
  arrange(newID) %>% 
  transmute(Ward = newID, this_status, that_status, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal)

# output tables for reuse
tab %>%   
  write_tsv(paste0("output/dAIC.tsv"))

tab_thisthat %>%   
  write_tsv(paste0("output/dAIC_thisthat.tsv"))

# output tables for presentation
tab %>%   
  mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 1), trim = T, nsmall=1)))) %>% 
  write_tsv(paste0("output/Tab2 - dAIC_intensity.tsv"))

tab_thisthat %>% 
  filter(this_status %in% c("PA", "PE"), that_status == "PA") %>% 
  arrange(this_status) %>% 
  mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 1), trim = T, nsmall=1)))) %>% 
  write_tsv(paste0("output/Tab34 - dAIC_intensity_thisthat.tsv"))




# catHosp

ward_catHosp_aic = ward_intervals_catHosp_loglik %>% 
  ## sum likelihoods
  group_by(ward_id, newID, this_catHosp, that_catHosp) %>% 
  summarise(across(starts_with("loglik"), sum)
            # , eta = eta_median[1]
  ) %>% 
  pivot_longer(-c("ward_id", "newID", "this_catHosp", "that_catHosp"), names_sep = "_", names_to = c(".value", "Model")) %>% 
  mutate(K = case_when(Model == "constant" ~ 1
                       , Model == "linear" ~ 1
                       , Model == "nonlin" ~ 2
                       , Model == "cdiurnal" ~ 2
                       , Model == "ldiurnal" ~ 2
                       , Model == "ndiurnal" ~ 4
  )) %>% 
  mutate(AIC = 2*K - 2*loglik) %>% 
  group_by(ward_id, newID, this_catHosp, that_catHosp) %>% 
  mutate(dAIC = AIC - min(AIC)) %>% 
  pivot_wider(id_cols = c("ward_id", "newID", "this_catHosp", "that_catHosp"), names_from = "Model", values_from = "dAIC") %>% 
  ungroup

# ward_aic %>% 
#   filter(this_status == this_stat, that_status == that_stat)
# ward_nonlinear_main %>% 
#   filter(ward_id == 6) %>% 
#   filter(this_status == this_stat, that_status == that_stat)
# tab <- ward_aic %>% 
#   filter(ward_id != 15) %>% 
#   filter(this_status == "ALL", that_status == "ALL") %>% 
#   arrange(newID) %>% 
#   transmute(Ward = newID, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
#   mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))
# 
# tab_this <- ward_aic %>% 
#   filter(ward_id != 15) %>% 
#   filter(this_status != "ALL", that_status == "ALL") %>% 
#   arrange(newID) %>% 
#   transmute(Ward = newID, this_status, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
#   mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))
# 
# tab_that <- ward_aic %>% 
#   filter(ward_id != 15) %>% 
#   filter(this_status == "ALL", that_status != "ALL") %>% 
#   arrange(newID) %>% 
#   transmute(Ward = newID, that_status, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
#   mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))

tab_catHosp_allthisthat <- ward_catHosp_aic %>% 
  filter(ward_id != 15) %>%
  filter(!is.na(constant)) %>% 
  arrange(newID) %>% 
  transmute(Ward = newID, this_catHosp, that_catHosp, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal)

tab_catHosp_allthisthat %>% 
  write_tsv(paste0("output/dAIC_catHosp_allthisthat.tsv"))


### SENSITIVITY ANALYSIS ###

tabP = list()
this_boundary = 8 #default value
for(this_popcut in popcut_vec){
  ward_intervals_Pedit = ward_intervals %>% 
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
  
  
  
  
  ward_intervals_Ploglik <- ward_intervals_Pedit %>% 
    {rbind(., mutate(., ward_id = 0))} %>% 
    left_join(newID_ref) %>% 
    filter(!is.na(newID)) %>% 
    left_join(ward_nonlinear_sensAnal %>% filter(Period == "All", popcut == this_popcut) %>% transmute(ward_id, this_status, that_status, phi_nondiurnal = phi_median)) %>%
    left_join(ward_nonlinear_sensAnal %>% filter(Period != "All", popcut == this_popcut) %>% transmute(ward_id, Period, this_status, that_status, phi_diurnal = phi_median)) %>%
    mutate(phi_diurnal = ifelse(is.na(phi_diurnal), phi_nondiurnal, phi_diurnal)) %>% 
    # group_by(ward_id, this_status, that_status) %>% 
    # filter(form_minutes > 0) %>% 
    # mutate(var = var(contactIntensity)) %>% 
    ## calculate all the scaling parameters and corresponding rates for each model
    group_by(ward_id, this_status, that_status) %>%
    mutate(w_constant = sum(break_minutes)/sum(form_minutes)
           , w_linear = sum(break_minutes)/sum(form_minutes*Nthat)
           , w_nonlin = sum(break_minutes)/sum(form_minutes*Nthat^phi_nondiurnal)) %>% 
    group_by(ward_id, this_status, that_status, Period) %>% #group by day/night to calculate the day and night parameters
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
  
  
  ward_Paic = ward_intervals_Ploglik %>% 
    ## sum likelihoods
    group_by(ward_id, newID, this_status, that_status) %>% 
    summarise(across(starts_with("loglik"), sum)
              # , eta = eta_median[1]
    ) %>% 
    pivot_longer(-c("ward_id", "newID", "this_status", "that_status"), names_sep = "_", names_to = c(".value", "Model")) %>% 
    mutate(K = case_when(Model == "constant" ~ 1
                         , Model == "linear" ~ 1
                         , Model == "nonlin" ~ 2
                         , Model == "cdiurnal" ~ 2
                         , Model == "ldiurnal" ~ 2
                         , Model == "ndiurnal" ~ 4
    )) %>% 
    mutate(AIC = 2*K - 2*loglik) %>% 
    group_by(ward_id, newID, this_status, that_status) %>% 
    mutate(dAIC = AIC - min(AIC)) %>% 
    pivot_wider(id_cols = c("ward_id", "newID", "this_status", "that_status"), names_from = "Model", values_from = "dAIC") %>% 
    ungroup
  
  tabP[[paste0("popcut", this_popcut)]] <- ward_Paic %>% 
    filter(ward_id != 15) %>% 
    filter(this_status == "ALL", that_status == "ALL") %>% 
    arrange(newID) %>% 
    transmute(Ward = newID, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal)
  
}


tabPs <- rbind(tabP[[paste0("popcut", 0)]] %>% mutate(popcut = 0), 
               read_tsv("output/dAIC.tsv") %>% mutate(popcut = 10), 
               # tab[[paste0("popcut", 10)]] %>% mutate(popcut = 10), 
               tabP[[paste0("popcut", 20)]] %>% mutate(popcut = 20)) %>% 
  select(Ward, popcut, everything()) %>% 
  arrange(Ward, popcut)

tabPs %>% 
  write_tsv(paste0("output/dAIC_popcuts.tsv"))

tabPs %>% 
  mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 1), trim = T, nsmall=1)))) %>% 
  write_tsv(paste0("output/SuppTab1 - dAIC_popcuts.tsv"))



tabB = list()
this_popcut = 10 #default value
for(this_boundary in boundary_vec){
  
  ward_intervals_Bedit = ward_intervals %>% 
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
  
  
  
  
  ward_intervals_Bloglik <- ward_intervals_Bedit %>% 
    {rbind(., mutate(., ward_id = 0))} %>% 
    left_join(newID_ref) %>% 
    filter(!is.na(newID)) %>% 
    left_join(ward_nonlinear_sensAnal %>% filter(Period == "All", popcut == this_popcut) %>% transmute(ward_id, this_status, that_status, phi_nondiurnal = phi_median)) %>%
    left_join(ward_nonlinear_sensAnal %>% filter(Period != "All", popcut == this_popcut, boundary == this_boundary) %>% transmute(ward_id, Period, this_status, that_status, phi_diurnal = phi_median)) %>%
    mutate(phi_diurnal = ifelse(is.na(phi_diurnal), phi_nondiurnal, phi_diurnal)) %>% 
    # group_by(ward_id, this_status, that_status) %>% 
    # filter(form_minutes > 0) %>% 
    # mutate(var = var(contactIntensity)) %>% 
    ## calculate all the scaling parameters and corresponding rates for each model
    group_by(ward_id, this_status, that_status) %>%
    mutate(w_constant = sum(break_minutes)/sum(form_minutes)
           , w_linear = sum(break_minutes)/sum(form_minutes*Nthat)
           , w_nonlin = sum(break_minutes)/sum(form_minutes*Nthat^phi_nondiurnal)) %>% 
    group_by(ward_id, this_status, that_status, Period) %>% #group by day/night to calculate the day and night parameters
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
  
  
  ward_Baic = ward_intervals_Bloglik %>% 
    ## sum likelihoods
    group_by(ward_id, newID, this_status, that_status) %>% 
    summarise(across(starts_with("loglik"), sum)
              # , eta = eta_median[1]
    ) %>% 
    pivot_longer(-c("ward_id", "newID", "this_status", "that_status"), names_sep = "_", names_to = c(".value", "Model")) %>% 
    mutate(K = case_when(Model == "constant" ~ 1
                         , Model == "linear" ~ 1
                         , Model == "nonlin" ~ 2
                         , Model == "cdiurnal" ~ 2
                         , Model == "ldiurnal" ~ 2
                         , Model == "ndiurnal" ~ 4
    )) %>% 
    mutate(AIC = 2*K - 2*loglik) %>% 
    group_by(ward_id, newID, this_status, that_status) %>% 
    mutate(dAIC = AIC - min(AIC)) %>% 
    pivot_wider(id_cols = c("ward_id", "newID", "this_status", "that_status"), names_from = "Model", values_from = "dAIC") %>% 
    ungroup
  
  tabB[[paste0("boundary", this_boundary)]] <- ward_Baic %>% 
    filter(ward_id != 15) %>% 
    filter(this_status == "ALL", that_status == "ALL") %>% 
    arrange(newID) %>% 
    transmute(Ward = newID, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal)
  
}


tabBs <- rbind(tabB[[paste0("boundary", 4)]] %>% mutate(boundary = 4), 
               tabB[[paste0("boundary", 6)]] %>% mutate(boundary = 6), 
               # tab[[paste0("boundary", 8)]] %>% mutate(boundary = 8), 
               read_tsv("output/dAIC.tsv") %>% mutate(boundary = 8),
               tabB[[paste0("boundary", 10)]] %>% mutate(boundary = 10)) %>% 
  select(Ward, boundary, everything()) %>% 
  arrange(Ward, boundary)

tabBs %>% 
  write_tsv(paste0("output/dAIC_boundaries.tsv"))

tabBs %>% 
  mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 1), trim = T, nsmall=1)))) %>% 
  write_tsv(paste0("output/SuppTab2 - dAIC_boundaries.tsv"))

