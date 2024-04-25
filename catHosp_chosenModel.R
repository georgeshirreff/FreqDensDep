ward_rates_catHosp <- read_csv("FreqDens_ward_rates_catHosp.csv"
                       , locale = locale(tz = "CET"))

# ward_nonlinear_catHosp <- read_csv("output/ward_nonlinear_bayesIntensity_catHosp.csv")


ward_rates_catHosp_edit = ward_rates_catHosp %>% 
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


ward_rates_catHosp_loglik <- ward_rates_catHosp_edit %>% 
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



# ward_rates_loglik %>%
#   filter(this_status == "PA", that_status == "PA", ward_id %in% c(2, 5, 8, 13))

ward_catHosp_aic = ward_rates_catHosp_loglik %>% 
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
  transmute(Ward = newID, this_catHosp, that_catHosp, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
  mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))

tab_catHosp_allthisthat %>% 
  write_tsv(paste0("AICtables/dAIC_intensity_catHosp_allthisthat.tsv"))


chosen_model_catHosp_allthisthat <- read_tsv(paste0("AICtables/dAIC_intensity_catHosp_allthisthat.tsv")) %>% 
  pivot_longer(-c("Ward", "this_catHosp", "that_catHosp"), names_to = "ChosenModel", values_to = "dAIC") %>% 
  rename(newID = Ward) %>% 
  filter(dAIC == 0) %>% 
  select(-dAIC)

