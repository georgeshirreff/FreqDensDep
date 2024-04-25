library(tidyverse)
library(ggplot2)
# library(FME)
library(ggh4x)
library(modi)
library(scales)

#### make plots comparing models #####
weighted.var.se <- function(x, w, na.rm=FALSE)
  #  Computes the variance of a weighted mean following Cochran 1977 definition
{
  if (na.rm) { w <- w[i <- !is.na(x)]; x <- x[i] }
  n = length(w)
  xWbar = weighted.mean(x,w,na.rm=na.rm)
  wbar = mean(w)
  out = n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
  return(out)
}

ward_rates <- read_csv("FreqDens_ward_rates.csv"
                       , locale = locale(tz = "CET"))
newID_COVIDstat <- read_csv("../R_objects/typeID_COVIDstat.csv")

newID_ref <- newID_COVIDstat %>% 
  select(ward_id, newID) %>% 
  arrange(newID) %>% 
  rbind(tibble(ward_id = 0, newID = "All wards")) %>% 
  mutate(newID = factor(newID, levels = newID))


# read results of Bayesian analysis
ward_nonlinear_main <- read_csv(paste0("output/ward_nonlinear_bayesIntensity_multiTarget.csv"))

# ward_nonlinear_catHosp <- read_csv("output/ward_nonlinear_bayes_catHosp.csv") %>%
#   select(-any_of(c("phi", 'phi_opt_loglik', 'chi', 'chi_opt_loglik'))) %>% rename(phi_median = phi_bayes, chi_median = chi_bayes)


### calculate AIC tables ###

this_popcut = ward_nonlinear_main$popcut %>% unique
this_boundary = ward_nonlinear_main$this_boundary %>% unique %>% {.[. > 0]}


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
                            T ~ Period)) %>% # this should never apply in the post-processing
  filter(Nthat > 0, form_minutes > 0) # do the standard filtering for negligible time periods


ward_rates_loglik <- ward_rates_edit %>% 
  {rbind(., mutate(., ward_id = 0))} %>% 
  left_join(newID_ref) %>% 
  filter(!is.na(newID)) %>% 
  left_join(ward_nonlinear_main %>% filter(Period == "All") %>% transmute(ward_id, this_status, that_status, phi_nondiurnal = phi_median)) %>%
  left_join(ward_nonlinear_main %>% filter(Period != "All") %>% transmute(ward_id, Period, this_status, that_status, phi_diurnal = phi_median)) %>%
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
  


# ward_rates_loglik %>%
#   filter(this_status == "PA", that_status == "PA", ward_id %in% c(2, 5, 8, 13))

ward_aic = ward_rates_loglik %>% 
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
    
    # ward_aic %>% 
    #   filter(this_status == this_stat, that_status == that_stat)
    # ward_nonlinear_main %>% 
    #   filter(ward_id == 6) %>% 
    #   filter(this_status == this_stat, that_status == that_stat)
    tab <- ward_aic %>% 
      filter(ward_id != 15) %>% 
      filter(this_status == "ALL", that_status == "ALL") %>% 
      arrange(newID) %>% 
      transmute(Ward = newID, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
      mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))
    
    tab_this <- ward_aic %>% 
      filter(ward_id != 15) %>% 
      filter(this_status != "ALL", that_status == "ALL") %>% 
      arrange(newID) %>% 
      transmute(Ward = newID, this_status, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
      mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))
    
    tab_that <- ward_aic %>% 
      filter(ward_id != 15) %>% 
      filter(this_status == "ALL", that_status != "ALL") %>% 
      arrange(newID) %>% 
      transmute(Ward = newID, that_status, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
      mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))

    tab_thisthat <- ward_aic %>% 
      filter(ward_id != 15) %>% 
      filter(this_status != "ALL", that_status != "ALL") %>% 
      arrange(newID) %>% 
      transmute(Ward = newID, this_status, that_status, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
      mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))
    
    # tab_papepepa <- ward_aic %>% 
    #   filter(ward_id != 15) %>% 
    #   filter(this_status %in% c("PA", "PE"), that_status %in% c("PA", "PE")) %>% 
    #   arrange(newID) %>% 
    #   transmute(Ward = newID, this_status, that_status, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
    #   mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))
      
    tab %>%   
      write_tsv(paste0("AICtables/dAIC_intensity.tsv"))
      # write_tsv(paste0("AICtables/dAIC_intensity_popcut", this_popcut, "boundary", this_boundary, ".tsv"))
    
    # tab_this %>%   
    #   write_tsv(paste0("AICtables/dAIC_intensity_this_popcut", this_popcut, "boundary", this_boundary, ".tsv"))
    # 
    # tab_that %>%   
    #   write_tsv(paste0("AICtables/dAIC_intensity_that_popcut", this_popcut, "boundary", this_boundary, ".tsv"))
    
    tab_thisthat %>%   
      write_tsv(paste0("AICtables/dAIC_intensity_thisthat.tsv"))
    
    # tab_papepepa %>%   
    #   write_tsv(paste0("AICtables/dAIC_intensity_papepepa.tsv"))
    
  






ward_nonlinear_main



pair_cols = RColorBrewer::brewer.pal(6, name="Paired")

# rates_N_allModels <- ward_rates_edit %>% 
#   {rbind(., mutate(., ward_id = 0))} %>% 
#   left_join(newID_ref) %>% 
#   filter(!is.na(newID)) %>% 
#   filter(Nthat > 0) %>% 
#   left_join(ward_nonlinear_main %>% filter(Period == "All") %>% select(-Period, -ends_with("loglik"))) %>%
#   left_join(ward_nonlinear_main %>% filter(Period != "All") %>% transmute(ward_id, Period, this_status, that_status, eta_period = eta_median)) %>%
#   mutate(Nthat_cat = cut(Nthat, seq(0, by = 10, max(Nthat) + 10))) %>% 
#   # filter(this_status == "ALL", that_status == "ALL") %>% 
#   filter(form_minutes > 0) %>% 
#   # rename_func_bool(process) %>% 
#   # rename(event_rate = form_rate, number_events = number_forms, event_minutes = form_minutes, nonlin = phi, nonlin_period = phi_period) %>% mutate(nonlin_expr = paste0(Period, "\n(phi=", round(nonlin_period, 1), ")")) %>%
#   # rename(event_rate = break_rate, number_events = number_breaks, event_minutes = break_minutes, nonlin = chi, nonlin_period = chi_period) %>% mutate(nonlin_expr = paste0(Period, "\n(chi=", round(nonlin_period, 1), ")")) %>%
#   group_by(ward_id, Period, Nthat_cat) %>% ### PREPARE THE DATA ###
#   mutate(mean_Nthat = mean(Nthat)
#          , low_Nthat = min(Nthat)
#          , high_Nthat = max(Nthat)
#          , obs_intensity = sum(contactIntensity*form_minutes)/sum(form_minutes)
#          # , lo_event_rate = ifelse(sum(event_minutes) == 0, NaN, weighted.quantile(event_rate, event_minutes, prob = 0.25))
#          # , hi_event_rate = ifelse(sum(event_minutes) == 0, NaN, weighted.quantile(event_rate, event_minutes, prob = 0.75))
#          , var_intensity = ifelse(sum(contactIntensity) == 0, NaN, weighted.var(contactIntensity, form_minutes, na.rm = T))
#          # , mean_break_rate = sum(break_rate*break_minutes, na.rm = T)/sum(break_minutes)
#   ) %>% 
#   
#   group_by(ward_id, this_status, that_status) %>%
#   mutate(w_constant = sum(break_minutes)/sum(form_minutes)
#          , w_linear = sum(break_minutes)/sum(form_minutes*Nthat)
#          , w_nonlin = sum(break_minutes)/sum(form_minutes*Nthat^eta_median)) %>% 
#   group_by(ward_id, this_status, that_status, Period) %>% #group by day/night to calculate the day and night parameters
#   mutate(w_cdiurnal = sum(break_minutes)/sum(form_minutes)
#          , w_ldiurnal = sum(break_minutes)/sum(form_minutes*Nthat)
#          , w_ndiurnal = sum(break_minutes)/sum(form_minutes*Nthat^eta_period)) %>% 
#   mutate(rate_constant = w_constant #calculate the rates for both non-diurnal and diurnal models
#          , rate_linear = w_linear*Nthat
#          , rate_nonlin = w_nonlin*Nthat^eta_median
#          , rate_cdiurnal = w_cdiurnal
#          , rate_ldiurnal = w_ldiurnal*Nthat
#          , rate_ndiurnal = w_ndiurnal*Nthat^eta_period
#   )
  
  
plot_N_allModels <- ward_rates_loglik %>% 
  # plot_N_allModels <- rates_N_allModels %>% 
    # filter(plot_row == r) %>% 
    # mutate(newID = factor(newID, levels = levels(rates_N_allModels$newID))) %>% 
    # filter(Nthat > 0) %>% 
    # filter(form_minutes > 0) %>% 
    filter(this_status == "ALL", that_status == "ALL") %>% 
    ggplot(aes(x = Nthat)) + 
    geom_point(aes(y = contactIntensity, colour = ".data"), alpha = 0.2, size = 0.2) + 
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_event_rate, colour = "data"), size=1) + 
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_event_rate, colour = "data", size = weight_Nthat)) + 
    # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +  
    # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) + 
    geom_line(aes(y = rate_constant, colour = "Constant")) +
    geom_line(aes(y = rate_linear, colour = "Linear")) +
    geom_line(aes(y = rate_nonlin, colour = "Non-linear")) +
    geom_line(aes(y = rate_cdiurnal, colour = "Constant\ndiurnal")) +
    geom_line(aes(y = rate_ldiurnal, colour = "Linear\ndiurnal")) +
    geom_line(aes(y = rate_ndiurnal, colour = "Non-linear\ndiurnal")) +
    scale_colour_manual(values = c(`.data` = "blue"
                                   , `Constant` = pair_cols[1]
                                   , `Constant\ndiurnal` = pair_cols[2]
                                   , `Linear` = pair_cols[3]
                                   , `Linear\ndiurnal` = pair_cols[4]
                                   , `Non-linear` = pair_cols[5]
                                   , `Non-linear\ndiurnal` = pair_cols[6]
    ))+
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
    facet_nested_wrap(.~newID + Period, scales = "free_x",  nrow = 2) +
    # guides(colour = F) +
    scale_y_log10() +
    theme_bw() + 
    labs(y = "Contact minutes per person minute", x = "Persons present", colour = "Model") #+ 
    # coord_cartesian(ylim = c(1e-3, NA))
  
  plot_N_allModels
  
 

  greek = "phi"
  
  
  
  # plot_N_ndiurnal <- rates_N_allModels %>% 
  plot_N_ndiurnal <- ward_rates_loglik %>% 
    # filter(plot_row == r) %>% 
    # mutate(newID = factor(newID, levels = levels(rates_N_allModels$newID))) %>% 
    # filter(Nthat > 0) %>% 
    # filter(form_minutes > 0) %>% 
    filter(this_status == "ALL", that_status == "ALL") %>% 
    mutate(param_label = paste0(greek, "==", round(phi_diurnal, 1))) %>% 
    ggplot(aes(x = Nthat)) + 
    geom_point(aes(y = contactIntensity, colour = ".data"), alpha = 0.2, size = 0.2) + 
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_event_rate, colour = "data"), size=1) + 
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_event_rate, colour = "data", size = weight_Nthat)) + 
    # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +  
    # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) + 
    # geom_line(aes(y = rate_constant, colour = "Constant")) +
    # geom_line(aes(y = rate_linear, colour = "Linear")) +
    # geom_line(aes(y = rate_nonlin, colour = "Non-linear")) +
    # geom_line(aes(y = rate_cdiurnal, colour = "Constant\ndiurnal")) +
    # geom_line(aes(y = rate_ldiurnal, colour = "Linear\ndiurnal")) +
    geom_line(aes(y = rate_ndiurnal, colour = "Non-linear\ndiurnal")) +
    geom_text(aes(x = Inf, y = 3, label = param_label), parse = T, hjust = 1.0, vjust = 1) +
    scale_colour_manual(values = c(`.data` = "blue"
                                   , `Constant` = pair_cols[1]
                                   , `Constant\ndiurnal` = pair_cols[2]
                                   , `Linear` = pair_cols[3]
                                   , `Linear\ndiurnal` = pair_cols[4]
                                   , `Non-linear` = pair_cols[5]
                                   , `Non-linear\ndiurnal` = pair_cols[6]
    ))+
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
    facet_nested_wrap(.~newID + Period, scales = "free_x",  nrow = 2) +
    # guides(colour = F) +
    scale_y_log10() +
    theme_bw() + 
    labs(y = "Contact minutes per person minute", x = "Persons present", colour = "Model")

  
  
  # chosen_model <- read_tsv(paste0("AICtables/dAIC_intensity_cutoff1.tsv")) %>% 
  chosen_model <- read_tsv(paste0("AICtables/dAIC_intensity.tsv")) %>% 
    pivot_longer(-Ward, names_to = "ChosenModel", values_to = "dAIC") %>% 
    rename(newID = Ward) %>% 
    filter(dAIC == 0) %>% 
    select(-dAIC)
  
  # chosen_model_this <- read_tsv(paste0("AICtables/dAIC_intensity_this_popcut", this_popcut, "boundary", this_boundary, ".tsv")) %>% 
  #   pivot_longer(-c("Ward", "this_status"), names_to = "ChosenModel", values_to = "dAIC") %>% 
  #   rename(newID = Ward) %>% 
  #   filter(dAIC == 0) %>% 
  #   select(-dAIC)
  
  # chosen_model_that <- read_tsv(paste0("AICtables/dAIC_intensity_that_popcut", this_popcut, "boundary", this_boundary, ".tsv")) %>% 
  #   pivot_longer(-c("Ward", "that_status"), names_to = "ChosenModel", values_to = "dAIC") %>% 
  #   rename(newID = Ward) %>% 
  #   filter(dAIC == 0) %>% 
  #   select(-dAIC)
  
  chosen_model_thisthat <- read_tsv(paste0("AICtables/dAIC_intensity_thisthat.tsv")) %>% 
    pivot_longer(-c("Ward", "this_status", "that_status"), names_to = "ChosenModel", values_to = "dAIC") %>% 
    rename(newID = Ward) %>% 
    filter(dAIC == 0) %>% 
    select(-dAIC)
  
  
   
  toPlot_N_chosenModels <- ward_rates_loglik %>%
  # plot_N_chosenModels <- ward_rates_loglik %>% 
  # plot_N_chosenModels <- rates_N_allModels %>% 
    left_join(chosen_model) %>%
    mutate(newID = factor(newID, levels = newID_ref$newID)) %>% 
    # mutate(newID = factor(newID, levels = levels(rates_N_allModels$newID))) %>% 
    # filter(plot_row == r) %>% 
    # filter(Nthat > 0) %>% 
    # filter(form_minutes > 0) %>% 
    # filter(this_status == "ALL", that_status == "ALL") %>% 
    mutate(rate_chosenmodel = case_when(ChosenModel == "Constant" ~ rate_constant
                                        , ChosenModel == "Linear" ~ rate_linear
                                        , ChosenModel == "Non-linear" ~ rate_nonlin
                                        , ChosenModel == "Constant diurnal" ~ rate_cdiurnal
                                        , ChosenModel == "Linear diurnal" ~ rate_ldiurnal
                                        , ChosenModel == "Non-linear diurnal" ~ rate_ndiurnal
    )) %>% 
    mutate(param_label = case_when(ChosenModel == "Non-linear" ~ paste0(greek, "==", round(phi_nondiurnal, 1))
                                   , ChosenModel == "Non-linear diurnal" ~ paste0(greek, "==", round(phi_diurnal, 1))
                                   , ChosenModel == "Constant" ~ paste0(greek, "==", 0)
                                   , ChosenModel == "Constant diurnal" ~ paste0(greek, "==", 0)
                                   , ChosenModel == "Linear" ~ paste0(greek, "==", 1)
                                   , ChosenModel == "Linear diurnal" ~ paste0(greek, "==", 1)
                                   , T ~ "")) %>% 
    mutate(model_label = case_when(ChosenModel == "Constant" ~ "Frequency dependent"
                                   , ChosenModel == "Constant diurnal" ~ "Frequency dependent - diurnal"
                                   , ChosenModel == "Linear" ~ "Linear density dependent"
                                   , ChosenModel == "Linear diurnal" ~ "Linear density dependent - diurnal"
                                   , ChosenModel == "Non-linear" ~ "Non-linear"
                                   , ChosenModel == "Non-linear diurnal" ~ "Non-linear - diurnal"
    )) %>% 
    mutate(model_label = factor(model_label, levels = c(".data"
                                                        , "Frequency dependent"
                                                        , "Frequency dependent - diurnal"
                                                        , "Linear density dependent"
                                                        , "Linear density dependent - diurnal"
                                                        , "Non-linear"
                                                        , "Non-linear - diurnal")))
  
  plot_N_chosenModels = toPlot_N_chosenModels %>% 
    filter(this_status == "ALL", that_status == "ALL") %>%
    ggplot(aes(x = Nthat)) + 
    geom_point(aes(y = contactIntensity, colour = ".data"), alpha = 0.3, size = 0.3) + 
    # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +
    # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) +
    geom_line(aes(x = Nthat, y = rate_chosenmodel, colour = model_label)) +
    geom_text(aes(x = Inf, y = 3, label = param_label), parse = T, hjust = 1.0, vjust = 1) +
    scale_colour_manual(values = c(`.data` = "blue"
                                   , `Frequency dependent` = pair_cols[1]
                                   , `Frequency dependent - diurnal` = pair_cols[2]
                                   , `Linear density dependent` = pair_cols[3]
                                   , `Linear density dependent - diurnal` = pair_cols[4]
                                   , `Non-linear` = pair_cols[5]
                                   , `Non-linear - diurnal` = pair_cols[6]
    ))+
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
    facet_nested_wrap(.~newID + Period
                      , scales = "free_x"
                      ,  nrow = 2) +
    # guides(colour = F) +
    scale_y_log10() +
    theme_bw() + 
    labs(y = "Contact rate", x = "Persons present", colour = "Model")

  
  
  ggsave(plot_N_chosenModels, filename = paste0("Graphs/rates_over_N_chosenModels_intensity.png")
         , width = 35, height = 20, units = "cm", dpi = 1000)
  
  
  ggsave(plot_N_allModels, filename = paste0("Graphs/rates_over_N_allModels_intensity.png"), width = 35, height = 20, units = "cm")
  
  ggsave(plot_N_ndiurnal, filename = paste0("Graphs/rates_over_N_ndiurnal_intensity.png"), width = 35, height = 20, units = "cm")
  
# this_stat = that_stat = "PA"
  
  model_levels = c(".data", 
                   "Non-linear (negative)", 
                   "Non-linear (negative) - diurnal", 
                   "Frequency dependence", 
                   "Frequency dependence - diurnal", 
                   "Non-linear (partial)", 
                   "Non-linear (partial) - diurnal", 
                   "Linear density dependence", 
                   "Linear density dependence - diurnal", 
                   "Non-linear (more than linear)", 
                   "Non-linear (more than linear) - diurnal")
  
  model_colours = c(`.data` = "black", 
    `Non-linear (negative)` = alpha("darkblue", 0.6), 
    `Non-linear (negative) - diurnal` = alpha("darkblue", 1), 
    `Frequency dependence` = alpha("blue", 0.6), 
    `Frequency dependence - diurnal` = alpha("blue", 1), 
    `Non-linear (partial)` = alpha("green", 0.6), 
    `Non-linear (partial) - diurnal` = alpha("green", 1), 
    `Linear density dependence` = alpha("gold", 0.6), 
    `Linear density dependence - diurnal` = alpha("gold", 1), 
    `Non-linear (more than linear)` = alpha("orange", 0.6), 
    `Non-linear (more than linear) - diurnal` = alpha("orange", 1))
    
  
  toPlot_N_chosenModels_Slope <- ward_rates_loglik %>%
    # plot_N_chosenModels <- ward_rates_loglik %>% 
    # plot_N_chosenModels <- rates_N_allModels %>% 
    left_join(chosen_model) %>%
    mutate(newID = factor(newID, levels = newID_ref$newID)) %>% 
    # mutate(newID = factor(newID, levels = levels(rates_N_allModels$newID))) %>% 
    # filter(plot_row == r) %>% 
    # filter(Nthat > 0) %>% 
    # filter(form_minutes > 0) %>% 
    # filter(this_status == "ALL", that_status == "ALL") %>% 
    mutate(rate_chosenmodel = case_when(ChosenModel == "Constant" ~ rate_constant
                                        , ChosenModel == "Linear" ~ rate_linear
                                        , ChosenModel == "Non-linear" ~ rate_nonlin
                                        , ChosenModel == "Constant diurnal" ~ rate_cdiurnal
                                        , ChosenModel == "Linear diurnal" ~ rate_ldiurnal
                                        , ChosenModel == "Non-linear diurnal" ~ rate_ndiurnal
    )) %>% 
    mutate(param_label = case_when(ChosenModel == "Non-linear" ~ paste0(greek, "==", round(phi_nondiurnal, 1))
                                   , ChosenModel == "Non-linear diurnal" ~ paste0(greek, "==", round(phi_diurnal, 1))
                                   , ChosenModel == "Constant" ~ paste0(greek, "==", 0)
                                   , ChosenModel == "Constant diurnal" ~ paste0(greek, "==", 0)
                                   , ChosenModel == "Linear" ~ paste0(greek, "==", 1)
                                   , ChosenModel == "Linear diurnal" ~ paste0(greek, "==", 1)
                                   , T ~ "")) %>% 
    mutate(model_label = case_when(ChosenModel == "Constant" ~ "Frequency dependence"
                                   , ChosenModel == "Constant diurnal" ~ "Frequency dependence - diurnal"
                                   , ChosenModel == "Linear" ~ "Linear density dependence"
                                   , ChosenModel == "Linear diurnal" ~ "Linear density dependence - diurnal"
                                   , ChosenModel == "Non-linear" & phi_nondiurnal < 0 ~ "Non-linear (negative)"
                                   , ChosenModel == "Non-linear" & phi_nondiurnal < 1 ~ "Non-linear (partial)"
                                   , ChosenModel == "Non-linear" ~ "Non-linear (more than linear)"
                                   , ChosenModel == "Non-linear diurnal" & phi_diurnal < 0 ~ "Non-linear (negative) - diurnal"
                                   , ChosenModel == "Non-linear diurnal" & phi_diurnal < 1 ~ "Non-linear (partial) - diurnal"
                                   , ChosenModel == "Non-linear diurnal" ~ "Non-linear (more than linear) - diurnal"
    )) %>% 
    mutate(model_label = factor(model_label, levels = model_levels))
  
  
  
  
  plot_N_chosenModels_Slope = toPlot_N_chosenModels_Slope %>% 
    filter(this_status == "ALL", that_status == "ALL") %>%
    ggplot(aes(x = Nthat)) + 
    geom_point(aes(y = contactIntensity, colour = factor(".data", levels = model_levels)), alpha = 0.3, size = 0.3) + 
    # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +
    # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) +
    geom_line(aes(y = rate_chosenmodel, colour = model_label), linewidth = 1.5) +
    geom_text(aes(x = Inf, y = 3, label = param_label), parse = T, hjust = 1.0, vjust = 1) +
    scale_colour_manual(drop = F, values = model_colours) +
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
    facet_nested_wrap(.~newID + Period
                      , scales = "free_x"
                      ,  nrow = 2) +
    # guides(colour = F) +
    scale_y_log10() +
    theme_bw() + 
    labs(y = "Contact rate", x = "Persons present", colour = "Model")
  
  ggsave(plot_N_chosenModels_Slope, filename = paste0("Graphs/rates_over_N_chosenModels_intensitySlope.png")
         , width = 35, height = 20, units = "cm", dpi = 1000)
  

for(this_stat in c("PA", "PE")){
  for(that_stat in c("PA", "PE")){
    
    plot_N_chosenModels_Slope_thisthat <- toPlot_N_chosenModels_Slope %>% 
      filter(this_status == this_stat, that_status == that_stat) %>%
      ggplot(aes(x = Nthat)) + 
      geom_point(aes(y = contactIntensity, colour = factor(".data", levels = model_levels)), alpha = 0.3, size = 0.3) + 
      # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +
      # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) +
      geom_line(aes(y = rate_chosenmodel, colour = model_label), linewidth = 1.5) +
      geom_text(aes(x = Inf, y = 3, label = param_label), parse = T, hjust = 1.0, vjust = 1) +
      scale_colour_manual(drop = F, values = model_colours) + 
                          # values = c(`.data` = "black", 
                          #            `Non-linear (negative)` = alpha("darkblue", 0.6), 
                          #            `Non-linear (negative) - diurnal` = alpha("darkblue", 1), 
                          #            `Frequency dependence` = alpha("blue", 0.6), 
                          #            `Frequency dependence - diurnal` = alpha("blue", 1), 
                          #            `Non-linear (partial)` = alpha("green", 0.6), 
                          #            `Non-linear (partial) - diurnal` = alpha("green", 1), 
                          #            `Linear density dependence` = alpha("gold", 0.6), 
                          #            `Linear density dependence - diurnal` = alpha("gold", 1), 
                          #            `Non-linear (more than linear)` = alpha("orange", 0.6), 
                          #            `Non-linear (more than linear) - diurnal` = alpha("orange", 1))
      # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
      facet_nested_wrap(.~newID + Period
                        , scales = "free_x"
                        ,  nrow = 2) +
      # guides(colour = F) +
      scale_y_log10() +
      theme_bw() + 
      labs(y = "Contact rate", x = "Persons present", colour = "Model")
    
    

    # plot_N_chosenModels_thisthat <- toPlot_N_chosenModels %>% 
    #   filter(this_status == this_stat, that_status == that_stat) %>%
    #   ggplot(aes(x = Nthat)) + 
    #   geom_point(aes(y = contactIntensity, colour = ".data"), alpha = 0.3, size = 0.3) + 
    #   # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +
    #   # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) +
    #   geom_line(aes(x = Nthat, y = rate_chosenmodel, colour = model_label)) +
    #   geom_text(aes(x = Inf, y = 3, label = param_label), parse = T, hjust = 1.0, vjust = 1) +
    #   scale_colour_manual(values = c(`.data` = "blue"
    #                                  , `Frequency dependent` = pair_cols[1]
    #                                  , `Frequency dependent - diurnal` = pair_cols[2]
    #                                  , `Linear density dependent` = pair_cols[3]
    #                                  , `Linear density dependent - diurnal` = pair_cols[4]
    #                                  , `Non-linear` = pair_cols[5]
    #                                  , `Non-linear - diurnal` = pair_cols[6]
    #   ))+
    #   # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
    #   facet_nested_wrap(.~newID + Period
    #                     , scales = "free_x"
    #                     ,  nrow = 2) +
    #   # guides(colour = F) +
    #   scale_y_log10() +
    #   theme_bw() + 
    #   labs(y = "Contact rate", x = "Persons present", colour = "Model")
    
    # 
    # plot_N_chosenModels_thisthat <- rates_N_allModels %>% 
    #   # filter(plot_row == r) %>% 
    #   filter(Nthat > 0) %>% 
    #   filter(form_minutes > 0) %>% 
    #   filter(this_status == this_stat, that_status == that_stat) %>% 
    #   left_join(chosen_model_thisthat) %>%
    #   mutate(newID = factor(newID, levels = levels(rates_N_allModels$newID))) %>% 
    #   mutate(rate_chosenmodel = case_when(ChosenModel == "Constant" ~ rate_constant
    #                                       , ChosenModel == "Linear" ~ rate_linear
    #                                       , ChosenModel == "Non-linear" ~ rate_nonlin
    #                                       , ChosenModel == "Constant diurnal" ~ rate_cdiurnal
    #                                       , ChosenModel == "Linear diurnal" ~ rate_ldiurnal
    #                                       , ChosenModel == "Non-linear diurnal" ~ rate_ndiurnal
    #   )) %>% 
    #   mutate(param_label = case_when(ChosenModel == "Non-linear" ~ paste0(greek, "==", round(eta_median, 1))
    #                                  , ChosenModel == "Non-linear diurnal" ~ paste0(greek, "==", round(eta_period, 1))
    #                                  , ChosenModel == "Constant" ~ paste0(greek, "==", 0)
    #                                  , ChosenModel == "Constant diurnal" ~ paste0(greek, "==", 0)
    #                                  , ChosenModel == "Linear" ~ paste0(greek, "==", 1)
    #                                  , ChosenModel == "Linear diurnal" ~ paste0(greek, "==", 1)
    #                                  , T ~ "")) %>% 
    #   mutate(model_label = case_when(ChosenModel == "Constant" ~ "Frequency dependent"
    #                                  , ChosenModel == "Constant diurnal" ~ "Frequency dependent - diurnal"
    #                                  , ChosenModel == "Linear" ~ "Linear density dependent"
    #                                  , ChosenModel == "Linear diurnal" ~ "Linear density dependent - diurnal"
    #                                  , ChosenModel == "Non-linear" ~ "Non-linear"
    #                                  , ChosenModel == "Non-linear diurnal" ~ "Non-linear - diurnal"
    #   )) %>% 
    #   mutate(model_label = factor(model_label, levels = c(".data"
    #                                                       , "Frequency dependent"
    #                                                       , "Frequency dependent - diurnal"
    #                                                       , "Linear density dependent"
    #                                                       , "Linear density dependent - diurnal"
    #                                                       , "Non-linear"
    #                                                       , "Non-linear - diurnal"))) %>% 
    #   ggplot(aes(x = Nthat)) + 
    #   geom_point(aes(y = contactIntensity, colour = ".data"), alpha = 0.3, size = 0.3) + 
    #   # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +
    #   # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) +
    #   geom_line(aes(x = Nthat, y = rate_chosenmodel, colour = model_label)) +
    #   geom_text(aes(x = Inf, y = 3, label = param_label), parse = T, hjust = 1.0, vjust = 1) +
    #   scale_colour_manual(values = c(`.data` = "blue"
    #                                  , `Frequency dependent` = pair_cols[1]
    #                                  , `Frequency dependent - diurnal` = pair_cols[2]
    #                                  , `Linear density dependent` = pair_cols[3]
    #                                  , `Linear density dependent - diurnal` = pair_cols[4]
    #                                  , `Non-linear` = pair_cols[5]
    #                                  , `Non-linear - diurnal` = pair_cols[6]
    #   ))+
    #   # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
    #   facet_nested_wrap(.~newID + Period
    #                     , scales = "free_x"
    #                     ,  nrow = 2) +
    #   # guides(colour = F) +
    #   scale_y_log10() +
    #   theme_bw() + 
    #   labs(y = "Contact rate", x = "Persons present", colour = "Model")
    
    ggsave(plot_N_chosenModels_Slope_thisthat, filename = paste0("Graphs/rates_over_N_chosenModels_intensity_", this_stat, that_stat, ".png")
           , width = 35, height = 20, units = "cm", dpi = 1000)
    
  }
}  
  

  plot_N_allModels_oneward <- ward_rates_loglik %>% 
  # plot_N_allModels_oneward <- rates_N_allModels %>% 
    filter(ward_id == 16) %>% 
    # filter(plot_row == r) %>% 
    # mutate(newID = factor(newID, levels = levels(rates_N_allModels$newID))) %>% 
    # filter(Nthat > 0) %>% 
    # filter(form_minutes > 0) %>% 
    filter(this_status == "ALL", that_status == "ALL") %>% 
    mutate(Period_label = ifelse(Period == "Day", "Day (08:00-19:59)", "Night (20:00-07:59)")) %>% 
    ggplot(aes(x = Nthat)) + 
    geom_point(aes(y = contactIntensity, colour = ".data"), alpha = 0.4, size = 0.5) + 
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_event_rate, colour = "data"), size=1) + 
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_event_rate, colour = "data", size = weight_Nthat)) + 
    # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +  
    # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) + 
    geom_line(aes(y = rate_constant, colour = "Frequency dependent")) +
    geom_line(aes(y = rate_linear, colour = "Linear density dependent")) +
    geom_line(aes(y = rate_nonlin, colour = "Non-linear")) +
    geom_line(aes(y = rate_cdiurnal, colour = "Frequency dependent - diurnal")) +
    geom_line(aes(y = rate_ldiurnal, colour = "Linear density dependent - diurnal")) +
    geom_line(aes(y = rate_ndiurnal, colour = "Non-linear - diurnal")) +
    scale_colour_manual(values = c(`.data` = "blue"
                                   , `Frequency dependent` = pair_cols[1]
                                   , `Frequency dependent - diurnal` = pair_cols[2]
                                   , `Linear density dependent` = pair_cols[3]
                                   , `Linear density dependent - diurnal` = pair_cols[4]
                                   , `Non-linear` = pair_cols[5]
                                   , `Non-linear - diurnal` = pair_cols[6]
    ))+
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
    facet_grid(.~Period_label) +
    # guides(colour = F) +
    # scale_y_log10() +
    theme_bw() + 
    labs(y = "Contact rate", x = "Persons present", colour = "Model") + 
    theme(legend.title = element_text(size = 15))
  # coord_cartesian(ylim = c(1e-3, NA))
  
  
  ggsave(plot_N_allModels_oneward, filename = paste0("Graphs/plot_N_allModels_oneward.png"), width = 20, height = 10, units = "cm", dpi = 1000)
  


### NOW PLOT THE ETA VALUES THEMSELVES ###

error_bar_width = 0.5
  
ward_nonlinear_plottable = ward_nonlinear_main %>% 
  left_join(newID_ref) %>% 
  filter(!is.na(newID)) %>% 
  select(ward_id, newID, this_status, that_status, Period, phi_median, phi_lo, phi_hi)



grade_plot = ward_nonlinear_plottable %>% 
  pivot_wider(id_cols = c("newID", "this_status", "that_status"), names_from = "Period", values_from = c("phi_median", "phi_lo", "phi_hi")) %>%
  left_join(rbind(chosen_model %>% mutate(this_status = "ALL", that_status = "ALL")
                  , chosen_model_thisthat)) %>% 
  mutate(newID = factor(newID, levels = newID_ref$newID)) %>% 
  mutate(Day_median = case_when(grepl("Constant", ChosenModel) ~ 0,
                                grepl("Linear", ChosenModel) ~ 1,
                                ChosenModel == "Non-linear diurnal" ~ phi_median_Day,
                                ChosenModel == "Non-linear" ~ phi_median_All, 
                                T ~ NA_real_), 
         Night_median = case_when(grepl("Constant", ChosenModel) ~ 0,
                                  grepl("Linear", ChosenModel) ~ 1,
                                  ChosenModel == "Non-linear diurnal" ~ phi_median_Night,
                                  ChosenModel == "Non-linear" ~ phi_median_All, 
                                  T ~ NA_real_),
         Day_lo = case_when(ChosenModel == "Non-linear diurnal" ~ phi_lo_Day,
                            ChosenModel == "Non-linear" ~ phi_lo_All, 
                            T ~ NA_real_),
         Day_hi = case_when(ChosenModel == "Non-linear diurnal" ~ phi_hi_Day,
                            ChosenModel == "Non-linear" ~ phi_hi_All, 
                            T ~ NA_real_),
         Night_lo = case_when(ChosenModel == "Non-linear diurnal" ~ phi_lo_Night,
                              ChosenModel == "Non-linear" ~ phi_lo_All, 
                              T ~ NA_real_),
         Night_hi = case_when(ChosenModel == "Non-linear diurnal" ~ phi_hi_Night,
                              ChosenModel == "Non-linear" ~ phi_hi_All, 
                              T ~ NA_real_),
         Day_colour = case_when(grepl("Constant", ChosenModel) ~ "Frequency dependence", 
                                grepl("Linear", ChosenModel) ~ "Linear density dependence", 
                                Day_median < 0 ~ "Non-linear (negative)", 
                                Day_median > 1 ~ "Non-linear (more than linear)",
                                Day_median > 0 & Day_median < 1 ~ "Non-linear (partial)",
                                T ~ NA_character_),
         Night_colour = case_when(grepl("Constant", ChosenModel) ~ "Frequency dependence", 
                                  grepl("Linear", ChosenModel) ~ "Linear density dependence", 
                                  Night_median < 0 ~ "Non-linear (negative)", 
                                  Night_median > 1 ~ "Non-linear (more than linear)",
                                  Night_median > 0 & Night_median < 1 ~ "Non-linear (partial)", 
                                  T ~ NA_character_)
  ) %>% 
  select(newID, this_status, that_status, ChosenModel, starts_with("Day"), starts_with("Night")) %>% 
  pivot_longer(-c("newID", "this_status", "that_status", "ChosenModel"), names_to = c("Period", ".value"), names_sep = "_") %>% 
  mutate(isDiurnal = ifelse(grepl("diurnal", ChosenModel), "Diurnal model", "Non-diurnal model")) %>% 
  mutate(colour = factor(colour, levels = c("Non-linear (negative)", 
                                            "Frequency dependence", 
                                            "Non-linear (partial)", 
                                            "Linear density dependence", 
                                            "Non-linear (more than linear)")))


grade_plot %>% 
  filter(this_status == "ALL", that_status == "ALL") %>% 
  mutate(colour = factor(colour, levels = c("Non-linear (negative)", 
                                            "Frequency dependence", 
                                            "Non-linear (partial)", 
                                            "Linear density dependence", 
                                            "Non-linear (more than linear)"))) %>% 
  ggplot(aes(y = newID)) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  geom_point(aes(x = median, colour = colour, shape = isDiurnal)) + 
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = colour), height=error_bar_width) +
  facet_nested(.~Period) + 
  theme_bw() + 
  labs(y = "", x = expression(varphi), colour = "Best fit model", shape = "") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(values = c(`Non-linear (negative)` = "darkblue", 
                                 `Frequency dependence` = "blue", 
                                 `Non-linear (partial)` = "green", 
                                 `Linear density dependence` = "gold", 
                                 `Non-linear (more than linear)` = "orange")) +
  scale_shape_manual(values = c(`Diurnal model` = "circle", `Non-diurnal model` = "square")) + 
  coord_cartesian(xlim = c(-1, 2.8))

ggsave(paste0("Graphs/phi_daynight.png"), width = 22, height = 12, units = "cm", dpi = 1000)


grade_plot %>% 
  filter(this_status == "ALL", that_status == "ALL") %>% 
  filter(Period == "Day") %>% 
  mutate(colour = factor(colour, levels = c("Non-linear (negative)", 
                                            "Frequency dependence", 
                                            "Non-linear (partial)", 
                                            "Linear density dependence", 
                                            "Non-linear (more than linear)"))) %>% 
  ggplot(aes(y = newID)) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  geom_point(aes(x = median, colour = colour, shape = isDiurnal)) + 
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = colour), height=error_bar_width) +
  facet_nested(.~Period) + 
  theme_bw() + 
  labs(y = "", x = expression(varphi), colour = "Best fit model", shape = "") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(drop = F, values = c(`Non-linear (negative)` = "darkblue", 
                                 `Frequency dependence` = "blue", 
                                 `Non-linear (partial)` = "green", 
                                 `Linear density dependence` = "gold", 
                                 `Non-linear (more than linear)` = "orange")) +
  scale_shape_manual(values = c(`Diurnal model` = "circle", `Non-diurnal model` = "square")) + 
  coord_cartesian(xlim = c(-1, 2))

ggsave(paste0("Graphs/phi_day.png"), width = 16, height = 12, units = "cm", dpi = 1000)




grade_plot %>% 
  mutate(thisstat_label = case_when(this_status == "PA" ~ "Patient"
                                    , this_status == "PE" ~ "HCW"
                                    , this_status == "V" ~ "Visitor"
                                    , this_status == "ALL" ~ "All")) %>% 
  mutate(thatstat_label = case_when(that_status == "PA" ~ "Patient"
                                    , that_status == "PE" ~ "HCW"
                                    , that_status == "V" ~ "Visitor"
                                    , that_status == "ALL" ~ "All")) %>% 
  mutate(thisthat_label = paste0(thisstat_label, ">", thatstat_label)) %>% 
  # filter(this_status %in% c("PA", "PE"), that_status %in% c("PA", "PE")) %>% 
  filter((this_status == "ALL" & that_status == "ALL") |(this_status %in% c("PA", "PE") & that_status %in% c("PA"))) %>% 
  ggplot(aes(y = newID)) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  geom_point(aes(x = median, colour = colour, shape = isDiurnal)) + 
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = colour), height=error_bar_width) +
  facet_nested(.~thisthat_label + Period) + 
  theme_bw() + 
  labs(y = "", x = expression(varphi), colour = "Best fit model", shape = "") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(values = c(`Non-linear (negative)` = "darkblue", 
                                 `Frequency dependence` = "blue", 
                                 `Non-linear (partial)` = "green", 
                                 `Linear density dependence` = "gold", 
                                 `Non-linear (more than linear)` = "orange")) +
  scale_shape_manual(values = c(`Diurnal model` = "circle", `Non-diurnal model` = "square")) + 
  geom_hline(aes(yintercept = 1.5), linetype = "solid")

ggsave(paste0("Graphs/phi_daynight_ALL_toPAchosen.png"), width = 30, height = 12, units = "cm")
# ggsave(paste0("Graphs/phi_daynight_toPAchosen.png"), width = 20, height = 12, units = "cm")
# ggsave(paste0("Graphs/phi_daynight_thisthatchosen.png"), width = 30, height = 12, units = "cm")

grade_plot %>% 
  mutate(thisstat_label = case_when(this_status == "PA" ~ "Patient"
                                    , this_status == "PE" ~ "HCW"
                                    , this_status == "V" ~ "Visitor"
                                    , this_status == "ALL" ~ "All")) %>% 
  mutate(thatstat_label = case_when(that_status == "PA" ~ "Patient"
                                    , that_status == "PE" ~ "HCW"
                                    , that_status == "V" ~ "Visitor"
                                    , that_status == "ALL" ~ "All")) %>% 
  mutate(thisthat_label = paste0(thisstat_label, ">", thatstat_label)) %>% 
  filter(this_status %in% c("PE"), that_status %in% c("PA")) %>%
  # filter(this_status %in% c("PA", "PE"), that_status %in% c("PA", "PE")) %>%
  # filter((this_status == "ALL" & that_status == "ALL") |(this_status %in% c("PA", "PE") & that_status %in% c("PA"))) %>% 
  ggplot(aes(y = newID)) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  geom_point(aes(x = median, colour = colour, shape = isDiurnal)) + 
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = colour), height=error_bar_width) +
  facet_nested(.~thisthat_label + Period) + 
  theme_bw() + 
  labs(y = "", x = expression(varphi), colour = "Best fit model", shape = "") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(values = c(`Non-linear (negative)` = "darkblue", 
                                 `Frequency dependence` = "blue", 
                                 `Non-linear (partial)` = "green", 
                                 `Linear density dependence` = "gold", 
                                 `Non-linear (more than linear)` = "orange")) +
  scale_shape_manual(values = c(`Diurnal model` = "circle", `Non-diurnal model` = "square")) + 
  geom_hline(aes(yintercept = 1.5), linetype = "solid") + 
  coord_cartesian(xlim = c(-1, 2.8))


# ggsave(paste0("Graphs/phi_daynight_ALL_toPAchosen.png"), width = 30, height = 12, units = "cm")
ggsave(paste0("Graphs/phi_daynight_PEtoPAchosen.png"), width = 16, height = 12, units = "cm")
# ggsave(paste0("Graphs/phi_daynight_thisthatchosen.png"), width = 30, height = 12, units = "cm")

grade_plot %>% 
  mutate(thisstat_label = case_when(this_status == "PA" ~ "Patient"
                                    , this_status == "PE" ~ "HCW"
                                    , this_status == "V" ~ "Visitor"
                                    , this_status == "ALL" ~ "All")) %>% 
  mutate(thatstat_label = case_when(that_status == "PA" ~ "Patient"
                                    , that_status == "PE" ~ "HCW"
                                    , that_status == "V" ~ "Visitor"
                                    , that_status == "ALL" ~ "All")) %>% 
  mutate(thisthat_label = paste0(thisstat_label, ">", thatstat_label)) %>% 
  filter(this_status %in% c("PA", "PE"), that_status %in% c("PA")) %>%
  # filter(this_status %in% c("PA", "PE"), that_status %in% c("PA", "PE")) %>%
  # filter((this_status == "ALL" & that_status == "ALL") |(this_status %in% c("PA", "PE") & that_status %in% c("PA"))) %>% 
  ggplot(aes(y = newID)) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  geom_point(aes(x = median, colour = colour, shape = isDiurnal)) + 
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = colour), height=error_bar_width) +
  facet_nested(.~thisthat_label + Period) + 
  theme_bw() + 
  labs(y = "", x = expression(varphi), colour = "Best fit model", shape = "") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(values = c(`Non-linear (negative)` = "darkblue", 
                                 `Frequency dependence` = "blue", 
                                 `Non-linear (partial)` = "green", 
                                 `Linear density dependence` = "gold", 
                                 `Non-linear (more than linear)` = "orange")) +
  scale_shape_manual(values = c(`Diurnal model` = "circle", `Non-diurnal model` = "square")) + 
  geom_hline(aes(yintercept = 1.5), linetype = "solid") + 
  coord_cartesian(xlim = c(-1, 2.8))

# ggsave(paste0("Graphs/phi_daynight_ALL_toPAchosen.png"), width = 30, height = 12, units = "cm")
ggsave(paste0("Graphs/phi_daynight_toPAchosen.png"), width = 22, height = 12, units = "cm")
# ggsave(paste0("Graphs/phi_daynight_thisthatchosen.png"), width = 30, height = 12, units = "cm")



grade_plot %>% 
  mutate(thisstat_label = case_when(this_status == "PA" ~ "Patient"
                                    , this_status == "PE" ~ "HCW"
                                    , this_status == "V" ~ "Visitor"
                                    , this_status == "ALL" ~ "All")) %>% 
  mutate(thatstat_label = case_when(that_status == "PA" ~ "Patient"
                                    , that_status == "PE" ~ "HCW"
                                    , that_status == "V" ~ "Visitor"
                                    , that_status == "ALL" ~ "All")) %>% 
  mutate(thisthat_label = paste0(thisstat_label, ">", thatstat_label)) %>% 
  # filter(this_status %in% c("PA", "PE"), that_status %in% c("PA", "PE")) %>% 
  filter((this_status == "ALL" & that_status == "ALL") |(this_status %in% c("PA", "PE") & that_status %in% c("PA"))) %>% 
  left_join(ward_rates_loglik %>% 
              group_by(newID, that_status) %>% 
              summarise(mean_Nthat = sum(Nthat*interval_durmin)/sum(interval_durmin))) %>% 
  mutate(mean_Nthat_110 = mean_Nthat*1.1, 
         # increase_110 = mean_Nthat_110^median/mean_Nthat^median - 1, 
         increase_110 = 1.1^median - 1, 
         ) %>% 
  mutate(increase_110 = increase_110 + ifelse(!is.na(colour) & median == 0, 0.01, 0)) %>% 
  mutate(newID = fct_rev(newID)) %>% 
  ggplot() + 
  # geom_vline(aes(xintercept = 0), linetype = "solid")+ 
  geom_bar(aes(x = newID, y = increase_110, fill = colour), stat = "identity") + 
  coord_flip() + 
  facet_nested(.~thisthat_label + Period) + 
  theme_bw() + 
  labs(x = "", y = "Increase in contact rate with 10% increase in population density", fill = "Best fit model") + 
  scale_fill_manual(values = c(`Non-linear (negative)` = "darkblue", 
                                 `Frequency dependence` = "blue", 
                                 `Non-linear (partial)` = "green", 
                                 `Linear density dependence` = "gold", 
                                 `Non-linear (more than linear)` = "orange")) +
  geom_vline(aes(xintercept = 1.5), linetype = "solid") + 
  scale_y_continuous(labels = scales::percent)
  
ggsave(paste0("Graphs/phi_daynight_ALL_toPAchosen_increase110.png"), width = 30, height = 12, units = "cm")




#### catHosp

ward_nonlinear_catHosp <- read_csv(paste0("output/ward_nonlinear_bayesIntensity_catHosp.csv"))

ward_nonlinear_catHosp_plottable = ward_nonlinear_catHosp %>% 
  left_join(newID_ref) %>% 
  filter(!is.na(newID)) %>% 
  select(ward_id, newID, this_catHosp, that_catHosp, Period, phi_median, phi_lo, phi_hi)

source("freqdens_catHosp_chosenModel.R")

grade_catHosp_plot = ward_nonlinear_catHosp_plottable %>% 
  pivot_wider(id_cols = c("newID", "this_catHosp", "that_catHosp"), names_from = "Period", values_from = c("phi_median", "phi_lo", "phi_hi")) %>%
  left_join(chosen_model_catHosp_allthisthat) %>% 
  mutate(newID = factor(newID, levels = newID_ref$newID)) %>% 
  mutate(Day_median = case_when(grepl("Constant", ChosenModel) ~ 0,
                                grepl("Linear", ChosenModel) ~ 1,
                                ChosenModel == "Non-linear diurnal" ~ phi_median_Day,
                                ChosenModel == "Non-linear" ~ phi_median_All, 
                                T ~ NA_real_), 
         Night_median = case_when(grepl("Constant", ChosenModel) ~ 0,
                                  grepl("Linear", ChosenModel) ~ 1,
                                  ChosenModel == "Non-linear diurnal" ~ phi_median_Night,
                                  ChosenModel == "Non-linear" ~ phi_median_All, 
                                  T ~ NA_real_),
         Day_lo = case_when(ChosenModel == "Non-linear diurnal" ~ phi_lo_Day,
                            ChosenModel == "Non-linear" ~ phi_lo_All, 
                            T ~ NA_real_),
         Day_hi = case_when(ChosenModel == "Non-linear diurnal" ~ phi_hi_Day,
                            ChosenModel == "Non-linear" ~ phi_hi_All, 
                            T ~ NA_real_),
         Night_lo = case_when(ChosenModel == "Non-linear diurnal" ~ phi_lo_Night,
                              ChosenModel == "Non-linear" ~ phi_lo_All, 
                              T ~ NA_real_),
         Night_hi = case_when(ChosenModel == "Non-linear diurnal" ~ phi_hi_Night,
                              ChosenModel == "Non-linear" ~ phi_hi_All, 
                              T ~ NA_real_),
         Day_colour = case_when(grepl("Constant", ChosenModel) ~ "Frequency dependence", 
                                grepl("Linear", ChosenModel) ~ "Linear density dependence", 
                                Day_median < 0 ~ "Non-linear (negative)", 
                                Day_median > 1 ~ "Non-linear (more than linear)",
                                Day_median > 0 & Day_median < 1 ~ "Non-linear (partial)",
                                T ~ NA_character_),
         Night_colour = case_when(grepl("Constant", ChosenModel) ~ "Frequency dependence", 
                                  grepl("Linear", ChosenModel) ~ "Linear density dependence", 
                                  Night_median < 0 ~ "Non-linear (negative)", 
                                  Night_median > 1 ~ "Non-linear (more than linear)",
                                  Night_median > 0 & Night_median < 1 ~ "Non-linear (partial)", 
                                  T ~ NA_character_)
  ) %>% 
  select(newID, this_catHosp, that_catHosp, ChosenModel, starts_with("Day"), starts_with("Night")) %>% 
  pivot_longer(-c("newID", "this_catHosp", "that_catHosp", "ChosenModel"), names_to = c("Period", ".value"), names_sep = "_") %>% 
  mutate(isDiurnal = ifelse(grepl("diurnal", ChosenModel), "Diurnal model", "Non-diurnal model")) %>% 
  mutate(colour = factor(colour, levels = c("Non-linear (negative)", 
                                            "Frequency dependence", 
                                            "Non-linear (partial)", 
                                            "Linear density dependence", 
                                            "Non-linear (more than linear)")))

grade_catHosp_plot$this_catHosp %>% unique
grade_catHosp_plot$that_catHosp %>% unique


grade_catHosp_plot %>% 
  mutate(thiscH_label = case_when(this_catHosp == "nurse" ~ "Nurse"
                                    , this_catHosp == "aux nurse" ~ "Auxiliary nurse"
                                    , this_catHosp == "physician" ~ "Physician"
                                    , this_catHosp == "administration" ~ "Administration"
                                  , T ~ NA)) %>% 
  mutate(thatcH_label = case_when(that_catHosp == "patient" ~ "Patients"
                                    # , that_catHosp == "ALL" ~ "All"
                                  )) %>% 
  mutate(thisthat_label = paste0(thiscH_label, ">", thatcH_label)) %>% 
  filter(!is.na(thiscH_label), 
         !is.na(thatcH_label)
        ) %>% 
  mutate(thisthat_label = factor(thisthat_label
                                 , levels = c("Nurse>Patients", 
                                              "Auxiliary nurse>Patients", 
                                              "Physician>Patients", 
                                              "Administration>Patients"))) %>% 
  ggplot(aes(y = newID)) + 
  geom_point(aes(x = median, colour = colour, shape = isDiurnal)) + 
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = colour), height=error_bar_width) +
  facet_nested(.~thisthat_label + Period) + 
  theme_bw() + 
  labs(y = "", x = expression(phi), colour = "Best fit model", shape = "") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(values = c(`Non-linear (negative)` = "darkblue", 
                                 `Frequency dependence` = "blue", 
                                 `Non-linear (partial)` = "green", 
                                 `Linear density dependence` = "gold", 
                                 `Non-linear (more than linear)` = "orange")) +
  scale_shape_manual(values = c(`Diurnal model` = "circle", `Non-diurnal model` = "square")) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  geom_hline(aes(yintercept = 1.5), linetype = "solid")
  

ggsave(paste0("Graphs/phi_daynight_catHosp.png"), width = 25, height = 12, units = "cm")



grade_catHosp_plot$this_catHosp %>% unique
grade_catHosp_plot$that_catHosp %>% unique


grade_catHosp_plot %>% 
  mutate(thiscH_label = case_when(this_catHosp == "nurse" ~ "Nurse"
                                  , this_catHosp == "aux nurse" ~ "Auxiliary nurse"
                                  # , this_catHosp == "physician" ~ "Physician"
                                  # , this_catHosp == "administration" ~ "Administration"
                                  , T ~ NA)) %>% 
  mutate(thatcH_label = case_when(that_catHosp == "patient" ~ "Patients"
                                  # , that_catHosp == "ALL" ~ "All"
  )) %>% 
  mutate(thisthat_label = paste0(thiscH_label, ">", thatcH_label)) %>% 
  filter(!is.na(thiscH_label), 
         !is.na(thatcH_label)
  ) %>% 
  mutate(thisthat_label = factor(thisthat_label
                                 , levels = c("Nurse>Patients", 
                                              "Auxiliary nurse>Patients", 
                                              "Physician>Patients", 
                                              "Administration>Patients"))) %>% 
  ggplot(aes(y = newID)) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  geom_point(aes(x = median, colour = colour, shape = isDiurnal)) + 
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = colour), height=error_bar_width) +
  facet_nested(.~thisthat_label + Period) + 
  theme_bw() + 
  labs(y = "", x = expression(varphi), colour = "Best fit model", shape = "") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(values = c(`Non-linear (negative)` = "darkblue", 
                                 `Frequency dependence` = "blue", 
                                 `Non-linear (partial)` = "green", 
                                 `Linear density dependence` = "gold", 
                                 `Non-linear (more than linear)` = "orange")) +
  scale_shape_manual(values = c(`Diurnal model` = "circle", `Non-diurnal model` = "square")) + 
  geom_hline(aes(yintercept = 1.5), linetype = "solid") + 
  coord_cartesian(xlim = c(-2, 3))


ggsave(paste0("Graphs/phi_daynight_catHosp_nurseaux.png"), width = 24, height = 12, units = "cm")



grade_catHosp_plot %>% 
  mutate(thiscH_label = case_when(#this_catHosp == "nurse" ~ "Nurse"
                                  #, this_catHosp == "aux nurse" ~ "Auxiliary nurse"
                                  this_catHosp == "physician" ~ "Physician"
                                  , this_catHosp == "administration" ~ "Administration"
                                  , T ~ NA)) %>% 
  mutate(thatcH_label = case_when(that_catHosp == "patient" ~ "Patients"
                                  # , that_catHosp == "ALL" ~ "All"
  )) %>% 
  mutate(thisthat_label = paste0(thiscH_label, ">", thatcH_label)) %>% 
  filter(!is.na(thiscH_label), 
         !is.na(thatcH_label)
  ) %>% 
  mutate(thisthat_label = factor(thisthat_label
                                 , levels = c("Nurse>Patients", 
                                              "Auxiliary nurse>Patients", 
                                              "Physician>Patients", 
                                              "Administration>Patients"))) %>% 
  ggplot(aes(y = newID)) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  geom_point(aes(x = median, colour = colour, shape = isDiurnal)) + 
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = colour), height=error_bar_width) +
  facet_nested(.~thisthat_label + Period) + 
  theme_bw() + 
  labs(y = "", x = expression(varphi), colour = "Best fit model", shape = "") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(values = c(`Non-linear (negative)` = "darkblue", 
                                 `Frequency dependence` = "blue", 
                                 `Non-linear (partial)` = "green", 
                                 `Linear density dependence` = "gold", 
                                 `Non-linear (more than linear)` = "orange")) +
  scale_shape_manual(values = c(`Diurnal model` = "circle", `Non-diurnal model` = "square")) + 
  geom_hline(aes(yintercept = 1.5), linetype = "solid") + 
  coord_cartesian(xlim = c(-2, 3))


ggsave(paste0("Graphs/phi_daynight_catHosp_physadmin.png"), width = 24, height = 12, units = "cm")


# chosen_model
# ward_nonlinear_plottable %>% 
# filter(this_status == "ALL", that_status == "ALL")
#   
#   ward_nonlinear_plottable %>% 
#   filter(this_status == "ALL", that_status == "ALL") %>% 
#   filter(Period != "All") %>%
#   left_join(chosen_model) %>% 
#   mutate(newID = factor(newID, levels = newID_ref$newID)) %>% 
#   mutate(FD = ifelse(ChosenModel == "Constant" | ChosenModel == "Constant diurnal", 0, NA_real_)
#          , LDD = ifelse(ChosenModel == "Linear" | ChosenModel == "Linear diurnal", 1, NA_real_)
#          , NL = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", phi_median, NA_real_)
#          , NL_lo = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", phi_lo, NA_real_)
#          , NL_hi = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", phi_hi, NA_real_)) %>% 
#   ggplot(aes(y = newID)) + 
#   # geom_point(aes(x = eta_median), colour = "darkgrey") + 
#   # geom_errorbarh(aes(xmin = eta_lo, xmax = eta_hi), height=error_bar_width, colour = "darkgrey") +
#   geom_point(aes(x = NL, colour = "Non-linear")) + 
#   geom_errorbarh(aes(xmin = NL_lo, xmax = NL_hi, colour = "Non-linear"), height=error_bar_width) +
#   geom_point(aes(x = FD, colour = "Frequency dependent")) + 
#   geom_point(aes(x = LDD, colour = "Linear density dependent")) + 
#   facet_nested(.~Period) + 
#   theme_bw() + 
#   labs(y = "", x = expression(phi), colour = "Best fit model") + 
#   scale_y_discrete(limits=rev) + 
#   scale_colour_manual(values = c(`Frequency dependent` = "blue", `Linear density dependent` = "green", `Non-linear` = "black")) + 
#   geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
#   geom_vline(aes(xintercept = 1), linetype = "dashed")
# 
# ggsave(paste0("Graphs/phi_daynight.png"), width = 20, height = 12, units = "cm", dpi = 1000)
# 
# 
# 
# 
# ward_nonlinear_plottable %>% 
#   filter(this_status %in% c("PA", "PE"), that_status %in% c("PA", "PE")) %>%
#   # filter((this_status %in% c("PA", "PE") & that_status %in% c("PA", "PE")) | (this_status == "ALL" & that_status == "ALL")) %>% 
#   filter(Period != "All") %>%
#   left_join(chosen_model_thisthat) %>%
#   # left_join(rbind(chosen_model %>% mutate(this_status = "ALL", that_status = "ALL")
#   #                 , chosen_model_thisthat)) %>% 
#   mutate(newID = factor(newID, levels = newID_ref$newID)) %>% 
#   mutate(FD = ifelse(ChosenModel == "Constant" | ChosenModel == "Constant diurnal", 0, NA_real_)
#          , DD = ifelse(ChosenModel == "Linear" | ChosenModel == "Linear diurnal", 1, NA_real_)
#          , NL = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", phi_median, NA_real_)
#          , NL_lo = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", phi_lo, NA_real_)
#          , NL_hi = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", phi_hi, NA_real_)) %>% 
#   mutate(thisstat_label = case_when(this_status == "PA" ~ "Patient"
#                                     , this_status == "PE" ~ "Staff"
#                                     , this_status == "V" ~ "Visitor"
#                                     , this_status == "ALL" ~ "All")) %>% 
#   mutate(thatstat_label = case_when(that_status == "PA" ~ "Patient"
#                                     , that_status == "PE" ~ "Staff"
#                                     , that_status == "V" ~ "Visitor"
#                                     , that_status == "ALL" ~ "All")) %>% 
#   mutate(thisthat_label = paste0(thisstat_label, ">", thatstat_label)) %>% 
#   ggplot(aes(y = newID)) + 
#   # geom_point(aes(x = eta_median), colour = "darkgrey") +
#   # geom_errorbarh(aes(xmin = eta_lo, xmax = eta_hi), height=error_bar_width, colour = "darkgrey") +
#   geom_point(aes(x = NL, colour = "Non-linear")) + 
#   geom_errorbarh(aes(xmin = NL_lo, xmax = NL_hi, colour = "Non-linear"), height=error_bar_width) +
#   geom_point(aes(x = NL, colour = "Non-linear")) + 
#   geom_point(aes(x = FD, colour = "Frequency dependent")) + 
#   geom_point(aes(x = DD, colour = "Linear density dependent")) + 
#   facet_nested(.~thisthat_label + Period) + 
#   theme_bw() + 
#   labs(y = "", x = expression(phi), colour = "Best fit model") + 
#   scale_y_discrete(limits=rev) + 
#   scale_colour_manual(values = c(`Frequency dependent` = "blue", `Linear density dependent` = "green", `Non-linear` = "black")) + 
#   geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
#   geom_vline(aes(xintercept = 1), linetype = "dashed")
# 
# ggsave(paste0("Graphs/eta_daynight_thisthatchosen.png"), width = 30, height = 12, units = "cm")
# 
# 
# ward_nonlinear_plottable %>% 
#   # filter(this_status %in% c("PA", "PE"), that_status %in% c("PA", "PE")) %>% 
#   filter((this_status %in% c("PA", "PE") & that_status %in% c("PA", "PE")) | (this_status == "ALL" & that_status == "ALL")) %>% 
#   filter(Period != "All") %>%
#   # left_join(chosen_model_thisthat) %>% 
#   left_join(rbind(chosen_model %>% mutate(this_status = "ALL", that_status = "ALL")
#         , chosen_model_thisthat)) %>% 
#   mutate(newID = factor(newID, levels = newID_ref$newID)) %>% 
#   mutate(FD = ifelse(ChosenModel == "Constant" | ChosenModel == "Constant diurnal", 0, NA_real_)
#          , DD = ifelse(ChosenModel == "Linear" | ChosenModel == "Linear diurnal", 1, NA_real_)
#          , NL = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", phi_median, NA_real_)
#          , NL_lo = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", phi_lo, NA_real_)
#          , NL_hi = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", phi_hi, NA_real_)) %>% 
#   mutate(thisstat_label = case_when(this_status == "PA" ~ "Patient"
#                                     , this_status == "PE" ~ "Staff"
#                                     , this_status == "V" ~ "Visitor"
#                                     , this_status == "ALL" ~ "All")) %>% 
#   mutate(thatstat_label = case_when(that_status == "PA" ~ "Patient"
#                                     , that_status == "PE" ~ "Staff"
#                                     , that_status == "V" ~ "Visitor"
#                                     , that_status == "ALL" ~ "All")) %>% 
#   mutate(thisthat_label = paste0(thisstat_label, ">", thatstat_label)) %>% 
#   ggplot(aes(y = newID)) + 
#   # geom_point(aes(x = eta_median), colour = "darkgrey") +
#   # geom_errorbarh(aes(xmin = eta_lo, xmax = eta_hi), height=error_bar_width, colour = "darkgrey") +
#   geom_point(aes(x = NL, colour = "Non-linear")) + 
#   geom_errorbarh(aes(xmin = NL_lo, xmax = NL_hi, colour = "Non-linear"), height=error_bar_width) +
#   geom_point(aes(x = NL, colour = "Non-linear")) + 
#   geom_point(aes(x = FD, colour = "Frequency dependent")) + 
#   geom_point(aes(x = DD, colour = "Linear density dependent")) + 
#   facet_nested(.~thisthat_label + Period) + 
#   theme_bw() + 
#   labs(y = "", x = expression(phi), colour = "Best fit model") + 
#   scale_y_discrete(limits=rev) + 
#   scale_colour_manual(values = c(`Frequency dependent` = "blue", `Linear density dependent` = "green", `Non-linear` = "black")) + 
#   geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
#   geom_vline(aes(xintercept = 1), linetype = "dashed")
# 
# ggsave(paste0("Graphs/eta_daynight_allthisthatchosen.png"), width = 35, height = 12, units = "cm")



### changes in peak intensity, what effect would they have?

plot_intensity_bar <- rates_N_allModels %>% 
  filter(this_status == "ALL", that_status == "ALL") %>% 
  group_by(newID, Period) %>% 
  left_join(chosen_model) %>% 
  mutate(chosen_eta = case_when(ChosenModel == "Constant" | ChosenModel == "Constant diurnal" ~ 0
                                , ChosenModel == "Linear" | ChosenModel == "Linear diurnal" ~ 1
                                , ChosenModel == "Non-linear"  ~ eta_median
                                , ChosenModel == "Non-linear diurnal" ~ eta_period)) %>% 
  mutate(w_chosenmodel = case_when(ChosenModel == "Constant" ~ w_constant
                                      , ChosenModel == "Linear" ~ w_linear
                                      , ChosenModel == "Non-linear" ~ w_nonlin
                                      , ChosenModel == "Constant diurnal" ~ w_cdiurnal
                                      , ChosenModel == "Linear diurnal" ~ w_ldiurnal
                                      , ChosenModel == "Non-linear diurnal" ~ w_ndiurnal)) %>% 
  # WHAT DOES THE W ACTUALLY MEAN? SCALE IN TERMS OF INTENSITY
  summarise(Nthat_peak = max(Nthat), w_chosenmodel = w_chosenmodel[1], chosen_eta = chosen_eta[1]) %>% 
  mutate(Nthat_110pct = Nthat_peak*1.1
            , Nthat_90pct = Nthat_peak*0.9
            , intensity_peak = w_chosenmodel * Nthat_peak^chosen_eta
            , intensity_110pct = w_chosenmodel * Nthat_110pct^chosen_eta
            , intensity_90pct = w_chosenmodel * Nthat_90pct^chosen_eta
            , change_110pct = intensity_110pct/intensity_peak - 1
            , change_90pct = intensity_90pct/intensity_peak - 1) %>% 
  View
  select(newID, Period, Nthat_peak, Nthat_90pct, Nthat_110pct, intensity_peak, intensity_90pct, change_90pct, intensity_110pct, change_110pct) %>% 
  pivot_longer(-c("newID", "Period"), names_sep = "_", names_to = c(".value", "Scenario")) %>% 
  mutate(Scenario = factor(Scenario, levels = c("90pct", "peak", "110pct"))
         , Scenario = fct_recode(Scenario, `90%\nof peak` = "90pct", `110%\nof peak` = "110pct")
         , change_label = ifelse(Scenario == "peak", "", paste0(round(change*100), "%"))) %>% 
  ggplot(aes(x = Scenario, y = intensity*60)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = change_label), size = 2, nudge_y = 1, colour = "red", size = 3) + 
  facet_nested_wrap(.~newID + Period, nrow = 2) + 
  theme_bw() + 
  labs(x = "", y = "Contact intensity (total cumulative contact minutes per hour)") + 
  theme(axis.text.x = element_text(size = 5))

ggsave(plot_intensity_bar, filename = paste0("Graphs/intensity.png")
       , width = 35, height = 20, units = "cm")



plot_intensity_bar <- rates_N_intensity %>% 
  left_join(newID_ref) %>% 
  group_by(newID, Period) %>% 
  filter(Nthat == max(Nthat)) %>% 
  transmute(newID, Nthat_peak = Nthat
            , w_form_ndiurnal, phi_period, rate_form_ndiurnal
            , w_break_ndiurnal, chi_period, rate_break_ndiurnal
            , intensity_peak = intensity_ndiurnal) %>% 
  mutate(Nthat_110pct = round(Nthat_peak*1.1)
         , rate_form_ndiurnal_110pct = w_form_ndiurnal*Nthat_110pct^phi_period
         , rate_break_ndiurnal_110pct = w_break_ndiurnal*Nthat_110pct^chi_period
         , intensity_110pct = rate_form_ndiurnal_110pct/rate_break_ndiurnal_110pct
         , change_110pct = intensity_110pct/intensity_peak
         , Nthat_90pct = round(Nthat_peak*0.9)
         , rate_form_ndiurnal_90pct = w_form_ndiurnal*Nthat_90pct^phi_period
         , rate_break_ndiurnal_90pct = w_break_ndiurnal*Nthat_90pct^chi_period
         , intensity_90pct = rate_form_ndiurnal_90pct/rate_break_ndiurnal_90pct
         , change_90pct = intensity_90pct/intensity_peak
  ) %>% 
  select(newID, Period, Nthat_peak, Nthat_90pct, Nthat_110pct, intensity_peak, intensity_90pct, change_90pct, intensity_110pct, change_110pct) %>% 
  pivot_longer(-c("newID", "Period"), names_sep = "_", names_to = c(".value", "Scenario")) %>% 
  mutate(Scenario = factor(Scenario, levels = c("90pct", "peak", "110pct"))
         , Scenario = fct_recode(Scenario, `90%\nof peak` = "90pct", `110%\nof peak` = "110pct")
         , change_label = ifelse(Scenario == "peak", "", paste0(round(change*100), "%"))) %>% 
  ggplot(aes(x = Scenario, y = intensity*60)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = change_label), nudge_y = 1, colour = "red", size = 3) + 
  facet_nested_wrap(.~newID + Period, nrow = 2) + 
  theme_bw() + 
  labs(x = "", y = "Contact intensity (total cumulative contact minutes per hour)") + 
  theme(axis.text.x = element_text(size = 5))

ggsave(plot_intensity_bar, filename = paste0("Graphs/intensity.png")
       , width = 35, height = 20, units = "cm")



# now plot the values of phi and chi by catHosp

ward_nonlinear_catHosp_plottable = ward_nonlinear_catHosp %>% 
  left_join(newID_ref) %>% 
  filter(!is.na(newID)) %>% 
  select(ward_id, newID, this_catHosp, that_catHosp, Period, lines, forms, breaks, phi_median, phi_lo, phi_hi, chi_median, chi_lo, chi_hi)

error_bar_width = 0.5



ward_nonlinear_catHosp_plottable %>% 
  filter(this_catHosp != "ALL", that_catHosp == "ALL") %>%
  # filter(that_catHosp == "ALL") %>%
  # filter(that_status == "ALL") %>% 
  mutate(thiscatHosp_label = factor(this_catHosp, levels = unique(this_catHosp))) %>% 
  filter(Period != "All") %>%
  ggplot(aes(y = newID)) + 
  geom_point(aes(x = phi_median)) + 
  geom_errorbarh(aes(xmin = phi_lo, xmax = phi_hi), height=error_bar_width) + 
  facet_nested_wrap(.~thiscatHosp_label + Period, nrow = 3) + 
  theme_bw() + 
  labs(y = "", x = expression(phi), colour = "") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(values = c(Patient = "red", HCW = "blue", Visitor = "black", All = "grey")) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed")

ggsave("Graphs/phi_daynight_thiscatHosp.png", width = 24, height = 24, units = "cm")


ward_nonlinear_catHosp_plottable %>% 
  filter(this_catHosp != "ALL", that_catHosp == "ALL") %>%
  # filter(that_catHosp == "ALL") %>%
  # filter(that_status == "ALL") %>% 
  mutate(thiscatHosp_label = factor(this_catHosp, levels = unique(this_catHosp))) %>% 
  filter(Period != "All") %>%
  ggplot(aes(y = newID)) + 
  geom_point(aes(x = chi_median)) + 
  geom_errorbarh(aes(xmin = chi_lo, xmax = chi_hi), height=error_bar_width) + 
  facet_nested_wrap(.~thiscatHosp_label + Period, nrow = 3) + 
  theme_bw() + 
  labs(y = "", x = expression(chi), colour = "") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(values = c(Patient = "red", HCW = "blue", Visitor = "black", All = "grey")) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed")

ggsave("Graphs/chi_daynight_thiscatHosp.png", width = 24, height = 24, units = "cm")




allwards_types <- bind_rows(
  ward_nonlinear_plottable %>% 
  filter(that_status == "ALL") %>%
  filter(ward_id == 0) %>% 
  # filter(that_catHosp == "ALL") %>%
  # filter(that_status == "ALL") %>% 
  filter(Period != "All") %>% 
  mutate(this_type = case_when(this_status == "PA" ~ "Patient"
                               , this_status == "PE" ~ "All HCW"
                               , this_status == "V" ~ "Visitor"
                               , this_status == "ALL" ~ "All"))
,
ward_nonlinear_catHosp_plottable %>% 
  filter(this_catHosp != "ALL", that_catHosp == "ALL") %>%
  filter(ward_id == 0) %>% 
  filter(Period != "All") %>% 
  # filter(that_catHosp == "ALL") %>%
  # filter(that_status == "ALL") %>% 
  mutate(this_type = this_catHosp) 
, 
) %>% 
  mutate(this_type = factor(this_type, levels = unique(this_type)))
  


phi_roles <- allwards_types %>% 
  mutate(process = "Formation") %>% 
  ggplot(aes(y = this_type)) + 
  geom_point(aes(x = phi_median)) + 
  geom_errorbarh(aes(xmin = phi_lo, xmax = phi_hi), height=error_bar_width) + 
  facet_nested(.~process + Period) + 
  theme_bw() + 
  labs(y = "", x = expression(chi)) + 
  scale_y_discrete(limits=rev) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed")

chi_roles <- allwards_types %>% 
  mutate(process = "Breaking") %>% 
  ggplot(aes(y = this_type)) + 
  geom_point(aes(x = chi_median)) + 
  geom_errorbarh(aes(xmin = chi_lo, xmax = chi_hi), height=error_bar_width) + 
  facet_nested(.~process + Period) + 
  theme_bw() + 
  labs(y = "", x = expression(chi)) + 
  scale_y_discrete(limits=rev) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  scale_x_continuous(breaks= pretty_breaks(n = 3))

phi_chi_roles <- arrangeGrob(phi_roles, chi_roles + theme(axis.text.y = element_blank()), ncol = 2, widths = c(5.5, 4))

ggsave(plot = phi_chi_roles, "Graphs/phi_chi_roles.png", width = 20, height = 12, units = "cm")


