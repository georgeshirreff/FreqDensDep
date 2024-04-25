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


ward_rates_catHosp <- read_csv("FreqDens_ward_rates_catHosp.csv"
                       , locale = locale(tz = "CET"))
newID_COVIDstat <- read_csv("../R_objects/typeID_COVIDstat.csv")

newID_ref <- newID_COVIDstat %>% 
  select(ward_id, newID) %>% 
  arrange(newID) %>% 
  rbind(tibble(ward_id = 0, newID = "All wards")) %>% 
  mutate(newID = factor(newID, levels = newID))


# ward_nonlinear <- read_csv("output/ward_nonlinear_bayes.csv")


# likelihood_distribution = "norm"
likelihood_distribution = "exp"
# likelihood_distribution = "gamma"


# ward_nonlinear_cutoff0 <- read_csv(paste0("output/ward_nonlinear_bayesIntensity_", likelihood_distribution, ".csv")) %>%
#   filter(cutoff == 0) %>% 
#   # filter(this_status == "ALL", that_status == "ALL") %>% 
#   rename(eta_median = eta_bayes)
# 
# ward_nonlinear_cutoff1 <- read_csv(paste0("output/ward_nonlinear_bayesIntensity_", likelihood_distribution, ".csv")) %>%
#   filter(cutoff == 1) %>% 
#   # filter(this_status == "ALL", that_status == "ALL") %>% 
#   rename(eta_median = eta_bayes)
# 
# ward_nonlinear_cutoff2 <- read_csv(paste0("output/ward_nonlinear_bayesIntensity_", likelihood_distribution, ".csv")) %>%
#   filter(cutoff == 2) %>% 
#   # filter(this_status == "ALL", that_status == "ALL") %>% 
#   rename(eta_median = eta_bayes)
# 
# 
# # ward_nonlinear_cutoff1 <- read_csv("output/ward_nonlinear_bayes_main.csv")  %>%
# #   select(-any_of(c("phi", 'phi_opt_loglik', 'chi', 'chi_opt_loglik'))) %>% rename(phi_median = phi_bayes, chi_median = chi_bayes)
# # ward_nonlinear_cutoff2 <- read_csv("output/ward_nonlinear_bayes_cutoff.csv") %>%
# #   filter(cutoff == 2) %>% 
# #   select(-any_of(c("phi", 'phi_opt_loglik', 'chi', 'chi_opt_loglik'))) %>% rename(phi_median = phi_bayes, chi_median = chi_bayes)
# 
# 
# ward_nonlinear_main <- read_csv(paste0("output/ward_nonlinear_bayesIntensity_", likelihood_distribution, ".csv")) %>%
#   filter(cutoff == 1) %>% 
#   rename(eta_median = eta_bayes)
# 
# # ward_nonlinear_main <- read_csv("output/ward_nonlinear_bayes_main.csv") %>% 
# #   select(-any_of(c("phi", 'phi_opt_loglik', 'chi', 'chi_opt_loglik'))) %>% rename(phi_median = phi_bayes, chi_median = chi_bayes)

ward_nonlinear_catHosp <- read_csv("output/ward_nonlinear_bayesIntensity_exp_catHosp.csv") %>%
    rename(eta_median = eta_bayes)
  

### calculate AIC tables ###



# for(process in c("form", "break")){


      cutoff = 1
      print(c(cutoff))
        
    ward_rates_edit = ward_rates_catHosp %>% 
      group_by(ward_id) %>% 
      filter(as.numeric(difftime(time, min(time), units = "hours")) >= cutoff) %>% 
      ungroup
    
    ward_nonlinear_edit = ward_nonlinear_catHosp

    
    # total breakdown across different types of contact 
    
    # ward_rates_loglik_form <- ward_rates %>%
    # ward_rates_loglik_form_cutoffone <- ward_rates %>% group_by(ward_id) %>% filter(as.numeric(difftime(time, min(time), units = "hours")) > 1) %>% ungroup %>% 
    # ward_rates_loglik_form_cutofftwo <- ward_rates %>% group_by(ward_id) %>% filter(as.numeric(difftime(time, min(time), units = "hours")) > 2) %>% ungroup %>% 
    # ward_rates_loglik_break <- ward_rates %>%
    # ward_rates_loglik_break_cutoffone <- ward_rates %>% ward_rates %>% group_by(ward_id) %>% filter(as.numeric(difftime(time, min(time), units = "hours")) > 1) %>% ungroup %>% 
    # ward_rates_loglik_break_cutofftwo <- ward_rates %>% group_by(ward_id) %>% filter(as.numeric(difftime(time, min(time), units = "hours")) > 2) %>% ungroup %>% 

    
    ward_rates_loglik <- ward_rates_edit %>% 
      # filter(this_status == "ALL", that_status == "ALL") %>% 
      filter(Nthat > 0) %>% 
      mutate(Period = ifelse(hour(time) %in% 6:17, "Day", "Night")) %>% 
      {rbind(., mutate(., ward_id = 0))} %>% 
      left_join(ward_nonlinear_edit %>% filter(Period == "All") %>% select(-Period, -ends_with("loglik"))) %>%
      left_join(ward_nonlinear_edit %>% filter(Period != "All") %>% transmute(ward_id, Period, this_catHosp, that_catHosp, eta_period = eta_median)) %>%
 
      group_by(ward_id, this_catHosp, that_catHosp) %>% 
      filter(form_minutes > 0) %>% 
      mutate(var = var(contactIntensity)) %>% 
      group_by(ward_id, this_catHosp, that_catHosp) %>%
      mutate(w_constant = sum(break_minutes)/sum(form_minutes)
             , w_linear = sum(break_minutes)/sum(form_minutes*Nthat)
             , w_nonlin = sum(break_minutes)/sum(form_minutes*Nthat^eta_median)) %>% 
      group_by(ward_id, this_catHosp, that_catHosp, Period) %>% #group by day/night to calculate the day and night parameters
      mutate(w_cdiurnal = sum(break_minutes)/sum(form_minutes)
             , w_ldiurnal = sum(break_minutes)/sum(form_minutes*Nthat)
             , w_ndiurnal = sum(break_minutes)/sum(form_minutes*Nthat^eta_period)) %>% 
      mutate(rate_constant = w_constant #calculate the rates for both non-diurnal and diurnal models
             , rate_linear = w_linear*Nthat
             , rate_nonlin = w_nonlin*Nthat^eta_median
             , rate_cdiurnal = w_cdiurnal
             , rate_ldiurnal = w_ldiurnal*Nthat
             , rate_ndiurnal = w_ndiurnal*Nthat^eta_period
      ) %>%
      mutate(across(starts_with("rate")
                    , .fns = list(loglik = function(r) dnorm(x = contactIntensity, rate = 1/r, log = T))
                    # , .fns = list(loglik = function(r) dnorm(x = contactIntensity, mean = r, sd = sqrt(var), log = T))
                    , .names = "{.fn}_{.col}")) %>% 
      rename_all(function(x) gsub("loglik_rate_", "loglik_", x)) %>% 
      group_by(ward_id, this_catHosp, that_catHosp) %>% 
      summarise(across(starts_with("loglik"), sum)
                # , eta = eta_median[1]
                ) %>% 
      pivot_longer(-c("ward_id", "this_catHosp", "that_catHosp"), names_sep = "_", names_to = c(".value", "Model")) %>% 
      mutate(K = case_when(Model == "constant" ~ 1
                           , Model == "linear" ~ 1
                           , Model == "nonlin" ~ 2
                           , Model == "cdiurnal" ~ 2
                           , Model == "ldiurnal" ~ 2
                           , Model == "ndiurnal" ~ 4
      ))
    
    # # in some instances the non-linear model has poorer likelihood than the constant model
    # # but this is probably due to rounding errors when phi is near zero
    # ward_rates_loglik_break %>% 
    #   group_by(ward_id, this_status, that_status) %>% 
    #   mutate(report = case_when(loglik[Model == "nonlin"] < loglik[Model == "constant"] ~ "Nonlin worse than const"
    #                             , loglik[Model == "nonlin"] < loglik[Model == "linear"] ~ "Nonlin worse than linear"
    #                             , loglik[Model == "ndiurnal"] < loglik[Model == "cdiurnal"] ~ "Nonlin diurnal worse than constant diurnal"
    #                             , loglik[Model == "ndiurnal"] < loglik[Model == "ldiurnal"] ~ "Nonlin diurnal worse than linear diurnal"
    #                             , min(loglik[Model %in% c("diurnal", "ldiurnal", "ndiurnal")]) < min(loglik[Model %in% c("constant", "linear", "nonlin")]) ~ "Diurnal worse than non-diurnal"
    #                             , T ~ "Fine"
    #   )) %>%
    #   filter(report != "Fine") %>% 
    #   filter(report != "Nonlin diurnal worse than constant diurnal") %>% 
    #   filter(Model %in% c("nonlin", "constant")) %>% 
    #   print(n = 1000)
    
    
    ward_aic = ward_rates_loglik %>% 
      filter(that_catHosp %in% c("ALL", "patient")) %>% 
      # ward_aic_form = ward_rates_loglik_form %>%
      # ward_aic_break = ward_rates_loglik_break %>%
      # ward_aic_form_cutoff = ward_rates_loglik_form_cutoff %>%
      # ward_aic_break_cutoff = ward_rates_loglik_break_cutoff %>%
      mutate(AIC = 2*K - 2*loglik) %>% 
      # filter(this_status == "ALL", that_status == "ALL") %>%
      # filter(this_status != "ALL", that_status == "ALL") %>%
      # filter(this_status == "ALL", that_status != "ALL") %>%
      # filter(this_status != "ALL", that_status != "ALL") %>%
      group_by(ward_id, this_catHosp, that_catHosp) %>% 
      mutate(dAIC = AIC - min(AIC)) %>% 
      pivot_wider(id_cols = c("ward_id", "this_catHosp", "that_catHosp"), names_from = "Model", values_from = "dAIC") %>% 
      left_join(newID_ref) %>% 
      ungroup
    
    tab_catHosp <- ward_aic %>% 
      filter(ward_id != 15) %>% 
      filter(this_catHosp != "ALL") %>% 
      arrange(newID) %>% 
      transmute(Ward = newID, this_catHosp, that_catHosp, Constant = constant, Linear = linear, `Non-linear` = nonlin, `Constant diurnal` = cdiurnal, `Linear diurnal` = ldiurnal, `Non-linear diurnal` = ndiurnal) %>% 
      mutate_if(is.numeric, function(x) paste0(ifelse(x == 0, "0", format(round(x, 2), trim = T, nsmall=2))))
    
 
    
      
    tab_catHosp %>%   
      write_tsv(paste0("AICtables/dAIC_intensity_catHosp_cutoff", cutoff, ".tsv"))
    
    
  






ward_nonlinear_main

cutoff = 1
ward_rates_main = ward_rates %>% 
  group_by(ward_id) %>% 
  filter(as.numeric(difftime(time, min(time), units = "hours")) >= cutoff) %>% 
  ungroup

pair_cols = RColorBrewer::brewer.pal(6, name="Paired")

# process = "form"
rates_N_allModels <- ward_rates_catHosp %>% 
  mutate(Period = ifelse(hour(time) %in% 6:17, "Day", "Night")) %>% 
  {rbind(., mutate(., ward_id = 0))} %>% 
  left_join(newID_ref) %>% 
  filter(!is.na(newID)) %>% 
  filter(Nthat > 0) %>% 
  left_join(ward_nonlinear_edit %>% filter(Period == "All") %>% select(-Period, -ends_with("loglik"))) %>%
  left_join(ward_nonlinear_edit %>% filter(Period != "All") %>% transmute(ward_id, Period, this_catHosp, that_catHosp, eta_period = eta_median)) %>%
  mutate(Nthat_cat = cut(Nthat, seq(0, by = 10, max(Nthat) + 10))) %>% 
  # filter(this_status == "ALL", that_status == "ALL") %>% 
  filter(form_minutes > 0) %>% 
  # rename_func_bool(process) %>% 
  # rename(event_rate = form_rate, number_events = number_forms, event_minutes = form_minutes, nonlin = phi, nonlin_period = phi_period) %>% mutate(nonlin_expr = paste0(Period, "\n(phi=", round(nonlin_period, 1), ")")) %>%
  # rename(event_rate = break_rate, number_events = number_breaks, event_minutes = break_minutes, nonlin = chi, nonlin_period = chi_period) %>% mutate(nonlin_expr = paste0(Period, "\n(chi=", round(nonlin_period, 1), ")")) %>%
  group_by(ward_id, Period, Nthat_cat) %>% ### PREPARE THE DATA ###
  mutate(mean_Nthat = mean(Nthat)
         , low_Nthat = min(Nthat)
         , high_Nthat = max(Nthat)
         , obs_intensity = sum(contactIntensity*form_minutes)/sum(form_minutes)
         # , lo_event_rate = ifelse(sum(event_minutes) == 0, NaN, weighted.quantile(event_rate, event_minutes, prob = 0.25))
         # , hi_event_rate = ifelse(sum(event_minutes) == 0, NaN, weighted.quantile(event_rate, event_minutes, prob = 0.75))
         , var_intensity = ifelse(sum(contactIntensity) == 0, NaN, weighted.var(contactIntensity, form_minutes, na.rm = T))
         # , mean_break_rate = sum(break_rate*break_minutes, na.rm = T)/sum(break_minutes)
  ) %>% 
  
  group_by(ward_id, this_catHosp, that_catHosp) %>%
  mutate(w_constant = sum(break_minutes)/sum(form_minutes)
         , w_linear = sum(break_minutes)/sum(form_minutes*Nthat)
         , w_nonlin = sum(break_minutes)/sum(form_minutes*Nthat^eta_median)) %>% 
  group_by(ward_id, this_catHosp, that_catHosp, Period) %>% #group by day/night to calculate the day and night parameters
  mutate(w_cdiurnal = sum(break_minutes)/sum(form_minutes)
         , w_ldiurnal = sum(break_minutes)/sum(form_minutes*Nthat)
         , w_ndiurnal = sum(break_minutes)/sum(form_minutes*Nthat^eta_period)) %>% 
  mutate(rate_constant = w_constant #calculate the rates for both non-diurnal and diurnal models
         , rate_linear = w_linear*Nthat
         , rate_nonlin = w_nonlin*Nthat^eta_median
         , rate_cdiurnal = w_cdiurnal
         , rate_ldiurnal = w_ldiurnal*Nthat
         , rate_ndiurnal = w_ndiurnal*Nthat^eta_period
  )
  
  
  
  plot_N_allModels <- rates_N_allModels %>% 
    # filter(plot_row == r) %>% 
    mutate(newID = factor(newID, levels = levels(rates_N_allModels$newID))) %>% 
    filter(Nthat > 0) %>% 
    filter(form_minutes > 0) %>% 
    filter(this_catHosp == "nurse", that_catHosp == "ALL") %>% 
    ggplot(aes(x = Nthat)) + 
    geom_point(aes(y = contactIntensity, colour = ".data"), alpha = 0.2, size = 0.2) + 
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_event_rate, colour = "data"), size=1) + 
    # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_event_rate, colour = "data", size = weight_Nthat)) + 
    geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +  
    geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) + 
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
  
 
  
  # plot_N <- rates_N_allModels %>% 
  #   # filter(plot_row == r) %>% 
  #   filter(Nthat > 0) %>% 
  #   mutate(nonlin_expr = paste0(Period, "\n(", ifelse(process == "form", "phi", "chi"), "=", round(nonlin_period, 1), ")")) %>% 
  #   ggplot(aes(x = Nthat)) + 
  #   geom_point(aes(y = event_rate, colour = ".data"), alpha = 0.2, size = 0.2) + 
  #   geom_point(aes(x = mean_Nthat, y = mean_event_rate, colour = ".data"), alpha = 1) +  
  #   geom_errorbar(aes(x = mean_Nthat, ymin = mean_event_rate - sqrt(var_event_rate), ymax = mean_event_rate + sqrt(var_event_rate), colour = ".data"), alpha = 1) + 
  #   # geom_errorbar(aes(x = mean_Nthat, ymin = lo_event_rate, ymax = hi_event_rate, colour = "data"), alpha = 1) + 
  #   # geom_point(aes(y = break_rate, colour = "break"), alpha = 0.1) + 
  #   # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_event_rate, colour = "data"), size=1) + 
  #   geom_line(aes(y = rate_ndiurnal, colour = "Non-linear\ndiurnal")) + 
  #   scale_colour_manual(values = c(`.data` = "blue"
  #                                  , `Constant` = pair_cols[1]
  #                                  , `Constant\ndiurnal` = pair_cols[2]
  #                                  , `Linear` = pair_cols[3]
  #                                  , `Linear\ndiurnal` = pair_cols[4]
  #                                  , `Non-linear` = pair_cols[5]
  #                                  , `Non-linear\ndiurnal` = pair_cols[6]
  #   ))+
  #   # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
  #   facet_nested_wrap(.~newID + nonlin_expr, scales = "free_x", nrow = 2) +
  #   # guides(colour = F) +
  #   scale_y_log10() +
  #   theme_bw() + 
  #   labs(y = "Rate per person minute", x = "Persons present", colour = "Contact process") + 
  #   coord_cartesian(ylim = c(1e-3, NA))
  greek = "phi"
  
  
  
  chosen_model_catHosp <- read_tsv(paste0("AICtables/dAIC_intensity_catHosp_cutoff1.tsv")) %>% 
    pivot_longer(-c("Ward", "this_catHosp", "that_catHosp"), names_to = "ChosenModel", values_to = "dAIC") %>% 
    rename(newID = Ward) %>% 
    filter(dAIC == 0) %>% 
    select(-dAIC)
  
  
  # 
  # plot_N_chosenModels <- rates_N_allModels %>% 
  #   left_join(chosen_model) %>%
  #   mutate(newID = factor(newID, levels = levels(rates_N_allModels$newID))) %>% 
  #   # filter(plot_row == r) %>% 
  #   filter(Nthat > 0) %>% 
  #   filter(form_minutes > 0) %>% 
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
  #   ggplot(aes(x = Nthat)) + 
  #   geom_point(aes(y = contactIntensity, colour = ".data"), alpha = 0.3, size = 0.3) + 
  #   # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +
  #   # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) +
  #   geom_line(aes(x = Nthat, y = rate_chosenmodel, colour = ChosenModel)) +
  #   geom_text(aes(x = Inf, y = 3, label = param_label), parse = T, hjust = 1.0, vjust = 1) +
  #   scale_colour_manual(values = c(`.data` = "blue"
  #                                  , `Constant` = pair_cols[1]
  #                                  , `Constant diurnal` = pair_cols[2]
  #                                  , `Linear` = pair_cols[3]
  #                                  , `Linear diurnal` = pair_cols[4]
  #                                  , `Non-linear` = pair_cols[5]
  #                                  , `Non-linear diurnal` = pair_cols[6]
  #   ))+
  #   # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
  #   facet_nested_wrap(.~newID + Period
  #                     , scales = "free_x"
  #                     ,  nrow = 2) +
  #   # guides(colour = F) +
  #   scale_y_log10() +
  #   theme_bw() + 
  #   labs(y = "Contact minutes per person minute", x = "Persons present", colour = "Model")
  
  
  plot_N_chosenModels <- rates_N_allModels %>% 
    left_join(chosen_model) %>%
    mutate(newID = factor(newID, levels = levels(rates_N_allModels$newID))) %>% 
    # filter(plot_row == r) %>% 
    filter(Nthat > 0) %>% 
    filter(form_minutes > 0) %>% 
    filter(this_catHosp == "nurse", that_catHosp == "ALL") %>% 
    mutate(rate_chosenmodel = case_when(ChosenModel == "Constant" ~ rate_constant
                                        , ChosenModel == "Linear" ~ rate_linear
                                        , ChosenModel == "Non-linear" ~ rate_nonlin
                                        , ChosenModel == "Constant diurnal" ~ rate_cdiurnal
                                        , ChosenModel == "Non-linear" ~ rate_ldiurnal
                                        , ChosenModel == "Non-linear diurnal" ~ rate_ndiurnal
    )) %>% 
    mutate(param_label = case_when(ChosenModel == "Non-linear" ~ paste0(greek, "==", round(eta_median, 1))
                                   , ChosenModel == "Non-linear diurnal" ~ paste0(greek, "==", round(eta_period, 1))
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
                                                        , "Non-linear - diurnal"))) %>% 
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
  
  
  



### NOW PLOT THE ETA VALUES THEMSELVES ###

ward_nonlinear_main
ward_nonlinear_plottable = ward_nonlinear_catHosp %>% 
  left_join(newID_ref) %>% 
  filter(!is.na(newID)) %>% 
  select(ward_id, newID, this_catHosp, that_catHosp, Period, eta_median, eta_lo, eta_hi)

error_bar_width = 0.5


ward_nonlinear_plottable %>% 
  # filter(this_status == "ALL", that_status == "ALL") %>% 
  # filter(newID == "All wards") %>% 
  filter(this_catHosp != "ALL") %>% 
  # filter(that_catHosp == "ALL") %>% 
  filter(this_catHosp != "investigation") %>%
  filter(this_catHosp %in% c("nurse", "aux nurse", "student nurse", "physician", "patient")) %>% 
  
  filter(Period != "All") %>%
  left_join(chosen_model_catHosp) %>% 
  mutate(newID = factor(newID, levels = newID_ref$newID)) %>% 
  mutate(FD = ifelse(ChosenModel == "Constant" | ChosenModel == "Constant diurnal", 0, NA_real_)
         , DD = ifelse(ChosenModel == "Linear" | ChosenModel == "Linear diurnal", 1, NA_real_)
         , NL = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", eta_median, NA_real_)
         , NL_lo = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", eta_lo, NA_real_)
         , NL_hi = ifelse(ChosenModel == "Non-linear" | ChosenModel == "Non-linear diurnal", eta_hi, NA_real_)) %>% 
  ggplot(aes(y = newID)) + 
  geom_point(aes(x = eta_median), colour = "darkgrey") + 
  geom_errorbarh(aes(xmin = eta_lo, xmax = eta_hi), height=error_bar_width, colour = "darkgrey") +
  geom_point(aes(x = NL, colour = "Non-linear")) + 
  geom_errorbarh(aes(xmin = NL_lo, xmax = NL_hi, colour = "Non-linear"), height=error_bar_width) +
  geom_point(aes(x = NL, colour = "Non-linear")) + 
  geom_point(aes(x = FD, colour = "Frequency dependent")) + 
  geom_point(aes(x = DD, colour = "Linear density dependent")) + 
  facet_grid(this_catHosp~that_catHosp + Period) + 
  theme_bw() + 
  labs(y = "", x = expression(phi), colour = "Best fit model") + 
  scale_y_discrete(limits=rev) + 
  scale_colour_manual(values = c(`Frequency dependent` = "blue", `Linear density dependent` = "green", `Non-linear` = "black")) + 
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed")

ggsave(paste0("Graphs/eta_daynight_catHosp.png"), width = 30, height = 30, units = "cm", dpi = 1000)

