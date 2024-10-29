library(tidyverse)
library(ggplot2)
# library(FME)
library(ggh4x)
library(modi)
library(scales)

source("2.1_load_analysis.R")


# read in tables produced by make_aicTables.R

chosen_model <- read_tsv(paste0("output/dAIC.tsv")) %>% 
  pivot_longer(-Ward, names_to = "ChosenModel", values_to = "dAIC") %>% 
  rename(newID = Ward) %>% 
  filter(dAIC == 0) %>% 
  select(-dAIC) %>% 
  group_by(newID) %>% #selects only the simplest model, in case of a tie
  slice(1)

chosen_model_thisthat <- read_tsv(paste0("output/dAIC_thisthat.tsv")) %>% 
  pivot_longer(-c("Ward", "this_status", "that_status"), names_to = "ChosenModel", values_to = "dAIC") %>% 
  rename(newID = Ward) %>% 
  filter(dAIC == 0) %>% 
  select(-dAIC) %>% 
  group_by(newID, this_status, that_status) %>% #selects only the simplest model, in case of a tie
  slice(1)


chosen_model_catHosp_allthisthat <- read_tsv(paste0("output/dAIC_catHosp_allthisthat.tsv")) %>% 
  pivot_longer(-c("Ward", "this_catHosp", "that_catHosp"), names_to = "ChosenModel", values_to = "dAIC") %>% 
  rename(newID = Ward) %>% 
  filter(dAIC == 0) %>% 
  select(-dAIC) %>% 
  group_by(newID, this_catHosp, that_catHosp) %>% #selects only the simplest model, in case of a tie
  slice(1)


### PLOT GRAPHS SHOWING THE BEST MODEL SLOPES ###

# slope_model_colours = c(`.data` = "black", 
#                   `Non-linear (negative)` = alpha("darkblue", 0.6), 
#                   `Non-linear (negative) - diurnal` = alpha("darkblue", 1), 
#                   `Frequency dependence` = alpha("blue", 0.6), 
#                   `Frequency dependence - diurnal` = alpha("blue", 1), 
#                   `Non-linear (partial)` = alpha("green", 0.6), 
#                   `Non-linear (partial) - diurnal` = alpha("green", 1), 
#                   `Linear density dependence` = alpha("gold", 0.6), 
#                   `Linear density dependence - diurnal` = alpha("gold", 1), 
#                   `Non-linear (more than linear)` = alpha("orange", 0.6), 
#                   `Non-linear (more than linear) - diurnal` = alpha("orange", 1))


slope_model_colours = c(`.data` = "black", 
                        `Non-linear (negative)` = "blue4", 
                        `Frequency dependence` = "blue", 
                        `Non-linear (partial)` = "green", 
                        `Linear density dependence` = "gold", 
                        `Non-linear (more than linear)` = "orange") 

greek = "phi"

toPlot_N_chosenModels_Slope <- ward_intervals_loglik %>% 
  filter((this_status == "ALL" & that_status == "ALL") | 
           ((this_status %in% c("PA", "PE")) & (that_status == "PA"))) %>% 
  ungroup %>% 
  select(newID, this_status, that_status, Period, starts_with("phi"), starts_with("w_"), starts_with("sd_")) %>% unique %>% 
  cross_join(tibble(Nthat = 1:100)) %>% 
  mutate(rate_constant = w_constant #calculate the rates for both non-diurnal and diurnal models
         , rate_linear = w_linear*Nthat
         , rate_nonlin = w_nonlin*Nthat^phi_nondiurnal
         , rate_cdiurnal = w_cdiurnal
         , rate_ldiurnal = w_ldiurnal*Nthat
         , rate_ndiurnal = w_ndiurnal*Nthat^phi_diurnal
         
         , rate_nonlin_lo = w_nonlin_lo*Nthat^phi_nondiurnal_lo
         , rate_nonlin_hi = w_nonlin_hi*Nthat^phi_nondiurnal_hi
         
         , rate_ndiurnal_lo = w_ndiurnal_lo*Nthat^phi_diurnal_lo
         , rate_ndiurnal_hi = w_ndiurnal_hi*Nthat^phi_diurnal_hi
         
         # , rate_constant_lo = w_constant - sd_constant #calculate the rates for both non-diurnal and diurnal models
         # , rate_linear_lo = (w_linear - sd_linear)*Nthat
         # , rate_nonlin_lo = (w_nonlin_lo - sd_nonlin_lo)*Nthat^phi_nondiurnal_lo
         # , rate_cdiurnal_lo = w_cdiurnal - sd_cdiurnal
         # , rate_ldiurnal_lo = (w_ldiurnal - sd_ldiurnal)*Nthat
         # , rate_ndiurnal_lo = (w_ndiurnal_lo - sd_ndiurnal_lo)*Nthat^phi_diurnal_lo
         # 
         # , rate_constant_hi = w_constant + sd_constant #calculate the rates for both non-diurnal and diurnal models
         # , rate_linear_hi = (w_linear + sd_linear)*Nthat
         # , rate_nonlin_hi = (w_nonlin_hi + sd_nonlin_hi)*Nthat^phi_nondiurnal_hi
         # , rate_cdiurnal_hi = w_cdiurnal + sd_cdiurnal
         # , rate_ldiurnal_hi = (w_ldiurnal + sd_ldiurnal)*Nthat
         # , rate_ndiurnal_hi = (w_ndiurnal_hi + sd_ndiurnal_hi)*Nthat^phi_diurnal_hi
  ) %>% 
  # left_join(chosen_model) %>%
  left_join(rbind(chosen_model %>% mutate(this_status = "ALL", that_status = "ALL")
                  , chosen_model_thisthat)) %>%
  mutate(newID = factor(newID, levels = newID_ref$newID)) %>% 
  mutate(rate_chosenmodel = case_when(ChosenModel == "Constant" ~ rate_constant
                                      , ChosenModel == "Linear" ~ rate_linear
                                      , ChosenModel == "Non-linear" ~ rate_nonlin
                                      , ChosenModel == "Constant diurnal" ~ rate_cdiurnal
                                      , ChosenModel == "Linear diurnal" ~ rate_ldiurnal
                                      , ChosenModel == "Non-linear diurnal" ~ rate_ndiurnal)) %>% 
  mutate(rate_chosenmodel_lo = case_when(ChosenModel == "Constant" ~ rate_constant
                                      , ChosenModel == "Linear" ~ rate_linear
                                      , ChosenModel == "Non-linear" ~ rate_nonlin_lo
                                      , ChosenModel == "Constant diurnal" ~ rate_cdiurnal
                                      , ChosenModel == "Linear diurnal" ~ rate_ldiurnal
                                      , ChosenModel == "Non-linear diurnal" ~ rate_ndiurnal_lo)) %>% 
  mutate(rate_chosenmodel_hi = case_when(ChosenModel == "Constant" ~ rate_constant
                                      , ChosenModel == "Linear" ~ rate_linear
                                      , ChosenModel == "Non-linear" ~ rate_nonlin_hi
                                      , ChosenModel == "Constant diurnal" ~ rate_cdiurnal
                                      , ChosenModel == "Linear diurnal" ~ rate_ldiurnal
                                      , ChosenModel == "Non-linear diurnal" ~ rate_ndiurnal_hi)) %>% 
  mutate(param_label = case_when(ChosenModel == "Non-linear" ~ paste0(greek, "==", round(phi_nondiurnal, 1))
                                 , ChosenModel == "Non-linear diurnal" ~ paste0(greek, "==", round(phi_diurnal, 1))
                                 , ChosenModel == "Constant" ~ paste0(greek, "==", 0)
                                 , ChosenModel == "Constant diurnal" ~ paste0(greek, "==", 0)
                                 , ChosenModel == "Linear" ~ paste0(greek, "==", 1)
                                 , ChosenModel == "Linear diurnal" ~ paste0(greek, "==", 1)
                                 , T ~ "")) %>% 
  mutate(chosen_model_slope = case_when(ChosenModel == "Constant" ~ "Frequency dependence"
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
  mutate(isDiurnal = ifelse(grepl("diurnal", ChosenModel), "Diurnal", "Non-diurnal"), 
         model_label = gsub(" - diurnal", "", chosen_model_slope)) %>% 
  mutate(model_label = factor(model_label, levels = names(slope_model_colours))) %>% 
  select(newID, this_status, that_status, Period, Nthat, ChosenModel, model_label, isDiurnal, param_label, rate_chosenmodel, rate_chosenmodel_lo, rate_chosenmodel_hi) %>% 
  filter(rate_chosenmodel <= 1)



plot_N_chosenModels_Slope = toPlot_N_chosenModels_Slope %>% 
  filter(this_status == "ALL", that_status == "ALL") %>%
  ggplot(aes(x = Nthat)) + 
  geom_line(aes(y = rate_chosenmodel, colour = model_label, linetype = isDiurnal), linewidth = 0.5) +
  geom_ribbon(aes(ymin = rate_chosenmodel_lo, ymax = rate_chosenmodel_hi, fill = model_label), alpha = 0.3, show.legend = F) +
  # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +
  # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) +
  # geom_line(aes(y = rate_chosenmodel, colour = model_label), linewidth = 1.0, show.legend = T) +
  # geom_line(aes(y = rate_chosenmodel, colour = model_label, linetype = diurnality), linewidth = 1.0, show.legend = T) +
  geom_point(data = ward_intervals_loglik %>% 
               filter(this_status == "ALL", that_status == "ALL"), aes(y = contactIntensity, colour = factor(".data", levels = names(slope_model_colours))), alpha = 0.3, size = 0.3) + 
  geom_text(data = toPlot_N_chosenModels_Slope %>% 
              filter(this_status == "ALL", that_status == "ALL") %>% 
              select(newID, Period, param_label) %>% unique, 
            aes(x = Inf, y = 3, label = param_label), parse = T, hjust = 1.0, vjust = 1) +
  scale_colour_manual(drop = F, values = slope_model_colours) +
  scale_fill_manual(values = slope_model_colours) +
  scale_linetype_manual(values = c(Diurnal = "solid", `Non-diurnal` = "dashed")) + 
  # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
  facet_nested_wrap(.~newID + Period
                    # , scales = "free_x"
                    ,  nrow = 2) +
  # guides(colour = F) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = c(0, 25, "", 75, ""), minor_breaks = NULL) + 
  theme_bw() + 
  labs(y = "Contact rate", x = "Persons present", colour = "Model", linetype = "") + 
  coord_cartesian(ylim = c(NA, 4))


ggsave(plot_N_chosenModels_Slope, filename = paste0("output/Fig2 - rates_over_N_chosenModels_intensitySlope.png")
       , width = 35, height = 20, units = "cm", dpi = 1000)
ggsave(plot_N_chosenModels_Slope, filename = paste0("output/Fig2 - rates_over_N_chosenModels_intensitySlope.pdf")
       , width = 35, height = 20, units = "cm", device = "pdf")  


for(this_stat in c("PA", "PE")){
  for(that_stat in c("PA")){
    
    plot_N_chosenModels_Slope_thisthat <- toPlot_N_chosenModels_Slope %>% 
      filter(this_status == this_stat, that_status == that_stat) %>%
      ggplot(aes(x = Nthat)) + 
      geom_line(aes(x = Nthat, y = rate_chosenmodel, colour = model_label, linetype = isDiurnal), linewidth = 0.5) +
      # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +
      # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) +
      # geom_line(aes(y = rate_chosenmodel, colour = model_label), linewidth = 1.0, show.legend = T) +
      # geom_line(aes(y = rate_chosenmodel, colour = model_label, linetype = diurnality), linewidth = 1.0, show.legend = T) +
      geom_point(data = ward_intervals_loglik %>% 
                   filter(this_status == this_stat, that_status == that_stat), 
                 aes(y = contactIntensity, colour = factor(".data", levels = names(slope_model_colours))), alpha = 0.3, size = 0.3) + 
      geom_text(data = toPlot_N_chosenModels_Slope %>% 
                  filter(this_status == this_stat, that_status == that_stat) %>% 
                  select(newID, Period, param_label) %>% unique, 
                aes(x = Inf, y = 3, label = param_label), parse = T, hjust = 1.0, vjust = 1) +
      scale_colour_manual(drop = F, values = slope_model_colours) +
      scale_linetype_manual(values = c(Diurnal = "solid", `Non-diurnal` = "dashed")) + 
      # geom_crossbar(aes(x = mean_Nthat, xmin = low_Nthat, xmax = high_Nthat, y = mean_break_rate, colour = "break"), size=1) + 
      facet_nested_wrap(.~newID + Period
                        # , scales = "free_x"
                        ,  nrow = 2) +
      # guides(colour = F) +
      scale_y_log10() +
      scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = c(0, 25, "", 75, ""), minor_breaks = NULL) + 
      theme_bw() + 
      labs(y = "Contact rate", x = "Persons present", colour = "Model", linetype = "")
    
    # 
    # plot_N_chosenModels_Slope_thisthat <- toPlot_N_chosenModels_Slope %>% 
    #   filter(this_status == this_stat, that_status == that_stat) %>%
    #   ggplot(aes(x = Nthat)) + 
    #   # geom_point(aes(x = mean_Nthat, y = obs_intensity, colour = ".data"), alpha = 1) +
    #   # geom_errorbar(aes(x = mean_Nthat, ymin = obs_intensity - sqrt(var_intensity), ymax = obs_intensity + sqrt(var_intensity), colour = ".data"), alpha = 1) +
    #   geom_line(aes(y = rate_chosenmodel, colour = model_label), linewidth = 1.0, show.legend = T) +
    #   geom_point(aes(y = contactIntensity, colour = factor(".data", levels = names(slope_model_colours))), alpha = 0.3, size = 0.3) + 
    #   geom_text(aes(x = Inf, y = 3, label = param_label), parse = T, hjust = 1.0, vjust = 1) +
    #   scale_colour_manual(drop = F, values = slope_model_colours) + 
    #   facet_nested_wrap(.~newID + Period
    #                     , scales = "free_x"
    #                     ,  nrow = 2) +
    #   # guides(colour = F) +
    #   scale_y_log10() +
    #   theme_bw() + 
    #   labs(y = "Contact rate", x = "Persons present", colour = "Model")
    

    ggsave(plot_N_chosenModels_Slope_thisthat, filename = paste0("output/SuppFig12 - rates_over_N_chosenModels_intensity_", this_stat, that_stat, ".png")
           , width = 35, height = 20, units = "cm", dpi = 1000)    
  }
}  



### NOW PLOT THE PHI VALUES THEMSELVES ###

model_cols = c(`Non-linear (negative)` = "blue4", 
               `Frequency dependence` = "blue", 
               `Non-linear (partial)` = "green", 
               `Linear density dependence` = "gold", 
               `Non-linear (more than linear)` = "orange")


error_bar_width = 0.5

grade_plot = ward_nonlinear_plottable %>% 
  pivot_wider(id_cols = c("newID", "this_status", "that_status"), names_from = "Period", values_from = c("phi_median", "phi_lo", "phi_hi")) %>%
  left_join(rbind(chosen_model %>% mutate(this_status = "ALL", that_status = "ALL")
                  , chosen_model_thisthat)) %>% 
  left_join(newID_ref %>% select(newID, newID_hospAbb)) %>% 
  select(-newID) %>% 
  rename(newID = newID_hospAbb) %>% 
  mutate(newID = factor(newID, levels = newID_ref$newID_hospAbb)) %>% 
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
  mutate(colour = factor(colour, levels = names(model_cols)))


grade_catHosp_plot = ward_nonlinear_catHosp_plottable %>% 
  pivot_wider(id_cols = c("newID", "this_catHosp", "that_catHosp"), names_from = "Period", values_from = c("phi_median", "phi_lo", "phi_hi")) %>%
  left_join(chosen_model_catHosp_allthisthat) %>% 
  left_join(newID_ref %>% select(newID, newID_hospAbb)) %>% 
  select(-newID) %>% 
  rename(newID = newID_hospAbb) %>% 
  mutate(newID = factor(newID, levels = newID_ref$newID_hospAbb)) %>% 
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
  mutate(colour = factor(colour, levels = names(model_cols)))




### PLOT THE PHI VALUES WITH THEIR CIs ###

grade_plot %>% 
  filter(!is.na(median)) %>% 
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
  scale_colour_manual(values = model_cols) +
  scale_shape_manual(values = c(`Diurnal model` = "circle", `Non-diurnal model` = "square")) + 
  geom_hline(aes(yintercept = 1.5), linetype = "solid")

ggsave(paste0("output/Fig1 - phi_daynight_ALL_toPAchosen.png"), width = 30, height = 12, units = "cm")
ggsave(paste0("output/Fig1 - phi_daynight_ALL_toPAchosen.pdf"), width = 30, height = 12, units = "cm", device = "pdf")



grade_catHosp_plot %>% 
  filter(!is.na(median)) %>% 
  
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
  geom_vline(aes(xintercept = 0), linetype = "dotted")+ 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  geom_point(aes(x = median, colour = colour, shape = isDiurnal)) + 
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = colour), height=error_bar_width) +
  facet_nested(.~thisthat_label + Period) + 
  theme_bw() + 
  labs(y = "", x = expression(phi), colour = "Best fit model", shape = "") + 
  scale_y_discrete(limits=rev) + 
  theme(axis.text.x = element_text(size = 7)) + 
  scale_colour_manual(values = model_cols) +
  scale_shape_manual(values = c(`Diurnal model` = "circle", `Non-diurnal model` = "square")) + 
  geom_hline(aes(yintercept = 1.5), linetype = "solid")

ggsave(paste0("output/Fig3 - phi_daynight_catHosp.png"), width = 25, height = 12, units = "cm")
ggsave(paste0("output/Fig3 - phi_daynight_catHosp.pdf"), width = 25, height = 12, units = "cm", device = "pdf")





### SIMULATION OF INCREASING POPULATION DENSITY ###

grade_plot %>% 
  filter(!is.na(median)) %>% 
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
  left_join(ward_intervals_loglik %>% 
              group_by(newID, that_status) %>% 
              summarise(mean_Nthat = sum(Nthat*interval_durmin)/sum(interval_durmin))) %>% 
  mutate(mean_Nthat_110 = mean_Nthat*1.1, 
         # increase_110 = mean_Nthat_110^median/mean_Nthat^median - 1, 
         increase_110 = 1.1^median - 1, 
         increase_110_lo = 1.1^lo - 1, 
         increase_110_hi = 1.1^hi - 1, 
  ) %>% 
  mutate(increase_110 = increase_110 + ifelse(!is.na(colour) & median == 0, 0.01, 0)) %>% 
  mutate(newID = fct_rev(newID)) %>% 
  ggplot() + 
  # geom_vline(aes(xintercept = 0), linetype = "solid")+ 
  geom_bar(aes(y = newID, x = increase_110, fill = colour), stat = "identity") + 
  geom_errorbarh(aes(y = newID, xmin = increase_110_lo, xmax = increase_110_hi), colour = "black", height = 0.5) +
  geom_errorbarh(aes(y = newID, xmin = increase_110_lo, xmax = increase_110_hi, colour = colour), alpha = 0.8, height = 0.5) +
  # coord_flip() + 
  facet_nested(.~thisthat_label + Period) + 
  theme_bw() + 
  labs(y = "", x = "Increase in contact rate with 10% increase in population density", fill = "Best fit model") + 
  scale_fill_manual(values = model_cols) +
  scale_colour_manual(values = model_cols) +
  geom_hline(aes(yintercept = 1.5), linetype = "solid") + 
  scale_x_continuous(labels = scales::percent) + 
  guides(colour = F)

ggsave(paste0("output/Fig4 - phi_daynight_ALL_toPAchosen_increase110.png"), width = 30, height = 12, units = "cm")
ggsave(paste0("output/Fig4 - phi_daynight_ALL_toPAchosen_increase110.pdf"), width = 30, height = 12, units = "cm", device = "pdf")






