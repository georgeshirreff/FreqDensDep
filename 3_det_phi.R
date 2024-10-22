library(deSolve)
library(tidyverse)
library(ggplot2)


det_phi <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    
    if(t < transfer_switch){
      transfer_in = transfer1
      transfer_out = transfer1  
    } else {
      transfer_in = transfer2
      transfer_out = transfer1  
      
    }
    
    
    beta = b/S0^phi
    N = S + I + R
    pI = I/N
    pI_index = 1/N
    FOI = (beta * pI * N^phi)
    FOI_index = (beta * pI_index * N^phi)
    
    dS = S0*transfer_in - S*FOI - S*transfer_out
    dI = S*FOI - gamma*I - I*transfer_out
    dR = gamma*I - R*transfer_out
  
    list(c(dS, dI, dR))
    
    })}
  
max_time = 100
pars = list(S0 = 1000
            , transfer1 = 0.1
            , transfer2 = 0.2
            , transfer_switch = 25
            , b = 2/5, phi = 0, gamma = 1/5)
inits <- c(S = pars$S0, I = 0, R = 0 )

phi_vec = seq(-0.5, 1.5, by = 0.5)
transfer2_vec = c(0.09, 0.1, 0.11)
# this_phi = phi_vec[1]
this_out = NULL
for(this_transfer2 in transfer2_vec){
  for(this_phi in phi_vec){
    this_out_piece <- ode(y = inits, times = 0:max_time, func = det_phi
                          , parms = pars %>% 
                            replace(c("phi", "transfer2"), c(this_phi, this_transfer2))) %>% 
      as.data.frame %>% as_tibble %>% 
      mutate(phi = this_phi, transfer2 = this_transfer2)
    
    if(is.null(this_out)){
      this_out <- this_out_piece
    } else {
      this_out <- rbind(this_out, this_out_piece)
    }
  }
}

par_labeller <- function(variable,value){
  if (variable=="name") {
    return(list(
      N="N",
      R0=expression(R[0])
    )[value])
  } else {
    return(value)
  }
}

model_colours = c(N = "red", 
                  `\"Non-linear (negative):\"~varphi==-0.5` = alpha("darkblue", 1), 
                  # `Non-linear (negative) - diurnal` = alpha("darkblue", 1), 
                  `\"Frequency dependence:\"~varphi==0` = alpha("blue", 1), 
                  # `Frequency dependence - diurnal` = alpha("blue", 1), 
                  `\"Non-linear (partial):\"~varphi==0.5` = alpha("green", 1), 
                  # `Non-linear (partial) - diurnal` = alpha("green", 1), 
                  `\"Linear density dependence:\"~varphi==1` = alpha("gold", 1), 
                  # `Linear density dependence - diurnal` = alpha("gold", 1), 
                  `\"Non-linear (more than linear):\"~varphi==1.5` = alpha("orange", 1) 
                  # `Non-linear (more than linear) - diurnal` = alpha("orange", 1)
                  )


this_out %>% 
  mutate(N = S + I + R
         , beta = pars$b/pars$S0^phi
         , R0 = beta*N^phi/pars$gamma
         , phi = factor(phi, levels = phi_vec)
  ) %>% 
  pivot_longer(-c("time", "phi", "transfer2")) %>% 
  mutate(demog_label = c("Population decreasing", "Population constant", "Population increasing")[match(transfer2, transfer2_vec)]) %>% 
  filter(demog_label != "Population constant") %>% 
  filter(name %in% c("N", "R0")) %>% 
  mutate(label = case_when(name == "N" ~ "N", 
                           phi == -0.5 ~ "\"Non-linear (negative):\"~varphi==-0.5", 
                           phi == 0 ~ "\"Frequency dependence:\"~varphi==0", 
                           phi == 0.5 ~ "\"Non-linear (partial):\"~varphi==0.5", 
                           phi == 1 ~ "\"Linear density dependence:\"~varphi==1", 
                           phi == 1.5 ~ "\"Non-linear (more than linear):\"~varphi==1.5", 
  )) %>% 
  mutate(label = factor(label, levels = names(model_colours))) %>% 
  # mutate(label = ifelse(name == "N", "N", as.character(phi)) 
  #        %>% factor(levels = c(phi_vec, "N"))
  #        ) %>% 
  # mutate(name_label = fct_relabel(name_label, levels = c("N", expression(R_0)))) %>% 
  # mutate(name_label = ifelse(name == "N", expression(N), expression(R_0))) %>% 
  ggplot(aes(x = time, y = value, colour = label)) + 
  geom_line() + 
  facet_grid(name~demog_label, scales = "free_y", labeller = par_labeller, switch = "y") + 
  theme_bw() + 
  # expand_limits(y = 0) + 
  # labs(y = "", colour = expression(italic(varphi))) + 
  labs(y = "", colour = "") + 
  # scale_colour_manual(values = c(`-0.5`="purple", `0`="deepskyblue4", `0.5`="darkturquoise", `1`="green3", `1.5`="gold", N="red"))+
  scale_colour_manual(values = model_colours, labels = lapply(names(model_colours), \(x) parse(text = x))) +
  theme(strip.text.y.left = element_text(angle = 0), strip.placement = "outside", strip.background = element_blank())


ggsave(filename = "output/Box1_det_phi.png", units = "cm", width = 18, height = 6, dpi = 600)
