
freq_phi_nonlinear <- function(pars, this_tib){
  phi = pars[1]
  
  loglik <- this_tib %>%
    mutate(w = sum(break_minutes)/sum(form_minutes*Nthat^phi)
           , modelIntensity  = w * Nthat^phi
           , loglik_phi = dexp(x = contactIntensity, rate = 1/modelIntensity, log = T)) %>%
    pull(loglik_phi) %>% sum
  
  if(is.na(loglik)) stop("loglik is invalid")
  
  -2*loglik
}





freq_dens_prior <- function(pars, prior_mean, prior_sd){
  phi = pars[1]
  
  loglik = dnorm(x = phi, mean = prior_mean, sd = prior_sd, log = T)
  
  -2*loglik
}


n_iterations = 2000
burnin = 200
max_iterations = 100000
effSize_threshold = 1000
prior_mean = 0.5
prior_sd = 1
param_jump = 0.2
param_lower = -3
param_upper = 3

if(total_formmin > 0 & sum(this_ward_intervals$break_minutes) > 0){
  effSize = 0
  current_phi = runif(1, 0, 1)
  phi_post = numeric(0)
  lik_post = numeric(0)
  while(effSize < effSize_threshold & length(phi_post) < max_iterations){
    mod_phi_post = modMCMC(f = freq_phi_nonlinear
                           , this_tib = this_ward_intervals
                           , prior = function(par) freq_dens_prior(par, prior_mean = prior_mean, prior_sd = prior_sd)
                           , p = c(phi = current_phi), jump = param_jump, lower = param_lower, upper = param_upper, niter = n_iterations)
    
    phi_post = c(phi_post, mod_phi_post$pars %>% unlist)
    lik_post = c(lik_post, mod_phi_post$prior + mod_phi_post$SS)
    
    if(effSize == 0){ #remove burnin
      phi_post = phi_post[-(1:burnin)]
      lik_post = lik_post[-(1:burnin)]
    }
    effSize = effectiveSize(phi_post)
    current_phi = last(phi_post)
    print(paste("effSize:", round(effSize)))
  }            
  
  
  phi_median = median(phi_post)
  phi_best = (phi_post[which.min(lik_post)]) %>% median #in case of multiple max values
  
  phi_lo = quantile(phi_post, 0.025)
  phi_hi = quantile(phi_post, 0.975)
  
  D = -2*lik_post
  DIC = 0.5*var(D) + mean(D)

  # random sample
  phi_vec100 = sample(phi_post, size = 100, replace = F)
  # regular sample
  # phi_vec100 = phi_post[seq(1, length(phi_post), length.out = 100) %>% round]
  

} else {
  phi_median = NA_real_
  phi_best = NA_real_
  
  phi_lo = NA_real_
  phi_hi = NA_real_
  
  phi_vec100 = rep(NA_real_, 100)
  
  DIC = NA_real_
}


