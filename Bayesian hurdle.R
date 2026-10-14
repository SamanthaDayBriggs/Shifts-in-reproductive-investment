library(rstanarm); library(parallel); library(loo); library(ggplot2) ; library(ggeffects) ;
library(car) ; library(tidybayes)


  full_dat <- read.csv("~/fulldat.csv")
  dat <- full_dat %>%
    group_by(pot_id) %>%
    mutate(total_flowers = sum(flower, na.rm = TRUE)) %>%
    slice_max(doy, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  dat$flowered <- ifelse(dat$total_flowers > 0, 1, 0)
  dat_live <- subset(dat, scam > 0)
  dat_pos <- subset(dat_live, total_flowers > 0)
  
  dat_live$level <- factor(dat_live$level)
  dat_pos$level  <- factor(dat_pos$level)
  
  #BAYESIAN HURDLE MODEL
  
  ##PROBABILITY OF PRODUCING FLOWERS
  ## Base model
  p <- stan_glmer(flowered ~ scam + level + salinity + co2 + age + comp+ location + (1|genotype)+(1|frame) + (1|position),
                  data = dat_live,
                  family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  # Environmental interactions
  p1 <- stan_glmer(flowered ~ scam + level*salinity + co2 + age + comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_live,
                   family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p2 <- stan_glmer(flowered ~ scam + level*co2 + level*salinity + age + comp+ location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_live,
                   family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p3 <- stan_glmer(flowered ~ scam + level*salinity*co2 + co2 + age + comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_live,
                   family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p4 <- stan_glmer(flowered ~ scam + level*salinity*comp+ age + co2 + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_live,
                   family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  ## Environmental and biomass
  p5 <- stan_glmer(flowered ~ scam * level +salinity + age + co2 +comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_live,
                   family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p6 <- stan_glmer(flowered ~ scam * salinity * level + age + co2 +comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_live,
                   family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  #Evolution of plasticity
  p7 <- stan_glmer(flowered ~ scam * age + level+ salinity + co2+ comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_live,
                   family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p8 <- stan_glmer(flowered ~ scam * level*age + salinity + co2+ comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_live,
                   family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p9 <- stan_glmer(flowered ~ scam * level + salinity*age + co2+ comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_live,
                   family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p10 <- stan_glmer(flowered ~ scam * level + salinity + co2*age + comp + location + (1|genotype)+(1|frame) + (1|position),
                    data = dat_live,
                    family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p11 <- stan_glmer(flowered ~ scam * level + salinity + co2 + comp*age + location + (1|genotype)+(1|frame) + (1|position),
                    data = dat_live,
                    family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  #GxE
  p12 <- stan_glmer(flowered ~ scam + level*salinity + age *co2 + comp + (1+ comp|genotype)+location +(1|frame) + (1|position),
                    data = dat_live,
                    family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p13 <- stan_glmer(flowered ~ scam + level +age* co2 + comp + (1+ salinity |genotype)+ location +(1|frame) + (1|position),
                    data = dat_live,
                    family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p14 <- stan_glmer(flowered ~ scam + level*salinity + age + comp + (1+ co2 |genotype)+ location +(1|frame) + (1|position),
                    data = dat_live,
                    family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  p15 <- stan_glmer(flowered ~ scam + salinity + age* co2 + comp + (1+ level |genotype)+ location +(1|frame) + (1|position),
                    data = dat_live,
                    family = binomial, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  

  ##NUMBER OF FLOWERS PRODUCED
  # Base model
  n <- stan_glmer(total_flowers ~ scam + level + salinity + co2 + age + comp+ location + (1|genotype)+(1|frame) + (1|position),
                  data = dat_pos,
                  family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  # Environmental interactions
  n1 <- stan_glmer(total_flowers ~ scam + level*salinity+ co2 + age + comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_pos,
                   family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  n2 <- stan_glmer(total_flowers ~ scam + level*co2+ level*salinity + age + comp+ location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_pos,
                   family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  n3 <- stan_glmer(total_flowers ~ scam + level*salinity*co2 + co2 + age + comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_pos,
                   family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  n4 <- stan_glmer(total_flowers ~ scam + level*salinity*comp + age + co2 + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_pos,
                   family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  ## Inundation and biomass
  n5 <- stan_glmer(total_flowers ~ scam * level +salinity + age * co2 +comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_pos,
                   family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  
  n6 <- stan_glmer(total_flowers ~ scam * salinity *level + age * co2 +comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_pos,
                   family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  #Evolution of plasticity
  n7 <- stan_glmer(total_flowers ~ scam * age + level+ salinity + co2+ comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_pos,
                   family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  n8 <- stan_glmer(total_flowers ~ scam * level*age + salinity + co2+ comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_pos,
                   family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  n9 <- stan_glmer(total_flowers ~ scam * level + salinity*age + age*co2+ comp + location + (1|genotype)+(1|frame) + (1|position),
                   data = dat_pos,
                   family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  n10 <- stan_glmer(total_flowers ~ scam * level + salinity + co2*age + comp + location + (1|genotype)+(1|frame) + (1|position),
                    data = dat_pos,
                    family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  n11 <- stan_glmer(total_flowers ~ scam * level + salinity + co2 + comp*age + location + (1|genotype)+(1|frame) + (1|position),
                    data = dat_pos,
                    family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  #GxE
  n12 <- stan_glmer(total_flowers ~ scam + level*salinity + location + co2*age + comp + (1+ comp|genotype) +(1|frame) + (1|position),
                    data = dat_pos,
                    family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  n13 <- stan_glmer(total_flowers ~ scam * level +salinity + location + co2*age + comp + (1+ salinity |genotype) +(1|frame) + (1|position),
                    data = dat_pos,
                    family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  n14 <- stan_glmer(total_flowers ~ scam + level*salinity + location +co2 + age + comp + (1+ co2 |genotype) +(1|frame) + (1|position),
                    data = dat_pos,
                    family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  n15 <- stan_glmer(total_flowers ~ scam + level*salinity + location + co2*age + comp + (1+ level |genotype) +(1|frame) + (1|position),
                    data = dat_pos,
                    family = neg_binomial_2, chains = 4, cores = 5, iter = 3000, warmup = 1000)
  
  prob_models <- list(p, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15)
  
  num_models <- list(n, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15)
  
  loos <- lapply(models, loo)
  
  loo_compare(loos)
  
  probresults <-summary(p1, probs = c(0.025, 0.5, 0.975))
  numresults2 <-summary(n13, probs = c(0.025, 0.5, 0.975))
  
  
  # Combined prediction function
  {
    mean_scam <- mean(c(dat_live$scam, dat_pos$scam))
    max(c(dat_live$scam, dat_pos$scam))
    scam_values <- seq(
      min(c(dat_live$scam, dat_pos$scam), na.rm = TRUE),
      max(c(dat_live$scam, dat_pos$scam), na.rm = TRUE),
      length.out = 5
    )
  
    # Create an expanded reference grid of values to predict on
    ref_grid_prob <- expand.grid(
      level = unique(dat_live$level),
      salinity = unique(dat_live$salinity),
      scam = scam_values,
      location = unique(dat_live$location),
      co2 = unique(dat_live$co2),
      age = unique(dat_live$age),
      comp = unique(dat_live$comp)
    )
    
    ref_grid_num <- expand.grid(
      level = unique(dat_pos$level),  # this excludes 5 and 6 because none of them reproduced
      salinity = unique(dat_pos$salinity),
      scam = scam_values,
      location = unique(dat_pos$location),
      co2 = unique(dat_pos$co2),
      age = unique(dat_pos$age),
      comp = unique(dat_pos$comp))
    
    prob_draws <- add_epred_draws(
      p1,
      newdata = ref_grid_prob,
      re_formula = NA)
    
    num_draws <- add_epred_draws(
      n13,
      newdata = ref_grid_num,
      re_formula = NA)
    
    compute_repro_effort <- function(prob_draws, num_draws, group_vars) {
      
      prob_summary <- prob_draws %>%
        group_by(across(all_of(c(".draw", group_vars)))) %>%
        summarise(prob = mean(.epred), .groups = "drop")
      
      num_summary <- num_draws %>%
        group_by(across(all_of(c(".draw", group_vars)))) %>%
        summarise(num = mean(.epred), .groups = "drop")
      
      combined <- left_join(prob_summary, num_summary,
                            by = c(".draw", group_vars)) %>%
        mutate(num = ifelse(is.na(num), 0, num),
               output = prob * num)
      
      combined %>%
        group_by(across(all_of(group_vars))) %>%
        summarise(
          prob_mean  = mean(prob),
          prob_lower = quantile(prob, 0.025),
          prob_upper = quantile(prob, 0.975),
          num_mean   = mean(num),
          num_lower  = quantile(num, 0.025),
          num_upper  = quantile(num, 0.975),
          out_mean   = mean(output),
          out_lower  = quantile(output, 0.025),
          out_upper  = quantile(output, 0.975),
          .groups = "drop")
    }
    summary_biomass_level <- compute_repro_effort(prob_draws, num_draws, 
                                                  c("level", "scam"))
    summary_level_sal <- compute_repro_effort(prob_draws, num_draws, 
                                              c("level", "salinity"))
    
    summary_age_co2 <- compute_repro_effort(prob_draws, num_draws, 
                                            c("co2", "age"))
    
    summary_age <- compute_repro_effort(prob_draws, num_draws, c("age"))
    
    summary_all <- compute_repro_effort(prob_draws, num_draws, 
                                        c("scam", "level", "salinity", "co2", "age", "comp"))
    
  
  