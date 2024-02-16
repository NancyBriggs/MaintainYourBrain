library(lme4)
library(emmeans)
library(tidyverse)

options(warn = 0)


cond =resid(., type = "pearson") ~ predict(., re.form = NA)
plannedtimes = factor(c(12,1),levels = c(1,12))
## fit all models

MYB_fit <- function(response, predictors, covariates, id, dat, times){
  
  covar = paste(covariates, collapse = " + ")
  ran_eff = paste0("(1| ",id,")")
  main = paste(predictors, collapse = "+")
  inter = paste(predictors, collapse = "*")
  
  dat <-  dat %>% 
    mutate_at(.vars = predictors[1], .funs = as.factor)
  
  unadj_mod <- lmer(paste(response, " ~ ", inter,  "+", ran_eff ), dat)
  unadj_null <- lmer(paste(response, " ~ ", main,  "+", ran_eff ), dat)
  
  adj_mod <- lmer(paste(response, " ~ ", covar , "+" , inter,  "+", ran_eff ), dat)
  adj_null <- lmer(paste(response, " ~ ", covar , "+" , main,  "+", ran_eff ), dat)
  
  
  
  return(list(adj_mod = adj_mod,
              unadj_null = unadj_null,
              unadj_mod = unadj_mod,
              adj_null = adj_null))
}


# inference on one model

inf_single <- function(mod_int, mod_main){
  
  anova_call <- anova(mod_main,mod_int)
  likrat <-  data.frame(Chisq = paste0(round(anova_call$Chisq[2],3),"(",
                                       anova_call$Df[2],")"),
                        pvalue = anova_call$`Pr(>Chisq)`[2])
  em <- emmeans(mod_int,  ~ rc_trial_time*group_type_id , type = "response",  
                at = list(rc_trial_time = plannedtimes) )
  cont <- contrast(em, interaction = c("consec", "consec"))
  CI <- confint(cont) 
  Su <- summary(cont)
  p_0_4 <- data.frame(zstat = round(Su$z.ratio,3),
                      pz = ifelse (Su$p.value < 0.001, "<0.001",paste("= ",round(Su$p.value,3))))
  cont <- cont %>% as.data.frame()
  coef <- try(cont$estimate)
  if(is.null(coef)){
    coef <- cont$ratio
  }
  confint <- try(CI[,c("asymp.LCL", "asymp.UCL")])
  if(class(confint) == "try-error"){
    confint <- try(CI[,c("lower.CL", "upper.CL")])
    names(confint) = c("asymp.LCL", "asymp.UCL")
  }
  
  
  cbind(coef = coef, confint, likrat,p_0_4)

}





## calculate all inference

MYB_inf <- function(fitted_mods, response){
  
  rbind(inf_single(fitted_mods$unadj_mod, fitted_mods$unadj_null),
        inf_single(fitted_mods$adj_mod, fitted_mods$adj_null)) %>% 
    as.data.frame() %>% 
    remove_rownames %>% 
    mutate(response = response,
           adjusted = c("No","Yes")) %>% 
    transmute(response, adjusted,
              planned = coef, asymp.LCL = asymp.LCL, asymp.UCL = asymp.UCL, 
              zstat, pz,
              Chisq, p.value = pvalue)
    
}

