## ----setup, include=FALSE---------------------------------------------------------------
require(knitr)
knitr::opts_chunk$set(echo = TRUE, cache=TRUE, message = F)
r <- getOption("repos")
r["CRAN"] <- "https://ftp.osuosl.org/pub/cran/"
options(repos = r)


## ----eval=F, include=F------------------------------------------------------------------
## devtools::install_github("glmmTMB/glmmTMB/glmmTMB")
## devtools::install_github("bbolker/broom.mixed", type ="source")
## install.packages("broom.mixed")


## ----packages---------------------------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#load or install these packages:
packages <- c("glmmTMB", "tidyverse", "broom.mixed", "bbmle", "sjPlot", "GGally", "lme4")
#run function to install packages
ipak(packages)


## ----load_prepare_data------------------------------------------------------------------
elk <- read.table("Data/lab7_elk_migrant.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
elk$elkuidF <- as.factor(elk$elkuid)
elk2 <- elk[complete.cases(elk[30:31]), ]
elk2$ctotrisk[elk2$ctotrisk>1]=1
elk2$totalherb_sc <- as.numeric(scale(elk2$totalherb))
elk2$ctotrisk_sc <- as.numeric(scale(elk2$ctotrisk))
elk2$ctotrisk2_sc <- as.numeric(scale(elk2$ctotrisk2))
elk2$riskforage_sc <- as.numeric(scale(elk2$riskforage))
elk2$for2_sc <- as.numeric(scale(elk2$for2))
elk2$risk2_sc <- as.numeric(scale(elk2$risk2))


## ---------------------------------------------------------------------------------------
forrisk_sc = glm(used~totalherb_sc+ctotrisk_sc+ctotrisk_sc*totalherb_sc, data=elk2,family=binomial(link="logit"))
summary(forrisk_sc)
ggcoef(forrisk_sc, exclude_intercept = TRUE)
elk2$naive.pred <-predict(forrisk_sc, type = "response")


## ---------------------------------------------------------------------------------------
fr.rc = glmer(used~totalherb_sc+ctotrisk_sc+totalherb_sc*ctotrisk_sc+(ctotrisk_sc|elkuid), data=elk2,family=binomial(link="logit"), verbose=FALSE)
summary(fr.rc)

fixef(fr.rc) # This is the fixed effects coefficients
ranef(fr.rc) # These are the random effects, which in this model is just (1|elkuid), so, one coefficient for each individual elk
elk2$fr.rc.pred <- predict(fr.rc, type = "response")
hist(elk2$fr.rc.pred)


## ---------------------------------------------------------------------------------------
elk2$fr.rc.pred2 <- predict(fr.rc, re.form = NA, type = "response")
summary(elk2$fr.rc.pred2)


## ---------------------------------------------------------------------------------------
elk2$fr.rc.pred3 <- predict(fr.rc, re.form = ~(1|elkuid) , type = "response")
summary(elk2$fr.rc.pred3)
hist(elk2$fr.rc.pred3)


## ----add_weights------------------------------------------------------------------------
elk2 <-
  elk2 %>% 
  as_tibble() %>% 
  mutate(w=if_else(used==0, 5000,1),
         elkuid = as.factor(elkuid),
         used = as.factor(used),
         log_risk = log(ctotrisk),
         log_risk_sc = as.numeric(scale(log_risk))) %>% 
  rename(totalherb2_sc = for2_sc)


## ----intercept_only_model---------------------------------------------------------------
system.time(
  forage_risk_r_int <- glmmTMB(used ~ totalherb_sc + totalherb2_sc + log_risk_sc + 
                               (1|elkuid),
                             weights = w, data=elk2, family=binomial,
                             map=list(theta=factor(NA)),
                             start=list(theta=log(1e3)))
)

summary(forage_risk_r_int)


## ----model_without_fixed_random_int-----------------------------------------------------
system.time(
  forage_risk_r_int_free <- glmmTMB(used ~ totalherb_sc + totalherb2_sc + log_risk_sc + 
                               (1|elkuid),
                             weights = w, data=elk2, family=binomial)
)

summary(forage_risk_r_int_free)
plot_model(forage_risk_r_int_free, transform=NULL)


## ----random_slope_risk_model------------------------------------------------------------
system.time(
  forage_risk_slope_risk <- glmmTMB(used ~ totalherb_sc + totalherb2_sc + log_risk_sc + 
                                    (1|elkuid) + (0 + log_risk_sc|elkuid), 
                                  data=elk2, family=binomial, weights = w, 
                                  map=list(theta=factor(c(NA, 1))),
                                  start=list(theta=c(log(1e3), 0)))
)

summary(forage_risk_slope_risk)
plot_model(forage_risk_slope_risk, transform=NULL)


## ----random_slope_both_model------------------------------------------------------------
system.time(forage_risk_slopes_both <- glmmTMB(used ~ totalherb_sc + totalherb2_sc + log_risk_sc + 
                                     (1|elkuid) + (0 + totalherb_sc|elkuid) + (0 + totalherb2_sc|elkuid) + 
                                     (0 + log_risk_sc|elkuid), 
                                  data=elk2, family=binomial, weights = w, 
                                  map=list(theta=factor(c(NA, 1:3))),
                                  start=list(theta=c(log(1e3), rep(0,3))))
)
            
summary(forage_risk_slopes_both)
sjPlot::plot_model(forage_risk_slopes_both, transform=NULL)


## ---------------------------------------------------------------------------------------
bbmle::AICtab(forage_risk_r_int, forage_risk_slope_risk, forage_risk_slopes_both)    


## ----ran_coefs_using_coef---------------------------------------------------------------
ranef(forage_risk_slopes_both)


## ---------------------------------------------------------------------------------------
coef(forage_risk_slopes_both)

as_tibble(rownames_to_column(coef(forage_risk_slopes_both)$cond$elkuid, "elkuid")) %>% 
  dplyr::select(-"(Intercept)")


## ---------------------------------------------------------------------------------------
as_tibble(rownames_to_column(coef(forage_risk_slopes_both)$cond$elkuid, "elkuid")) %>% 
  dplyr::select(totalherb_sc) %>% 
  ggplot(., aes(x=totalherb_sc)) +
  geom_histogram() +
  geom_vline(xintercept = fixef(forage_risk_slopes_both)$cond["totalherb_sc"], 
             color = "orange", linetype="dashed", size = 1) +
  theme_classic()


## ---------------------------------------------------------------------------------------
broom.mixed::tidy(forage_risk_slopes_both, effects = "ran_vals")


## ----ran_coefs_using_tidy---------------------------------------------------------------
forage_ran_coefs <- broom.mixed::tidy(forage_risk_slopes_both, effects = "ran_vals") %>%
  filter(term=="totalherb_sc") %>% 
  dplyr::select(elkuid=level, estimate, std.error) %>% 
  mutate(forage_coef = estimate + fixef(forage_risk_slopes_both)$cond["totalherb_sc"],
         conf_low = forage_coef - std.error*1.96,
         conf_high = forage_coef + std.error*1.96) 

forage_ran_coefs


## ----plot_random_coefs------------------------------------------------------------------
ggplot(forage_ran_coefs, aes(x=elkuid, y=forage_coef)) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_pointrange(aes(ymin = conf_low,
                        ymax = conf_high),
                    size=1) +
    xlab("Elk ID") +
    theme_bw(base_size = 15)


## ----plot_with_random_and_fixed_coefs---------------------------------------------------
fixed_coef <-
  broom.mixed::tidy(forage_risk_slopes_both, effects="fixed", conf.int = T) %>% 
  filter(term == "totalherb_sc")

ggplot(
  forage_ran_coefs, aes(x=elkuid, y=forage_coef)) +
  coord_flip() +
  geom_rect(ymin=fixed_coef$conf.low, ymax=fixed_coef$conf.high,
            xmin=-Inf,xmax=Inf, fill="red", alpha=0.01) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_hline(yintercept = fixef(forage_risk_slopes_both)$cond["totalherb_sc"],
             linetype="dashed", color="red", size=1) +
  geom_pointrange(aes(ymin = conf_low,
                      ymax = conf_high),
                  size=1) +
  xlab("Elk ID") +
  theme_bw(base_size = 15)


## ----run_unfixed_model------------------------------------------------------------------
system.time(forage_risk_slopes_both_UNFIXED <- glmmTMB(used ~ totalherb_sc + totalherb2_sc + log_risk_sc + 
                                     (1|elkuid) + (0 + totalherb_sc|elkuid) + (0 + totalherb2_sc|elkuid) + 
                                     (0 + log_risk_sc|elkuid), 
                                  data=elk2, family=binomial, weights = w)
)


## ----create_df_for_fixed_preds----------------------------------------------------------
forage_for_predict_population <-
  tibble(elkuid = NA,
         totalherb = seq(min(elk2$totalherb), max(elk2$totalherb), len=100),
         totalherb2 = totalherb^2,
         totalherb_sc = as.numeric(scale(totalherb)),
         totalherb2_sc = as.numeric(scale(totalherb2)),
         log_risk_sc = 0,
         w = NA) 

forage_for_predict_population


## ----make_fixed_preds-------------------------------------------------------------------
pop_pred_unfixed <-
  forage_for_predict_population %>% 
  mutate(pred_LP = predict(forage_risk_slopes_both_UNFIXED, ., re.form=NA),
         pred_real = exp(pred_LP - fixef(forage_risk_slopes_both_UNFIXED)$cond["(Intercept)"]))


## ----plot_fixed_preds-------------------------------------------------------------------
pop_pred_unfixed %>% 
  ggplot(., aes(x=totalherb, y=pred_real)) +
  geom_line(size=1) +
  theme_classic(base_size=15)

pop_pred_unfixed %>% 
  filter(pred_real == max(pred_real)) %>% 
  select(totalherb)


## ----echo=FALSE-------------------------------------------------------------------------
knitr::include_graphics("/Users/mark.hebblewhite/Box Sync/Teaching/UofMcourses/WILD562/Spring2021/Labs/Lab8/Figures/Avgar_figure4.PNG")


## ----plogis_example---------------------------------------------------------------------
  forage_for_predict_population %>% 
  mutate(pred_LP = predict(forage_risk_slopes_both_UNFIXED, ., re.form=NA),
         pred_real = plogis(pred_LP - fixef(forage_risk_slopes_both_UNFIXED)$cond["(Intercept)"])) %>% 
  ggplot(., aes(x=totalherb, y=pred_real)) +
  geom_line(size=1) +
  theme_classic(base_size=15)


## ----fixed_predict_with_predict---------------------------------------------------------
pop_pred_unfixed %>% 
  mutate(pred_01 = (pred_real - min(pred_real))/diff(range(pred_real))) %>% 
  ggplot(., aes(x=totalherb, y=pred_01)) +
  geom_line(size=1) +
  theme_classic(base_size=15)


## ----manual_fixed_effects_pred----------------------------------------------------------
coefs_fixed <-
  fixef(forage_risk_slopes_both_UNFIXED)$cond

pop_pred_unfixed_manual <- 
  forage_for_predict_population %>% 
  mutate(pred_real = exp(totalherb_sc*coefs_fixed["totalherb_sc"] +
                           totalherb2_sc*coefs_fixed["totalherb2_sc"]))

ggplot(pop_pred_unfixed_manual, aes(x=totalherb, y=pred_real)) +
  geom_line(size=1) + 
  theme_classic(base_size=15)
  
pop_pred_unfixed_manual %>% 
  filter(pred_real == max(pred_real)) %>% 
  select(totalherb)


## ----custom_predict_CI_function_for_fixed-----------------------------------------------
pred_CI <- function(model, newdata=NULL, alpha=0.05) {
  pred0 <- predict(model, re.form=NA, newdata=newdata) - fixef(model)$cond["(Intercept)"]
  X <- model.matrix(formula(model, fixed.only=TRUE)[-2], newdata)[,-1]
  V <- vcov(model)$cond[-1,-1]     
  pred_se <- sqrt(diag(X %*% V %*% t(X))) 
  crit <- -qnorm(alpha/2)
  pred_df <- as_tibble(exp(cbind(pred=pred0, conf.low=pred0-crit*pred_se,
                      conf.high=pred0+crit*pred_se)))
  bind_cols(newdata, pred_df)
}

pred_CI(forage_risk_slopes_both_UNFIXED, forage_for_predict_population)

pred_CI(forage_risk_slopes_both_UNFIXED, forage_for_predict_population) %>% 
  ggplot(., aes(x=totalherb, y=pred)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill="lightgray") +
  geom_line(size=1) + 
  theme_classic(base_size=15) 


## ----create_ind_df_for_predict----------------------------------------------------------
for_predict_ind <-
  elk2 %>% 
  select(elkuid, totalherb) %>%
  nest(data = c(totalherb)) %>%
  mutate(totalherb = map(data, ~seq(min(.), max(.), len=100))) %>%
  unnest(totalherb) %>% 
  mutate(totalherb2 = totalherb^2,
         totalherb_sc = as.numeric(scale(totalherb)),
         totalherb2_sc = as.numeric(scale(totalherb2)),
         log_risk = 0,
         w=NA) %>% 
  select(-data) 

for_predict_ind


## ----get_ind_coefs_for_predict----------------------------------------------------------
coefs_ind <-   
  rownames_to_column(coef(forage_risk_slopes_both)$cond$elkuid, "elkuid") %>% 
  rename(ran_int = `(Intercept)`,
         forage_coef = totalherb_sc,
         forage2_coef = totalherb2_sc) 

for_predict_ind <-
  for_predict_ind %>% 
  inner_join(coefs_ind)


## ---------------------------------------------------------------------------------------
for_predict_ind %>%
  mutate(pred = exp(totalherb_sc*forage_coef + totalherb2_sc*forage2_coef)) %>% 
  ggplot(., aes(x=totalherb, y=pred, color=elkuid)) +
  geom_line(size=1) +
  theme_classic(base_size=15)


## ---------------------------------------------------------------------------------------
for_predict_ind %>%
  mutate(pred = exp(totalherb_sc*forage_coef + totalherb2_sc*forage2_coef)) %>% 
  ggplot(., aes(x=totalherb, y=pred, color=elkuid)) +
  geom_line(size=1) +
  theme_classic(base_size=15) +
  coord_cartesian(xlim= c(0,200), ylim=c(0,10))


## ---------------------------------------------------------------------------------------
fixed_response <- pred_CI(forage_risk_slopes_both, forage_for_predict_population)

for_predict_ind %>%
  mutate(pred = exp(totalherb_sc*forage_coef + totalherb2_sc*forage2_coef)) %>% 
  ggplot(., aes(x=totalherb, y=pred, color=elkuid)) +
    geom_ribbon(data=fixed_response, aes(ymin = conf.low, ymax = conf.high), fill="lightgray", alpha = .4, color=NA) +
    geom_line(size=1) +
    geom_line(data = fixed_response, size = 2, color="black") + 
    theme_classic(base_size=15) +
    coord_cartesian(xlim= c(0,200), ylim=c(0,10))


## ---------------------------------------------------------------------------------------
for_predict_ind %>%
  mutate(pred = exp(predict(forage_risk_slopes_both, .) - ran_int)) %>% 
  ggplot(., aes(x=totalherb, y=pred, color=elkuid)) +
    geom_ribbon(data=fixed_response, aes(ymin = conf.low, ymax = conf.high), fill="lightgray", alpha = .4, color=NA) +
    geom_line(size=1) +
    geom_line(data = fixed_response, size = 2, color="black") + 
    theme_classic(base_size=15) +
    coord_cartesian(xlim= c(0,200), ylim=c(0,10))


## ---------------------------------------------------------------------------------------
elk2 %>%
  filter(used==0) %>%
  mutate(elkuid=NA,
         predict = exp(predict(forage_risk_slopes_both, ., re.form = NA) - fixef(forage_risk_slopes_both)$cond["(Intercept)"])) %>%
  ggplot(., aes(x = totalherb, y = predict)) +
  geom_smooth(alpha=0.2, size=1.6) +
  theme_classic(base_size=20)


## ---------------------------------------------------------------------------------------
elk2 <-  
  elk2 %>% 
    group_by(elkuid) %>% 
    mutate(mean_risk = mean(log_risk[used==0]),
           mean_herb = mean(totalherb[used==0]),
           mean_herb2 = mean(for2[used==0])) %>% 
    ungroup() %>% 
    mutate(mean_risk_sc = as.numeric(scale(mean_risk)),
           mean_herb_sc = as.numeric(scale(mean_herb)),
           mean_herb2_sc = as.numeric(scale(mean_herb2)))


## ---------------------------------------------------------------------------------------
hist(elk2$mean_risk_sc)
hist(elk2$mean_herb2)


## ----functional_response_model----------------------------------------------------------
forage_risk_slopes_both_FR <- glmmTMB(used ~ totalherb_sc + totalherb2_sc + log_risk_sc + 
                                        log_risk_sc:mean_risk_sc + totalherb_sc:mean_herb_sc + 
                                        totalherb2_sc:mean_herb2_sc + 
                                        (1|elkuid) + (0 + totalherb_sc|elkuid) + (0 + totalherb2_sc|elkuid) +
                                        (0 + log_risk_sc|elkuid), 
                                   data=elk2, family=binomial, weights = w, 
                                   map=list(theta=factor(c(NA, 1:3))),
                                   start=list(theta=c(log(1e3), rep(0,3))))

plot_model(forage_risk_slopes_both_FR, transform = NULL)


## ----AIC_model_comparison---------------------------------------------------------------
bbmle::AICtab(forage_risk_slopes_both, forage_risk_slopes_both_FR)

