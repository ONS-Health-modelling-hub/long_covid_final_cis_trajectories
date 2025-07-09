library(dplyr)
library(splines)
library(lme4)
library(survey)
library(car)
library(emmeans)
library(ggplot2)

source("filepath/emmeans_collate.r")

############################### PREPARE DATA ###############################

### load dataset
load("filepath/lc_trajectories_dataset.RData")

### set up factors
dat$age10_visit0 <- as.factor(dat$age10_visit0)

dat$sex <- as.factor(dat$sex)
levels(dat$sex) <- c("Male", "Female")

dat$white <- as.factor(dat$white)
levels(dat$white) <- c("Non-White", "White")

dat$gor9d <- as.factor(dat$gor9d)
levels(dat$gor9d) <- c("North East", "North West", "Yorkshire and the Humber",
                       "East Midlands", "West Midlands", "East of England",
                       "London", "South East", "South West",
                       "Northern Ireland", "Scotland", "Wales")

dat$imd_quintile <- as.factor(dat$imd_quintile)

dat$health_status_visit0[dat$health_status_visit0==3] <- 2
dat$health_status_visit0 <- as.factor(dat$health_status_visit0)
levels(dat$health_status_visit0) <- c("No health conditions",
                                      "Health conditions without activity-limitation",
                                      "Health conditions with activity-limitation")

dat$infection_period <- ifelse(dat$infection_date<='2020-12-11', "Pre-Alpha", "Alpha")
dat$infection_period <- as.factor(dat$infection_period)
dat$infection_period <- relevel(dat$infection_period, ref="Pre-Alpha")

### derive calendar day of infection (day 1 = first day of study period)
dat$infection_day <- as.numeric(dat$infection_date) - as.numeric(min(dat$infection_date)) + 1

### rescale sampling weights so that they sum to sample size
dat$svywgt <- dat$svywgt * (sum(!duplicated(dat$participant_id)) /
                            sum(dat$svywgt[!duplicated(dat$participant_id)]))

### derive weights for marginal probabilities: covariate distribution of participants
### at enrolment
emmeans_wgts <- dat %>%
  distinct(participant_id, .keep_all=TRUE) %>%
  group_by(health_status_visit0, imd_quintile, gor9d,
           white, sex, age10_visit0, .drop=FALSE) %>%
  summarise(n=sum(svywgt), .groups="drop")

### derive study period mid-point for marginal probabilities
mid_study <- median(c(min(dat$infection_day),
                      max(dat$infection_day)))

### restrict dataset to visits >= 12 weeks after infection
dat <- dat[dat$time_since_infection >= 84,]

### restrict dataset to visit <= 2.5 years after infection
dat <- dat[dat$time_since_infection <= 365.25*2.5,]

### restrict dataset to visits when participants reported their LC status
dat <- dat[!is.na(dat$lc_any),]

############################### OPTIMISE SPLINES ###############################

### grid search across internal and boundary knots for time since infection, according to BIC
### note: not worrying about clustering when fitting model as not interested in SEs
### just perform the grid search for the primary outcome (LC of any severity)

df_list <- 1:4
boundary_list <- list(c(0.00,1.00), c(0.01,0.99), c(0.05,0.95), c(0.10,0.90))
bic_out <- as.data.frame(matrix(NA, nrow=length(df_list), ncol=length(boundary_list)))

for(j in 1:length(boundary_list)) {

  bic_list <- as.list(NULL)
  
  for(i in 1:length(df_list)) {

    mod_spline <- glm(lc_any ~ ns(time_since_infection, df=df_list[i], Boundary.knots=quantile(time_since_infection, probs=boundary_list[[j]])),
                      data=dat, family=binomial, weights=svywgt)
  
    bic_list[[i]] <- BIC(mod_spline)

  }

  bic_out[,j] <- unlist(bic_list)
  
}

bic_out <- cbind(df_list-1, bic_out)
colnames(bic_out) <- c("n_internal_knots", paste0("boundary_knots_", boundary_list))
write.csv(bic_out, file="filepath/knots_grid_search_time_since_infection.csv", row.names=FALSE)

### pick out internal and boundary knots that minimise BIC
min_bic_row <- which(bic_out[,-1] == min(bic_out[,-1]), arr.ind=TRUE)[1,1]
min_bic_col <- which(bic_out[,-1] == min(bic_out[,-1]), arr.ind=TRUE)[1,2]
df_select1 <- df_list[min_bic_row]
boundary_select1 <- boundary_list[[min_bic_col]]

### repeat grid search, this for calendar time of infection

df_list <- 1:4
boundary_list <- list(c(0.00,1.00), c(0.01,0.99), c(0.05,0.95), c(0.10,0.90))
bic_out <- as.data.frame(matrix(NA, nrow=length(df_list), ncol=length(boundary_list)))

for(j in 1:length(boundary_list)) {

  bic_list <- as.list(NULL)
  
  for(i in 1:length(df_list)) {

    mod_spline <- glm(lc_any ~ ns(infection_day, df=df_list[i], Boundary.knots=quantile(infection_day, probs=boundary_list[[j]])),
                      data=dat, family=binomial, weights=svywgt)
  
    bic_list[[i]] <- BIC(mod_spline)

  }

  bic_out[,j] <- unlist(bic_list)
  
}

bic_out <- cbind(df_list-1, bic_out)
colnames(bic_out) <- c("n_internal_knots", paste0("boundary_knots_", boundary_list))
write.csv(bic_out, file="filepath/knots_grid_search_infection_day.csv", row.names=FALSE)

### pick out internal and boundary knots that minimise BIC
min_bic_row <- which(bic_out[,-1] == min(bic_out[,-1]), arr.ind=TRUE)[1,1]
min_bic_col <- which(bic_out[,-1] == min(bic_out[,-1]), arr.ind=TRUE)[1,2]
df_select2 <- df_list[min_bic_row]
boundary_select2 <- boundary_list[[min_bic_col]]

### print optimal knot positions
### NOTE: these are hard-coded in splines below as macro variables cannot be
### incorporated into survey design object
print(df_select1)
print(boundary_select1)
print(df_select2)
print(boundary_select2)

############################### CHECK FOR CLUSTERING ###############################

### frequency distribution of participants-within-households cluster size
cluster_sizes_hh <- dat %>%
  distinct(participant_id, .keep_all=TRUE) %>%
  group_by(hh_id_fake) %>%
  summarise(n_participants=n()) %>%
  group_by(n_participants) %>%
  summarise(freq_hhs=n()) %>%
  mutate(pct_hhs = freq_hhs / sum(freq_hhs) * 100,
         pct_participants = n_participants*freq_hhs / sum(n_participants*freq_hhs) * 100)

write.csv(cluster_sizes_hh, file="filepath/cluster_sizes_hh.csv", row.names=FALSE)

### frequency distribution of observations-within-participants cluster size
cluster_sizes_participant <- dat %>%
  group_by(participant_id) %>%
  summarise(n_observations=n()) %>%
  group_by(n_observations) %>%
  summarise(freq_participants=n()) %>%
  mutate(pct_participants = freq_participants / sum(freq_participants) * 100,
         pct_observations = n_observations*freq_participants / sum(n_observations*freq_participants) * 100)

write.csv(cluster_sizes_participant, file="filepath/cluster_sizes_participant.csv", row.names=FALSE)

### calculate participant- and household-level contributions to variance in LC
### note: using linear mixed-effect models (ie linear probability model) rather than
### logistic mixed-effect models due to processing constraints
lmer_mod1 <- lmer(lc_any ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))) +
                    ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                    age10_visit0 + sex + white + gor9d +
                    imd_quintile + health_status_visit0 +
                    (1|hh_id_fake) + (1|hh_id_fake:participant_id),
                  data=dat, weights=svywgt)

var_comps <- as.data.frame(VarCorr(lmer_mod1, comp="Variance"))
var_comps$contribution <- var_comps$vcov / sum(var_comps$vcov) * 100
write.csv(var_comps, file="LC trajectories/lmer_variance_components.csv", row.names=FALSE)

### test signifiance of participant random intercepts vs. no random effects, and
### participant-within-household random intercepts vs. participant-only random intercepts
lmer_mod2 <- lmer(lc_any ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))) +
                    ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                    age10_visit0 + sex + white + gor9d +
                    imd_quintile + health_status_visit0 +
                    (1|participant_id),
                  data=dat, weights=svywgt)

lmer_mod3 <- lm(lc_any ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))) +
                  ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                  age10_visit0 + sex + white + gor9d +
                  imd_quintile + health_status_visit0,
                data=dat, weights=svywgt)

lr_participant_re <- as.data.frame(anova(lmer_mod2, lmer_mod3))
lr_hh_re <- as.data.frame(anova(lmer_mod1, lmer_mod2))

write.csv(lr_participant_re, file="filepath/lmer_lr_participant_re.csv")
write.csv(lr_hh_re, file="filepath/lmer_lr_hh_re.csv")

############################### MAIN ANALYSIS ###############################

### set up survey design object
svy_obj <- svydesign(ids=~participant_id, strata=~cis20_samp, weights=~svywgt, data=dat)

### fit model - LC of any severity

mod_any <- svyglm(lc_any ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))) +
                    ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                    age10_visit0 + sex + white + gor9d +
                    imd_quintile + health_status_visit0,
                  family=binomial, design=svy_obj)

vcov_any <- vcov(mod_any)

coeff_any <- as.data.frame(summary(mod_any)$coefficients)
write.csv(coeff_any, file="filepath/coeff_any.csv")

lr_any <- as.data.frame(Anova(mod_any, test.statistic="Chisq"))
write.csv(lr_any, file="filepath/lr_any.csv")

### fit model - activity-limiting LC

mod_lim <- svyglm(lc_lim ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))) +
                    ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                    age10_visit0 + sex + white + gor9d +
                    imd_quintile + health_status_visit0,
                  family=binomial, design=svy_obj)

vcov_lim <- vcov(mod_lim)

coeff_lim <- as.data.frame(summary(mod_lim)$coefficients)
write.csv(coeff_lim, file="filepath/coeff_lim.csv")

lr_lim <- as.data.frame(Anova(mod_lim, test.statistic="Chisq"))
write.csv(lr_lim, file="filepath/lr_lim.csv")

### marginal probabilities

emmeans_any <- emmeans.collate(mod_any, vcov_any, emmeans_wgts$n, mid_study)
emmeans_lim <- emmeans.collate(mod_lim, vcov_lim, emmeans_wgts$n, mid_study)

write.csv(emmeans_any, file="filepath/emmeans_any.csv", row.names=FALSE)
write.csv(emmeans_lim, file="filepath/emmeans_lim.csv", row.names=FALSE)

emmeans_any$outcome <- "Any severity"
emmeans_lim$outcome <- "Activity limiting"
emmeans_all <- rbind(emmeans_any, emmeans_lim)
emmeans_all$outcome <- factor(emmeans_all$outcome, levels=c("Any severity", "Activity limiting"))
emmeans_all$time <- emmeans_all$time_since_infection / (365.25/12)
emmeans_all$prevalence <- emmeans_all$prob * 100
emmeans_all$lcl <- emmeans_all$asymp.LCL * 100
emmeans_all$ucl <- emmeans_all$asymp.UCL * 100

ggplot(emmeans_all, aes(x=time, y=prevalence)) +
  geom_line(linewidth=0.8, colour="#206095") +
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2, fill="#206095") +
  facet_wrap(~outcome) +
  xlab("Time since infection (months)") +
  ylab("Prevalence (percent)") +
  scale_x_continuous(breaks=seq(0,30,6), limits=c(0,30), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,NA), expand=expansion(mult=c(0,0.05))) +
  theme(
    axis.text = element_text(colour="black", size=11),
    axis.title = element_text(colour="black", size=11),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(colour="grey", linewidth=0.2),
    strip.text = element_text(colour="black", size=11, face="bold"),
    strip.background = element_blank(),
    plot.margin = margin(0.2,0.5,0.2,0.2,"cm"),
    panel.spacing = unit(0.75,"cm")
  )

ggsave(filename="filepath/emmeans.jpg", width=15, height=10, units="cm", dpi=300)

### pairwise comparison of LC prevalence at 12 weeks vs. 2.5 years

pairs_any <- as.data.frame(pairs(emmeans(
  object = mod_any,
  vcov. = vcov_any,
  specs = ~time_since_infection,
  weights = emmeans_wgts$n,
  at = list(time_since_infection = c(max(dat$time_since_infection), min(dat$time_since_infection)),
            infection_day = mid_study),
  type = "response",
  params = c("df_select1", "boundary_select1", "df_select2", "boundary_select2"),
  rg.limit = 10000
)))

pairs_lim <- as.data.frame(pairs(emmeans(
  object = mod_lim,
  vcov. = vcov_lim,
  specs = ~time_since_infection,
  weights = emmeans_wgts$n,
  at = list(time_since_infection = c(max(dat$time_since_infection), min(dat$time_since_infection)),
            infection_day = mid_study),
  type = "response",
  params = c("df_select1", "boundary_select1", "df_select2", "boundary_select2"),
  rg.limit = 10000
)))

pairs_any$outcome <- "Any severity"
pairs_lim$outcome <- "Activity limiting"
pairs_all <- rbind(pairs_any, pairs_lim)
pairs_all <- pairs_all[,c(8,1:3,6:7)]
write.csv(pairs_all, file="filepath/pairs.csv", row.names=FALSE)

############################ UNWEIGHTED & UNADJUSTED MODELS ############################

### fit model - LC of any severity

mod_any_unadj <- svyglm(lc_any ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))),
                        family=binomial, design=svy_obj)

vcov_any_unadj <- vcov(mod_any_unadj)

coeff_any_unadj <- as.data.frame(summary(mod_any_unadj)$coefficients)
write.csv(coeff_any_unadj, file="filepath/coeff_any_unadj.csv")

lr_any_unadj <- as.data.frame(Anova(mod_any_unadj, test.statistic="Chisq"))
write.csv(lr_any_unadj, file="filepath/lr_any_unadj.csv")

### fit model - activity-limiting LC

mod_lim_unadj <- svyglm(lc_lim ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))),
                        family=binomial, design=svy_obj)

vcov_lim_unadj <- vcov(mod_lim_unadj)

coeff_lim_unadj <- as.data.frame(summary(mod_lim_unadj)$coefficients)
write.csv(coeff_lim_unadj, file="filepath/coeff_lim_unadj.csv")

lr_lim_unadj <- as.data.frame(Anova(mod_lim_unadj, test.statistic="Chisq"))
write.csv(lr_lim_unadj, file="filepath/lr_lim_unadj.csv")

### marginal probabilities

emmeans_any_unadj <- as.data.frame(emmeans(
  object = mod_any_unadj,
  vcov. = vcov_any_unadj,
  specs = ~time_since_infection,
  at = list(time_since_infection = 84:max(dat$time_since_infection)),
  type = "response",
  params = c("df_select1", "boundary_select1", "df_select2", "boundary_select2")
))

write.csv(emmeans_any_unadj, file="filepath/emmeans_any_unadj.csv", row.names=FALSE)

emmeans_lim_unadj <- as.data.frame(emmeans(
  object = mod_lim_unadj,
  vcov. = vcov_lim_unadj,
  specs = ~time_since_infection,
  at = list(time_since_infection = 84:max(dat$time_since_infection)),
  type = "response",
  params = c("df_select1", "boundary_select1", "df_select2", "boundary_select2")
))

write.csv(emmeans_lim_unadj, file="filepath/emmeans_lim_unadj.csv", row.names=FALSE)

emmeans_any$model <- "Adjusted"
emmeans_any$outcome <- "Any severity"

emmeans_lim$model <- "Adjusted"
emmeans_lim$outcome <- "Activity limiting"

emmeans_any_unadj$model <- "Unadjusted"
emmeans_any_unadj$outcome <- "Any severity"

emmeans_lim_unadj$model <- "Unadjusted"
emmeans_lim_unadj$outcome <- "Activity limiting"

emmeans_all_unadj <- rbind(emmeans_any, emmeans_lim, emmeans_any_unadj, emmeans_lim_unadj)
emmeans_all_unadj$model <- factor(emmeans_all_unadj$model, levels=c("Unadjusted", "Adjusted"))
emmeans_all_unadj$outcome <- factor(emmeans_all_unadj$outcome, levels=c("Any severity", "Activity limiting"))
emmeans_all_unadj$time <- emmeans_all_unadj$time_since_infection / (365.25/12)
emmeans_all_unadj$prevalence <- emmeans_all_unadj$prob * 100
emmeans_all_unadj$lcl <- emmeans_all_unadj$asymp.LCL * 100
emmeans_all_unadj$ucl <- emmeans_all_unadj$asymp.UCL * 100

ggplot(emmeans_all_unadj, aes(x=time, y=prevalence, group=model, colour=model, fill=model)) +
  geom_line(linewidth=0.8) +
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2, linewidth=NA) +
  facet_wrap(~outcome) +
  xlab("Time since infection (months)") +
  ylab("Prevalence (percent)") +
  scale_x_continuous(breaks=seq(0,30,6), limits=c(0,30), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,NA), expand=expansion(mult=c(0,0.05))) +
  scale_colour_manual(values=c("#206095","#f39431")) +
  scale_fill_manual(values=c("#206095","#f39431")) +
  theme(
    axis.text = element_text(colour="black", size=11),
    axis.title = element_text(colour="black", size=11),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(colour="grey", linewidth=0.2),
    strip.text = element_text(colour="black", size=11, face="bold"),
    strip.background = element_blank(),
    plot.margin = margin(0.2,0.5,0.2,0.2,"cm"),
    panel.spacing = unit(0.75,"cm")
  )

ggsave(filename="filepath/emmeans_unadj.jpg", width=15, height=10, units="cm", dpi=300)

### pairwise comparison of LC prevalence at 12 weeks vs. 2.5 years

pairs_any_unadj <- as.data.frame(pairs(emmeans(
  object = mod_any_unadj,
  vcov. = vcov_any_unadj,
  specs = ~time_since_infection,
  at = list(time_since_infection = c(max(dat$time_since_infection), min(dat$time_since_infection))),
  type = "response",
  params = c("df_select1", "boundary_select1", "df_select2", "boundary_select2")
)))

pairs_lim_unadj <- as.data.frame(pairs(emmeans(
  object = mod_lim_unadj,
  vcov. = vcov_lim_unadj,
  specs = ~time_since_infection,
  at = list(time_since_infection = c(max(dat$time_since_infection), min(dat$time_since_infection))),
  type = "response",
  params = c("df_select1", "boundary_select1", "df_select2", "boundary_select2")
)))

pairs_any_unadj$outcome <- "Any severity"
pairs_lim_unadj$outcome <- "Activity limiting"
pairs_all_unadj <- rbind(pairs_any_unadj, pairs_lim_unadj)
pairs_all_unadj <- pairs_all_unadj[,c(8,1:3,6:7)]
write.csv(pairs_all_unadj, file="filepath/pairs_unadj.csv", row.names=FALSE)

####################### EFFECT MODIFICATION: SOCIO-DEMOGRAPHICS #######################

### list model variables
mod_vars <- c("age10_visit0", "sex", "white", "gor9d",
              "imd_quintile", "health_status_visit0",
              "ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9)))")

### loop over model variables
for(i in 1:(length(mod_vars)-1)) {

  ### select effect-modifier of interest
  modifier = mod_vars[i]
  
  ### list all other model-variabels as non-modifiers
  non_modifiers <- mod_vars[mod_vars!=modifier]
  
  ### construct the model formula
  spl <- paste0("ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1)))")
  int <- paste0(spl, "*", modifier)
  covs <- paste0(non_modifiers, collapse="+")
  rhs <- paste0(c(int, covs), collapse="+")
  mod_formula <- as.formula(paste0("lc_any ~ ", rhs))
  
  ### fit the model
  mod_het <- svyglm(formula=mod_formula, family=binomial, design=svy_obj)
  
  vcov_het <- vcov(mod_het)
  
  coeff_het <- as.data.frame(summary(mod_het)$coefficients)
  write.csv(coeff_het, file=paste0("filepath/coeff_", modifier, ".csv"))
  
  lr_het <- as.data.frame(Anova(mod_het, test.statistic="Chisq"))
  write.csv(lr_het, file=paste0("filepath/lr_", modifier, ".csv"))
  
  ### derive weights for marginal probabilities
  group_by_vars <- rev(mod_vars)
  group_by_vars <- group_by_vars[group_by_vars!=modifier]
  group_by_vars <- group_by_vars[group_by_vars!="ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9)))"]
  
  emmeans_wgts_het <- emmeans_wgts %>%
    group_by_at(group_by_vars, .drop=FALSE) %>%
    summarise(n=sum(n), .groups="drop")
  
  ### marginal probabilities
  emmeans_het <- as.data.frame(emmeans(
    object = mod_het,
    vcov. = vcov_het,
    specs = as.formula(paste0("~time_since_infection | ", modifier)),
    weights = emmeans_wgts_het$n,
    at = list(time_since_infection = c(84, seq(0.5,2.5,0.5)*365.25),
              infection_day = mid_study),
    type = "response",
    params = c("df_select1", "boundary_select1",
               "df_select2", "boundary_select2",
               "mid_study"),
    rg.limit = 100000
  ))
  
  write.csv(emmeans_het, file=paste0("filepath/emmeans_", modifier, ".csv"), row.names=FALSE)
  
}

####################### EFFECT MODIFICATION: VARIANT PERIOD #######################

### fit the model
mod_period <- svyglm(lc_any ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1)))*infection_period +
                       age10_visit0 + sex + white + gor9d +
                       imd_quintile + health_status_visit0,
                     family=binomial, design=svy_obj)

vcov_period <- vcov(mod_period)

coeff_period <- as.data.frame(summary(mod_period)$coefficients)
write.csv(coeff_period, file="filepath/coeff_infection_period.csv")

lr_period <- as.data.frame(Anova(mod_period, test.statistic="Chisq"))
write.csv(lr_period, file="filepath/lr_infection_period.csv")

### marginal probabilities
emmeans_period <- as.data.frame(emmeans(
  object = mod_period,
  vcov. = vcov_period,
  specs = ~time_since_infection | infection_period,
  weights = emmeans_wgts$n,
  at = list(time_since_infection = c(84, seq(0.5,2.5,0.5)*365.25)),
  type = "response",
  params = c("df_select1", "boundary_select1"),
  rg.limit = 100000
))

write.csv(emmeans_period, file=paste0("filepath/emmeans_infection_period.csv"), row.names=FALSE)

############################### SENSITIVITY ANALYSES ###############################

### SA1: censor follow-up at reinfection

dat_sa1 <- dat[dat$reinfected_to_date==0,]
svy_obj_sa1 <- svydesign(ids=~participant_id, strata=~cis20_samp, weights=~svywgt, data=dat_sa1)

mod_sa1_any <- svyglm(lc_any ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))) +
                          ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                          age10_visit0 + sex + white + gor9d +
                          imd_quintile + health_status_visit0,
                      family=binomial, design=svy_obj_sa1)

vcov_sa1_any <- vcov(mod_sa1_any)

coeff_sa1_any <- as.data.frame(summary(mod_sa1_any)$coefficients)
write.csv(coeff_sa1_any, file="filepath/coeff_sa1_any.csv")

lr_sa1_any <- as.data.frame(Anova(mod_sa1_any, test.statistic="Chisq"))
write.csv(lr_sa1_any, file="filepath/lr_sa1_any.csv")

mod_sa1_lim <- svyglm(lc_lim ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))) +
                        ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                        age10_visit0 + sex + white + gor9d +
                        imd_quintile + health_status_visit0,
                      family=binomial, design=svy_obj_sa1)

vcov_sa1_lim <- vcov(mod_sa1_lim)

coeff_sa1_lim <- as.data.frame(summary(mod_sa1_lim)$coefficients)
write.csv(coeff_sa1_lim, file="filepath/coeff_sa1_lim.csv")

lr_sa1_lim <- as.data.frame(Anova(mod_sa1_lim, test.statistic="Chisq"))
write.csv(lr_sa1_lim, file="filepath/lr_sa1_lim.csv")

emmeans_sa1_any <- emmeans.collate(mod_sa1_any, vcov_sa1_any, emmeans_wgts$n, mid_study)
emmeans_sa1_any <- emmeans_sa1_any[emmeans_sa1_any$time_since_infection <= max(dat$time_since_infection[dat$reinfected_to_date==0]),]
write.csv(emmeans_sa1_any, file="filepath/emmeans_sa1_any.csv", row.names=FALSE)

emmeans_sa1_lim <- emmeans.collate(mod_sa1_lim, vcov_sa1_lim, emmeans_wgts$n, mid_study)
emmeans_sa1_lim <- emmeans_sa1_lim[emmeans_sa1_lim$time_since_infection <= max(dat$time_since_infection[dat$reinfected_to_date==0]),]
write.csv(emmeans_sa1_lim, file="filepath/emmeans_sa1_lim.csv", row.names=FALSE)

emmeans_any$analysis <- "Main analysis"
emmeans_any$outcome <- "Any severity"

emmeans_lim$analysis <- "Main analysis"
emmeans_lim$outcome <- "Activity limiting"

emmeans_sa1_any$analysis <- "Censor at reinfection"
emmeans_sa1_any$outcome <- "Any severity"

emmeans_sa1_lim$analysis <- "Censor at reinfection"
emmeans_sa1_lim$outcome <- "Activity limiting"

emmeans_sa1 <- rbind(emmeans_any, emmeans_lim, emmeans_sa1_any, emmeans_sa1_lim)
emmeans_sa1$analysis <- factor(emmeans_sa1$analysis, levels=c("Censor at reinfection", "Main analysis"))
emmeans_sa1$outcome <- factor(emmeans_sa1$outcome, levels=c("Any severity", "Activity limiting"))
emmeans_sa1$time <- emmeans_sa1$time_since_infection / (365.25/12)
emmeans_sa1$prevalence <- emmeans_sa1$prob * 100
emmeans_sa1$lcl <- emmeans_sa1$asymp.LCL * 100
emmeans_sa1$ucl <- emmeans_sa1$asymp.UCL * 100

ggplot(emmeans_sa1, aes(x=time, y=prevalence, group=analysis, colour=analysis, fill=analysis)) +
  geom_line(linewidth=0.8) +
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2, linewidth=NA) +
  facet_wrap(~outcome) +
  xlab("Time since infection (months)") +
  ylab("Prevalence (percent)") +
  scale_x_continuous(breaks=seq(0,30,6), limits=c(0,30), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,NA), expand=expansion(mult=c(0,0.05))) +
  scale_colour_manual(values=c("#206095","#f39431")) +
  scale_fill_manual(values=c("#206095","#f39431")) +
  theme(
    axis.text = element_text(colour="black", size=11),
    axis.title = element_text(colour="black", size=11),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(colour="grey", linewidth=0.2),
    strip.text = element_text(colour="black", size=11, face="bold"),
    strip.background = element_blank(),
    plot.margin = margin(0.2,0.5,0.2,0.2,"cm"),
    panel.spacing = unit(0.75,"cm")
  )

ggsave(filename="filepath/emmeans_sa1.jpg", width=15, height=10, units="cm", dpi=300)

### SA2: keep only remote assessments

dat_sa2 <- dat[dat$remote_collection==1,]
svy_obj_sa2 <- svydesign(ids=~participant_id, strata=~cis20_samp, weights=~svywgt, data=dat_sa2)

mod_sa2_any <- svyglm(lc_any ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))) +
                        ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                        age10_visit0 + sex + white + gor9d +
                        imd_quintile + health_status_visit0,
                      family=binomial, design=svy_obj_sa2)

vcov_sa2_any <- vcov(mod_sa2_any)

coeff_sa2_any <- as.data.frame(summary(mod_sa2_any)$coefficients)
write.csv(coeff_sa2_any, file="filepath/coeff_sa2_any.csv")

lr_sa2_any <- as.data.frame(Anova(mod_sa2_any, test.statistic="Chisq"))
write.csv(lr_sa2_any, file="filepath/lr_sa2_any.csv")

mod_sa2_lim <- svyglm(lc_lim ~ ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))) +
                        ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                        age10_visit0 + sex + white + gor9d +
                        imd_quintile + health_status_visit0,
                      family=binomial, design=svy_obj_sa2)

vcov_sa2_lim <- vcov(mod_sa2_lim)

coeff_sa2_lim <- as.data.frame(summary(mod_sa2_lim)$coefficients)
write.csv(coeff_sa2_lim, file="filepath/coeff_sa2_lim.csv")

lr_sa2_lim <- as.data.frame(Anova(mod_sa2_lim, test.statistic="Chisq"))
write.csv(lr_sa2_lim, file="filepath/lr_sa2_lim.csv")

emmeans_sa2_any <- emmeans.collate(mod_sa2_any, vcov_sa2_any, emmeans_wgts$n, mid_study)
emmeans_sa2_any <- emmeans_sa2_any[emmeans_sa2_any$time_since_infection >= min(dat$time_since_infection[dat$remote_collection==1]),]
write.csv(emmeans_sa2_any, file="filepath/emmeans_sa2_any.csv", row.names=FALSE)

emmeans_sa2_lim <- emmeans.collate(mod_sa2_lim, vcov_sa2_lim, emmeans_wgts$n, mid_study)
emmeans_sa2_lim <- emmeans_sa2_lim[emmeans_sa2_lim$time_since_infection >= min(dat$time_since_infection[dat$remote_collection==1]),]
write.csv(emmeans_sa2_lim, file="filepath/emmeans_sa2_lim.csv", row.names=FALSE)

emmeans_any$analysis <- "Main analysis"
emmeans_any$outcome <- "Any severity"

emmeans_lim$analysis <- "Main analysis"
emmeans_lim$outcome <- "Activity limiting"

emmeans_sa2_any$analysis <- "Remote study assessments only"
emmeans_sa2_any$outcome <- "Any severity"

emmeans_sa2_lim$analysis <- "Remote study assessments only"
emmeans_sa2_lim$outcome <- "Activity limiting"

emmeans_sa2 <- rbind(emmeans_any, emmeans_lim, emmeans_sa2_any, emmeans_sa2_lim)
emmeans_sa2$analysis <- factor(emmeans_sa2$analysis, levels=c("Remote study assessments only", "Main analysis"))
emmeans_sa2$outcome <- factor(emmeans_sa2$outcome, levels=c("Any severity", "Activity limiting"))
emmeans_sa2$time <- emmeans_sa2$time_since_infection / (365.25/12)
emmeans_sa2$prevalence <- emmeans_sa2$prob * 100
emmeans_sa2$lcl <- emmeans_sa2$asymp.LCL * 100
emmeans_sa2$ucl <- emmeans_sa2$asymp.UCL * 100

ggplot(emmeans_sa2, aes(x=time, y=prevalence, group=analysis, colour=analysis, fill=analysis)) +
  geom_line(linewidth=0.8) +
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2, linewidth=NA) +
  facet_wrap(~outcome) +
  xlab("Time since infection (months)") +
  ylab("Prevalence (percent)") +
  scale_x_continuous(breaks=seq(0,30,6), limits=c(0,30), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,NA), expand=expansion(mult=c(0,0.05))) +
  scale_colour_manual(values=c("#206095","#f39431")) +
  scale_fill_manual(values=c("#206095","#f39431")) +
  theme(
    axis.text = element_text(colour="black", size=11),
    axis.title = element_text(colour="black", size=11),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(colour="grey", linewidth=0.2),
    strip.text = element_text(colour="black", size=11, face="bold"),
    strip.background = element_blank(),
    plot.margin = margin(0.2,0.5,0.2,0.2,"cm"),
    panel.spacing = unit(0.75,"cm")
  )

ggsave(filename="filepath/emmeans_sa2.jpg", width=15, height=10, units="cm", dpi=300)
