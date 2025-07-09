library(sparklyr)
library(tidyverse)
library(dbplyr)
library(dplyr)
library(splines)
library(survey)
library(car)
library(emmeans)
library(ggplot2)

source("filepath/emmeans_collate.r")

############################### PREPARE DATA ###############################

### set up the Spark connection
xl_config <- spark_config()
xl_config$spark.executor.memory <- "20g"
xl_config$spark.yarn.executor.memoryOverhead <- "2g"
xl_config$spark.executor.cores <- 5
xl_config$spark.dynamicAllocation.enabled <- "true"
xl_config$spark.dynamicAllocation.maxExecutors <- 12
xl_config$spark.sql.shuffle.partitions <- 240
#xl_config$spark.shuffle.service.enabled <- "true"
xl_config$spark.sql.legacy.timeParserPolicy = "LEGACY"
#xl_config$spark.sql.session.timeZone = "GMT"
xl_config$spark.sql.session.timeZone = "UTC+01:00"
#xl_config$spark.sql.parquet.datetimeRebaseModeInRead = "CORRECTED"
xl_config$spark.sql.parquet.int96RebaseModeInRead = "CORRECTED"
xl_config$spark.sql.analyzer.maxIterations = 1000

sc <- spark_connect(master = "yarn-client",
                    app_name = "ONS_session",
                    config = xl_config)

### read in study dataset
load("filepath/lc_trajectories_dataset.RData")

### read in symptoms from full visit-level CIS dataset
sympt <- spark_read_csv(sc, name="data_participant_clean", path="filepath/data_participant_clean_20230331.csv", 
                        header=TRUE, infer_schema=TRUE)

vars_of_interest <- c(
  "visit_id",
  "long_covid_fever",
  "long_covid_weakness_tiredness",
  "long_covid_diarrhoea",
  "long_covid_loss_of_smell",
  "long_covid_shortness_of_breath",
  "long_covid_vertigo_dizziness",
  "long_covid_trouble_sleeping",
  "long_covid_headache",
  "long_covid_nausea_vomiting",
  "long_covid_loss_of_appetite",
  "long_covid_sore_throat",
  "long_covid_chest_pain",
  "long_covid_worry_anxiety",
  "long_covid_memory_loss_confusion",
  "long_covid_muscle_ache",
  "long_covid_abdominal_pain",
  "long_covid_loss_of_taste",
  "long_covid_cough",
  "long_covid_palpitations",
  "long_covid_low_mood_not_enjoying",
  "long_covid_difficult_concentrate"
)

sympt <- sympt %>%
  select(all_of(vars_of_interest))

### keep symptoms for visits in study dataset
ids <- dat[,"visit_id",drop=FALSE]
ids <- copy_to(sc, ids, overwrite=TRUE)

ids <- ids %>%
  left_join(sympt, by=join_by(visit_id==visit_id)) %>%
  collect()

### join symptoms to study dataset
dat <- merge(x=dat, y=ids, by="visit_id", all.x=TRUE, all.y=FALSE)

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
### who ever reported having LC during follow-up
dat_ever_lc <- dat %>%
  filter(!is.na(lc_any) &
         time_since_infection >= 84 &
         time_since_infection <= 365.25*2.5) %>%
  group_by(participant_id) %>%
  summarise(ever_lc=max(lc_any))

dat <- merge(x=dat, y=dat_ever_lc, by="participant_id", all.x=TRUE, all.y=TRUE)

emmeans_wgts <- dat %>%
  filter(!is.na(ever_lc) & ever_lc==1) %>%
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

### restrict dataset to visits when participants reported that they had LC
dat <- dat[!is.na(dat$lc_any) & dat$lc_any==1,]

### list LC symptoms to loop over
sympts <- c(
  "long_covid_fever",
  "long_covid_weakness_tiredness",
  "long_covid_diarrhoea",
  "long_covid_loss_of_smell",
  "long_covid_shortness_of_breath",
  "long_covid_vertigo_dizziness",
  "long_covid_trouble_sleeping",
  "long_covid_headache",
  "long_covid_nausea_vomiting",
  "long_covid_loss_of_appetite",
  "long_covid_sore_throat",
  "long_covid_chest_pain",
  "long_covid_worry_anxiety",
  "long_covid_memory_loss_confusion",
  "long_covid_muscle_ache",
  "long_covid_abdominal_pain",
  "long_covid_loss_of_taste",
  "long_covid_cough",
  "long_covid_palpitations",
  "long_covid_low_mood_not_enjoying",
  "long_covid_difficult_concentrate"
)

### set up empty lists to store outputs
lr_list <- as.list(NULL)
emmeans_list <- as.list(NULL)
pairs_list <- as.list(NULL)

### loop over symptoms
for(i in 1:length(sympts)) {

  ### select symptom of interest
  sympt <- sympts[i]
  
  ### take a copy of the dataset
  dat2 <- dat
  
  ### clean the symptom variable
  dat2[[sympt]] <- ifelse(is.na(dat2[[sympt]]), 0, dat2[[sympt]])
  dat2[[sympt]] <- ifelse(dat2$lc_any==0 & dat2[[sympt]]==1, 0, dat2[[sympt]])
  
  ### set up survey design object
  svy_obj <- svydesign(ids=~participant_id, strata=~cis20_samp, weights=~svywgt, data=dat2)
  
  ### construct the model formula
  covs <- "ns(time_since_infection, df=4, Boundary.knots=quantile(time_since_infection, probs=c(0,1))) +
           ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
           age10_visit0 + sex + white + gor9d + imd_quintile + health_status_visit0"
  
  mod_formula <- as.formula(paste0(sympt, "~", covs))
  
  ### fit model
  mod_any <- svyglm(mod_formula, family=binomial, design=svy_obj)
  vcov_any <- vcov(mod_any)
  
  ### LR test of time trajectory
  lr_any <- as.data.frame(Anova(mod_any, test.statistic="Chisq"))
  lr_any <- lr_any[1,]
  lr_any$Symptom <- sympt
  lr_any <- lr_any[,c(4,1:3)]
  rownames(lr_any) <- NULL
  lr_list[[i]] <- lr_any
  
  ### marginal probabilities
  emmeans_any <- emmeans.collate(mod_any, vcov_any, emmeans_wgts$n, mid_study)
  emmeans_any <- emmeans_any[emmeans_any$time_since_infection >= min(dat2$time_since_infection) &
                             emmeans_any$time_since_infection <= max(dat2$time_since_infection),]
  
  emmeans_any$symptom <- sympt
  emmeans_any$time <- emmeans_any$time_since_infection / (365.25/12)
  emmeans_any$prevalence <- emmeans_any$prob * 100
  emmeans_any$lcl <- emmeans_any$asymp.LCL * 100
  emmeans_any$ucl <- emmeans_any$asymp.UCL * 100
  emmeans_any <- emmeans_any[,7:11]
  emmeans_list[[i]] <- emmeans_any
  
  ### pairwise comparison of LC prevalence at 12 weeks vs. 2.5 years
  pairs_any <- as.data.frame(pairs(emmeans(
    object = mod_any,
    vcov. = vcov_any,
    specs = ~time_since_infection,
    weights = emmeans_wgts$n,
    at = list(time_since_infection = c(max(dat2$time_since_infection), min(dat2$time_since_infection)),
              infection_day = mid_study),
    type = "response",
    rg.limit = 10000
  )))
  
  pairs_any$symptom <- sympt
  pairs_any <- pairs_any[,c(8,1:3,6:7)]
  
}

### collate outputs
lr_out <- lr_list[[1]]
emmeans_out <- emmeans_list[[1]]
pairs_out <- pairs_list[[1]]

if(length(sympts)>1) {
  for(i in 2:length(sympts)) {
    lr_out <- rbind(lr_out, lr_list[[i]])
    emmeans_out <- rbind(emmeans_out, emmeans_list[[i]])
    pairs_out <- rbind(pairs_out, pairs_list[[i]])
  }
}

write.csv(lr_out, file="filepath/lr_symptoms_denom_lc.csv", row.names=FALSE)
write.csv(emmeans_out, file="filepath/emmeans_symptoms_denom_lc.csv", row.names=FALSE)
write.csv(pairs_out, file="filepath/pairs_symptoms_denom_lc.csv", row.names=FALSE)

### plot heatmap
emmeans_out$symptom[emmeans_out$symptom=="long_covid_fever"] <- "Fever"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_weakness_tiredness"] <- "Weakness/tiredness"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_diarrhoea"] <- "Diarrhoea"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_loss_of_smell"] <- "Loss of smell"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_shortness_of_breath"] <- "Shortness of breath"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_vertigo_dizziness"] <- "Vertigo/dizziness"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_trouble_sleeping"] <- "Trouble sleeping"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_headache"] <- "Headache"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_nausea_vomiting"] <- "Nausea/vomiting"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_loss_of_appetite"] <- "Loss of appetite"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_sore_throat"] <- "Sore throat"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_chest_pain"] <- "Chest pain"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_worry_anxiety"] <- "Worry/anxiety"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_memory_loss_confusion"] <- "Memory loss or confusion"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_muscle_ache"] <- "Muscle ache"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_abdominal_pain"] <- "Abdominal pain"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_loss_of_taste"] <- "Loss of taste"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_cough"] <- "Cough"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_palpitations"] <- "Palpitations"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_low_mood_not_enjoying"] <- "Low mood/not enjoying anything"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_difficult_concentrate"] <- "Difficulty concentrating"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_runny_nose_sneezing"] <- "Runny nose or sneezing"
emmeans_out$symptom[emmeans_out$symptom=="long_covid_noisy_breathing"] <- "Noisy breathing"

prevalence_18w <- emmeans_out$prevalence[emmeans_out$time==max(emmeans_out$time)]
names(prevalence_18w) <- emmeans_out$symptom[emmeans_out$time==max(emmeans_out$time)]
sympt_ordered <- names(sort(prevalence_18w))
emmeans_out$symptom <- factor(emmeans_out$symptom, levels=sympt_ordered)

ggplot(emmeans_out, aes(x=time, y=symptom, fill=prevalence)) +
  geom_tile() +
  scale_fill_distiller(palette="Blues", direction=1) +
  scale_x_continuous(breaks=seq(6,30,6), expand=c(0,0)) +
  xlab("Time since infection (months)") +
  theme(
    axis.title.x=element_text(size=11, colour="black"),
    axis.title.y=element_blank(),
    axis.text = element_text(size=11, colour="black"),
    axis.ticks.x=element_line(size=0.5, colour="black"),
    axis.ticks.y=element_blank(),
    axis.line=element_blank(),
    panel.border=element_rect(size=0.5, colour="black", fill=NA),
    panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    legend.direction="horizontal",
    legend.justification="center",
    legend.text=element_text(size=11, colour="black", face="plain")
  )

ggsave(filename="filepath/emmeans_symptoms_denom_lc.jpg", width=15, height=15, units="cm", dpi=300)
