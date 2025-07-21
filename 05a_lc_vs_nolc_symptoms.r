library(sparklyr)
library(tidyverse)
library(dbplyr)
library(dplyr)
library(splines)
library(survey)
library(ggplot2)

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
  "sympt_now_fever",
  "sympt_now_muscle_ache_myalgia",
  "sympt_now_fatigue_weakness",
  "sympt_now_sore_throat",
  "sympt_now_cough",
  "sympt_now_shortness_of_breath",
  "sympt_now_headache",
  "sympt_now_nausea_vomiting",
  "sympt_now_abdominal_pain",
  "sympt_now_diarrhoea",
  "sympt_now_loss_of_taste",
  "sympt_now_loss_of_smell",
  "sympt_now_more_trouble_sleeping",
  "sympt_now_loss_of_appetite_or_e",
  "sympt_now_runny_nose_sneezing",
  "sympt_now_noisy_breathing_wheez",
  "sympt_now_chest_pain",
  "sympt_now_palpitations",
  "sympt_now_vertigo_dizziness",
  "sympt_now_worry_anxiety",
  "sympt_now_low_mood_not_enjoying",
  "sympt_now_memory_loss_confusion",
  "sympt_now_difficulty_concentrate"
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

### restrict dataset to visits >= 12 weeks after infection
dat <- dat[dat$time_since_infection >= 84,]

### restrict dataset to visit <= 2.5 years after infection
dat <- dat[dat$time_since_infection <= 365.25*2.5,]

### restrict dataset to visits when participants reported their LC status
dat <- dat[!is.na(dat$lc_any),]

### binary flag indicating whether participant ever reported LC during follow-up
dat_ever_lc <- dat %>%
  group_by(participant_id) %>%
  summarise(ever_lc=max(lc_any))

dat <- merge(x=dat, y=dat_ever_lc, by="participant_id", all.x=TRUE, all.y=TRUE)

### for participants who ever reported having LC, keep only visits when they had LC
dat <- dat[dat$ever_lc==0 | dat$lc_any==1,]

############################### MODELLING ###############################

### list symptoms to loop over
sympts <- c(
  "sympt_now_fever",
  "sympt_now_muscle_ache_myalgia",
  "sympt_now_fatigue_weakness",
  "sympt_now_sore_throat",
  "sympt_now_cough",
  "sympt_now_shortness_of_breath",
  "sympt_now_headache",
  "sympt_now_nausea_vomiting",
  "sympt_now_abdominal_pain",
  "sympt_now_diarrhoea",
  "sympt_now_loss_of_taste",
  "sympt_now_loss_of_smell",
  "sympt_now_more_trouble_sleeping",
  "sympt_now_loss_of_appetite_or_e",
  "sympt_now_runny_nose_sneezing",
  "sympt_now_noisy_breathing_wheez",
  "sympt_now_chest_pain",
  "sympt_now_palpitations",
  "sympt_now_vertigo_dizziness",
  "sympt_now_worry_anxiety",
  "sympt_now_low_mood_not_enjoying",
  "sympt_now_memory_loss_confusion",
  "sympt_now_difficulty_concentrate"
)

### set up empty list to store outputs
out_list <- as.list(NULL)

### loop over symptoms
for(i in 1:length(sympts)) {
  
  ### take a copy of the dataset
  dat2 <- dat
  
  ### select symptom of interest
  sympt <- sympts[i]
  dat2$sympt <- dat2[[sympt]]
  
  ### clean the symptom variable
  dat2$sympt <- ifelse(is.na(dat2$sympt), 0, dat2$sympt)
  
  ### drop visits before symptom question was introduced
  if(sympt %in% c("sympt_now_more_trouble_sleeping",
                  "sympt_now_loss_of_appetite_or_e",
                  "sympt_now_runny_nose_sneezing",
                  "sympt_now_noisy_breathing_wheez")) {
    dat2 <- dat2[dat2$visit_date >= "2021-10-01",]
  }
  
  if(sympt %in% c("sympt_now_chest_pain",
                  "sympt_now_palpitations",
                  "sympt_now_vertigo_dizziness",
                  "sympt_now_worry_anxiety",
                  "sympt_now_low_mood_not_enjoying",
                  "sympt_now_memory_loss_confusion",
                  "sympt_now_difficulty_concentrate")) {
    dat2 <- dat2[dat2$visit_date >= "2022-02-01",]
  }
  
  ### aggregate dataset to participant-level
  dat2 <- dat2 %>%
    group_by(participant_id, cis20_samp, svywgt, ever_lc, infection_day,
             age10_visit0, sex, white, gor9d, imd_quintile, health_status_visit0) %>%
    summarise(n_visits=n(), n_sympt=sum(sympt), .groups="drop")
  
  ### set up survey design object
  svy_obj <- svydesign(ids=~participant_id, strata=~cis20_samp, weights=~svywgt, data=dat2)
  
  ### descriptive stats
  desc <- dat2 %>%
    group_by(ever_lc) %>%
    summarise(n_participants = n(),
              person_symptoms = sum(n_sympt),
              person_visits = sum(n_visits)) %>%
    mutate(rate = person_symptoms / person_visits * 1000)
  
  ### unadjusted IRRs
  mod_unadj <- svyglm(n_sympt ~ ever_lc + offset(log(n_visits)),
                      family=poisson, design=svy_obj)
  
  coeffs_unadj <- as.data.frame(summary(mod_unadj)$coefficients)
  coeffs_unadj <- coeffs_unadj["ever_lc",,drop=FALSE]
  rownames(coeffs_unadj) <- NULL
  coeffs_unadj$irr_unadj <- exp(coeffs_unadj$Estimate)
  coeffs_unadj$irr_unadj_lcl <- exp(coeffs_unadj$Estimate - 1.96 * coeffs_unadj$`Std. Error`)
  coeffs_unadj$irr_unadj_ucl <- exp(coeffs_unadj$Estimate + 1.96 * coeffs_unadj$`Std. Error`)
  coeffs_unadj$irr_unadj_pvalue <- coeffs_unadj$`Pr(>|t|)`
  coeffs_unadj <- coeffs_unadj[,c("irr_unadj", "irr_unadj_lcl", "irr_unadj_ucl", "irr_unadj_pvalue")]
  coeffs_unadj <- rbind(NA, coeffs_unadj)
  
  ### adjusted IRRs
  mod_adj <- svyglm(n_sympt ~ ever_lc +
                      ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                      age10_visit0 + sex + white + gor9d + imd_quintile + health_status_visit0 +
                      offset(log(n_visits)),
                    family=poisson, design=svy_obj)
  
  coeffs_adj <- as.data.frame(summary(mod_adj)$coefficients)
  coeffs_adj <- coeffs_adj["ever_lc",,drop=FALSE]
  rownames(coeffs_adj) <- NULL
  coeffs_adj$irr_adj <- exp(coeffs_adj$Estimate)
  coeffs_adj$irr_adj_lcl <- exp(coeffs_adj$Estimate - 1.96 * coeffs_adj$`Std. Error`)
  coeffs_adj$irr_adj_ucl <- exp(coeffs_adj$Estimate + 1.96 * coeffs_adj$`Std. Error`)
  coeffs_adj$irr_adj_pvalue <- coeffs_adj$`Pr(>|t|)`
  coeffs_adj <- coeffs_adj[,c("irr_adj", "irr_adj_lcl", "irr_adj_ucl", "irr_adj_pvalue")]
  coeffs_adj <- rbind(NA, coeffs_adj)
  
  ### combine outputs
  comb <- cbind(desc, coeffs_unadj, coeffs_adj)
  comb$symptom <- sympt
  comb <- comb[,c(ncol(comb),1:(ncol(comb)-1))]
  out_list[[i]] <- comb
  
}

### collate outputs
out_df <- out_list[[1]]

if(length(sympts)>1) {
  for(i in 2:length(sympts)) {out_df <- rbind(out_df, out_list[[i]])}
}

write.csv(out_df, file="filepath/irr_symptoms.csv", row.names=FALSE)

### prepare data for plotting
plot_df <- out_df[out_df$ever_lc==1, c("symptom","irr_adj","irr_adj_lcl","irr_adj_ucl")]

plot_df$symptom[plot_df$symptom=="sympt_now_fever"] <- "Fever"
plot_df$symptom[plot_df$symptom=="sympt_now_muscle_ache_myalgia"] <- "Muscle ache"
plot_df$symptom[plot_df$symptom=="sympt_now_fatigue_weakness"] <- "Weakness/tiredness"
plot_df$symptom[plot_df$symptom=="sympt_now_sore_throat"] <- "Sore throat"
plot_df$symptom[plot_df$symptom=="sympt_now_cough"] <- "Cough"
plot_df$symptom[plot_df$symptom=="sympt_now_shortness_of_breath"] <- "Shortness of breath"
plot_df$symptom[plot_df$symptom=="sympt_now_headache"] <- "Headache"
plot_df$symptom[plot_df$symptom=="sympt_now_nausea_vomiting"] <- "Nausea/vomiting"
plot_df$symptom[plot_df$symptom=="sympt_now_abdominal_pain"] <- "Abdominal pain"
plot_df$symptom[plot_df$symptom=="sympt_now_diarrhoea"] <- "Diarrhoea"
plot_df$symptom[plot_df$symptom=="sympt_now_loss_of_taste"] <- "Loss of taste"
plot_df$symptom[plot_df$symptom=="sympt_now_loss_of_smell"] <- "Loss of smell"
plot_df$symptom[plot_df$symptom=="sympt_now_more_trouble_sleeping"] <- "Trouble sleeping"
plot_df$symptom[plot_df$symptom=="sympt_now_loss_of_appetite_or_e"] <- "Loss of appetite"
plot_df$symptom[plot_df$symptom=="sympt_now_runny_nose_sneezing"] <- "Runny nose or sneezing"
plot_df$symptom[plot_df$symptom=="sympt_now_noisy_breathing_wheez"] <- "Noisy breathing"
plot_df$symptom[plot_df$symptom=="sympt_now_chest_pain"] <- "Chest pain"
plot_df$symptom[plot_df$symptom=="sympt_now_palpitations"] <- "Palpitations"
plot_df$symptom[plot_df$symptom=="sympt_now_vertigo_dizziness"] <- "Vertigo/dizziness"
plot_df$symptom[plot_df$symptom=="sympt_now_worry_anxiety"] <- "Worry/anxiety"
plot_df$symptom[plot_df$symptom=="sympt_now_low_mood_not_enjoying"] <- "Low mood/not enjoying anything"
plot_df$symptom[plot_df$symptom=="sympt_now_memory_loss_confusion"] <- "Memory loss or confusion"
plot_df$symptom[plot_df$symptom=="sympt_now_difficulty_concentrate"] <- "Difficulty concentrating"

plot_df <- plot_df[order(plot_df$irr_adj),]
plot_df$symptom <- factor(plot_df$symptom, levels=plot_df$symptom)

### plot IRRs
ggplot(data=plot_df, aes(x=irr_adj, y=symptom)) +
  geom_vline(xintercept=1, linewidth=0.2, linetype="longdash", colour="grey30") +
  geom_point(size=0.8) +
  geom_errorbar(aes(xmin=irr_adj_lcl, xmax=irr_adj_ucl), width=0) +
  scale_x_continuous(breaks=c(1,2,4,8,16)) +
  coord_trans(x="log") +
  xlab("Incident rate ratio (log scale)") +
  ylab("") +
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
    panel.grid.minor=element_blank()
  )

ggsave(filename="filepath/irr_symptoms.jpg", width=15, height=15, units="cm", dpi=300)

############################### SENSITIVITY ANALYSES ###############################

### SA1: exclude participants who were infected before 11/11/2020 (12 weeks before the
### LC question was implemented on the CIS), as they may have previously had LC but
### recovered before the first time they answered the survey question, hence they would
### be misclassified as being in the non-LC exposure group in the main analysis

### set up empty list to store outputs
out_list_sa1 <- as.list(NULL)

### loop over symptoms
for(i in 1:length(sympts)) {
  
  ### take a copy of the dataset
  dat2 <- dat
  
  ### filter to infections from 11/11/2020
  dat2 <- dat2[dat2$infection_date >= '2020-11-11',]
  
  ### select symptom of interest
  sympt <- sympts[i]
  dat2$sympt <- dat2[[sympt]]
  
  ### clean the symptom variable
  dat2$sympt <- ifelse(is.na(dat2$sympt), 0, dat2$sympt)
  
  ### drop visits before symptom question was introduced
  if(sympt %in% c("sympt_now_more_trouble_sleeping",
                  "sympt_now_loss_of_appetite_or_e",
                  "sympt_now_runny_nose_sneezing",
                  "sympt_now_noisy_breathing_wheez")) {
    dat2 <- dat2[dat2$visit_date >= "2021-10-01",]
  }
  
  if(sympt %in% c("sympt_now_chest_pain",
                  "sympt_now_palpitations",
                  "sympt_now_vertigo_dizziness",
                  "sympt_now_worry_anxiety",
                  "sympt_now_low_mood_not_enjoying",
                  "sympt_now_memory_loss_confusion",
                  "sympt_now_difficulty_concentrate")) {
    dat2 <- dat2[dat2$visit_date >= "2022-02-01",]
  }
  
  ### aggregate dataset to participant-level
  dat2 <- dat2 %>%
    group_by(participant_id, cis20_samp, svywgt, ever_lc, infection_day,
             age10_visit0, sex, white, gor9d, imd_quintile, health_status_visit0) %>%
    summarise(n_visits=n(), n_sympt=sum(sympt), .groups="drop")
  
  ### set up survey design object
  svy_obj <- svydesign(ids=~participant_id, strata=~cis20_samp, weights=~svywgt, data=dat2)
  
  ### descriptive stats
  desc <- dat2 %>%
    group_by(ever_lc) %>%
    summarise(n_participants = n(),
              person_symptoms = sum(n_sympt),
              person_visits = sum(n_visits)) %>%
    mutate(rate = person_symptoms / person_visits * 1000)
  
  ### unadjusted IRRs
  mod_unadj <- svyglm(n_sympt ~ ever_lc + offset(log(n_visits)),
                      family=poisson, design=svy_obj)
  
  coeffs_unadj <- as.data.frame(summary(mod_unadj)$coefficients)
  coeffs_unadj <- coeffs_unadj["ever_lc",,drop=FALSE]
  rownames(coeffs_unadj) <- NULL
  coeffs_unadj$irr_unadj <- exp(coeffs_unadj$Estimate)
  coeffs_unadj$irr_unadj_lcl <- exp(coeffs_unadj$Estimate - 1.96 * coeffs_unadj$`Std. Error`)
  coeffs_unadj$irr_unadj_ucl <- exp(coeffs_unadj$Estimate + 1.96 * coeffs_unadj$`Std. Error`)
  coeffs_unadj$irr_unadj_pvalue <- coeffs_unadj$`Pr(>|t|)`
  coeffs_unadj <- coeffs_unadj[,c("irr_unadj", "irr_unadj_lcl", "irr_unadj_ucl", "irr_unadj_pvalue")]
  coeffs_unadj <- rbind(NA, coeffs_unadj)
  
  ### adjusted IRRs
  mod_adj <- svyglm(n_sympt ~ ever_lc +
                      ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                      age10_visit0 + sex + white + gor9d + imd_quintile + health_status_visit0 +
                      offset(log(n_visits)),
                    family=poisson, design=svy_obj)
  
  coeffs_adj <- as.data.frame(summary(mod_adj)$coefficients)
  coeffs_adj <- coeffs_adj["ever_lc",,drop=FALSE]
  rownames(coeffs_adj) <- NULL
  coeffs_adj$irr_adj <- exp(coeffs_adj$Estimate)
  coeffs_adj$irr_adj_lcl <- exp(coeffs_adj$Estimate - 1.96 * coeffs_adj$`Std. Error`)
  coeffs_adj$irr_adj_ucl <- exp(coeffs_adj$Estimate + 1.96 * coeffs_adj$`Std. Error`)
  coeffs_adj$irr_adj_pvalue <- coeffs_adj$`Pr(>|t|)`
  coeffs_adj <- coeffs_adj[,c("irr_adj", "irr_adj_lcl", "irr_adj_ucl", "irr_adj_pvalue")]
  coeffs_adj <- rbind(NA, coeffs_adj)
  
  ### combine outputs
  comb <- cbind(desc, coeffs_unadj, coeffs_adj)
  comb$symptom <- sympt
  comb <- comb[,c(ncol(comb),1:(ncol(comb)-1))]
  out_list_sa1[[i]] <- comb
  
}

### collate outputs
out_df_sa1 <- out_list_sa1[[1]]

if(length(sympts)>1) {
  for(i in 2:length(sympts)) {out_df_sa1 <- rbind(out_df_sa1, out_list_sa1[[i]])}
}

write.csv(out_df_sa1, file="filepath/irr_symptoms_sa1.csv", row.names=FALSE)

### prepare data for plotting
plot_df_sa1 <- out_df_sa1[out_df_sa1$ever_lc==1, c("symptom","irr_adj","irr_adj_lcl","irr_adj_ucl")]

plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_fever"] <- "Fever"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_muscle_ache_myalgia"] <- "Muscle ache"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_fatigue_weakness"] <- "Weakness/tiredness"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_sore_throat"] <- "Sore throat"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_cough"] <- "Cough"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_shortness_of_breath"] <- "Shortness of breath"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_headache"] <- "Headache"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_nausea_vomiting"] <- "Nausea/vomiting"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_abdominal_pain"] <- "Abdominal pain"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_diarrhoea"] <- "Diarrhoea"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_loss_of_taste"] <- "Loss of taste"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_loss_of_smell"] <- "Loss of smell"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_more_trouble_sleeping"] <- "Trouble sleeping"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_loss_of_appetite_or_e"] <- "Loss of appetite"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_runny_nose_sneezing"] <- "Runny nose or sneezing"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_noisy_breathing_wheez"] <- "Noisy breathing"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_chest_pain"] <- "Chest pain"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_palpitations"] <- "Palpitations"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_vertigo_dizziness"] <- "Vertigo/dizziness"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_worry_anxiety"] <- "Worry/anxiety"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_low_mood_not_enjoying"] <- "Low mood/not enjoying anything"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_memory_loss_confusion"] <- "Memory loss or confusion"
plot_df_sa1$symptom[plot_df_sa1$symptom=="sympt_now_difficulty_concentrate"] <- "Difficulty concentrating"

plot_df$analysis <- "Main analysis"
plot_df_sa1$analysis <- "Infected from 11 November 2020"

plot_df_sa1$symptom <- factor(plot_df_sa1$symptom, levels=plot_df$symptom)

plot_df_sa1 <- rbind(plot_df, plot_df_sa1)

plot_df_sa1$analysis <- as.factor(plot_df_sa1$analysis)
plot_df_sa1$analysis <- relevel(plot_df_sa1$analysis, ref="Main analysis")

### plot IRRs
ggplot(data=plot_df_sa1, aes(x=irr_adj, y=symptom, group=analysis, colour=analysis, fill=analysis)) +
  geom_vline(xintercept=1, linewidth=0.2, linetype="longdash", colour="grey30") +
  geom_point(size=0.8, position=position_dodge(0.5)) +
  geom_errorbar(aes(xmin=irr_adj_lcl, xmax=irr_adj_ucl), width=0, position=position_dodge(0.5)) +
  scale_x_continuous(breaks=c(1,2,4,8,16)) +
  scale_colour_manual(values=c("#206095","#f39431")) +
  coord_trans(x="log") +
  xlab("Incident rate ratio (log scale)") +
  ylab("") +
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
    legend.direction="vertical",
    legend.justification="center",
    legend.text=element_text(size=11, colour="black", face="plain")
  )

ggsave(filename="filepath/irr_symptoms_sa1.jpg", width=15, height=15, units="cm", dpi=300)
