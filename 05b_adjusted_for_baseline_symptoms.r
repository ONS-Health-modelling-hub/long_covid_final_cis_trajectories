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
  "sympt_now_loss_of_smell"
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

### drop participants who didn't responded to the survey 30-90 days before infection
dat$baseline_visit <- ifelse(as.numeric(dat$infection_date - dat$visit_date >= 30) &
                             as.numeric(dat$infection_date - dat$visit_date <= 90), 1, 0)

dat_baseline <- dat %>%
  group_by(participant_id) %>%
  summarise(baseline_response=max(baseline_visit))

dat <- merge(x=dat, y=dat_baseline, by="participant_id", all.x=TRUE, all.y=TRUE)

dat <- dat[dat$baseline_response==1,]

### derive indicators for symptoms reported 30-90 days before infection
dat$sympt_now_fever_baseline <- ifelse(!is.na(dat$sympt_now_fever) & dat$sympt_now_fever==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_muscle_ache_myalgia_baseline <- ifelse(!is.na(dat$sympt_now_muscle_ache_myalgia) & dat$sympt_now_muscle_ache_myalgia==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_fatigue_weakness_baseline <- ifelse(!is.na(dat$sympt_now_fatigue_weakness) & dat$sympt_now_fatigue_weakness==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_sore_throat_baseline <- ifelse(!is.na(dat$sympt_now_sore_throat) & dat$sympt_now_sore_throat==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_cough_baseline <- ifelse(!is.na(dat$sympt_now_cough) & dat$sympt_now_cough==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_shortness_of_breath_baseline <- ifelse(!is.na(dat$sympt_now_shortness_of_breath) & dat$sympt_now_shortness_of_breath==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_headache_baseline <- ifelse(!is.na(dat$sympt_now_headache) & dat$sympt_now_headache==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_nausea_vomiting_baseline <- ifelse(!is.na(dat$sympt_now_nausea_vomiting) & dat$sympt_now_nausea_vomiting==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_abdominal_pain_baseline <- ifelse(!is.na(dat$sympt_now_abdominal_pain) & dat$sympt_now_abdominal_pain==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_diarrhoea_baseline <- ifelse(!is.na(dat$sympt_now_diarrhoea) & dat$sympt_now_diarrhoea==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_loss_of_taste_baseline <- ifelse(!is.na(dat$sympt_now_loss_of_taste) & dat$sympt_now_loss_of_taste==1 & dat$baseline_visit==1, 1, 0)
dat$sympt_now_loss_of_smell_baseline <- ifelse(!is.na(dat$sympt_now_loss_of_smell) & dat$sympt_now_loss_of_smell==1 & dat$baseline_visit==1, 1, 0)

dat_baseline_sympts <- dat %>%
  group_by(participant_id) %>%
  summarise(sympt_now_fever_baseline=max(sympt_now_fever_baseline),
            sympt_now_muscle_ache_myalgia_baseline=max(sympt_now_muscle_ache_myalgia_baseline),
            sympt_now_fatigue_weakness_baseline=max(sympt_now_fatigue_weakness_baseline),
            sympt_now_sore_throat_baseline=max(sympt_now_sore_throat_baseline),
            sympt_now_cough_baseline=max(sympt_now_cough_baseline),
            sympt_now_shortness_of_breath_baseline=max(sympt_now_shortness_of_breath_baseline),
            sympt_now_headache_baseline=max(sympt_now_headache_baseline),
            sympt_now_nausea_vomiting_baseline=max(sympt_now_nausea_vomiting_baseline),
            sympt_now_abdominal_pain_baseline=max(sympt_now_abdominal_pain_baseline),
            sympt_now_diarrhoea_baseline=max(sympt_now_diarrhoea_baseline),
            sympt_now_loss_of_taste_baseline=max(sympt_now_loss_of_taste_baseline),
            sympt_now_loss_of_smell_baseline=max(sympt_now_loss_of_smell_baseline))

dat$sympt_now_fever_baseline <- NULL
dat$sympt_now_muscle_ache_myalgia_baseline <- NULL
dat$sympt_now_fatigue_weakness_baseline <- NULL
dat$sympt_now_sore_throat_baseline <- NULL
dat$sympt_now_cough_baseline <- NULL
dat$sympt_now_shortness_of_breath_baseline <- NULL
dat$sympt_now_headache_baseline <- NULL
dat$sympt_now_nausea_vomiting_baseline <- NULL
dat$sympt_now_abdominal_pain_baseline <- NULL
dat$sympt_now_diarrhoea_baseline <- NULL
dat$sympt_now_loss_of_taste_baseline <- NULL
dat$sympt_now_loss_of_smell_baseline <- NULL

dat <- merge(x=dat, y=dat_baseline_sympts, by="participant_id", all.x=TRUE, all.y=TRUE)

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
  "sympt_now_loss_of_smell"
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
  
  baseline_sympt <- paste0(sympt, "_baseline")
  dat2$baseline_sympt <- dat2[[baseline_sympt]]
  
  ### clean the symptom variable
  dat2$sympt <- ifelse(is.na(dat2$sympt), 0, dat2$sympt)
  
  ### aggregate dataset to participant-level
  dat2 <- dat2 %>%
    group_by(participant_id, cis20_samp, svywgt, ever_lc, infection_day,
             age10_visit0, sex, white, gor9d, imd_quintile, health_status_visit0,
             baseline_sympt) %>%
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
  
  ### IRRs not adjusted for baseline symptoms
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
  
  ### IRRs adjusted for baseline symptoms
  mod_adj2 <- svyglm(n_sympt ~ ever_lc +
                       ns(infection_day, df=4, Boundary.knots=quantile(infection_day, probs=c(0.1,0.9))) +
                       age10_visit0 + sex + white + gor9d + imd_quintile + health_status_visit0 +
                       baseline_sympt +
                       offset(log(n_visits)),
                     family=poisson, design=svy_obj)
  
  coeffs_adj2 <- as.data.frame(summary(mod_adj2)$coefficients)
  coeffs_adj2 <- coeffs_adj2["ever_lc",,drop=FALSE]
  rownames(coeffs_adj2) <- NULL
  coeffs_adj2$irr_adj2 <- exp(coeffs_adj2$Estimate)
  coeffs_adj2$irr_adj2_lcl <- exp(coeffs_adj2$Estimate - 1.96 * coeffs_adj2$`Std. Error`)
  coeffs_adj2$irr_adj2_ucl <- exp(coeffs_adj2$Estimate + 1.96 * coeffs_adj2$`Std. Error`)
  coeffs_adj2$irr_adj2_pvalue <- coeffs_adj2$`Pr(>|t|)`
  coeffs_adj2 <- coeffs_adj2[,c("irr_adj2", "irr_adj2_lcl", "irr_adj2_ucl", "irr_adj2_pvalue")]
  coeffs_adj2 <- rbind(NA, coeffs_adj2)
  
  ### combine outputs
  comb <- cbind(desc, coeffs_adj, coeffs_adj2)
  comb$symptom <- sympt
  comb <- comb[,c(ncol(comb),1:(ncol(comb)-1))]
  out_list[[i]] <- comb
  
}

### collate outputs
out_df <- out_list[[1]]

if(length(sympts)>1) {
  for(i in 2:length(sympts)) {out_df <- rbind(out_df, out_list[[i]])}
}

write.csv(out_df, file="filepath/irr_symptoms_adjusted_for_baseline.csv", row.names=FALSE)

### prepare data for plotting
out_df_adj <- out_df[out_df$ever_lc==1, c("symptom","irr_adj","irr_adj_lcl","irr_adj_ucl")]
out_df_adj$model <- "Not adjusted for baseline symptoms"

out_df_adj$symptom[out_df_adj$symptom=="sympt_now_fever"] <- "Fever"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_muscle_ache_myalgia"] <- "Muscle ache"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_fatigue_weakness"] <- "Weakness/tiredness"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_sore_throat"] <- "Sore throat"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_cough"] <- "Cough"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_shortness_of_breath"] <- "Shortness of breath"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_headache"] <- "Headache"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_nausea_vomiting"] <- "Nausea/vomiting"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_abdominal_pain"] <- "Abdominal pain"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_diarrhoea"] <- "Diarrhoea"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_loss_of_taste"] <- "Loss of taste"
out_df_adj$symptom[out_df_adj$symptom=="sympt_now_loss_of_smell"] <- "Loss of smell"

out_df_adj2 <- out_df[out_df$ever_lc==1, c("symptom","irr_adj2","irr_adj2_lcl","irr_adj2_ucl")]
out_df_adj2$model <- "Adjusted for baseline symptoms"
colnames(out_df_adj2) <- colnames(out_df_adj)

out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_fever"] <- "Fever"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_muscle_ache_myalgia"] <- "Muscle ache"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_fatigue_weakness"] <- "Weakness/tiredness"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_sore_throat"] <- "Sore throat"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_cough"] <- "Cough"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_shortness_of_breath"] <- "Shortness of breath"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_headache"] <- "Headache"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_nausea_vomiting"] <- "Nausea/vomiting"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_abdominal_pain"] <- "Abdominal pain"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_diarrhoea"] <- "Diarrhoea"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_loss_of_taste"] <- "Loss of taste"
out_df_adj2$symptom[out_df_adj2$symptom=="sympt_now_loss_of_smell"] <- "Loss of smell"

plot_df <- rbind(out_df_adj, out_df_adj2)
rownames(plot_df) <- NULL

plot_df$model <- factor(plot_df$model, levels=c("Not adjusted for baseline symptoms",
                                                "Adjusted for baseline symptoms"))

out_df_adj <- out_df_adj[order(out_df_adj$irr_adj),]
plot_df$symptom <- factor(plot_df$symptom, levels=out_df_adj$symptom)

### plot IRRs
ggplot(data=plot_df, aes(x=irr_adj, y=symptom, group=model, colour=model, fill=model)) +
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

ggsave(filename="filepath/irr_symptoms_adjusted_for_baseline.jpg", width=15, height=13, units="cm", dpi=300)
