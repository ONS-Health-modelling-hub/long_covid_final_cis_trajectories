library(sparklyr)
library(tidyverse)
library(dbplyr)

############################### READ IN DATASETS ###############################

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

### read in CIS visit-level dataset
cis_visit <- spark_read_csv(sc, name="data_participant_clean", path="filepath/data_participant_clean_20230331.csv", 
                            header=TRUE, infer_schema=TRUE)

vars_of_interest <- c(
  "participant_id",
  "hh_id_fake",
  "visit_id",
  "visit_number",
  "visit_date",
  "dataset",
  "age_at_visit",
  "sex",
  "ethnicityg",
  "country",
  "gor9d",
  "cis20_samp",
  "imd_samp",
  "health_conditions",
  "health_conditions_impact",
  "covid_admitted",
  "long_covid_have_symptoms",
  "reduce_activities_long_covid"
)

dat <- cis_visit %>%
  select(all_of(vars_of_interest)) %>% 
  mutate(visit_date = to_date(visit_date, "ddMMMyyyy"))

### remove duplicate visits
dat <- dat %>%
  window_order(participant_id, visit_date, desc(visit_id)) %>%
  distinct(participant_id, visit_date, .keep_all=TRUE)

### join vaccination data
cis_vacc <- spark_read_csv(sc, name="data_participant_vaccination", path="filepath/data_participant_vaccination_20230331.csv", 
                           header=TRUE, infer_schema=TRUE)

cis_vacc <- cis_vacc %>%
  select(participant_id, covid_vaccine_date1) %>%
  mutate(covid_vaccine_date1 = to_date(covid_vaccine_date1, "ddMMMyyyy"))

dat <- dat %>%
  left_join(cis_vacc, by=join_by(participant_id==participant_id))

### read in swab positivity episodes dataset (including tests from CIS and T&T)
swab_pos <- spark_read_csv(sc, name="swab_clearance_pos", path="filepath/swab_clearance_pos_all_updt_20230330_new.csv", 
                           header=TRUE, infer_schema=TRUE)

swab_pos <- swab_pos %>%
  mutate(first_pos_date = to_date(first_pos_date, "ddMMMyyyy"),
         first_pos_date_cis = to_date(first_pos_date_cis, "ddMMMyyyy"))

### join index date of first episode, ignoring any T&T results before mass community
### testing was rolled out (likely to reflect high clinical need)
swab_pos_first <- swab_pos %>%
  filter(episode_updt==1) %>%
  mutate(test_source = ifelse(!is.na(first_pos_date_cis) &
                                first_pos_date==first_pos_date_cis, "CIS", "T&T"),
         first_pos_date = ifelse(test_source=="T&T" &
                                   first_pos_date < '2020-09-01', NA, first_pos_date),
         first_pos_date = ifelse(is.na(first_pos_date), first_pos_date_cis, first_pos_date)) %>%
  filter(!is.na(first_pos_date)) %>%
  select(participant_id, first_pos_date, test_source) %>%
  rename(infection_date = first_pos_date)

dat <- dat %>%
  left_join(swab_pos_first, by=join_by(participant_id==participant_id))

### join index date of second episode (ie reinfection)
swab_pos_second <- swab_pos %>%
  filter(episode_updt==2) %>%
  select(participant_id, first_pos_date) %>%
  rename(reinfection_date = first_pos_date)

dat <- dat %>%
  left_join(swab_pos_second, by=join_by(participant_id==participant_id))

### read in Census 2021 populations for weighting
pops <- read.csv(file="filepath/pop2021_uk.csv")
pops <- copy_to(sc, pops, overwrite=TRUE)

############################### EARLIEST DATES ###############################

### find dates of first and last CIS visits
dat_first_last_visit <- dat %>%
  group_by(participant_id) %>%
  summarise(first_visit_date = min(visit_date),
            last_visit_date = max(visit_date))

dat <- dat %>%
  left_join(dat_first_last_visit, by=join_by(participant_id==participant_id))

### set infection date to NA if after last CIS visit date
dat <- dat %>%
  mutate(infection_date = ifelse(!is.na(infection_date) & infection_date > last_visit_date, NA, infection_date),
         reinfection_date = ifelse(!is.na(reinfection_date) & reinfection_date > last_visit_date, NA, reinfection_date))

### set LC responses before 3 Feb 2021 (when the question was implemented) to NA
dat <- dat %>%
  mutate(long_covid_have_symptoms = ifelse(visit_date < "2021-02-03", NA, long_covid_have_symptoms))

### find dates of first and last response to LC question
dat_first_last_lc_response <- dat %>%
  filter(!is.na(long_covid_have_symptoms)) %>%
  group_by(participant_id) %>%
  summarise(first_lc_response_date = min(visit_date),
            last_lc_response_date = max(visit_date))

dat <- dat %>%
  left_join(dat_first_last_lc_response, by=join_by(participant_id==participant_id))

### find dates of first and last positive response to LC question
dat_first_last_lc_yes <- dat %>%
  filter(!is.na(long_covid_have_symptoms) & long_covid_have_symptoms==1) %>%
  group_by(participant_id) %>%
  summarise(first_lc_yes_date = min(visit_date),
            last_lc_yes_date = max(visit_date))

dat <- dat %>%
  left_join(dat_first_last_lc_yes, by=join_by(participant_id==participant_id))

############################## SOCIO-DEMOGRAPHICS ##############################

### retrieve charactersitics at enrolment visit
dat0 <- dat %>%
  filter(visit_date==first_visit_date) %>%
  select(participant_id, age_at_visit, health_conditions, health_conditions_impact) %>%
  rename(age_at_visit0 = age_at_visit,
         health_conditions_visit0 = health_conditions,
         health_conditions_impact_visit0 = health_conditions_impact)

dat <- dat %>%
  left_join(dat0, by=join_by(participant_id==participant_id))

### define 10-year age-band variables
dat <- dat %>%
  mutate(age10 = case_when(age_at_visit>=2 & age_at_visit<=17 ~ "02-17",
                           age_at_visit>=18 & age_at_visit<=29 ~ "18-29",
                           age_at_visit>=30 & age_at_visit<=39 ~ "30-39",
                           age_at_visit>=40 & age_at_visit<=49 ~ "40-49",
                           age_at_visit>=50 & age_at_visit<=59 ~ "50-59",
                           age_at_visit>=60 & age_at_visit<=69 ~ "60-69",
                           age_at_visit>=70 & age_at_visit<=79 ~ "70-79",
                           age_at_visit>=80 ~ "80+",
                           TRUE ~ "Missing"),
         age10_visit0 = case_when(age_at_visit0>=2 & age_at_visit0<=17 ~ "02-17",
                                  age_at_visit0>=18 & age_at_visit0<=29 ~ "18-29",
                                  age_at_visit0>=30 & age_at_visit0<=39 ~ "30-39",
                                  age_at_visit0>=40 & age_at_visit0<=49 ~ "40-49",
                                  age_at_visit0>=50 & age_at_visit0<=59 ~ "50-59",
                                  age_at_visit0>=60 & age_at_visit0<=69 ~ "60-69",
                                  age_at_visit0>=70 & age_at_visit0<=79 ~ "70-79",
                                  age_at_visit0>=80 ~ "80+",
                                  TRUE ~ "Missing"))

### impute missing health conditions and impact variables
dat <- dat %>%
  mutate(health_conditions = ifelse(is.na(health_conditions), 0, health_conditions),
         health_conditions_impact = ifelse(is.na(health_conditions_impact), 0, health_conditions_impact),
         health_conditions_visit0 = ifelse(is.na(health_conditions_visit0), 0, health_conditions_visit0),
         health_conditions_impact_visit0 = ifelse(is.na(health_conditions_impact_visit0), 0, health_conditions_impact_visit0))

### define health/disability status variables
dat <- dat %>%
  mutate(health_status = case_when(health_conditions==1 & health_conditions_impact==0 ~ 1,
                                   health_conditions==1 & health_conditions_impact==1 ~ 2,
                                   health_conditions==1 & health_conditions_impact==2 ~ 3,
                                   TRUE ~ 0),
         health_status_visit0 = case_when(health_conditions_visit0==1 & health_conditions_impact_visit0==0 ~ 1,
                                   health_conditions_visit0==1 & health_conditions_impact_visit0==1 ~ 2,
                                   health_conditions_visit0==1 & health_conditions_impact_visit0==2 ~ 3,
                                   TRUE ~ 0))

### define white/non-white variable
dat <- dat %>%
  mutate(white = ifelse(ethnicityg==1, 1, 0))

### define IMD quintile variable
dat <- dat %>%
  mutate(imd_quintile_eng = case_when(imd_samp>(0*32844/5) & imd_samp<=(1*32844/5) ~ 1,
                                      imd_samp>(1*32844/5) & imd_samp<=(2*32844/5) ~ 2,
                                      imd_samp>(2*32844/5) & imd_samp<=(3*32844/5) ~ 3,
                                      imd_samp>(3*32844/5) & imd_samp<=(4*32844/5) ~ 4,
                                      imd_samp>(4*32844/5) & imd_samp<=(5*32844/5) ~ 5,
                                      TRUE ~ -999),
         
         imd_quintile_wal = case_when(imd_samp>(0*1909/5) & imd_samp<=(1*1909/5) ~ 1,
                                      imd_samp>(1*1909/5) & imd_samp<=(2*1909/5) ~ 2,
                                      imd_samp>(2*1909/5) & imd_samp<=(3*1909/5) ~ 3,
                                      imd_samp>(3*1909/5) & imd_samp<=(4*1909/5) ~ 4,
                                      imd_samp>(4*1909/5) & imd_samp<=(5*1909/5) ~ 5,
                                      TRUE ~ -999),
         
         imd_quintile_sco = case_when(imd_samp>(0*6976/5) & imd_samp<=(1*6976/5) ~ 1,
                                      imd_samp>(1*6976/5) & imd_samp<=(2*6976/5) ~ 2,
                                      imd_samp>(2*6976/5) & imd_samp<=(3*6976/5) ~ 3,
                                      imd_samp>(3*6976/5) & imd_samp<=(4*6976/5) ~ 4,
                                      imd_samp>(4*6976/5) & imd_samp<=(5*6976/5) ~ 5,
                                      TRUE ~ -999),
         
         imd_quintile_ni = case_when(imd_samp>(0*890/5) & imd_samp<=(1*890/5) ~ 1,
                                     imd_samp>(1*890/5) & imd_samp<=(2*890/5) ~ 2,
                                     imd_samp>(2*890/5) & imd_samp<=(3*890/5) ~ 3,
                                     imd_samp>(3*890/5) & imd_samp<=(4*890/5) ~ 4,
                                     imd_samp>(4*890/5) & imd_samp<=(5*890/5) ~ 5,
                                     TRUE ~ -999),
         
         imd_quintile = case_when(country==0 ~ imd_quintile_eng,
                                  country==1 ~ imd_quintile_wal,
                                  country==3 ~ imd_quintile_sco,
                                  country==2 ~ imd_quintile_ni,
                                  TRUE ~ -999))

############################## TIME-VARYING VARIABLES ##############################

### derive LC flags
dat <- dat %>%
  mutate(long_covid_have_symptoms = ifelse(!is.na(long_covid_have_symptoms) & long_covid_have_symptoms==1 & (is.na(infection_date) | datediff(visit_date, infection_date) < 84), 0, long_covid_have_symptoms),
         reduce_activities_long_covid = ifelse(is.na(long_covid_have_symptoms), NA, reduce_activities_long_covid),
         reduce_activities_long_covid = ifelse(!is.na(long_covid_have_symptoms) & long_covid_have_symptoms==0, 0, reduce_activities_long_covid),
         reduce_activities_long_covid = ifelse(is.na(reduce_activities_long_covid) & !is.na(long_covid_have_symptoms), long_covid_have_symptoms, reduce_activities_long_covid),
         lc_any = long_covid_have_symptoms,
         lc_lim = case_when(reduce_activities_long_covid %in% 0:1 ~ 0,
                            reduce_activities_long_covid %in% 2:3 ~ 1,
                            TRUE ~ NA))

### derive other time-varying variables
dat <- dat %>%
  mutate(time_since_infection = datediff(visit_date, infection_date),
         infected_12w_ago = ifelse(!is.na(infection_date) & datediff(visit_date, infection_date) >= 84, 1, 0),
         infected_to_date = ifelse(!is.na(infection_date) & visit_date >= infection_date, 1, 0),
         reinfected_to_date = ifelse(!is.na(reinfection_date) & visit_date >= reinfection_date, 1, 0),
         lc_to_date = ifelse(!is.na(first_lc_yes_date) & visit_date >= first_lc_yes_date, 1, 0),
         post_lc = ifelse(!is.na(last_lc_yes_date) & visit_date > last_lc_yes_date, 1, 0),
         remote_collection = ifelse(!is.na(dataset) & dataset==3, 1, 0),
         covid_admitted = ifelse(is.na(covid_admitted), 0, covid_admitted))

############################## SURVEY WEIGHTS ##############################

### count adult sample participants by age group, sex and geography, join
### population totals, and derive survey weights
svywgts <- dat %>%
  distinct(participant_id, .keep_all=TRUE) %>%
  filter(age_at_visit0 >= 18) %>%
  group_by(age10_visit0, sex, gor9d) %>%
  summarise(sample = n(), .groups="drop") %>%
  left_join(pops, by=join_by(age10_visit0==age10, sex==sex, gor9d==gor9d)) %>%
  mutate(svywgt = population / sample)

### join survey weights to visit-level dataset
dat <- dat %>%
  left_join(svywgts, by=join_by(age10_visit0==age10_visit0, sex==sex, gor9d==gor9d))

############################## FILTERING ##############################

### initial sample waterfall
n_people <- c(unlist(dat %>% distinct(participant_id) %>% summarise(n()) %>% collect()))
names(n_people) <- "All CIS participants"

### filter to participants aged 18+ at enrolment, and update sample waterfall
dat <- dat %>%
  filter(age_at_visit0 >= 18)

n_people <- c(n_people, unlist(dat %>% distinct(participant_id) %>% summarise(n()) %>% collect()))
names(n_people)[length(n_people)] <- "Aged 18+ years at enrolment"

### filter to participants with first positive swab before 1 May 2021 (when Delta became
### dominant), and update sample waterfall
dat <- dat %>%
  filter(!is.na(infection_date) & infection_date < '2021-05-01')

n_people <- c(n_people, unlist(dat %>% distinct(participant_id) %>% summarise(n()) %>% collect()))
names(n_people)[length(n_people)] <- "First positive swab before 1 May 2021 (start of Delta period)"

### filter to participants who were not vaccinated against COVID-19 at first positive swab,
### and update sample waterfall
dat <- dat %>%
  filter(is.na(covid_vaccine_date1) | covid_vaccine_date1 > infection_date)

n_people <- c(n_people, unlist(dat %>% distinct(participant_id) %>% summarise(n()) %>% collect()))
names(n_people)[length(n_people)] <- "Not vaccinated against COVID-19 at first positive swab"

### filter to participants who never reported being hospitalised with acute COVID-19
dat_ever_admitted <- dat %>%
  group_by(participant_id) %>%
  summarise(ever_admitted = max(covid_admitted))

dat <- dat %>%
  left_join(dat_ever_admitted, by=join_by(participant_id==participant_id)) %>%
  filter(ever_admitted==0)

n_people <- c(n_people, unlist(dat %>% distinct(participant_id) %>% summarise(n()) %>% collect()))
names(n_people)[length(n_people)] <- "Never reported being hospitalised with acute COVID-19"

### write out sample waterfall
sample_waterfall <- as.data.frame(n_people)
colnames(sample_waterfall) <- "Participants"
write.csv(sample_waterfall, file="filepath/sample_waterfall.csv")

############################## FINALISE DATASET ##############################

### restrict dataset to variables of interest
dropvars <- c(
  "covid_vaccine_date1",
  "long_covid_have_symptoms",
  "reduce_activities_long_covid",
  "imd_samp",
  "imd_quintile_eng",
  "imd_quintile_wal",
  "imd_quintile_sco",
  "imd_quintile_ni",
  "covid_admitted",
  "dataset",
  "sample",
  "population"
)

dat <- dat %>%
  select(-all_of(dropvars))

### save dataset locally
dat <- collect(dat)
dat <- as.data.frame(dat)
save(dat, file="filepath/lc_trajectories_dataset.RData")
