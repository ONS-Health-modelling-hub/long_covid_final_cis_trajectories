library(dplyr)
library(ggplot2)

source("filepath/cov_dist_cont.r")
source("filepath/cov_dist_cat.r")

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

### calculate date of last response to LC question within follow-up period
dat_last_lc_response_within_fu <- dat %>%
  filter(!is.na(lc_any) &
         time_since_infection >= 84 &
         time_since_infection <= 365.25*2.5) %>%
  group_by(participant_id) %>%
  summarise(last_lc_response_date_within_fu = max(visit_date))

dat <- dat %>%
  left_join(dat_last_lc_response_within_fu, by=join_by(participant_id==participant_id))

### flag if participant completed follow-up to end of CIS (i.e. last response to the LC
### question was within 6 weeks of the CIS ending)
dat$cis_completed <- ifelse(!is.na(dat$last_lc_response_date) &
                              dat$last_lc_response_date >= '2023-01-30', 1, 0)

### flag if participant dropped out within 2.5 years of infection, for reasons other than
### end of the CIS
dat$dropped_out <- ifelse(!is.na(dat$last_lc_response_date) &
                            as.numeric(dat$last_lc_response_date - dat$infection_date) < 365.25*2.5 &
                            dat$last_lc_response_date < '2023-01-30', 1, 0)

### create person-level version of the dataset
dat_pers <- dat[!duplicated(dat$participant_id),]

############################### SAMPLE COUNTS OVER TIME ###############################

### create new dataset
dat_counts <- dat

### restrict dataset to visits >= 12 weeks after infection
dat_counts <- dat_counts[dat_counts$time_since_infection >= 84,]

### restrict dataset to visits when participants reported their LC status
dat_counts <- dat_counts[!is.na(dat_counts$lc_any),]

### calculate months since infection
dat_counts$months_since_infection <- floor(dat_counts$time_since_infection / (365.25/12))
dat_counts$months_since_infection <- ifelse(dat_counts$months_since_infection %in% 2:3, 2.3, dat_counts$months_since_infection)

### count number of participants and number with LC by month since infection
sample_counts_over_time <- dat_counts %>%
  group_by(participant_id, months_since_infection) %>%
  summarise(lc_any=max(lc_any), .groups="drop") %>%
  group_by(months_since_infection) %>%
  summarise(n=n(), lc_any=sum(lc_any))

write.csv(sample_counts_over_time, file="filepath/sample_counts_over_time.csv", row.names=FALSE)

######################### COVARIATE DISTRIBUTIONS AT ENROLMENT #########################

### create new dataset
dat_covs <- dat

### derive indicator for 'all participants'
dat_covs$all_participants <- 1

### derive indicator if participant ever reported LC during follow-up
dat_ever_lc <- dat_covs %>%
  filter(!is.na(lc_any) &
         time_since_infection >= 84 &
         time_since_infection <= 365.25*2.5) %>%
  group_by(participant_id) %>%
  summarise(ever_lc=max(lc_any))

dat_covs <- merge(x=dat_covs, y=dat_ever_lc, by="participant_id",
                  all.x=TRUE, all.y=TRUE)

### keep enrolment visit
dat_covs <- dat_covs[dat_covs$visit_date==dat_covs$first_visit_date,]

### covariate distributions for all participants - continuous variables
cov_dist_cont_all <- cov.dist.cont(
  vars = c("age_at_visit0"),
  dataset = dat_covs,
  exposure = "all_participants"
)

write.csv(cov_dist_cont_all, file="filepath/cov_dist_cont_all.csv", row.names=FALSE)

### covariate distributions for all participants - categorical variables
cov_dist_cat_all <- cov.dist.cat(
  vars = c("age10_visit0", "sex", "white", "gor9d",
           "imd_quintile", "health_status_visit0"),
  dataset = dat_covs,
  exposure = "all_participants"
)

write.csv(cov_dist_cat_all, file="filepath/cov_dist_cat_all.csv", row.names=FALSE)

### covariate distributions by ever-LC status - continuous variables
cov_dist_cont_lc <- cov.dist.cont(
  vars = c("age_at_visit0"),
  dataset = dat_covs[!is.na(dat_covs$ever_lc),],
  exposure = "ever_lc"
)

write.csv(cov_dist_cont_lc, file="filepath/cov_dist_cont_lc.csv", row.names=FALSE)

### covariate distributions by ever-LC status - categorical variables
cov_dist_cat_lc <- cov.dist.cat(
  vars = c("age10_visit0", "sex", "white", "gor9d",
           "imd_quintile", "health_status_visit0"),
  dataset = dat_covs[!is.na(dat_covs$ever_lc),],
  exposure = "ever_lc"
)

write.csv(cov_dist_cat_lc, file="filepath/cov_dist_cat_lc.csv", row.names=FALSE)

############################### FOLLOW-UP TIME STATS ###############################

### create new dataset
dat_fu <- dat

### make dataset one row per participant
dat_fu <- dat_fu[!duplicated(dat_fu$participant_id),]

### calculate follow-up time from infection date to date of last response to LC question
dat_fu$futime <- as.numeric(dat_fu$last_lc_response_date_within_fu - dat_fu$infection_date)

### summary stats for follow-up time
futime_summary_stats <- as.data.frame(as.matrix(summary(dat_fu$futime)))
colnames(futime_summary_stats) <- "FU_time"

write.csv(futime_summary_stats, file="filepath/futime_summary_stats.csv")

### plot distributon of follow-up time
ggplot(dat_fu, aes(x=futime)) +
  geom_density(colour="black", fill=NA, linewidth=1) +
  xlab("Follow-up (days)") +
  ylab("Density") +
  scale_x_continuous(expand=expansion(mult=c(0,0))) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme(
    axis.text = element_text(colour="black", size=11),
    axis.title = element_text(colour="black", face="bold", size=11),
    panel.border = element_rect(colour="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(filename="filepath/futime_density.jpg", width=15, height=10, units="cm", dpi=300)

### calculate percentage of participants with complete follow-up to end of CIS
cis_completed <- sum(dat_pers$cis_completed) / nrow(dat_pers) * 100
cis_completed <- as.data.frame(cis_completed)
colnames(cis_completed) <- "Pct_completed_CIS"
rownames(cis_completed) <- c("All participants")
write.csv(cis_completed, file="filepath/cis_completed.csv")

### calculate percentage of participants who dropped out within 2.5 years of infection,
### for reasons other than end of the CIS
dropped_out <- sum(dat_pers$dropped_out) / nrow(dat_pers) * 100
dropped_out <- as.data.frame(dropped_out)
colnames(dropped_out) <- "Pct_dropped_out"
rownames(dropped_out) <- c("All participants")
write.csv(dropped_out, file="filepath/dropped_out.csv")

############ FOLLOW-UP TIME STATS - BY LC STATUS 12-26 WEEKS POST-INFECTION ############

### create new dataset
dat_fu_lc <- dat

### restrict dataset to visits 12-26 weeks after infection
dat_fu_lc <- dat_fu_lc[dat_fu_lc$time_since_infection >= 84 &
                         dat_fu_lc$time_since_infection <= 182,]

### restrict dataset to visits when participants reported their LC status
dat_fu_lc <- dat_fu_lc[!is.na(dat_fu_lc$lc_any),]

### keep each participant's first visit within the 12-26 weeks interval (index visit)
dat_fu_lc <- dat_fu_lc[order(dat_fu_lc$participant_id, dat_fu_lc$visit_date),]
dat_fu_lc <- dat_fu_lc[!duplicated(dat_fu_lc$participant_id),]

### calculate follow-up time from the index visit
dat_fu_lc$futime_post_lc_response <- as.numeric(dat_fu_lc$last_lc_response_date_within_fu - dat_fu_lc$visit_date)

### summary stats for follow-up time
futime_summary_stats_lc <- data.frame(
  Without_LC = as.matrix(summary(dat_fu_lc$futime_post_lc_response[dat_fu_lc$lc_any==0])),
  With_LC = as.matrix(summary(dat_fu_lc$futime_post_lc_response[dat_fu_lc$lc_any==1]))
)

write.csv(futime_summary_stats_lc, file="filepath/futime_summary_stats_lc.csv")

### plot distributon of follow-up time
dat_fu_lc$lc_label <- as.factor(dat_fu_lc$lc_any)
levels(dat_fu_lc$lc_label) <- c("Without Long Covid", "With Long Covid")

ggplot(dat_fu_lc, aes(x=futime_post_lc_response, group=lc_label, fill=lc_label, colour=lc_label)) +
  geom_density(fill=NA, linewidth=1) +
  xlab("Follow-up (days)") +
  ylab("Density") +
  scale_colour_manual(values=c("blue","orange")) +
  scale_x_continuous(expand=expansion(mult=c(0,0))) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme(
    axis.text = element_text(colour="black", size=11),
    axis.title = element_text(colour="black", face="bold", size=11),
    panel.border = element_rect(colour="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(colour="black", size=11)
  )

ggsave(filename="filepath/futime_density_lc.jpg", width=15, height=10, units="cm", dpi=300)

#### calculate percentage of participants with complete follow-up to end of CIS
cis_completed_lc <- c(sum(dat_fu_lc$cis_completed[dat_fu_lc$lc_any==0]) / nrow(dat_fu_lc[dat_fu_lc$lc_any==0,]) * 100,
                      sum(dat_fu_lc$cis_completed[dat_fu_lc$lc_any==1]) / nrow(dat_fu_lc[dat_fu_lc$lc_any==1,]) * 100)

cis_completed_lc <- as.data.frame(cis_completed_lc)
colnames(cis_completed_lc) <- "Pct_completed_CIS"
rownames(cis_completed_lc) <- c("Without Long Covid", "With Long Covid")
write.csv(cis_completed_lc, file="filepath/cis_completed_lc.csv")

### calculate percentage of participants who dropped out within 2.5 years of infection,
### for reasons other than end of the CIS
dropped_out_lc <- c(sum(dat_fu_lc$dropped_out[dat_fu_lc$lc_any==0]) / nrow(dat_fu_lc[dat_fu_lc$lc_any==0,]) * 100,
                    sum(dat_fu_lc$dropped_out[dat_fu_lc$lc_any==1]) / nrow(dat_fu_lc[dat_fu_lc$lc_any==1,]) * 100)

dropped_out_lc <- as.data.frame(dropped_out_lc)
colnames(dropped_out_lc) <- "Pct_dropped_out"
rownames(dropped_out_lc) <- c("Without Long Covid", "With Long Covid")
write.csv(dropped_out_lc, file="filepath/dropped_out_lc.csv")
