emmeans.collate <- function(mod_obj, vcov_obj, emmeans_weights, emmeans_infection_day) {
  
  start_day <- seq(84,max(dat$time_since_infection),10)
  end_day <- seq(93,max(dat$time_since_infection),10)
  if(length(end_day) < length(start_day)) {end_day <- c(end_day, max(dat$time_since_infection))}
  
  emmeans_out_list <- as.list(NULL)
  
  for(i in 1:length(start_day)) {
  
    emmeans_out_list[[i]] <- as.data.frame(emmeans(
      object = mod_obj,
      vcov. = vcov_obj,
      specs = ~time_since_infection,
      weights = emmeans_weights,
      at = list(time_since_infection = start_day[i]:end_day[i],
                infection_day = emmeans_infection_day),
      type = "response",
      params = c("df_select1", "boundary_select1",
                 "df_select2", "boundary_select2"),
      rg.limit = 100000
    ))
  
  }
  
  emmeans_out <- emmeans_out_list[[1]]
  if(length(emmeans_out_list)>1) {
    for(i in 2:length(emmeans_out_list)) {
      emmeans_out <- rbind(emmeans_out, emmeans_out_list[[i]])
    }
  }
  
  return(emmeans_out)
  
}
