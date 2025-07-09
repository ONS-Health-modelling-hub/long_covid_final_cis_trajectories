cov.dist.cont <- function(vars, dataset, exposure) {

out_list1 <- as.list(NULL)

for(i in 1:length(vars)) {
    
    ds0 <- dataset[dataset[[exposure]]==0,]
    ds1 <- dataset[dataset[[exposure]]==1,]
    
    mu_all <- mean(dataset[[vars[i]]])
    mu0 <- mean(ds0[[vars[i]]])
    mu1 <- mean(ds1[[vars[i]]])
    
    sd_all <- sd(dataset[[vars[i]]])
    sd0 <- sd(ds0[[vars[i]]])
    sd1 <- sd(ds1[[vars[i]]])
    
    var_all <- var(dataset[[vars[i]]])
    var0 <- var(ds0[[vars[i]]])
    var1 <- var(ds1[[vars[i]]])
    
    abs_std_diff <- abs((mu1 - mu0) / sqrt((var1 + var0) / 2))
    
    out_list1[[i]] <- as.data.frame(cbind(vars[i], "mean",
                                          mu_all, sd_all, mu0, sd0, mu1, sd1, abs_std_diff))
  
}

out_df1 <- out_list1[[1]]
if(length(vars)>1) {
  for(i in 2:length(vars)) {out_df1 <- rbind(out_df1, out_list1[[i]])}
}

rownames(out_df1) <- NULL
colnames(out_df1)[1:2] <- c("characteristic", "level")

return(out_df1)

}
