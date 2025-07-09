cov.dist.cat <- function(vars, dataset, exposure) {
  
  out_list1 <- as.list(NULL)
  
  for(i in 1:length(vars)) {
    
    dataset[[vars[i]]] <- as.factor(dataset[[vars[i]]])
    dataset[[vars[i]]] <- droplevels(dataset[[vars[i]]])
    
    ds0 <- dataset[dataset[[exposure]]==0,]
    ds1 <- dataset[dataset[[exposure]]==1,]
    
    n_all <- as.numeric(table(dataset[[vars[i]]]))
    n0 <- as.numeric(table(ds0[[vars[i]]]))
    n1 <- as.numeric(table(ds1[[vars[i]]]))
    
    p_all <- n_all / nrow(dataset)
    p0 <- n0 / nrow(ds0)
    p1 <- n1 / nrow(ds1)
    
    var_all <- p_all*(1-p_all)
    var0 <- p0*(1-p0)
    var1 <- p1*(1-p1)
    
    K = nlevels(dataset[[vars[i]]])
    
    if(K==1) {abs_std_diff <- 0}
    
    if(K>1) {
      diff <- p1[-1] - p0[-1]
      k=rep(2:K,times=K-1)
      l=rep(2:K,each=K-1)
      s <- ifelse(k==l, (var0[k] + var1[l]) / 2, (p0[k]*p0[l] + p1[k]*p1[l]) / 2)
      s <- matrix(s, nrow=K-1, byrow=FALSE)
      abs_std_diff <- abs(as.numeric(sqrt(diff %*% solve(s) %*% diff)))
      abs_std_diff <- c(abs_std_diff, rep("",K-1))
    }
    
    out_list1[[i]] <- as.data.frame(cbind(c(vars[i], rep("", nlevels(dataset[[vars[i]]])-1)),
                                          levels(dataset[[vars[i]]]),
                                          n_all, p_all, n0, p0, n1, p1, abs_std_diff))
    
  }
  
  out_df1 <- out_list1[[1]]
  if(length(vars)>1) {
    for(i in 2:length(vars)) {out_df1 <- rbind(out_df1, out_list1[[i]])}
  }
  
  rownames(out_df1) <- NULL
  colnames(out_df1)[1:2] <- c("characteristic", "level")
  
  return(out_df1)
  
}
