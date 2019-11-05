## For each gene and treatment group, saves model coefficients and predicted y-values for ## each x in each permutation (1000 permutations)
##
## permutation_env.RData includes: transcript counts normalized using fpm() command on dds 
## object in DESeq2 and dds object, which includes sample information (can be generated
## using deseq2.R
##
## all_vs_13.csv includes family IDs and hazard ratio information, which can be calculated
## using survival_analysis.R
##
## permutation loop written by J. Willoughby

setwd("/scratch/snyder/h/harder/ASperm_permutations/log_haz_values_incl_y_val_cals/")
load("permutation_env.RData")

fam.data <- read.csv("all_vs_13.csv")
full.names <- deg.names

OUTPUT1 <- NULL
for (i in full.names) {
  temp <- fpm.dds[which(rownames(fpm.dds)==i),]
  temp <- as.data.frame(temp)
  temp <- cbind(temp, dds.w.lengths@colData@listData[["treatment"]])
  temp <- cbind(temp, dds.w.lengths@colData@listData[["family"]])
  temp <- temp[order(temp[,3], temp[,2]),]
  temp <- cbind(temp, rownames(temp))
  colnames(temp) <- c("fpm","treatment","family","fry.id")
  temp <- merge(x=temp, y=fam.data[,c(2,6)], by="family")
  # temp$fry.id <- as.character(temp$fry.id)
  temp$replic <- c(rep(seq(1,2,1), 9))
  
  title <- do.call(rbind, strsplit(i,"[|]"))
  
  ## run permutations and save results
  OUTPUT <- NULL
  perms <- 1000
  
  #permutation loop
  for(p in 1:perms){ 
    #find the wanted records, one rep per treatment/hazard (family) combination
    startpos = seq(1, nrow(temp), 2) #this requires replicates to be ordered one after each other
    temp2 = NULL
    for(r in 1:length(startpos)){
      grab = sample(c(startpos[r], (startpos[r]+1)), 1) #this is the row number to select
      temp2 = rbind(temp2, temp[grab,])
    }#r
    print(temp2)
    
    #run regressions
    t0 = temp2[temp2$treatment=="untreated", ]
    t1 = temp2[temp2$treatment=="treated", ]
    
    regout = lm(fpm~log.haz.ratios, data=t0)
    newdata <- expand.grid(log.haz.ratios=seq(from=min(t0$log.haz.ratios), to=max(t0$log.haz.ratios), length=9))
    pm <- predict(regout, newdata, interval="confidence")
    newdata$p <- pm[,1]
    newdata$pmin <- pm[,2]
    newdata$pmax <- pm[,3]
    tosave = c(summary(regout)$coefficients[1,1], summary(regout)$coefficients[2,1], newdata$p, summary(regout)$coefficients[,4][2], summary(regout)$adj.r.squared)
    
    regout = lm(fpm~log.haz.ratios, data=t1)
    newdata <- expand.grid(log.haz.ratios=seq(from=min(t1$log.haz.ratios), to=max(t1$log.haz.ratios), length=9))
    pm <- predict(regout, newdata, interval="confidence")
    newdata$p <- pm[,1]
    newdata$pmin <- pm[,2]
    newdata$pmax <- pm[,3]
    tosave = c(tosave, summary(regout)$coefficients[1,1], summary(regout)$coefficients[2,1], newdata$p, summary(regout)$coefficients[,4][2], summary(regout)$adj.r.squared)
    OUTPUT = rbind(OUTPUT, tosave)
    rownames(OUTPUT)[length(rownames(OUTPUT))] <- i
  }#p
  
  OUTPUT1 <- rbind(OUTPUT1, OUTPUT)
  colnames(OUTPUT1)[1:13] <- c("unt.int","unt.slope","unt.1","unt.2","unt.3","unt.4","unt.5","unt.6","unt.7","unt.8","unt.9","unt.p","unt.adj.r2")
  colnames(OUTPUT1)[14:26] <- c("tre.int","tre.slope","tre.1","tre.2","tre.3","tre.4","tre.5","tre.6","tre.7","tre.8","tre.9","tre.p","tre.adj.r2")
}

## columns:
## untreated / treated
## intercept / haz. rank / intercept / haz. rank
# OUTPUT1

write.csv(OUTPUT1, "permutation_output_coeffs_and_pred_ys.csv")
