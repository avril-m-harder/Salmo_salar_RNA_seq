## Survival analyses for S. salar - 2017 spawning
## Some code written by C. Searle

library(OIsurv);library(survival);library(KMsurv);library(survsim);library(broom);library(plyr);library(dplyr);library(ggfortify);library(ggplot2);library(gridExtra);library(survminer)

setwd("/Users/Avril/Documents/atlantic_salmon/2017_spawning/2017_survival_analysis/")

## load data
data <- read.table("2017_R_t1_surv_and_thiamine_FINAL.csv", header=T, as.is=TRUE, sep=",")
attach(data)

#--------- Survival analysis for all families X treatments
sur.obj2 <- Surv(time = time , event = delta) ## create survival object, event = status (died=1 or censored=0) 
sur.obj2

srvfit2 = survfit(sur.obj2 ~ group, data = data, conf.type = "plain") ## compare survival distributions of populations
str(srvfit2)
summary(srvfit2)

#--------- Post-hoc comparisons to determine which groups are significantly different
post.diff <- pairwise_survdiff(Surv(time,delta)~group,data = data,rho = 0, p.adjust.method = "BH")
post.diff
summary(post.diff)
str(post.diff)

# -------- Analyze families independently and plot w/ ggplot
data$group <- as.factor(data$group)
i = 1
j = 1
no.fams <- length(levels(data$group))/2
result.plots <- vector(mode="list", length=no.fams)

for (j in 1:no.fams) {
  sub.group <- subset(data, group %in% c(levels(data$group)[i],levels(data$group)[i+1]))
  names(sub.group) <- names(data)
  attach(sub.group)
  sur.obj2 <- Surv(time = time , event = delta)
  srvfit2 = survfit(sur.obj2 ~ group, data = sub.group, conf.type = "plain")
  srvfit2
  print(srvfit2)
  fam = gsub('.{1}$',"",group[i])
  result.plots[[j]] <- ggsurvplot(srvfit2, size=1, palette=c("firebrick3", "gold3"), conf.int=T, pval=T, legend.labs=c("Untreated", "Treated"), ncensor.plot=F, xlab="Days post-fertilization", ylab="Proportion surviving", title=substitute(paste("Survival: Family", f), list(f=fam)))
  detach(sub.group)
  i <- i+2
  j <- j+1
}

pdf("/Users/Avril/Desktop/take_1_survivorship_plots.pdf", width=6, height=6)
result.plots
dev.off()

## mean ages at death for each family X treatment combo
tapply(time[delta==1],group[delta==1],mean)

# -------- Creating a composite plot of all families
sur.obj2 <- Surv(time = time , event = delta) ## create survival object, event = status (died=1 or censored=0) 
sur.obj2

srvfit2 = survfit(sur.obj2 ~ treatment, data = data, conf.type = "plain") ## compare survival distributions of populations
str(srvfit2)
summary(srvfit2)

pdf("/Users/Avril/Desktop/surv_all_families.pdf", width=8, height=6)
  ggsurvplot(srvfit2, data, size=1, palette=c("firebrick3", "gold3"), conf.int=T, pval=T, legend.labs=c("Untreated", "Treated"), ncensor.plot=F, xlab="Days post-fertilization", ylab="Proportion surviving", font.x=c(20), font.y=c(20), font.tickslab=c(20))
dev.off()

# -------- Fit Cox proportional hazards regression model to get family hazard ratio values
untreated.data <- subset(data, data$treatment=="-" )
untreated.data$group <- factor(untreated.data$group, ordered=FALSE)
untreated.data$group <- relevel(untreated.data$group, ref="13-")
coxfit2 <- coxph(Surv(untreated.data$time, untreated.data$delta) ~ untreated.data$group, data = untreated.data)
summary(coxfit2)
haz.ratios <- summary(coxfit2)[["coefficients"]][,2]
low.95 <- exp(confint(coxfit2))[,1]
up.95 <- exp(confint(coxfit2))[,2]

all.vs.13unt <- cbind(haz.ratios, low.95, up.95)
write.csv(as.data.frame(all.vs.13unt), file="/Users/Avril/Desktop/all_vs_13_untreated.csv")
