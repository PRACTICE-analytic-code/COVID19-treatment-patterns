# Aim 1 Tables and Figures
# Author: Daniel Vader
# Date created: 12/01/2020

library(tidyverse)

adat <- readRDS("/data/analytic/PRACTICE-analytic-data_2021-08-23.rds") %>%
  mutate(ecog.c = ifelse(ecog_b30a14 >= 2, 1, 0))

################################################################################
# ECOG NA assessment
################################################################################
ecogvars.c <- c("ecog_b14a7", "ecog_b28a7", "ecog_b56a7", "ecog_b28a14", 
                "ecog_b28a21", "ecog_b28a28")
adat.t <- adat %>% filter(studyperiod != "Washout 2020", studyperiod != "Washout 2019") %>%
  as.data.frame()

etab <- tibble(name = character(), nmiss = numeric(), ntot = numeric(), propmiss = numeric())

for(i in 1:length(ecogvars.c)){
  tn <- tn <- table(is.na(dplyr::pull(adat.t, ecogvars.c[i])))
  tp <- prop.table(tn)
  
  etab <- etab %>% add_row(name = ecogvars.c[i], nmiss=tn[2], ntot=nrow(adat.t), propmiss = tp[2])
}

write.csv(etab, "tables/intermediary/ecogmiss.csv", na = "")

################################################################################
# Descriptive tables (todo add categorical time variable) 
################################################################################
library(stddiff)

tabvars.c <- c("gend", "reth", "insclass", "cancer", "practicetype", "ecog.c", 
               "denovo_met", "liverf_adv", "kidneyf_adv",
               "steroid_adv", "opioid_adv")#, "myelo.yes") 
# and comorbidities (liver, renal, steroid use, opioid use) once coded

# tabvars.c <- c("gend", "reth", "insclass", "cancer", "practicetype", "ecog.c", 
#                "liverf_adv", "kidneyf_adv",
#                "steroid_adv", "opioid_adv")

tabvars.n <- c("age")

# Use NA in tables?
tabna <- "no"

adat.t <- adat %>% 
  filter(studyperiod != "Washout 2020", studyperiod != "Washout 2019") %>%
  mutate(postcovid = ifelse(precovid == "pre", 0, 1)) %>%
  as.data.frame()

# adat.denov.y <- adat.t %>% filter(denovo_met == 1)
# adat.denov.n <- adat.t %>% filter(denovo_met == 0)
# adat.denov.na <- adat.t %>% filter(is.na(denovo_met))

tab <- tibble(name = character(), n2019pre = numeric(), prop2019pre = numeric(), 
              prop2019preq = numeric(),
              
              n2019post = numeric(), prop2019post = numeric(), 
              prop2019postq = numeric(), break1 = character(),
              
              n2020pre = numeric(), prop2020pre = numeric(),
              prop2020preq = numeric(),
              
              n2020post = numeric(), prop2020post = numeric(), 
              prop2020postq = numeric(), break2 = character(),
              
              ntot = numeric(), proptot = numeric(), proptotq = numeric(), 
              break3 = character(),
              p=numeric(),
              
              smd2019 = numeric(), smd2020 = numeric(), 
              smdpre = numeric(), smdpost= numeric(),
              smddd = numeric())

# Function for calculating covariance matrix for categorical smds
smdcov <- function(ptab){
  tat <- ptab[,1] %>% as.matrix()
  tac <- ptab[,2] %>% as.matrix()
  
  s.mat <- matrix(nrow=nrow(ptab)-1, ncol=nrow(ptab)-1)
  for(k in 2:nrow(ptab)){
    for(l in 2:nrow(ptab)){
      if(k == l){
        s.mat[k-1, l-1] <- (tat[k]*(1-tat[k])+tac[k]*(1-tac[k]))/2
      } else{
        s.mat[k-1, l-1] <- (tat[k]*tac[l]+tac[k]*tac[l])/2
      }
    }
  }
  return(s.mat)
}

for(i in 1:length(tabvars.c)){
  tn <- table(dplyr::pull(adat.t, tabvars.c[i]), dplyr::pull(adat.t, studyperiod), useNA=tabna)
  tnm <- rownames(tn)
  tp <- prop.table(tn, margin=2)
  pval <- chisq.test(tn)
  pval <- pval$p.value
  
  tnt <- table(dplyr::pull(adat.t, tabvars.c[i]),useNA=tabna)
  tpt <- prop.table(tnt)
  
  # Calculate and extract standardized mean difference by diagnosis year
  smd19 <- stddiff.category(adat.t[adat.t$diagnosisyear == "2019",], 
                            gcol=which((colnames(adat.t) == "postcovid")), 
                            vcol=which((colnames(adat.t) == tabvars.c[i])))[1,5]
  smd20 <- stddiff.category(adat.t[adat.t$diagnosisyear == "2020",], 
                            gcol=which((colnames(adat.t) == "postcovid")), 
                            vcol=which((colnames(adat.t) == tabvars.c[i])))[1,5]
  smdpre <- stddiff.category(adat.t[adat.t$postcovid == 0,], 
                            gcol=which((colnames(adat.t) == "diagnosisyear")), 
                            vcol=which((colnames(adat.t) == tabvars.c[i])))[1,5]
  smdpost <- stddiff.category(adat.t[adat.t$postcovid == 1,], 
                            gcol=which((colnames(adat.t) == "diagnosisyear")), 
                            vcol=which((colnames(adat.t) == tabvars.c[i])))[1,5]
  
  # SMD DD for reviewer (based on: https://stat.ethz.ch/R-manual/R-devel/library/base/html/t.html)
  ta2019 <- prop.table(table(dplyr::pull(adat.t[adat.t$diagnosisyear==2019,], tabvars.c[i]), 
                    dplyr::pull(adat.t[adat.t$diagnosisyear==2019,], postcovid), 
                    useNA="no"), margin=2)
  ta2020 <- prop.table(table(dplyr::pull(adat.t[adat.t$diagnosisyear==2020,], tabvars.c[i]), 
                             dplyr::pull(adat.t[adat.t$diagnosisyear==2020,], postcovid), 
                             useNA="no"), margin=2)
  
  # Calculate dd covariance matrix for categorical vars
  s.mat <- (smdcov(ta2019)*2 + smdcov(ta2020)*2)/2
  smdddr <- (ta2020[-1,1]-ta2020[-1,2]) - (ta2019[-1,1]-ta2019[-1,2])
  
  # Calculate SMD DD
  smddd <- sqrt(t(smdddr)%*%solve(s.mat)%*%smdddr)
  
  for(j in 1:length(tnm)){
    tab <- tab %>% add_row(
      name = paste(tabvars.c[i], ": ", tnm[[j]], sep=""),
      n2019pre = tn[[j, 3]],
      prop2019pre = tp[[j, 3]],
      n2019post = tn[[j, 1]],
      prop2019post = tp[[j,1]],
      n2020pre = tn[[j, 4]],
      prop2020pre = tp[[j,4]],
      n2020post = tn[[j,2]],
      prop2020post = tp[[j,2]],
      ntot = tnt[[j]],
      proptot = tpt[[j]],
      p=pval,
      smd2019=smd19,
      smd2020=smd20,
      smdpre=smdpre,
      smdpost=smdpost,
      smddd=smddd
    )
  }
}

for(i in 1:length(tabvars.n)){
  #tm <- aggregate(. ~ adat.t$studyperiod, adat.t[,tabvars.n[i]], mean, na.rm=T) 
  #tsd <- aggregate(. ~ adat.t$studyperiod, adat.t[,tabvars.n[i]], sd, na.rm=T) 
  tm <- aggregate(x=adat.t[,tabvars.n[i]], by=list(adat.t$studyperiod), median, na.rm=T) 
  tsd25 <- aggregate(x=adat.t[,tabvars.n[i]], by=list(adat.t$studyperiod), quantile, .25, na.rm=T)
  tsd75 <- aggregate(x=adat.t[,tabvars.n[i]], by=list(adat.t$studyperiod), quantile, .75, na.rm=T)
  pval <- aov(dplyr::pull(adat.t, tabvars.n[i]) ~ dplyr::pull(adat.t, studyperiod))
  pval <- summary(pval)[[1]][["Pr(>F)"]][[1]]
  
  #om <- mean(pull(adat.t, tabvars.n[i]), na.rm=T)
  #osd <- sd(pull(adat.t, tabvars.n[i]), na.rm=T)
  om <- median(pull(adat.t, tabvars.n[i]), na.rm=T)
  osd25 <- quantile(pull(adat.t, tabvars.n[i]), .25, na.rm=T)
  osd75 <- quantile(pull(adat.t, tabvars.n[i]), .75, na.rm=T)
  
  # Calculate and extract standardized mean difference by diagnosis year
  smd19 <- stddiff.numeric(adat.t[adat.t$diagnosisyear == "2019",], 
                            gcol=which((colnames(adat.t) == "postcovid")), 
                            vcol=which((colnames(adat.t) == tabvars.n[i])))[7]
  smd20 <- stddiff.numeric(adat.t[adat.t$diagnosisyear == "2020",], 
                            gcol=which((colnames(adat.t) == "postcovid")), 
                            vcol=which((colnames(adat.t) == tabvars.n[i])))[7]
  smdpre <- stddiff.category(adat.t[adat.t$postcovid == 0,], 
                             gcol=which((colnames(adat.t) == "diagnosisyear")), 
                             vcol=which((colnames(adat.t) == tabvars.n[i])))[7]
  smdpost <- stddiff.category(adat.t[adat.t$postcovid == 1,], 
                              gcol=which((colnames(adat.t) == "diagnosisyear")), 
                              vcol=which((colnames(adat.t) == tabvars.n[i])))[7]
  
  tme <- aggregate(x=adat.t[,tabvars.n[i]], by=list(adat.t$studyperiod), mean, na.rm=T) 
  varms <- aggregate(x=adat.t[,tabvars.n[i]], by=list(adat.t$studyperiod), var, na.rm=T) 
  smddd <- (tme[1,2]-tme[3,2]+tme[2,2]-tme[4,2])/sqrt((varms[1,2]+varms[3,2]+varms[2,2]+varms[4,2])/4)
  
  
  tab <- tab %>% add_row(
    name = paste(tabvars.n[i], "(mean, SD)"),
    n2019pre = tm[3, 2],
    prop2019pre = tsd25[3,2],
    prop2019preq = tsd75[3,2],
    n2019post = tm[1,2],
    prop2019post = tsd25[1,2],
    prop2019postq = tsd75[1,2],
    n2020pre = tm[4,2],
    prop2020pre = tsd25[4,2],
    prop2020preq = tsd75[4,2],
    n2020post = tm[2,2],
    prop2020post = tsd25[2,2],
    prop2020postq = tsd75[2,2],
    ntot = om,
    proptot = osd25,
    proptotq = osd75,
    p = pval,
    smd2019=smd19,
    smd2020=smd20,
    smdpre=smdpre,
    smdpost=smdpost,
    smddd=smddd
  )
}

#write.csv(tab, "tables/intermediary/baseline-table_2021-08-23.csv", na = "")
write.csv(tab, "tables/intermediary/baseline-table-na_2021-10-19.csv", na = "")
#write.csv(tab, "tables/intermediary/baseline-table-denov-y.csv", na="")
#write.csv(tab, "tables/intermediary/baseline-table-denov-n.csv", na="")
#write.csv(tab, "tables/intermediary/baseline-table-denov-na.csv", na="")

# State Frequency Table ########################################################
statetab <- ifelse(adat.t$practicetype == "ACADEMIC", "Academic Practice",
                   ifelse(is.na(adat.t$state), "Unknown",
                          adat.t$state)) %>% table()
write.csv(statetab, "tables/intermediary/statetab.csv", na = "")

################################################################################
### KM Curves
################################################################################
library(survival)
library(survminer)
library(viridis)
library(ggquickeda)
library(RColorBrewer)
library(plotrix)

sumprob <- function(fit, s=T){
  if(s){
    fitsum <- cbind(fit$strata,
                    1 -fit$surv, 
                    1-fit$upper, 
                    1-fit$lower) %>% as.data.frame()
    colnames(fitsum) <- c("var", "est", "lci", "uci")
  } else {
    fitsum <- cbind(1 -fit$surv, 
                    1-fit$upper, 
                    1-fit$lower) %>% as.data.frame()
    colnames(fitsum) <- c("est", "lci", "uci")
  }
  
  return(fitsum)
}

adat.nowash <- adat %>% 
  filter(!(adat$studyperiod %in% c("Washout 2019", "Washout 2020"))) %>%
  mutate(osurv.ind5 = ifelse(osurv.ind == 1 & osurv.months > 5, 0, osurv.ind),
         osurv.months5 = ifelse(osurv.months > 5, 5, osurv.months),
         pfsurv.ind90 = ifelse(pfsurv.ind == 1 & pfsurv.days > 90, 0, pfsurv.ind),
         pfsurv.days90 = ifelse(pfsurv.days >90, 90, pfsurv.days),
         pfsurv.trt.ind90 = ifelse(pfsurv.trt.ind == 1 & pfsurv.trt.days > 90, 0, pfsurv.trt.ind),
         pfsurv.trt.days90 = ifelse(pfsurv.trt.days >90, 90, pfsurv.trt.days),
         sp2 = factor(studyperiod, levels=c("Jan 1-Mar 8, 2019", 
                                            "Apr 8-Jul 31, 2019",
                                            "Jan 1-Mar 8, 2020",
                                            "Apr 8-Jul 31, 2020")),
         surv.ind.f = factor(surv.ind, levels=c(0,1), labels=c("Censored", "Treated")))


### Some stats for the paper text ###

# Overall TTI
tti90o.fit <- survfit(Surv(adat.nowash$surv.days90, adat.nowash$surv.ind90)~1)
tti90o.fit.s <- summary(tti90o.fit, time=30) 
tti90o.fit # print median survival time
sumprob(tti90o.fit.s, s=F) %>% apply(MARGIN=2, FUN=round, digits=3) # print 30 day survival

# Unadjusted TTI by cancer
tti90canc.fit <- survfit(Surv(adat.nowash$surv.days90, adat.nowash$surv.ind90)~adat.nowash$cancer)
tticanc.fit.s <- summary(tti90canc.fit, time=30) 
sumprob(tticanc.fit.s, s=T) %>% apply(MARGIN=2, FUN=round, digits=3)

# Unadjusted TTI by denovo metastatic status (note: some missingness)
tti90denovo.fit <- survfit(Surv(adat.nowash$surv.days90, adat.nowash$surv.ind90)~adat.nowash$denovo_met)
tti90denovo.fit.s <- summary(tti90denovo.fit, time=30) 
sumprob(tti90denovo.fit.s, s=T) %>% apply(MARGIN=2, FUN=round, digits=3)


### Period*Year and cancer-specific period*year survival curves. Fit and plot. ###

# TTI right censored at 90 days
tti90.fit <- survfit(Surv(adat.nowash$surv.days90, adat.nowash$surv.ind90) ~ adat.nowash$studyperiod)

# Extract unadjusted PP at 30 days by time period and overall
tti.fit.s <- summary(tti90.fit, time=30) 
tti.fit.sp <- cbind(tti.fit.s$strata,
                    1 -tti.fit.s$surv, 
                    1-tti.fit.s$upper, 
                    1-tti.fit.s$lower) %>% as.data.frame() %>%
  mutate(V1 = factor(V1, levels=c(1,2,3,4), labels=c("AprJul2019", 
                                                     "AprJul2020",
                                                     "JanMar2019",
                                                     "JanMar2020")))
colnames(tti.fit.sp) <- c("period", "est", "lci", "uci")
write.csv(tti.fit.sp, "tables/intermediary/unadjusted-30daysurv-probs.csv", na = "")  

splot90 <- ggsurvplot(tti90.fit, data=adat.nowash,
                      fun="pct",
                    conf.int = T,
                    risk.table = T,
                    risk.table.height = 0.25,
                    #legend.labs = c("[2019] Apr 8-Jul 31", "[2020] Apr 8-Jul 31", 
                    #                "[2019] Jan 1-Mar 8", "[2020] Jan 1-Mar 8"),
                    ggtheme = theme_bw(),
                    xlab = "Time (days)",
                    title = "Time to treatment initiation",
                    xlim = c(0,90),
                    ylab = "Probability not treated (%)")


tiff("figures/PRACTICE_aim1_survival90.tiff", res=300, width = 8, height = 7, units = 'in')
splot90
dev.off()

# TTI
splot.surv2 <- ggplot(adat.nowash, aes(time=surv.days90, 
                                            status = surv.ind90, 
                                            fill=sp2, 
                                            color=sp2,
                                       group=sp2)) +
  geom_km() + 
  geom_kmband() +
  ylim(0,1) +
  xlab("Time (days)") +
  ylab("Probability not treated")

splot.surv2$labels$group <- "Study Period"
splot.surv2$labels$colour <- "Study Period"
splot.surv2$labels$fill <- "Study Period"

tiff("figures/PRACTICE_aim1_survival90.tiff", res=300, width = 6, height = 4, units = 'in')
splot.surv2
dev.off()

# TTI no 90 day censoring
splot.surv3 <- ggplot(adat.nowash, aes(time=surv.days, 
                                       status = surv.ind, 
                                       #fill=sp2, 
                                       color=sp2,
                                       group=sp2)) +
  geom_km() + 
  ylim(0,1) +
  xlab("Time (days)") +
  ylab("Probability not treated") +
  geom_vline(xintercept=90) +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(0), max(690), by=90))

splot.surv3$labels$group <- "Study Period"
splot.surv3$labels$colour <- "Study Period"
tiff("figures/RR/rr-km-surv.tiff", res=150, width = 8, height = 7, units = 'in')
splot.surv3
dev.off()


tiff("figures/RR/rr-surv-hist.tiff", res=150, width = 8, height = 7, units = 'in')
ggplot(adat.nowash, aes(x=surv.days, fill=surv.ind.f)) +
  geom_histogram() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(0), max(690), by=90)) +
  xlab("Days") +
  ylab("Count") +
  labs(fill="Survival\nIndicator")
dev.off()

# TTI by cancer
splot.surv.facet <- ggplot(adat.nowash, aes(time=surv.days90, 
                                            status = surv.ind90, 
                                            fill=sp2, 
                                            color=sp2,
                                            group=sp2)) +
  geom_km() + 
  geom_kmband() +
  #geom_kmticks() +
  facet_wrap(~cancer) + 
  ylim(0,1) +
  xlab("Time (days)") +
  ylab("Probability not treated")

splot.surv.facet$labels$group <- "Study Period"
splot.surv.facet$labels$colour <- "Study Period"
splot.surv.facet$labels$fill <- "Study Period"

tiff("figures/PRACTICE_aim1_survival90_facet.tiff", res=300, width = 8, height = 7, units = 'in')
splot.surv.facet
dev.off()

# Stats overall survival #######################################################
osurv.fit <- survfit(Surv(osurv.months5, osurv.ind5) ~ studyperiod, data=adat.nowash)

splot.osurv <- ggsurvplot(osurv.fit, data=adat.nowash,
                    conf.int = T,
                    risk.table = T,
                    risk.table.height = 0.25,
                    legend.labs = c("[2019] Apr 8-Jul 31", "[2020] Apr 8-Jul 31", 
                                    "[2019] Jan 1-Mar 8", "[2020] Jan 1-Mar 8"),
                    xlab = "Time(months) - Events at T = 0 ocurred the same month as diagnosis",
                    title = "Overall survival, subjects right censored at 5 months",
                    ggtheme = theme_bw())

tiff("figures/PRACTICE_aim1_overallsurvival5.tiff", res=300, width = 8, height = 7, units = 'in')
splot.osurv
dev.off()

# Facet by cancer
splot.osurv.facet <- ggplot(adat.nowash, aes(time=osurv.months5, status = osurv.ind5, 
                        fill=studyperiod, color=studyperiod,group=studyperiod)) +
  geom_km() + 
  geom_kmband() +
  #geom_kmticks() +
  facet_wrap(~cancer) + 
   ylim(0,1) +
  xlab("time(months)") +
  ylab("overall survival")

tiff("figures/PRACTICE_aim1_overallsurvival5_facet.tiff", res=300, width = 8, height = 7, units = 'in')
splot.osurv.facet
dev.off()

# Stats progression free survival: breast, NSCLC and RCC only
adat.pfsurv <- adat.nowash %>% filter(cancer %in% c("bre", "nsc", "rcc"))
pfsurv.fit <- survfit(Surv(adat.pfsurv$pfsurv.days90, adat.pfsurv$pfsurv.ind90) ~ adat.pfsurv$studyperiod)

splot.pfsurv <- ggsurvplot(pfsurv.fit, data=adat.pfsurv,
                          conf.int = T,
                          risk.table = T,
                          risk.table.height = 0.25,
                          legend.labs = c("[2019] Apr 8-Jul 31", "[2020] Apr 8-Jul 31", 
                                          "[2019] Jan 1-Mar 8", "[2020] Jan 1-Mar 8"),
                          xlab = "Time(days)",
                          title = "Progression free survival (Breast, NSCLC, and RCC only), subjects right censored at 90 days",
                          ggtheme = theme_bw())

tiff("figures/PRACTICE_aim1_progfreesurvival90.tiff", res=300, width = 8, height = 7, units = 'in')
splot.pfsurv
dev.off()

# Facet by cancer
splot.pfsurv.facet <- ggplot(adat.pfsurv, aes(time=pfsurv.days90, status = pfsurv.ind90, 
                                            fill=studyperiod, color=studyperiod,group=studyperiod)) +
  geom_km() + 
  geom_kmband() +
  #geom_kmticks() +
  facet_wrap(~cancer) + 
  ylim(0,1) +
  xlab("time(days)") +
  ylab("progression free survival from diagnosis")

tiff("figures/PRACTICE_aim1_progfreesurvival90_facet.tiff", res=300, width = 8, height = 2.5, units = 'in')
splot.pfsurv.facet
dev.off()

# Stats progression free survival: breast, NSCLC and RCC only. Start time at treatment
adat.pfsurv.trt <- adat.nowash %>% filter(cancer %in% c("bre", "nsc", "rcc"), pfsurv.trt.days90 >= 0)
pfsurv.trt.fit <- survfit(Surv(adat.pfsurv.trt$pfsurv.trt.days90, 
                               adat.pfsurv.trt$pfsurv.trt.ind90) ~ 
                            adat.pfsurv.trt$studyperiod)

splot.pfsurv.trt <- ggsurvplot(pfsurv.trt.fit, data=adat.pfsurv.trt,
                           conf.int = T,
                           risk.table = T,
                           risk.table.height = 0.25,
                           legend.labs = c("[2019] Apr 8-Jul 31", "[2020] Apr 8-Jul 31", 
                                           "[2019] Jan 1-Mar 8", "[2020] Jan 1-Mar 8"),
                           ylab = "PFS",
                           xlab = "Time(days), T=0 at treatment start",
                           title = "PFS (Breast, NSCLC, and RCC only) subjects right censored at 90 days.",
                           ggtheme = theme_bw())

tiff("figures/PRACTICE_aim1_progfreesurvival90_trt.tiff", res=300, width = 8, height = 7, units = 'in')
splot.pfsurv.trt
dev.off()

# Facet by cancer
splot.pfsurv.trt.facet <- ggplot(adat.pfsurv.trt, aes(time=pfsurv.trt.days90, status = pfsurv.trt.ind90, 
                                              fill=studyperiod, color=studyperiod,group=studyperiod)) +
  geom_km() + 
  geom_kmband() +
  #geom_kmticks() +
  facet_wrap(~cancer) + 
  ylim(0,1) +
  xlab("time(days), T=0 at treatment start") +
  ylab("PFS ")

tiff("figures/PRACTICE_aim1_progfreesurvival90_trt_facet.tiff", res=300, width = 8, height = 2.5, units = 'in')
splot.pfsurv.trt.facet
dev.off()


# Therapy classifier tables ####################################################
trt.classify <- 
  rbind(read_csv("/data/therapy-classifiers/bre_first-line_therapy_2021-08-18.csv") %>%
          mutate(cancer="Breast"), #%>% rename(targeted.yes = 'targeted or non-HR positive.yes'),
        read_csv("/data/therapy-classifiers/nsc_first-line_therapy_2021-08-18.csv") %>% 
          mutate(cancer="NSCLC"),
        read_csv("/data/therapy-classifiers/ucc_first-line_therapy_2021-08-18.csv") %>% 
          mutate(cancer="UCC"),
        read_csv("/data/therapy-classifiers/pro_first-line_therapy_2021-08-18.csv") %>% 
          mutate(cancer="Prostate"),
        read_csv("/data/therapy-classifiers/rcc_first-line_therapy_2021-03-30.csv") %>% 
          mutate(cancer="RCC"),
        read_csv("/data/therapy-classifiers/crc_first-line_therapy_2021-08-18.csv") %>% 
          mutate(cancer="CRC"),
        read_csv("/data/therapy-classifiers/pan_first-line_therapy_2021-03-30.csv") %>% 
          mutate(cancer="Pancreatic"))

trt.classify2 <- trt.classify %>%
  mutate(clinical = ifelse(is.na(cs.drug) | cs.drug == 0, "No", "Yes"), 
         nccn = ifelse(clinical == "No", 
                       ifelse(is.na(nccn.no) | nccn.no == 0, "Yes", "No"), NA),
         targeted = ifelse(nccn == "Yes",
                           ifelse(is.na(targeted.yes) | targeted.yes == 0, "No", "Yes"), NA),
         myelo = ifelse(targeted == "No",
                        ifelse(is.na(myelo.yes) | myelo.yes == 0, "No", "Yes"), NA)
         ) %>% select(cancer, drug.name, clinical, nccn, targeted, myelo)

write.csv(trt.classify2, "tables/intermediary/trt-classify-tab.csv", na = "")

# Effect forest plots ##########################################################
a1.canc <- read_csv("tables/intermediary/aim1_pp_cancer-interaction-full-p.csv") %>%
  mutate(intpar = ifelse(is.na(ccrc), "Combined", 
                               ifelse(ccrc == 1, "Colorectal",
                               ifelse(cnsc == 1, "NSCLC",
                               ifelse(cpan == 1, "Pancreatic",
                               ifelse(cpro == 1, "Prostate",
                               ifelse(crcc==1, "RCC",
                               ifelse(cucc==1, "UCC", "Breast"))))))),
         intpar = factor(intpar, levels=rev(c("Combined",
                                          "Breast",
                                          "Colorectal",
                                          "NSCLC",
                                          "Pancreatic",
                                          "Prostate",
                                          "RCC",
                                          "UCC"))),
         model = "Cancer"
         ) %>%
  select(model, intpar, est.dd, lci.dd, uci.dd)

a1.age <- read_csv("tables/intermediary/aim1_pp_age-interaction-full.csv") %>%
  mutate(intpar = ifelse(is.na(agegt75), "Combined",
                      ifelse(agegt75 == 1, "Age > 75", "Age <= 75")),
         intpar = factor(intpar, levels=rev(c("Combined",
                                              "Age <= 75",
                                              "Age > 75"))),
         model = "Age"
         ) %>%
  filter(intpar != "Combined") %>%
  select(model, intpar, est.dd, lci.dd, uci.dd)

a1.reth <- read_csv("tables/intermediary/aim1_pp_cancer-interaction-reth-full.csv") %>%
  mutate(intpar = ifelse(is.na(rblack), "Combined",
                      ifelse(rblack == 1, "Black",
                      ifelse(rhisp == 1, "Hispanic",
                      ifelse(rother == 1, "Other", "White")))),
         intpar = factor(intpar, levels=rev(c("Combined",
                                        "White",
                                        "Black",
                                        "Hispanic",
                                        "Other"))),
         model = "Race/Ethnicity"
         ) %>%
  filter(intpar != "Combined") %>%
  select(model, intpar, est.dd, lci.dd, uci.dd)

a1.all <- rbind(a1.canc, a1.reth, a1.age) %>%
  add_row(model = "break", intpar = " ") %>%
  add_row(model = "break", intpar = "  ") %>%
  add_row(model = "break", intpar = "   ") %>%
  mutate(
    intpar = ifelse(intpar == "Combined", 1,
                    ifelse(intpar == " ", 2,
                    ifelse(intpar == "Breast", 3,
                    ifelse(intpar == "Colorectal", 4,
                    ifelse(intpar == "NSCLC", 5,
                    ifelse(intpar == "Pancreatic", 6,
                    ifelse(intpar == "Prostate", 7,
                    ifelse(intpar == "RCC", 8,
                    ifelse(intpar == "UCC", 9,
                    ifelse(intpar == "  ", 10,
                    ifelse(intpar == "White", 14,
                    ifelse(intpar == "Black", 11,
                    ifelse(intpar == "Hispanic", 12,
                    ifelse(intpar == "Other", 13,
                    ifelse(intpar == "   ", 15,
                    ifelse(intpar == "Age <= 75", 16,
                    ifelse(intpar == "Age > 75", 17, NA
                           ))))))))))))))))),
      intpar = factor(intpar, levels=rev(1:17), labels=rev(c("Combined (n=14136)",
                                                    " ",
                                                    "Breast (n=1646)",
                                                    "Colorectal (n=2603)",
                                                    "NSCLC (n=5841)",
                                                    "Pancreatic (n=1268)",
                                                    "Prostate (n=1226)",
                                                    "RCC (n=732)",
                                                    "UCC (n=820)",
                                                    "  ",
                                                    "Black (n=1325)",
                                                    "Hispanic (n=732)",
                                                    "Other (n=1963)",
                                                    "White (n=8478)",
                                                    "   ",
                                                    "Age \u2264 75 (n=9857)",
                                                    "Age > 75 (n=4279)"))),
    model = ifelse(intpar == "Combined", "Combined", model)
  ) 

# Plot (No divider lines)
ylabel <- "Adjusted probability TTI <30 days, \n difference-in-differences (percentage points)"
p.all <- 
  ggplot(data=a1.all,
               aes(x=intpar, 
                   y=est.dd, 
                   ymin=lci.dd, 
                   ymax=uci.dd)) +
  geom_hline(yintercept=0, 
             color="#BBBBBB") +
  #scale_x_discrete(breaks=c("Combined", "Breast")) +
 # geom_vline(xintercept=5,
  #           color="black") +
  geom_point(shape=18, 
             size=2.5, na.rm=T) +
  geom_errorbar(aes(width=.3)) +
  
  # geom_vline(xintercept = 16,
  #            color="#000000",
  #            linetype=3
  #            ) +
  # geom_vline(xintercept=3,
  #            color="#000000",
  #            linetype=3) +
  # geom_vline(xintercept = 8,
  #            color="#000000",
  #            linetype=3) +
  
  ylab(ylabel) +
  xlab("") +
  scale_y_continuous(labels = scales::percent) +
  #scale_color_viridis(discrete=T) +
  theme_classic() + 
  labs(title = "Figure 1") +
  theme(legend.position="none", 
        axis.ticks.y = element_blank(),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot") +
  coord_flip()
  
tiff("figures/a1_dd-forest-plots.tiff", res=300, units="in", width = 5, height = 4)
p.all
dev.off()

# Plot (with divider lines)
p.all2 <- 
  ggplot(data=a1.all[a1.all$intpar != " ",],
         aes(x=intpar, 
             y=est.dd, 
             ymin=lci.dd, 
             ymax=uci.dd)) +
  geom_hline(yintercept=0, 
             color="#BBBBBB") +
  #scale_x_discrete(breaks=c("Combined", "Breast")) +
  # geom_vline(xintercept=5,
  #           color="black") +
  geom_point(shape=18, 
             size=2.5, na.rm=T) +
  geom_errorbar(aes(width=.3)) +
  
  geom_vline(xintercept=3,
             color="#000000",
             linetype=3) +
  geom_vline(xintercept = 8,
             color="#000000",
             linetype=3) +
  
  ylab(ylabel) +
  xlab("") +
  scale_y_continuous(labels = scales::percent) +
  #scale_color_viridis(discrete=T) +
  theme_classic() +
  theme(legend.position="none", axis.ticks.y = element_blank()) +
  coord_flip() 

tiff("figures/a1_dd-forest-plots2.tiff", res=300, units="in", width = 5.5, height = 4)
p.all2
dev.off()



# Old plotting approach 

p.canc <- 
  ggplot(data=a1.canc,
         aes(x=intpar, y=est.dd, ymin=lci.dd, ymax=uci.dd)) +
  geom_pointrange(aes(group=intpar)) +
  #geom_hline(aes(fill=cancer), linetype=2) +
  geom_errorbar(aes(ymin=lci.dd, ymax=uci.dd, group=intpar, width=.3)) +
  ylab("") +
  xlab("") +
  ylim(-.6, .6) +
  theme_light() +
  scale_color_viridis(discrete=T) +
  theme(legend.position="none")+
  coord_flip()

p.reth <- 
  ggplot(data=a1.reth,
         aes(x=intpar, y=est.dd, ymin=lci.dd, ymax=uci.dd)) +
  geom_pointrange(aes(group=intpar)) +
  #geom_hline(aes(fill=cancer), linetype=2) +
  geom_errorbar(aes(ymin=lci.dd, ymax=uci.dd, group=intpar, width=.3)) +
  ylab("") +
  xlab("") +
  ylim(-.6, .6) +
  theme_light() +
  scale_color_viridis(discrete=T) +
  theme_minimal()
  theme(legend.position="none")+
  coord_flip()

p.age <- 
  ggplot(data=a1.age,
         aes(x=intpar, y=est.dd, ymin=lci.dd, ymax=uci.dd)) +
  geom_pointrange(aes(group=intpar)) +
  #geom_hline(aes(fill=cancer), linetype=2) +
  geom_errorbar(aes(ymin=lci.dd, ymax=uci.dd, group=intpar, width=.3)) +
  ylab("Estimated difference in differences") +
  xlab("") +
  ylim(-.6, .6) +
  theme_minimal() +
  scale_color_viridis(discrete=T) +
  theme(legend.position="none")+
  coord_flip()
  
p.all <- ggarrange(p.canc, p.reth, p.age,
            heights = c(1,6/8,4/8),
            nrow=3)
tiff("figures/a1_dd-forest-plots.tiff", res=300, units="in", width = 4, height = 7)
p.all
dev.off()

# Center tables ################################################################
library(janitor)
sitetab <- tabyl(adat.nowash, practiceid, studyperiod, sort=T) %>%
  adorn_percentages("col") %>% adorn_pct_formatting(digits=5) %>% 
  adorn_ns(position = "front")
  

ct <- table(adat.nowash$practiceid, adat.nowash$studyperiod) %>% as.matrix.data.frame()
ctt <- table(adat.nowash$practiceid) %>% as.vector()
cp <- prop.table(ct, margin = 2) %>% as.matrix.data.frame()

sitetab <- cbind(ct[,1], cp[,1], ct[,2], cp[,2], ct[,3], cp[,3], ct[,4], cp[,4]) %>% as.data.frame()
sitetab$total <- ctt
write.csv(sitetab, "tables/sitetab.csv")
#names(sitetab) <- c("Apr 8-Jul 31, 2019", "Apr 8-Jul 31, 2020", "Jan 1-Mar 8, 2019", "Jan 1-Mar 8, 2020")


# TTI histogram? ###############################################################


