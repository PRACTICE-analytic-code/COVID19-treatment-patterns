# Aim 1 multiple imputation code
# Author: Daniel Vader

# Set path to local library for running extra packages on the cluster
# .libPaths(c(.libPaths(), "rlibs"))

# Load packages
library(tidyverse)
library(mice)
library(micemd)

# Set small number for ln transformation of positive continuous variables.
# td <- .0001

# Read in analytic data set
adat <- readRDS("../../data/analytic/PRACTICE-analytic-data_2021-08-23.rds") %>% 
  #filter(!(studyperiod %in% c("Washout 2019", "Washout 2020"))) %>%
  mutate(patientid = as.factor(patientid),
         practiceid = as.numeric(factor(practiceid)),
         practicetype = as.factor(practicetype),
         reth = as.factor(reth),
         insclass = as.factor(insclass),
         cancer = as.factor(cancer),
         denovo_met = factor(denovo_met, levels=c(0,1), labels=c("no", "yes")),
         #groupstage = as.factor(groupstage),
         studyperiod = as.factor(studyperiod),
         year2020 = ifelse(diagnosisyear == 2020, 1, 0),
         postcovid = ifelse(precovid == "post", 1, 0),
         ecog = ecog_b30a14,
         weight = weight_diag_enh,
         height = height_diag_enh,
         glom = result_glom_diag_adv,
         bili = result_bili_diag_adv,
         alan = result_alan_diag_adv,
         albu = result_albu_diag_adv,
         aspa = result_aspa_diag_adv,
         calc = result_calc_diag_adv,
         crea = result_crea_diag_adv) %>%
  select(patientid, practiceid, practicetype, gend, reth, age, insclass,
         cancer, ecog, denovo_met,
         weight, height, 
         glom, bili, alan, 
         albu, aspa, calc,
         crea,
         steroid_adv, opioid_adv, 
         year2020, postcovid, studyperiod, diagnosisday,
         surv.ind90, surv.days90, osurv.months, osurv.ind) 

# Set seed
#sample(1:100000, 1) # 49647
set.seed(49647)

################################################################################
# Attach the subject-specific Nelson-Aaalen estimator of the cumulative hazard rate
# Re: White & Royston (2009) doi: 10.1002/sim.3618 
################################################################################
#  Time to treatment initiation (TTI)
adat$hrna.tti <- nelsonaalen(as.data.frame(adat), timevar=surv.days90, statusvar=surv.ind90)

# Overall survival (OS)
adat$hrna.os <- nelsonaalen(as.data.frame(adat), timevar=osurv.months, statusvar=osurv.ind)

# Remove vars we won't use anymore
adat2 <- adat %>% select(-osurv.months) 

################################################################################
# Initialize imputation methods vector
################################################################################
# set cluster index (practiceid)
ic <- 2
# get initial cluster vector
meth <- micemd::find.defaultMethod(adat2, ind.clust=ic)

# Replace method for continuous variables with 2-level approach: 2l.norm
meth[c("ecog", "weight", "height", "glom",
       "bili", "alan", "albu",
       "aspa", "calc", 
       "crea")] <- "2l.norm"

# Replace jomo with polytomous regression (not multi-level)
meth[c("reth", "gend")] <- "polyreg"

# Replace 2l.glm.norm with logistic regression
meth["denovo_met"]  <- "logreg"

################################################################################
# Create predictor matrix
################################################################################
# Initialize matrix
pm <- make.predictorMatrix(adat2) 

# Set cluster variable for continuous variables.
pm[ic,ic] <- 0
pm[-ic,ic] <- -2
pm[c("reth", "gend", "denovo_met"),ic] <- 0

# Zero out variable columns to ignore in regression equations
pm[, c("patientid", "surv.days90", "year2020", "postcovid")] <- 0

# Zero out non imputed variable rows in matrix (these variables all have
   # complete data). Note, this 0-ing everything but the variables
   # we want to impute.
pm[-which(rownames(pm) %in% c("reth", "gend", "insclass", "ecog", "denovo_met",
                              "weight", "height",
                              "glom", "bili", "alan", 
                              "albu", "aspa", "calc",
                              "crea")),] <- 0

################################################################################
# Test imputation (for debugging)
################################################################################
#ini <- mice(adat2, predictorMatrix=pm, method = meth, maxit=1, m=1)
#meth <- ini$method
#t <- mice::complete(ini, "long", include = F)


################################################################################
# Impute missing data and save
################################################################################
adat.mi <- mice(adat2, method=meth, predictorMatrix=pm, m=10, maxit=5)
saveRDS(adat.mi, paste("/data/analytic/PRACTICE-imputed-data_", Sys.Date(), ".rds", sep=""))


################################################################################
# Create derived variables and save
################################################################################
# Re-read data after crash
#adat.mi <- readRDS("/data/analytic/PRACTICE-imputed-data_2021-03-31.rds")

adat.mi.l <- mice::complete(adat.mi, "long", include = T) %>%
  mutate(
    agegt75 = ifelse(age > 75, 1, 0),
    rblack = ifelse(reth == "Black", 1, 0),
    rhisp = ifelse(reth == "Hispanic", 1, 0),
    rother = ifelse(reth == "Other", 1, 0),
    rwhite = ifelse(reth == "White", 1, 0),
    gendm = ifelse(gend == "Male", 1, 0),
    icomm = ifelse(insclass == "commercial", 1, 0),
    igov = ifelse(insclass == "gov", 1, 0),
    iother = ifelse(insclass == "other", 1, 0),
    ecog = squeeze(ecog, c(0,5)),
    weight = squeeze(weight, c(min(adat$weight, na.rm=T), max(adat$weight, na.rm=T))),
    height = squeeze(height, c(min(adat$height, na.rm=T), max(adat$height, na.rm=T))),
    glom = squeeze(glom, c(min(adat$glom, na.rm=T), max(adat$glom, na.rm=T))),
    bili = squeeze(bili, c(min(adat$bili, na.rm=T), max(adat$bili, na.rm=T))),
    alan = squeeze(alan, c(min(adat$alan, na.rm=T), max(adat$alan, na.rm=T))),
    albu = squeeze(albu, c(min(adat$albu, na.rm=T), max(adat$albu, na.rm=T))),
    aspa = squeeze(aspa, c(min(adat$aspa, na.rm=T), max(adat$aspa, na.rm=T))),
    calc = squeeze(calc, c(min(adat$calc, na.rm=T), max(adat$calc, na.rm=T))),
    crea = squeeze(crea, c(min(adat$crea, na.rm=T), max(adat$crea, na.rm=T))),
    cbre = ifelse(cancer == "bre", 1, 0),
    ccrc = ifelse(cancer == "crc", 1, 0),
    cnsc = ifelse(cancer == "nsc", 1, 0),
    cpan = ifelse(cancer == "pan", 1, 0),
    cpro = ifelse(cancer == "pro", 1, 0),
    crcc = ifelse(cancer == "rcc", 1, 0),
    cucc = ifelse(cancer == "ucc", 1, 0),
    ecog.c = ifelse(ecog >= 2, 1, 0),
    liverf = ifelse(glom < 15, 1, 0),
    kidneyf = ifelse(bili > 2, 1, 0),
    studypf = ifelse(studyperiod == "Jan 1-Mar 8, 2019", 1, 
                     ifelse(studyperiod == "Apr 8-Jul 31, 2019", 2,
                            ifelse(studyperiod == "Jan 1-Mar 8, 2020", 3, 4))),
    studypf = factor(studypf, levels = c(1,2,3,4), labels=c("janmar2019", "aprjul2019", "janmar2020", "aprjul2020"))
  )

#adat.mi.l$ecog.c <- with(adat.mi.l, ifelse(ecog >= 2, 1, 0))
#adat.mi.l$bmi <- with(adat.mi.l, weight_diag_enh / height_diag_enh / height_diag_enh)
#adat.mi.l$liverf <- with(adat.mi.l, ifelse(glom < 15, 1, 0))
#adat.mi.l$kidneyf <- with(adat.mi.l, ifelse(bili > 2, 1, 0))
#adat.mi.l$insclass2 <- fct_collapse(adat.mi.l$insclass,
#                                    commercial = "commercial",
#                                    gov = "gov",
#                                    other = c("other", "selfpay"))

adat.mi2 <- as.mids(adat.mi.l)

# Save imputed data with derived variables
saveRDS(adat.mi2, paste("/data/analytic/PRACTICE-imputed-data-derv_", Sys.Date(), ".rds", sep=""))

adat.mi3 <- adat.mi.l %>% 
  filter(studyperiod != "Washout 2020", studyperiod != "Washout 2019") %>%
  as.mids()

# Save imputed data with derived variables
saveRDS(adat.mi3, paste("/data/analytic/A1_PRACTICE-imputed-data-derv_", Sys.Date(), ".rds", sep=""))


################################################################################
# Tabulate categorical vars
################################################################################

adat.mi.cont <- adat.mi.l %>% select(.imp, .id, ecog, weight, height, glom, bili, alan, albu,
                                     aspa, calc, crea) %>% as.mids()

tc1_0 <- table(adat.mi2$data$ecog.c)
tc1_0.p <- prop.table(tc1_0)
tc1 <- with(adat.mi2, table(ecog.c))$analyses
tc1.p <- with(adat.mi2, prop.table(table(ecog.c)))$analyses

tc2_0 <- table(adat.mi2$data$liverf)
tc2_0.p <- prop.table(tc2_0)
tc2 <- with(adat.mi2, table(liverf))$analyses
tc2.p <- with(adat.mi2, prop.table(table(liverf)))$analyses

tc3_0 <- table(adat.mi2$data$kidneyf)
tc3_0.p <- prop.table(tc3_0)
tc3 <- with(adat.mi2, table(kidneyf))$analyses
tc3.p <- with(adat.mi2, prop.table(table(kidneyf)))$analyses

tab <- tibble(names = character(), imp = numeric(),
              no = numeric(), nop=numeric(), yes = numeric(), yesp = numeric())

tab <- tab %>% add_row(names = "ecog >= 2", imp = 0,
                   no = tc1_0[1], nop = tc1_0.p[1],
                   yes = tc1_0[2], yesp = tc1_0.p[2])
for(i in 1:10){
  tab <- tab %>% add_row(names = "ecog >= 2", imp = i,
                         no = tc1[[i]][1], nop = tc1.p[[i]][1],
                         yes = tc1[[i]][2], yesp = tc1.p[[i]][2])
}

tab <- tab %>% add_row(names = "liverf", imp = 0,
                       no = tc2_0[1], nop = tc2_0.p[1],
                       yes = tc2_0[2], yesp = tc2_0.p[2])
for(i in 1:10){
  tab <- tab %>% add_row(names = "liverf", imp = i,
                         no = tc2[[i]][1], nop = tc2.p[[i]][1],
                         yes = tc2[[i]][2], yesp = tc2.p[[i]][2])
}

tab <- tab %>% add_row(names = "kidneyf", imp = 0,
                       no = tc3_0[1], nop = tc3_0.p[1],
                       yes = tc3_0[2], yesp = tc3_0.p[2])
for(i in 1:10){
  tab <- tab %>% add_row(names = "kidneyf", imp = i,
                         no = tc3[[i]][1], nop = tc3.p[[i]][1],
                         yes = tc3[[i]][2], yesp = tc3.p[[i]][2])
}

write.csv(tab, "tables/intermediary/PRACTICE_catimp-summary.csv")

################################################################################
# Create diagnositc plots/tables
################################################################################
tiff("figures/PRACTICE_imp-boxwisk-cont.tiff", res=300, width = 6, height = 5, units = 'in')
mice::bwplot(adat.mi.cont, main="Complete case data (all data) vs imputed (NAs only)")
dev.off()

tiff("figures/PRACTICE_imp-dens-cont.tiff", res=300, width = 6, height = 5, units = 'in')
densityplot(adat.mi.cont)
dev.off()
