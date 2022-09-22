# Aim 1 analytic code
# Author: Daniel Vader

#.libPaths(c(.libPaths(), "rlibs"))

library(tidyverse)
library(mice)
library(survival)
library(splines)
library(ggplot2)
library(survminer)

# Read in analytic data set
adat <- readRDS("../../data/analytic/A1_PRACTICE-imputed-data-derv_2021-08-23.rds")
  

################################################################################
# Fit models
################################################################################
# Unadjusted model
# Fit model w/ no cancer*period interaction
umod <- with(adat, 
               survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                 studypf, 
                               robust=T, cluster=practiceid,
                               ties = "breslow",
                               x = T
               )
)

saveRDS(umod, "models/aim1_unadjusted.rds")

# Fit model w/ no cancer*period interaction
m.nint <- with(adat, 
          survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                              rblack + rhisp + rother + # rwhite ref
                              age + gendm + 
                              igov + iother + # commercial ref 
                              ecog.c + opioid_adv + 
                              ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                              diagnosisday + studypf, 
                            robust=T, cluster=practiceid,
                            ties = "breslow",
                            x = T
                          )
          )

saveRDS(m.nint, "models/aim1_no-cancer-interaction.rds")

# Plot martingale residuals for continuous data
ra <- c(sapply(m.nint[[4]], FUN=residuals, type="martingale", simplify = T))
adat.l <- mice::complete(adat, "long", include=F)
adat.l$m_resid <- ra

ra2 <- c(sapply(m.nint[[4]], FUN=residuals, simplify = F))
adat.l$df_resid <- ra2

resplot.age <- ggplot(data=adat.l, aes(x=age, y=m_resid)) +
  geom_point(shape = ".") +
  #geom_smooth(method="loess") + 
  xlab("Age") +
  ylab("Martingale Residuals") +
  theme_bw()

resplot.day <- ggplot(data=adat.l, aes(x=diagnosisday, y=m_resid)) +
  geom_point(shape = ".") +
  #geom_smooth(method="loess") + 
  xlab("Diagnosis Day") +
  ylab("Martingale Residuals") + 
  theme_bw()

ggcoxzph(cox.zph(m.nint[[4]][[1]]))

resplot.all <- ggplot(data=adat.l, aes(y=df_resid)) +
  geom_point(shape = ".") +
  #geom_smooth(method="loess") + 
  xlab("Diagnosis Day") +
  ylab("Martingale Residuals") + 
  theme_bw()


tiff("/programs/dan/figures/a1_mart-resid-plots.tiff", res=300, width = 8, height = 4, units = 'in')
ggarrange(resplot.age, resplot.day, nrow=1, ncol=2)
dev.off()

par(mfrow=c(5,5))
plot(cox.zph(m.nint[[4]][[1]]), pch=".")

for(i in 1:10){
  tiff(paste0("/programs/dan/figures/a1_cat-resid-plots_", i, ".tiff"),
       res=300, 
       width = 8, height = 12, units = 'in')
  par(mfrow=c(5,4))
  plot(cox.zph(m.nint[[4]][[i]]), var=1, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=2, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=3, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=4, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=5, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=6, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=7, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=8, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=9, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=10, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=11, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=12, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=13, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=14, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=15, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=16, pch=".")
  plot(cox.zph(m.nint[[4]][[i]]), var=17, pch=".")
  dev.off()
}


# Fit model w/ cancer*period interaction
m.wint <- with(adat, 
                 survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                     rblack + rhisp + rother + # rwhite ref
                                     age + gendm + 
                                     igov + iother + # commercial ref 
                                     ecog.c + opioid_adv + 
                                     ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                                     diagnosisday + studypf +
                                   
                                     ccrc*studypf + 
                                     cnsc*studypf + 
                                     cpan*studypf +
                                     cpro*studypf + 
                                     crcc*studypf + 
                                     cucc*studypf
                                 , 
                                 robust=T, cluster=practiceid,
                                 ties = "breslow",
                                 x = T
                 )
)
saveRDS(m.wint, "models/aim1_cancer-interaction.rds")

# Compare
pool.compare(m.wint, m.nint)

for(i in 1:10){
  tiff(paste0("/programs/dan/figures/a1_cat-resid-plots-int_", i, ".tiff"),
       res=300, 
       width = 8, height = 12, units = 'in')
  par(mfrow=c(4,4))
  plot(cox.zph(m.wint[[4]][[i]]), var=9, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=10, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=11, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=12, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=13, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=14, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=15, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=17, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=18, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=19, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=20, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=21, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=22, pch=".")
  plot(cox.zph(m.wint[[4]][[i]]), var=23, pch=".")
  dev.off()
}

# Fit model w/ raceeth*period interaction
m.wint.reth <- with(adat, 
               survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                 rblack + rhisp + rother + # rwhite ref
                                 age + gendm + 
                                 igov + iother + # commercial ref 
                                 ecog.c + opioid_adv + 
                                 ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                                 diagnosisday + studypf +

                                 rblack*studypf +
                                 rhisp*studypf +
                                 rother*studypf
                               , 
                               robust=T, cluster=practiceid,
                               ties = "breslow",
                               x = T
               )
)
saveRDS(m.wint.reth, "models/aim1_cancer-interaction-reth.rds")

# Compare
pool.compare(m.wint.reth, m.nint)

# Fit model w/ denovo metastitic*period interaction
m.wint.denov <- with(adat, 
                    survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                      rblack + rhisp + rother + # rwhite ref
                                      age + gendm + 
                                      igov + iother + # commercial ref 
                                      ecog.c + opioid_adv + 
                                      ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                                      diagnosisday + denovo_met + studypf + 
                                      
                                      denovo_met*studypf
                                    , 
                                    robust=T, cluster=practiceid,
                                    ties = "breslow",
                                    x = T
                    )
)
saveRDS(m.wint.denov, "models/aim1_cancer-interaction-denov.rds")

# Fit model w/ insurance*period interaction
m.wint.insur <- with(adat, 
                    survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                      rblack + rhisp + rother + # rwhite ref
                                      age + gendm + 
                                      igov + iother + # commercial ref 
                                      ecog.c + opioid_adv + 
                                      ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                                      diagnosisday + studypf + denovo_met +
                                      
                                      igov*studypf + iother*studypf
                                    , 
                                    robust=T, cluster=practiceid,
                                    ties = "breslow",
                                    x = T
                    )
)
saveRDS(m.wint.insur, "models/aim1_cancer-interaction-insur.rds")

# Compare
pool.compare(m.wint.insur, m.nint)

# Fit model age*period interaction
m.wint.age <- with(adat, 
                     survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                       rblack + rhisp + rother + # rwhite ref
                                       age + gendm + 
                                       igov + iother + # commercial ref 
                                       ecog.c + opioid_adv + 
                                       ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                                       diagnosisday + studypf +
                                       
                                       age*studypf
                                     , 
                                     robust=T, cluster=practiceid,
                                     ties = "breslow",
                                     x = T
                     )
)
saveRDS(m.wint.age, "models/aim1_age-interaction.rds")

# Compare
pool.compare(m.wint.age, m.nint)


################################################################################
# Export model fit statistics
################################################################################
m.nint <- readRDS("models/aim1_no-cancer-interaction.rds")
m.wint <- readRDS("models/aim1_cancer-interaction.rds")

# Get and store fit statistics function
g.aicbic <- function(m){
  a <- m$analyses
  t <- tibble(aic=numeric(), bic=numeric())
  
  for(i in 1:length(a)){
    aic <- AIC(a[[i]])
    bic <- BIC(a[[i]])
    t <- t %>% add_row(aic=aic, bic=bic)
  }
  return(t)
}

f.wint <- g.aicbic(m.wint)
f.nint <- g.aicbic(m.nint)

ftable <- bind_cols(f.nint, f.wint, f.nint-f.wint)
colnames(ftable) <- c("AIC:no cancer interaction", "BIC:no cancer interaction",
                        "AIC:cancer interaction", "BIC:cancer interaction",
                        "AIC:noi - AIC:yesi", "BIC:noi - BIC:yesi")

write.csv(ftable, "tables/aim1_fitstats.csv")



################################################################################
# Tabulate estimates
################################################################################
# Table that reports the hazard ratios for time periods compared to ref

# Unadjusted model
p.umod <- pool(umod) # pool results
est.umod <- summary(p.umod, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.nint.fm, "tables/intermediary/aim1_unadjusted_2021-10-11.csv")

# No cancer interaction model
p.nint <- pool(m.nint) # pool results
est.nint.fm <- summary(p.nint, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.nint.fm, "tables/intermediary/aim1_no-cancer-interaction_fm_2021-08-23.csv")

# Cancer interaction model
p.wint <- pool(m.wint) # pool results
est.wint.fm <- summary(p.wint, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.wint.fm, "tables/intermediary/aim1_cancer-interaction_fm_2021-08-23.csv")

# Race/ethnicity interaction model
p.wint.reth <- pool(m.wint.reth) # pool results
est.wint.reth.fm <- summary(p.wint.reth, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.wint.reth.fm, "tables/intermediary/aim1_reth-interaction_fm_2021-08-23.csv")

# Denovo metastatic interaction model
p.wint.denov <- pool(m.wint.denov) # pool results
est.wint.denov.fm <- summary(p.wint.denov, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.wint.denov.fm, "tables/intermediary/aim1_denov-interaction_fm_2021-08-23.csv")

# Insurance interaction model
p.wint.insur <- pool(m.wint.insur) # pool results
est.wint.insur.fm <- summary(p.wint.insur, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.wint.insur.fm, "tables/intermediary/aim1_insur-interaction_fm_2021-08-23.csv")

# Age interaction model
p.wint.age <- pool(m.wint.age) # pool results
est.wint.age.fm <- summary(p.wint.age, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.wint.age.fm, "tables/intermediary/aim1_age-interaction_fm_2021-08-23.csv")



################################################################################
# Grab predicted probabilities of survival at 30 days by: 
   # year, period, and cancer (only in cancer interaction model)
################################################################################
m.umod <- readRDS("models/aim1_unadjusted.rds")
m.nint <- readRDS("models/aim1_no-cancer-interaction.rds")
m.wint <- readRDS("models/aim1_cancer-interaction.rds")
m.wint.reth <- readRDS("models/aim1_cancer-interaction-reth.rds")
m.wint.denov <- readRDS("models/aim1_cancer-interaction-denov.rds")
m.wint.insur <- readRDS("models/aim1_cancer-interaction-insur.rds")
m.wint.age <- readRDS("models/aim1_age-interaction.rds")

adat <- readRDS("../../data/analytic/A1_PRACTICE-imputed-data-derv_2021-08-23.rds") 

source("1-5_functions.R")

# No interaction model estimates ###############################################
pp.umod <- survpp(m.umod, adat, cp=F)
#saveRDS(pp.nint, "models/aim1_pp_no-cancer-interaction-marginal.rds")
write.csv(pp.umod, "tables/intermediary/aim1_pp_unadjusted.csv")

#pp.nint.c <- survpp(m.nint, adat, cp = T)
#write.csv(pp.nint.c, "tables/intermediary/aim1_pp_no-cancer-interaction-conditional.csv")

# Unadjusted model estimates
pp.wint <- survpp(umod, adat)

# Primary model estimates (cancer interaction) #################################
pp.wint <- survpp(m.wint, adat) # Marginal
#saveRDS(pp.wint, "models/aim1_pp_cancer-interaction-marginal2.rds")
write.csv(pp.wint, "tables/intermediary/aim1_pp_cancer-interaction-marginal2.csv")

pp.wint.c <- survpp(m.wint, adat, cp = T) # Conditional on cancer
#saveRDS(pp.wint.c, "models/aim1_pp_cancer-interaction-conditional2.rds")
write.csv(pp.wint.c, "tables/intermediary/aim1_pp_cancer-interaction-conditional2.csv")

pp.wint2 <- bind_rows(pp.wint, pp.wint.c) # Join marginal and conditional
#saveRDS(pp.wint2, "models/aim1_pp_cancer-interaction-full.rds")
write.csv(pp.wint2, "tables/intermediary/aim1_pp_cancer-interaction-full-p.csv")


# Race / ethnicity estimates ###################################################
pp.wint.reth <- survpp(m.wint.reth, adat) # Marginal
pp.wint.reth.c <- survpp(m.wint.reth, adat, cp = T, cvar="reth") # Conditional on race/eth

pp.wint.reth2 <- bind_rows(pp.wint.reth, pp.wint.reth.c) # Join marginal and conditional
write.csv(pp.wint.reth2, "tables/intermediary/aim1_pp_cancer-interaction-reth-full.csv")


# Denovo metastatic estimates ##################################################
pp.wint.denov <- survpp(m.wint.denov, adat) # Marginal
pp.wint.denov.c <- survpp(m.wint.denov, adat, cp = T, cvar="denov") # Conditional on denovo/met

pp.wint.denov2 <- bind_rows(pp.wint.denov, pp.wint.denov.c) # Join marginal and conditional
write.csv(pp.wint.denov2, "tables/intermediary/aim1_pp_cancer-interaction-denov-full.csv")


# Insurance estimates ##########################################################
pp.wint.insur <- survpp(m.wint.insur, adat) # Marginal
pp.wint.insur.c <- survpp(m.wint.insur, adat, cp = T, cvar="insur") # Conditional on insurance

pp.wint.insur2 <- bind_rows(pp.wint.insur, pp.wint.insur.c) # Join marginal and conditional
write.csv(pp.wint.insur2, "tables/intermediary/aim1_pp_cancer-interaction-insur-full.csv")


# Age estimates ################################################################
pp.wint.age <- survpp(m.wint.age, adat) # Marginal
pp.wint.age.c <- survpp(m.wint.age, adat, cp = T, cvar="age") # Conditional on age

pp.wint.age2 <- bind_rows(pp.wint.age, pp.wint.age.c) # Join marginal and conditional
write.csv(pp.wint.age2, "tables/intermediary/aim1_pp_age-interaction-full.csv")









################################################################################

pp.wint.plotd <- tribble(
  ~year, ~period, ~est, ~lci, ~uci,
  0, 0, pp.wint$est.janmar2019, pp.wint$lci.janmar2019, pp.wint$uci.janmar2019,
  0, 1, pp.wint$est.aprjul2019, pp.wint$lci.aprjul2019, pp.wint$uci.aprjul2019,
  1, 0, pp.wint$est.janmar2020, pp.wint$lci.janmar2020, pp.wint$uci.janmar2020,
  1, 1, pp.wint$est.aprjul2020, pp.wint$lci.aprjul2020, pp.wint$uci.aprjul2020
  ) %>% 
  mutate(Year = factor(year, levels=c(0,1), labels=c("2019", "2020")),
         Period = factor(period, levels=c(0,1), labels=c("Jan:Mar", "Apr:Jul")))



p1 <- ggplot(data=pp.wint.plotd, aes(x=Period, y = est, ymin = lci, ymax = uci)) +
  geom_pointrange(aes(col=Period)) +
  geom_errorbar(aes(ymin=lci, ymax=uci,col=Period),width=0.2,cex=.5)+
  facet_wrap(~Year, strip.position="left",nrow=2,scales = "free_y")  +
  ylab("Probability") +
  xlab("") +
  ggtitle("Probability of treatment within 30 days") +
  #ylim(.5,.7) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle=180)) +
  coord_flip()

tiff("figures/PRACTICE_aim1_cprob.tiff",res = 300, width=5, height=4, units="in")
p1
dev.off()



