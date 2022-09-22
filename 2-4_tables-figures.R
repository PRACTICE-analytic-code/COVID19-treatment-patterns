# Aim 2 Tables and Figures
# Author: Daniel Vader
# Date created: 02/16/2021

library(tidyverse)
library(stddiff)

adat <- readRDS("/data/analytic/A2_PRACTICE-analytic-data_2021-08-23.rds") %>%
  filter(studyperiod != "Washout 2020", 
         studyperiod != "Washout 2019"
         ) %>%
  mutate(postcovid = ifelse(precovid == "pre", 0, 1),
         ecog.c = ifelse(ecog_b30a14 >= 2, 1, 0)) %>% 
  as.data.frame()

################################################################################
# Descriptive tables (todo add categorical time variable) 
################################################################################

# function for building a table across study time periods
  # dat: dataset
  # cvars: vector of categorical variable names
  # nvars: vector of numeric variable names
  # una: option for tabulating with or without NAs. Uses table() options.
basetab <- function(dat, cvars, nvars, una="ifany"){
  
  # Create tibble to store tabulations
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
                
                smd2019 = numeric(), smd2020 = numeric())
  
  # Add categorical variables
  for(i in 1:length(tabvars.c)){
    tn <- table(dplyr::pull(dat, cvars[i]), dplyr::pull(dat, studyperiod), useNA=una)
    tnm <- rownames(tn)
    tp <- prop.table(tn, margin=2)
    pval <- chisq.test(tn)
    pval <- pval$p.value
    
    tnt <- table(dplyr::pull(dat, cvars[i]),useNA=una)
    tpt <- prop.table(tnt)
    
    # Calculate and extract standardized mean difference by diagnosis year
    smd19 <- stddiff.category(dat[dat$diagnosisyear == "2019",], 
                              gcol=which((colnames(dat) == "postcovid")), 
                              vcol=which((colnames(dat) == cvars[i])))[1,5]
    smd20 <- stddiff.category(dat[dat$diagnosisyear == "2020",], 
                              gcol=which((colnames(dat) == "postcovid")), 
                              vcol=which((colnames(dat) == cvars[i])))[1,5]
    
    for(j in 1:length(tnm)){
      tab <- tab %>% add_row(
        name = paste(cvars[i], ": ", tnm[[j]], sep=""),
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
        smd2020=smd20
      )
    }
  }
  
  # Add numeric variables
  for(i in 1:length(tabvars.n)){
    tm <- aggregate(x=dat[,nvars[i]], by=list(dat$studyperiod), median, na.rm=T) 
    tsd25 <- aggregate(x=dat[,nvars[i]], by=list(dat$studyperiod), quantile, .25, na.rm=T)
    tsd75 <- aggregate(x=dat[,nvars[i]], by=list(dat$studyperiod), quantile, .75, na.rm=T) 
    pval <- aov(dplyr::pull(dat, nvars[i]) ~ dplyr::pull(dat, studyperiod))
    pval <- summary(pval)[[1]][["Pr(>F)"]][[1]]
    
    om <- median(pull(dat, nvars[i]), na.rm=T)
    osd25 <- quantile(pull(dat, nvars[i]), .25, na.rm=T)
    osd75 <- quantile(pull(dat, nvars[i]), .75, na.rm=T)
    
    # Calculate and extract standardized mean difference by diagnosis year
    smd19 <- stddiff.numeric(dat[dat$diagnosisyear == "2019",], 
                             gcol=which((colnames(dat) == "postcovid")), 
                             vcol=which((colnames(dat) == nvars[i])))[7]
    smd20 <- stddiff.numeric(dat[dat$diagnosisyear == "2020",], 
                             gcol=which((colnames(dat) == "postcovid")), 
                             vcol=which((colnames(dat) == nvars)))[7]
    
    tab <- tab %>% add_row(
      name = paste(nvars[i], "(mean, SD)"),
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
      smd2020=smd20
    )
  }
  # Calc total row
  tn <- table(dplyr::pull(adat, studyperiod), useNA=una)
  tnm <- "Total"
  tp <- prop.table(tn)
  pval <- NA
  
  tab <- tab %>% add_row(
    name = tnm,
    n2019pre = tn[3],
    prop2019pre = tp[3],
    n2019post = tn[1],
    prop2019post = tp[1],
    n2020pre = tn[4],
    prop2020pre = tp[4],
    n2020post = tn[2],
    prop2020post = tp[2],
    ntot = nrow(dat)
  )
  return(tab)
}

tabvars.c <- c("gend", "reth", "insclass", "cancer", "practicetype", "ecog.c", 
               "denovo_met", "liverf_adv", "kidneyf_adv",
               "steroid_adv", "opioid_adv")#, "myelo.yes") 
# and comorbidities (liver, renal, steroid use, opioid use) once coded
tabvars.n <- c("age")

tab <- basetab(adat, tabvars.c, tabvars.n)

tab.rmna <- basetab(adat, tabvars.c, tabvars.n, "no")

write.csv(tab, "tables/intermediary/2-1_baseline-table.csv", na = "")
write.csv(tab.rmna, "tables/intermediary/2-2_baseline-table_nona.csv", na = "")


# Effect forest plots ##########################################################
a2.canc <- read_csv("tables/intermediary/aim2_pp_cancer-interaction-full_2021-08-23.csv") %>%
  mutate(intpar = ifelse(is.na(cnsc), "Combined", 
                                ifelse(cnsc == 1, "NSCLC",
                                ifelse(cpro == 1, "Prostate",
                                ifelse(cucc==1, "UCC", "Breast")))),
         intpar = factor(intpar, levels=rev(c("Combined",
                                              "Breast",
                                              "NSCLC",
                                              "Prostate",
                                              "UCC"))),
         model = "Cancer"
  ) %>%
  select(model, intpar, est.dd, lci.dd, uci.dd)

a2.age <- read_csv("tables/intermediary/aim2_pp_age-interaction-full_2021-08-23.csv") %>%
  mutate(intpar = ifelse(is.na(agegt75), "Combined",
                         ifelse(agegt75 == 1, "Age > 75", "Age <= 75")),
         intpar = factor(intpar, levels=rev(c("Combined",
                                              "Age <= 75",
                                              "Age > 75"))),
         model = "Age"
  ) %>%
  filter(intpar != "Combined") %>%
  select(model, intpar, est.dd, lci.dd, uci.dd)

a2.reth <- read_csv("tables/intermediary/aim2_pp_reth-interaction-full_2021-08-23.csv") %>%
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

a2.all <- rbind(a2.canc, a2.reth, a2.age) %>%
  add_row(model = "break", intpar = " ") %>%
  add_row(model = "break", intpar = "  ") %>%
  add_row(model = "break", intpar = "   ") %>%
  mutate(
    intpar = ifelse(intpar == "Combined", 1,
                    ifelse(intpar == " ", 2,
                    ifelse(intpar == "Breast", 3,
                    ifelse(intpar == "NSCLC", 4,
                    ifelse(intpar == "Prostate", 5,
                    ifelse(intpar == "UCC", 6,
                    ifelse(intpar == "  ", 7,
                    ifelse(intpar == "White", 11,
                    ifelse(intpar == "Black", 8,
                    ifelse(intpar == "Hispanic", 9,
                    ifelse(intpar == "Other", 10,
                    ifelse(intpar == "   ", 12,
                    ifelse(intpar == "Age <= 75", 13,
                    ifelse(intpar == "Age > 75", 14, NA
                           )))))))))))))),
    
    intpar = factor(intpar, levels=rev(1:14), labels=rev(c("Combined (n=5962)",
                                                           " ",
                                                           "Breast (n=1101)",
                                                           "NSCLC (n=3518)",
                                                           "Prostate (n=822)",
                                                           "UCC (n=521)",
                                                           "  ",
                                                           "Black (n=532)",
                                                           "Hispanic (n=258)",
                                                           "Other (n=754)",
                                                           "White (n=3742)",
                                                           "   ",
                                                           "Age \u2264 75 (n=4053)",
                                                           "Age > 75 (n=1909)"))),
    model = ifelse(intpar == "Combined", "Combined", model)
  ) 

# Plot (no divider lines)
ylabel <- "Adjusted probability TTI <30 days, \n difference-in-differences (percentage points)"
p2.all <- 
  ggplot(data=a2.all,
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
  ylab("Estimated difference in differences \n (percentage points)") +
  xlab("") +
  scale_y_continuous(labels = scales::percent) +
  #scale_color_viridis(discrete=T) +
  theme_classic() + 
  labs(title = "Figure 2") +
  theme(legend.position="none", 
        axis.ticks.y = element_blank(),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot") +
  coord_flip() 

tiff("figures/a2_dd-forest-plots.tiff", res=300, units="in", width = 4, height = 3)
p2.all
dev.off()

# Plot (with divider lines)
p2.all2 <- 
  ggplot(data=a2.all[a2.all$intpar != " ",],
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
  ylab("Estimated difference in differences \n (percentage points)") +
  xlab("") +
  scale_y_continuous(labels = scales::percent) +
  #scale_color_viridis(discrete=T) +
  theme_classic() +
  theme(legend.position="none", axis.ticks.y = element_blank()) +
  coord_flip() 

tiff("figures/a2_dd-forest-plots2.tiff", res=300, units="in", width = 4, height = 3)
p2.all2
dev.off()


# Old plotting approach
p.canc2 <- 
  ggplot(data=a2.canc,
         aes(x=intpar, y=est.dd, ymin=lci.dd, ymax=uci.dd)) +
  geom_pointrange(aes(group=intpar)) +
  #geom_hline(aes(fill=cancer), linetype=2) +
  geom_errorbar(aes(ymin=lci.dd, ymax=uci.dd, group=intpar, width=.3)) +
  ylab("") +
  xlab("") +
  ylim(-.6, .6) +
  scale_color_viridis(discrete=T) +
  theme(legend.position="none", axis.ticks.y = element_blank())+
  coord_flip()

p.reth2 <- 
  ggplot(data=a2.reth,
         aes(x=intpar, y=est.dd, ymin=lci.dd, ymax=uci.dd)) +
  geom_pointrange(aes(group=intpar)) +
  #geom_hline(aes(fill=cancer), linetype=2) +
  geom_errorbar(aes(ymin=lci.dd, ymax=uci.dd, group=intpar, width=.3)) +
  ylab("") +
  xlab("") +
  ylim(-.6, .6) +
  scale_color_viridis(discrete=T) +
  theme(legend.position="none", axis.ticks.y = element_blank())+
  coord_flip()

p.age2 <- 
  ggplot(data=a2.age,
         aes(x=intpar, y=est.dd, ymin=lci.dd, ymax=uci.dd)) +
  geom_pointrange(aes(group=intpar)) +
  #geom_hline(aes(fill=cancer), linetype=2) +
  geom_errorbar(aes(ymin=lci.dd, ymax=uci.dd, group=intpar, width=.3)) +
  ylab("Estimated difference in differences") +
  xlab("") +
  ylim(-.6, .6) +
  scale_color_viridis(discrete=T) +
  theme(legend.position="none", axis.ticks.y = element_blank())+
  coord_flip()

p.all2 <- ggarrange(p.canc2, p.reth2, p.age2,
                   labels = c("A", 
                              "B",
                              "C"),
                   heights = c(1,1,4/5),
                   nrow=3)
tiff("figures/a2_dd-forest-plots.tiff", res=300, units="in", width = 4, height = 6)
p.all2
dev.off()
