# PRACTICE: Data step. Import cleaned data and merge.
# Author: Daniel Vader
# Date created: 11/19/2020

library(tidyverse)
library(lubridate)
library(zoo)

# Indicator for whether to run code that outputs tables of first line therapies
   # and drugs for external assessment.
assess.external <- F

# Debug mode off/on (determines whether to rm temp data from env)
data.debug <- T

statadir <- "/data/stata_files/Data022021/"
# Load primary data
ucc.b <- haven::read_dta(paste0(statadir, "Adv_UCC_022021_all_wide_20201216.dta"))
bre.b <- haven::read_dta(paste0(statadir, "Met_Breast_022021_all_wide_20201216.dta"))
crc.b <- haven::read_dta(paste0(statadir, "Met_CRC_022021_all_wide_20201216.dta"))
pan.b <- haven::read_dta(paste0(statadir, "Met_Pancreatic_022021_all_wide_20201216.dta"))
pro.b <- haven::read_dta(paste0(statadir, "Met_Prostate_022021_all_wide_20201216.dta"))
rcc.b <- haven::read_dta(paste0(statadir, "Met_RCC_022021_all_wide_20201216.dta"))
nsc.b <- haven::read_dta(paste0(statadir, "NSCLC_022021_all_wide_20201216.dta"))

# Variable selction function
var.select <- function(dat, cancer){
  # Rename columns. We're interested in metastatic/advanced diagnosis.
  dat <- if(cancer %in% c("nsc", "ucc")){
    dat %>% rename(diagnosisdateold = diagnosisdate,
                   diagnosisdate = advanceddiagnosisdate,
                   weight_diag_enh = weight_advdiag,
                   weightunit_diag_enh = unit_weight_advdiag,
                   height_diag_enh = height_advdiag,
                   heightunit_diag_enh = unit_height_advdiag)
  } else{
    dat %>% rename(diagnosisdateold = diagnosisdate,
                   diagnosisdate = metdiagnosisdate,
                   weight_diag_enh = weight_metdiag,
                   weightunit_diag_enh = unit_weight_metdiag,
                   height_diag_enh = height_metdiag,
                   heightunit_diag_enh = unit_height_metdiag) %>%
      # Replace _met suffix with _adv so that names match across data files.
      rename_with(function(x){
        stringr::str_replace_all(x, "_met", "_adv")
        })
  }
  # Select variables we will be using (not strictly necessary, but nice
     # to filter out the chaff)
  dat %>%
    dplyr::select(patientid, practiceid, 
                  #demo
                  practicetype, birthyear, gender,
                  race, ethnicity, state, 
                  
                  #insurance
                  payercategory_1, # there are 15 payer categories with start and
                  #end dates. May need to load them all for classification.
                  
                  #treatment
                  linename_1, starttherapydate_1, endtherapydate_1, 
                  ismaintenancetherapy_1, 
                  #diagnosis
                  diagnosisdate, diagnosisdateold, 
                  
                  #comorbidities - stage, ECOG
                  # Note: may need ECOG pre/post treatment
                  groupstage, 
                  ecogvalue_postadvdiag, ecogddif_postadvdiag, ecogdate_postadvdiag,
                  ecogvalue_preadvdiag, ecogddif_preadvdiag, ecogdate_preadvdiag,
                  
                  # mortality
                  deathdate, OSurv,
                  bmi_advdiag, bmi_line1,
                  weight_diag_enh, weightunit_diag_enh, height_diag_enh,
                  heightunit_diag_enh, weight_line1, unit_weight_line1,
                  height_line1, unit_height_line1,
                  
                  # censoring (last structured activity)
                  censdate,
                  
                  # Structured activity info for exclusion criteria
                  visits90, ttfirstvisit, firstvisitdate,
                  
                  # comorbidities - drug, liver failure, kidney failure
                  units_glom_diag_adv, result_glom_diag_adv, testdate_glom_diag_adv,
                  ddate_glom_diag_adv, units_glom_line1, testdate_glom_line1,
                  result_glom_line1, result_bili_diag_adv, testdate_bili_diag_adv, 
                  units_bili_diag_adv, result_bili_line1, 
                  testdate_bili_line1, units_bili_line1,
                  
                  # drugs
                  drgenct_painagent_diag, drgenct_steroid_diag, 
                  drgenct_painagent_line1, drgenct_steroid_line1,
                  # Pull specific drug info (not needed)
                  #matches("^drct.*diag$")
                  
                  # Other variables to be used in imputation
                  result_alan_diag_adv, testdate_alan_diag_adv, 
                  units_alan_diag, result_alan_line1,
                  testdate_alan_line1, units_alan_line1,
                  
                  result_albu_diag_adv, testdate_albu_diag_adv,
                  units_albu_diag_adv, result_albu_line1,
                  testdate_albu_line1, units_albu_line1,
                  
                  result_aspa_diag_adv, testdate_aspa_diag_adv,
                  units_aspa_diag_adv, result_aspa_line1,
                  testdate_aspa_line1, units_aspa_line1,
                  
                  result_calc_diag_adv, testdate_calc_diag_adv,
                  units_calc_diag_adv, result_calc_line1,
                  testdate_calc_line1, units_calc_line1,
                  
                  result_crea_diag_adv, testdate_crea_diag_adv,
                  units_crea_diag_adv, result_crea_line1,
                  testdate_crea_line1, units_crea_line1,
                  
                  # Pull insurance info
                  matches("^insure*"), matches("^payer*")
                ) %>%
    # create column indicating cancer type
    dplyr::mutate(cancer=cancer, 
                  censdate=ymd(censdate),
                  linename_1=ifelse(linename_1 == "", NA, linename_1))
}

# Drop unneccessary variables
ucc <- ucc.b %>% var.select("ucc")
bre <- bre.b %>% var.select("bre")
crc <- crc.b %>% var.select("crc")
pan <- pan.b %>% var.select("pan")
pro <- pro.b %>% var.select("pro")
rcc <- rcc.b %>% var.select("rcc")
nsc <- nsc.b %>% var.select("nsc")

# Grab progression free survival data
bre.pfs <- readr::read_csv("/data/original/022021/progression/enhanced_metbreastprogression.csv") %>%
  select(PatientID, ProgressionDate) %>%
  dplyr::rename_all(.funs = tolower) %>%
  group_by(patientid) %>% 
  slice_min(progressiondate, with_ties=F) 
rcc.pfs <- readr::read_csv("/data/original/022021/progression/enhanced_metrccprogression.csv") %>%
  select(PatientID, ProgressionDate) %>%
  dplyr::rename_all(.funs = tolower)%>%
  group_by(patientid) %>% 
  slice_min(progressiondate, with_ties=F) 
nsc.pfs <- readr::read_csv("/data/original/022021/progression/enhanced_advnsclcprogression.csv") %>%
  select(PatientID, ProgressionDate) %>% 
  dplyr::rename_all(.funs = tolower) %>%
  group_by(patientid) %>% 
  slice_min(progressiondate, with_ties=F) 

bre <- left_join(bre, bre.pfs, by="patientid")
rcc <- left_join(rcc, rcc.pfs, by="patientid")
nsc <- left_join(nsc, nsc.pfs, by="patientid")

# Grab Prostate ADT 
pro.adt <- read_csv("/data/original/022021/prostate/enhanced_metpc_adt.csv",
                    col_types = list(col_character(), 
                                     col_character(), 
                                     col_date(format = ""),
                                     col_date(format = ""),
                                     col_character())) %>%
  rename(patientid = PatientID, adt_start = StartDate) %>%
  left_join(pro[,c("patientid", "diagnosisdate")]) %>%
  filter(TreatmentStatus == "Continuing", adt_start >= diagnosisdate) %>% # We only care about records where patients continuing treatment
  group_by(patientid) %>% 
  slice_min(adt_start, with_ties = F) %>% # Keep only first record per subject
  select(patientid, adt_start)

pro <- left_join(pro, pro.adt, by="patientid")

# Merge data. n = 164604 (164188 unique)
adat <- bind_rows(ucc, bre, crc, pan, pro, rcc, nsc) 


# Remove initial data from environment
if(data.debug == F){
  rm(ucc.b, bre.b, crc.b, pan.b, pro.b, rcc.b, nsc.b,
     ucc, bre, crc, pan, pro, rcc, nsc,
     bre.pfs, rcc.pfs, nsc.pfs, pro.adt)
}

# Select cohort ################################################################
# count unique patients function
cpts <- function(){length(unique(adat$patientid))}

# Remove patients in multiple cohorts. N =  163773 N dropped = 415
#adat <- adat %>% dplyr::filter(!(patientid %in% patientid[duplicated(patientid)]))

# Save vector of duplicated ids across cohorts
dup <- adat$patientid[duplicated(adat$patientid)]

# Drop all subjects outside study timeframe. N = 20261; N dropped = 143,927
adat <- adat %>% 
  filter((diagnosisdate >= ymd("2019/01/01") & diagnosisdate <= ymd("2019/07/31")) |
           (diagnosisdate >= ymd("2020/01/01") & diagnosisdate <= ymd("2020/07/31")) )

# Drop all subjects with no follow-up visit within 90 days. N = 18630; 1631 subjects lost.
adat <- adat %>% filter(ttfirstvisit <= 90 )

# Drop all subjects with fewer than 2 follow-up visits within 90 days. 
# N = 17355, 1,275 subjects lost.
adat <- adat %>% filter(visits90 > 1)

# Exclude patients in multiple cohorts. N=17,289. 66 sub lost.
adat <- adat %>% filter(!(patientid %in% dup))

# Drop all subjects where treatment occurs before advanced diagnosis.
# N = 16607. Lost = 682.
adat <- adat %>% filter(diagnosisdate <= starttherapydate_1 | is.na(starttherapydate_1))

# Drop all subjects treated with non-NCCN compliant therapy or clinical study drug
   # N = 19715 N dropped = 471 
trt.classify <- 
  rbind(read_csv("/data/therapy-classifiers/bre_first-line_therapy_2021-08-18.csv") %>%
          mutate(cancer="bre"), #%>% rename(targeted.yes = 'targeted or non-HR positive.yes'),
        read_csv("/data/therapy-classifiers/nsc_first-line_therapy_2021-08-18.csv") %>% 
          mutate(cancer="nsc"),
        read_csv("/data/therapy-classifiers/ucc_first-line_therapy_2021-08-18.csv") %>% 
          mutate(cancer="ucc"),
        read_csv("/data/therapy-classifiers/pro_first-line_therapy_2021-08-18.csv") %>% 
          mutate(cancer="pro"),
        read_csv("/data/therapy-classifiers/rcc_first-line_therapy_2021-03-30.csv") %>% 
          mutate(cancer="rcc"),
        read_csv("/data/therapy-classifiers/crc_first-line_therapy_2021-08-18.csv") %>% 
          mutate(cancer="crc"),
        read_csv("/data/therapy-classifiers/pan_first-line_therapy_2021-03-30.csv") %>% 
          mutate(cancer="pan"))
  
#names(trt.classify) <- c("drug.name", "n", "cs.drug", "nccn.no", "targeted.yes", "myelo.yes",
#                         "cancer")

trt.classify.a1 <- trt.classify %>% 
  filter(cs.drug == 1 | nccn.no == 1) %>% 
  select(drug.name, cancer)

adat <- adat %>% 
  filter(
    !(cancer == "bre" & linename_1 %in% 
        trt.classify.a1$drug.name[trt.classify.a1$cancer == "bre"]),
    !(cancer == "nsc" & linename_1 %in% 
        trt.classify.a1$drug.name[trt.classify.a1$cancer == "nsc"]),
    !(cancer == "ucc" & linename_1 %in% 
        trt.classify.a1$drug.name[trt.classify.a1$cancer == "ucc"]),
    !(cancer == "pro" & linename_1 %in% 
        trt.classify.a1$drug.name[trt.classify.a1$cancer == "pro"]),
    !(cancer == "rcc" & linename_1 %in% 
        trt.classify.a1$drug.name[trt.classify.a1$cancer == "rcc"]),
    !(cancer == "crc" & linename_1 %in% 
        trt.classify.a1$drug.name[trt.classify.a1$cancer == "crc"]),
    !(cancer == "pan" & linename_1 %in% 
        trt.classify.a1$drug.name[trt.classify.a1$cancer == "pan"]))

# Grab myelosuppresive vs non-myelosuppresive classifier
adat <- trt.classify %>% 
  select(drug.name, myelo.yes, targeted.yes, cancer) %>% 
  rename(linename_1 = drug.name) %>%
  mutate(myelo.yes = ifelse(is.na(myelo.yes), 0, myelo.yes),
         targeted.yes = ifelse(is.na(targeted.yes), 0, targeted.yes)) %>%
  right_join(adat, by=c("linename_1", "cancer"))

#t <- adat %>% filter(grepl('clinical study drug', linename_1, ignore.case = T)) %>% select(cancer, linename_1) 


# Create tables for external classification and assessment
if(assess.external){
  # Create table of first line therapies #######################################
  cancers <- unique(adat$cancer)
  for(i in 1:length(cancers)){
    adat$linename_1[adat$cancer == cancers[i]] %>% table() %>% as_tibble() %>%
      readr::write_csv(paste("tables/intermediary/", cancers[i], "_first-line-therapy", 
                             "2021-03-28", ".csv", sep=""))
  }
  
  for(i in 1:length(cancers)){
    
    t <- read_csv(paste0("/data/therapy-classifiers/", cancers[i], "_first-line_therapy_2021-01-18.csv" )) 
    names(t) <- c("drug.name", "n", "cs.drug", "nccn.no", "targeted.yes", "myelo.yes", "new")
    t <- t %>% mutate(new = 0)
    n <- unique(adat$linename_1[adat$cancer  == cancers[i]])
    n0 <- n[!(n %in% t$drug.name)] 
    t1 <- t %>% add_row(drug.name = n0, n=NA, cs.drug=NA, nccn.no=NA, 
                        targeted.yes=NA, myelo.yes=NA, new=1) %>% filter(drug.name != "")
    t1 <- t1[with(t1, order(-new, drug.name)),]
    readr::write_csv(t1, paste0("tables/", cancers[i], "_first-line_therapy_2021-03-28.csv"), na="")
  }
  
  
  # Output drug names for assessment ###########################################
  adat %>% select(matches("^drct.*diag$")) %>% names() %>% substring(6) %>% 
    substr(1,nchar(.)-5) %>%
    write.csv("drugs.csv")
  
  # Merge changes Ravi + Sam made with new files 
  #pro.m <- read_csv("tables/pro_first-line-therapy2020-11-23.csv") %>% select(-n)
  #pro.sam <- read_csv("tables/pro_first-line-therapy_SUT.csv")
  #pro.sam2 <- left_join(pro.m, pro.sam, by=".")
  #readr::write_csv(pro.sam2, "tables/pro_first-line-therapy2020-11-23_SUT.csv", na="")
  
  #pan.m <- read_csv("tables/pan_first-line-therapy2020-11-23.csv") %>% select(-n)
  #pan.ravi <- read_csv("tables/pan_first-line-therapy_RP.csv")
  #pan.ravi2 <- left_join(pan.m, pan.ravi, by=".")
  #readr::write_csv(pan.ravi2, "tables/pan_first-line-therapy2020-11-23_RP.csv", na="")
  
} # end assess.external


# Create variables #############################################################
elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

adat <- adat %>% dplyr::mutate(
  # Demographics
    # Set age as difference between diagnosis date and birth year (Date set to Jan 1)
  age = as.numeric(format(diagnosisdate, '%Y')) - birthyear,
  reth = ifelse(ethnicity == "Hispanic or Latino" & !is.na(ethnicity), 
                        "Hispanic",
         ifelse(race == "", NA,
         ifelse(race == "Hispanic or Latino", "Hispanic",
         ifelse(race == "Black or African American", "Black",
         ifelse(race %in% c("Asian", "Other Race"), "Other", race))))),
  gend = ifelse(gender == "F", "Female",
                ifelse(gender == "M", "Male", NA)),
  state = ifelse(state == "", NA, state),
  
  # Comorbidities/drugs
  denovo_met = ifelse(groupstage %in% c("IV", "IVA", "IVB", "IVC", 
                                       "Stage IV", "Stage IVA", "Stage IVB"),
                     1, 
                     ifelse(groupstage %in% c("Unknown", 
                                                 "Unknown / Not documented",
                                                 "Unknown/not documented",
                                                 "Not documented",
                                                 "Group stage is not reported"),
                               NA, 0)),
  liverf_adv = ifelse(is.na(result_glom_diag_adv), NA,
                      ifelse(result_glom_diag_adv < 15, 1, 0)),
  kidneyf_adv = ifelse(is.na(result_bili_diag_adv), NA,
                       ifelse(result_bili_diag_adv > 2, 1, 0)),
  steroid_adv = ifelse(is.na(drgenct_steroid_diag), 0,
                       ifelse(drgenct_steroid_diag > 0, 1, 0)),
  opioid_adv = ifelse(is.na(drgenct_painagent_diag), 0,
                      ifelse(drgenct_painagent_diag > 0, 1, 0)),
  liverf_line1 = ifelse(is.na(result_glom_line1), NA,
                        ifelse(result_glom_line1 < 15, 1, 0)),
  kidneyf_line1 = ifelse(is.na(result_bili_line1), NA,
                         ifelse(result_bili_line1 > 2, 1, 0)),
  steroid_line1 = ifelse(is.na(drgenct_steroid_line1), 0,
                         ifelse(drgenct_steroid_line1 > 0, 1, 0)),
  opioid_line1 = ifelse(is.na(drgenct_painagent_line1), 0,
                        ifelse(drgenct_painagent_line1 > 0, 1, 0)),
  
  # Set adt treatment indicator for Aim 2 (and part of Aim 3)
     # Patients who receive any non-ADT therapy within 60 days of
     # diagnosis will be classified using that therapy. Patients with no non-ADT
     # therapy will be classified using ADT.
  therapydays = as.numeric(difftime(starttherapydate_1, 
                                               diagnosisdate,
                                               units = "days")),
  adtdays = as.numeric(difftime(adt_start, 
                                           diagnosisdate,
                                           units = "days")),
  trtadt = ifelse(is.na(therapydays) | therapydays > 60, 
                    # Main therapy not started in 60 days. 
                    # If ADT not started in 60 days, 0. IF ADT started, 1.
                    ifelse(is.na(adtdays) | adtdays > 60, 0, 1), 
                    # Main therapy started in 60 days, 0.
                    0
                  ),
  myelo.yes = ifelse(trtadt == 1, 0, myelo.yes),
  targeted.yes = ifelse(trtadt == 1, 0, targeted.yes),
  
  # Aim 2 time used for sample selection
  surv.days.a2 = ifelse(trtadt == 1, adtdays, therapydays),
  
  # Period
  diagnosisyear = ifelse(diagnosisdate >= ymd("2020/01/01") & diagnosisdate < ymd("2021/01/01"),
                 "2020", "2019"),
  diagnosisday = yday(diagnosisdate),
  precovid = ifelse(diagnosisdate %within% interval(ymd("2019/01/01"), ymd("2019/03/08")) |  
                      diagnosisdate %within% interval(ymd("2020/01/01"), ymd("2020/03/08")),
                   "pre",
             ifelse(diagnosisdate %within% interval(ymd("2019/04/08"), ymd("2019/07/31")) |  
                      diagnosisdate %within% interval(ymd("2020/04/08"), ymd("2020/07/31")),
                    "post", "washout")
             ),
  studyperiod = ifelse(diagnosisyear == "2019",
                ifelse(precovid == "pre", "Jan 1-Mar 8, 2019",
                       ifelse(precovid == "post", "Apr 8-Jul 31, 2019",
                              "Washout 2019")),
                ifelse(precovid == "pre", "Jan 1-Mar 8, 2020",
                       ifelse(precovid == "post", "Apr 8-Jul 31, 2020",
                              "Washout 2020"))
                ),
  starttherapydate_new = as.Date(ifelse(is.na(starttherapydate_1), 
                                ifelse(is.na(adt_start), NA, adt_start), # starttherapydate_1 is missing
                                ifelse(is.na(adt_start), starttherapydate_1, # not missing, adt missing
                                       ifelse(adt_start < starttherapydate_1, 
                                              adt_start, starttherapydate_1)
                                       )
                                )),
                                
  # Create time to treatment initiation survival variables
  surv.days = ifelse(!is.na(starttherapydate_new),
                     as.numeric(difftime(starttherapydate_new, 
                                         diagnosisdate,
                                         units = "days")),
                     as.numeric(difftime(censdate, 
                                         diagnosisdate,
                                         units = "days"))
                     ),
  surv.ind = ifelse(is.na(starttherapydate_new), 0, 1
                    ),
  # Add TTI survival indicator/time censored at 90 days after advanced diagnosis
  surv.ind90 = ifelse(surv.ind == 1 & surv.days > 90, 0, surv.ind),
  surv.days90 = ifelse(surv.days >90, 90, surv.days),
  
  # Create overall survival variables (time months )
  osurv.months = ifelse(!is.na(OSurv), 
                        elapsed_months(deathdate, diagnosisdate),
                        elapsed_months(censdate, diagnosisdate)
                        ),
  osurv.ind = ifelse(is.na(OSurv), 0, 1
                     ),
  # Create progression free survival variables for breast, rcc, and nsclc.
     # Other cancers do not have progression var.
  pfsurv.ind = ifelse(!(cancer %in% c("bre", "rcc", "nsc")), NA,
                      # If progression observed, set to 1, else 0.
                      ifelse(!is.na(progressiondate), 1, 0)
                      ),
  pfsurv.days = ifelse(is.na(pfsurv.ind), NA,
                       ifelse(pfsurv.ind == 1, 
                              as.numeric(difftime(progressiondate, 
                                                  diagnosisdate,
                                                  units = "days")),
                              # Set time for censored subjects
                              as.numeric(difftime(censdate, 
                                                  diagnosisdate,
                                                  units = "days")))
                              
                       ), #end pfsurv.days
  # Progression free survival after treatment
  pfsurv.trt.ind = ifelse(!(cancer %in% c("bre", "rcc", "nsc") | is.na(starttherapydate_new)), NA,
                      # If progression or death observed, set to 1, else 0.
                      ifelse(!is.na(progressiondate) | osurv.ind == 1, 1, 0)
  ),
  pfsurv.trt.days = ifelse(is.na(pfsurv.ind), NA,
                       ifelse(pfsurv.ind == 1,
                              # Set time for observed event subjects
                              # If no progression date, set to overall survival time
                              ifelse(is.na(progressiondate), OSurv,
                                     # If no death date use progression date to determine time
                                     ifelse(is.na(deathdate), 
                                            as.numeric(difftime(progressiondate, 
                                                                starttherapydate_new,
                                                                units = "days")),
                                            # If both, take the first occuring
                                            # (this should always be progression
                                            # but could be deathdate due to
                                            # how deaths were coded in these
                                            # data)
                                            ifelse(deathdate < progressiondate,
                                                   OSurv,
                                                   as.numeric(difftime(progressiondate, 
                                                                       starttherapydate_new,
                                                                       units = "days")))
                                     )
                              ),
                              # Set time for censored subjects
                              as.numeric(difftime(censdate, 
                                                  starttherapydate_new,
                                                  units = "days")))
  ) #end pfsurv.days
  
  
  ) #end mutate

ins.classify <- function(x){
  # Get number of observed payer category entries
  pcl <- x %>% select(matches("^payercategory*")) %>% length()
  icat <- NULL
  cv <- rep(NA, nrow(adat))
  
    for(i in 1:pcl){
      # Get entry i
      pc <- x %>% pull(paste("payercategory_", i, sep=""))
      isd <- x %>% pull(paste("insurestartdate_", i, sep=""))
      ied <- x %>% pull(paste("insureenddate_", i, sep=""))
      # if end date but not start date is missing, set end date to 1 year after
      # start date
      isd <- if_else(is.na(isd) & !is.na(ied), ied - years(1), isd)
      # if start date but not end date is missing, set start date to Jan 1, 2022
      ied <- if_else(is.na(ied) & !is.na(isd), ymd("2022/01/01"), ied)
      
      # Get start/end date interval
      iint <- interval(isd, ied)
      
      # Check if advanced diagnosis date in insured period. If yes, store.
      # If do categories (potentiall just table 1), use hierarchy, Commercial -> Medicare/medicaid ->
      icat <- ifelse(pc == "", "blank",
                     ifelse(pc == "Self Pay", "selfpay",
                         ifelse(pc %in% c("Medicare", "Medicaid", "Other Government Program"), "gov",
                                ifelse(pc %in% c("Other Payer - Type Unknown", 
                                                 "Patient Assistance Program",
                                                 "Workers Compensation"), "other",
                                       ifelse(pc == "Commercial Health Plan", "commercial", "drop catch"))))) %>%
        as.vector()
        
      cv <- ifelse(is.na(iint) | is.na(iint), cv,
              ifelse(!(x$diagnosisdate %within% iint), cv,
                     ifelse(is.na(cv), icat,
                       ifelse(icat == "commercial" | cv == "commercial", "commercial",
                              ifelse(icat == "gov" | cv == "gov", "gov",
                              ifelse(!is.na(icat) & icat != cv, "other", cv))))))

    } # end loop
  
  cv <- ifelse(is.na(cv) | cv == "selfpay", "other", cv)
    
  return(cv)
}

# Consolidate insurance: Medicaid, medicare, self-pay, commercial, other, unknown/missing
  # Is no insurance covering date equal to no insurance?
adat$insclass <- ins.classify(adat)

# Choose ECOG window
ecogcap <- function(dat, tb, ta){
  ad <- dat$ecogdate_postadvdiag
  bd <- dat$ecogdate_preadvdiag
  av <- dat$ecogvalue_postadvdiag
  bv <- dat$ecogvalue_preadvdiag
  dd <- dat$diagnosisdate
  
  # Get days from diagnosis
  adiff <- ad - dd 
  bdiff <- dd - bd 
  
  # Check if measurement falls within the specified range
  adiff <- ifelse(adiff <= ta, adiff, NA) 
  bdiff <- ifelse(bdiff <= tb, bdiff, NA) 
  
  # Return closest available value recorded in the specified range. If tied, take
     # average. 
  fv <- ifelse(is.na(bdiff),
               ifelse(is.na(adiff), NA, av), 
               ifelse(is.na(adiff), bv,
                      ifelse(adiff < bdiff, av,
                             ifelse(bdiff < adiff, bv, (bv + av)/2))))
}

# Create sample ecog ranges for evaluation
adat$ecog_b14a7 <- ecogcap(adat, 14, 7)
adat$ecog_b28a7 <- ecogcap(adat, 28, 7)
adat$ecog_b56a7 <- ecogcap(adat, 56, 7)
adat$ecog_b28a14 <- ecogcap(adat, 28, 14)
adat$ecog_b28a21 <- ecogcap(adat, 28, 21)
adat$ecog_b28a28 <- ecogcap(adat, 28, 28)
adat$ecog_b30a14 <- ecogcap(adat, 30, 14)

# Drop variables that we are not using any more
adat2 <- adat %>% select(-matches("^payercategory*"), -matches("^insure"),
                         -matches("^weightunit*"), -matches("^heightunit*"), 
                         -unit_height_line1, -unit_weight_line1, 
                         -matches("^testdate*"), -matches("^drgenct*"))

# Save data
saveRDS(adat2, paste("/data/analytic/PRACTICE-analytic-data_", Sys.Date(), ".rds", sep=""))
