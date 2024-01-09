## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo=TRUE,
  error=TRUE,
  warning=TRUE,
  message=FALSE
)

## ----setup--------------------------------------------------------------------
library(apmx)
library(dplyr)
library(tidyr)

EX <- as.data.frame(EX)
PC <- as.data.frame(PC)
DM <- as.data.frame(DM)
LB <- as.data.frame(LB)

## -----------------------------------------------------------------------------
ex <- EX %>%
  dplyr::mutate(CMT = 1) %>%
  dplyr::select(USUBJID, STUDYID, EXSTDTC, VISIT, EXSTDY, EXTPTNUM, EXDOSE,
                CMT, EXTRT, EXTPT, EXROUTE, EXDOSFRQ, EXDOSU)

## -----------------------------------------------------------------------------
pc <- PC %>%
  dplyr::filter(PCSTAT=="Y") %>%
  dplyr::mutate(CMT = 2,
                TPT = dplyr::case_when(PCTPT=="<1 hour Pre-dose" ~ 0,
                                       PCTPT=="30 minutes post-dose" ~ 0.5/24,
                                       PCTPT=="1 hour post-dose" ~ 1/24,
                                       PCTPT=="2 hours post-dose" ~ 2/24,
                                       PCTPT=="4 hours post-dose" ~ 4/24,
                                       PCTPT=="6 hours post-dose" ~ 6/24,
                                       PCTPT=="8 hours post-dose" ~ 8/24,
                                       PCTPT=="12 hours post-dose" ~ 12/24,
                                       PCTPT=="24 hours post-dose" ~ 24/24,
                                       PCTPT=="48 hours post-dose" ~ 48/24)) %>%
  dplyr::select(USUBJID, PCDTC, PCDY, VISIT, TPT, PCSTRESN,
                PCLLOQ, CMT, PCTEST, PCTPT, PCSTRESU)

## -----------------------------------------------------------------------------
df_simple <- apmx::pk_build(ex = ex, pc = pc)

## -----------------------------------------------------------------------------
df_simple <- apmx::pk_build(ex = ex, #dataframe of prepared dose events
                            pc = pc, #dataframe of prepared pc observation events
                            time.units = "days", #can be set to days or hours.
                            #NOTE: units of TPT in ex and pc should match this unit
                            cycle.length = NA, #must be in units of days, will reset NTLC to 0
                            na = -999, #replaces missing nominal times and covariates with a numeric value
                            time.rnd = NULL, #rounds all time values to x decimal places
                            amt.rnd = NULL, #rounds calculated dose values to x decimal places
                            dv.rnd = NULL, #rounds observation columns to x decimal places
                            impute = NA, #imputation method for missing times
                            sparse = 3) #threshold for calculating sparse/serial distinctions


## -----------------------------------------------------------------------------
df_simple <- apmx::pk_build(ex, pc, time.rnd = 3)

## -----------------------------------------------------------------------------
dm <- DM %>%
  dplyr::select(USUBJID, AGE, SEX, RACE, ETHNIC) %>%
  dplyr::mutate(AGEU = "years") #AGE is continuous and requires a unit

## -----------------------------------------------------------------------------
lb <- LB %>% #select the desired labs
  dplyr::filter(LBCOMPFL=="Y") %>%
  dplyr::filter(LBVST %in% c("Baseline (D1)", "Screening")) %>%
  dplyr::filter(LBPARAMCD %in% c("ALB", "AST", "ALT", "BILI", "CREAT")) %>%
  dplyr::mutate(LBORRES = as.numeric(LBORRES))

lb <- lb %>% #select the lab collected immediately prior to first dose
  dplyr::arrange(USUBJID, LBPARAMCD, LBDT) %>%
  dplyr::group_by(USUBJID, LBPARAMCD) %>%
  dplyr::filter(row_number()==max(row_number())) %>%
  dplyr::ungroup()

lb <- lb %>% #finish formatting and add units since all labs are continuous
  dplyr::select(USUBJID, LBPARAMCD, LBORRES) %>%
  tidyr::pivot_wider(names_from = "LBPARAMCD", values_from = "LBORRES") %>%
  dplyr::mutate(ALBU = "g/dL",
                ASTU = "IU/L",
                ALTU = "IU/L",
                BILIU = "mg/dL",
                CREATU = "mg/dL")

## -----------------------------------------------------------------------------
tast <- LB %>%
  dplyr::filter(LBCOMPFL=="Y") %>%
  dplyr::filter(LBPARAMCD=="AST") %>%
  dplyr::mutate(LBORRES = as.numeric(LBORRES)) %>%
  dplyr::select(USUBJID, DTIM = LBDT, AST = LBORRES) %>%
  dplyr::mutate(ASTU = "IU/L")

## -----------------------------------------------------------------------------
talt <- LB %>%
  dplyr::filter(LBCOMPFL=="Y") %>%
  dplyr::filter(LBPARAMCD=="ALT") %>%
  dplyr::mutate(LBORRES = as.numeric(LBORRES)) %>%
  dplyr::select(USUBJID, DTIM = LBDT, ALT = LBORRES) %>%
  dplyr::mutate(ALTU = "IU/L")

## -----------------------------------------------------------------------------
pd <- LB %>%
  dplyr::filter(LBCOMPFL=="Y") %>%
  dplyr::filter(LBPARAM=="glucose") %>%
  dplyr::mutate(DTIM = paste(LBDT, "00:00"),
                VISIT = LBVST,
                NDAY = case_when(VISIT=="Screening" ~ -15,
                                 VISIT=="Baseline (D1)" ~ 1,
                                 VISIT=="Visit 2 (D8)" ~ 8,
                                 VISIT=="Visit 3 (D15)" ~ 15,
                                 VISIT=="Visit 4 (D29)" ~ 29,
                                 VISIT=="End of Treatment" ~ 45),
                TPT = 0,
                TPTC = LBTPT,
                ODV = as.numeric(LBORRES),
                DVIDU = LBORRESU,
                LLOQ = NA,
                CMT = 3,
                DVID = LBPARAM) %>%
  dplyr::select(USUBJID, DTIM, NDAY, VISIT, TPT,
                ODV, LLOQ, CMT, DVID, TPTC, DVIDU)

## -----------------------------------------------------------------------------
df_full <- apmx::pk_build(ex = ex, pc = pc, pd = pd,
                          sl.cov = list(dm, lb),
                          tv.cov = list(tast, talt),
                          time.rnd = 3)

## -----------------------------------------------------------------------------
unique(df_simple$EVID)
unique(df_full$EVID)

## -----------------------------------------------------------------------------
unique(df_simple$DVID)
unique(df_full$DVID)

unique(df_simple$DVIDC)
unique(df_full$DVIDC)

## -----------------------------------------------------------------------------
apmx::cov_find(df_full, cov = "categorical", type = "numeric")
apmx::cov_find(df_full, cov = "categorical", type = "character")
apmx::cov_find(df_full, cov = "continuous", type = "numeric")
apmx::cov_find(df_full, cov = "units", type = "character")

## -----------------------------------------------------------------------------
df_full <- apmx::pk_build(ex = ex, pc = pc, pd = pd,
                          sl.cov = list(dm, lb),
                          tv.cov = list(tast, talt),
                          time.rnd = 3,
                          cov.rnd = NULL, #rounds observation columns to x decimal places
                          BDV = FALSE, #calculates baseline dependent variable for PD events
                          DDV = FALSE, #calculates change (delta) from baseline for PD events
                          PDV = FALSE, #calculates percent change from baseline for PD events
                          demo.map = TRUE, #adds specific numeric mapping for SEX, RACE, and ETHNIC variables
                          tv.cov.fill = "downup", #fill pattern for time-varying covariates
                          keep.other = TRUE) #keep or drop all EVID = 2 rows

## -----------------------------------------------------------------------------
df_full <- apmx::pk_build(ex = ex, pc = pc, pd = pd,
                          sl.cov = list(dm, lb),
                          tv.cov = list(tast, talt),
                          time.rnd = 3, dv.rnd = 3,
                          BDV = TRUE, DDV = TRUE, PDV = TRUE,
                          keep.other = FALSE)

## -----------------------------------------------------------------------------
tast <- LB %>%
  dplyr::filter(LBCOMPFL=="Y") %>%
  dplyr::filter(LBPARAMCD=="AST") %>%
  dplyr::mutate(NTFD = case_when(LBVST=="Screening" ~ -15, #calculate NTFD from visit code
                                 LBVST=="Baseline (D1)" ~ 1,
                                 LBVST=="Visit 2 (D8)" ~ 8,
                                 LBVST=="Visit 3 (D15)" ~ 15,
                                 LBVST=="Visit 4 (D29)" ~ 29,
                                 LBVST=="End of Treatment" ~ 45)) %>%
  dplyr::mutate(AST = as.numeric(LBORRES)) %>%
  dplyr::select(USUBJID, NTFD, AST, ASTU = LBORRESU)

df_cov_apply <- apmx::pk_build(ex = ex, pc = pc,
                               sl.cov = list(dm, lb),
                               time.rnd = 3, dv.rnd = 3,
                               BDV = TRUE, DDV = TRUE, PDV = TRUE,
                               keep.other = FALSE) %>%
  apmx::cov_apply(tast, time.by = "NTFD")

## -----------------------------------------------------------------------------
df_cov_apply <- apmx::pk_build(ex = ex, pc = pc,
                               time.rnd = 3, dv.rnd = 3,
                               BDV = TRUE, DDV = TRUE, PDV = TRUE,
                               keep.other = FALSE) %>%
  apmx::cov_apply(dm) %>%
  apmx::cov_apply(lb) %>%
  apmx::cov_apply(talt, time.by = "DTIM") %>%
  apmx::cov_apply(tast, time.by = "NTFD")

## -----------------------------------------------------------------------------
exposure <- data.frame(ID = 1:22, #exposure metrics
                       MAX = 1001:1022,
                       MIN = 101:122,
                       AVG = 501:522)

parameters <- data.frame(ID = 1:22, #individual clearance and central volume estimates
                         CL = seq(0.1, 2.2, 0.1),
                         VC = seq(1, 11.5, 0.5))

## -----------------------------------------------------------------------------
df_cov_apply <- apmx::pk_build(ex = ex, pc = pc,
                               time.rnd = 3, dv.rnd = 3,
                               BDV = TRUE, DDV = TRUE, PDV = TRUE,
                               keep.other = FALSE) %>%
  apmx::cov_apply(dm) %>%
  apmx::cov_apply(lb) %>%
  apmx::cov_apply(talt, time.by = "DTIM", keep.other = FALSE) %>%
  apmx::cov_apply(tast, time.by = "NTFD", keep.other = FALSE) %>%
  apmx::cov_apply(exposure, id.by = "ID", exp = TRUE) %>%
  apmx::cov_apply(parameters, id.by = "ID", ebe = TRUE)

## -----------------------------------------------------------------------------
apmx::cov_find(df_cov_apply, cov = "categorical", type = "numeric")
apmx::cov_find(df_cov_apply, cov = "categorical", type = "character")
apmx::cov_find(df_cov_apply, cov = "continuous", type = "numeric")
apmx::cov_find(df_cov_apply, cov = "units", type = "character")
apmx::cov_find(df_cov_apply, cov = "exposure", type = "numeric")
apmx::cov_find(df_cov_apply, cov = "empirical bayes estimate", type = "numeric")

## -----------------------------------------------------------------------------
warning <- df_full %>%
  dplyr::filter(USUBJID=="ABC102-01-005")

nrow(warning)
warning$DVIDC

## -----------------------------------------------------------------------------
warning$C
warning$TIMEF

## -----------------------------------------------------------------------------
ex_error <- ex[, -5]

apmx::pk_build(ex_error, pc)

## -----------------------------------------------------------------------------
ex_error <- ex
ex_error$USUBJID <- 1:42

apmx::pk_build(ex_error, pc)

## -----------------------------------------------------------------------------
ex_error <- ex
ex_error$USUBJID[5] <- NA

apmx::pk_build(ex_error, pc)

## -----------------------------------------------------------------------------
ex_error <- ex
ex_error$ADDL <- 1

apmx::pk_build(ex_error, pc)

## -----------------------------------------------------------------------------
ex_error <- ex
ex_error$EXSTDTC <- substr(ex_error$EXSTDTC, 1, 10)

apmx::pk_build(ex_error, pc)

## -----------------------------------------------------------------------------
ex_error <- ex
ex_error$EXSTDY <- 0

apmx::pk_build(ex_error, pc)

## -----------------------------------------------------------------------------
ex_error <- ex
ex_error$ADDL <- 1
ex_error$II <- c(rep(1, 41), NA)

apmx::pk_build(ex_error, pc)

## -----------------------------------------------------------------------------
apmx::pk_build(ex)

## -----------------------------------------------------------------------------
pc_error <- pc
pc_error$PCSTRESN[10] <- 0

apmx::pk_build(ex, pc_error)

## -----------------------------------------------------------------------------
ex_error <- ex %>%
  select(-STUDYID)

apmx::pk_build(ex_error, pc)

## -----------------------------------------------------------------------------
dm_error <- dm
dm_error$USUBJID[2] <- "ABC102-01-001"

apmx::pk_build(ex, pc, sl.cov=dm_error)

## -----------------------------------------------------------------------------
apmx::pk_build(ex, pc, time.units="minutes")

## -----------------------------------------------------------------------------
apmx::pk_build(ex, pc, pd, DDV=TRUE, PDV==TRUE)

## -----------------------------------------------------------------------------
ex_error <- ex
ex_error$NSEX <- 0

apmx::pk_build(ex_error, pc, sl.cov = dm)

## -----------------------------------------------------------------------------
dm_error <- dm %>%
  select(-AGEU)

apmx::pk_build(ex, pc, sl.cov = dm_error)

## -----------------------------------------------------------------------------
dm_warning <- dm
dm_warning <- dm_warning[1:4,]

df_warning <- apmx::pk_build(ex, pc, sl.cov=dm_warning)

## -----------------------------------------------------------------------------
df_warning <- apmx::pk_build(ex, pc, sl.cov = list(dm_warning, lb))

## -----------------------------------------------------------------------------
pd_warning <- pd
pd_warning <- pd[3:nrow(pd_warning), ]

df_warning <- apmx::pk_build(ex, pc, pd_warning, BDV=TRUE)

## -----------------------------------------------------------------------------
df_warning <- apmx::pk_build(ex, pc, pd_warning)

## -----------------------------------------------------------------------------
pc_warning <- pc
pc_warning$TPT[1] <- 0.07

df_warning <- apmx::pk_build(ex, pc_warning,
                             time.rnd = 3)

## -----------------------------------------------------------------------------
ex_warning <- ex
ex_warning$EXDOSE[1] <- NA

df_warning <- apmx::pk_build(ex_warning, pc,
                             time.rnd = 3)

## -----------------------------------------------------------------------------
pc_warning <- pc
pc_warning[2, ] <- pc_warning[1, ]
pc_warning$PCSTRESN[2] <- 1400

df_warning <- apmx::pk_build(ex, pc_warning,
                             time.rnd = 3)

## -----------------------------------------------------------------------------
dm_warning <- dm %>%
  rename(ETHNICITY = ETHNIC)

df_warning <- apmx::pk_build(ex, pc, sl.cov = dm_warning)

## -----------------------------------------------------------------------------
lb_warning <- lb
lb_warning$ALT[1] <- 31

df_warning <- apmx::pk_build(ex, pc, sl.cov = lb_warning, tv.cov = talt)

## -----------------------------------------------------------------------------
pc_impute <- pc
pc_impute$PCDTC[c(4, 39, 73, 128)] <- NA

df_impute <- apmx::pk_build(ex, pc_impute,
                            time.rnd = 3)

## -----------------------------------------------------------------------------
df_impute_1 <- apmx::pk_build(ex, pc_impute,
                              time.rnd = 3, impute = 1)

## -----------------------------------------------------------------------------
nrow(df_impute_1[is.na(df_impute_1$ATFD),]) #number of rows with missing ATFD

imputed_events_1 <- df_impute_1 %>%
  dplyr::filter(IMPDV==1 | IMPEX==1)

## -----------------------------------------------------------------------------
times_check_1 <- df_impute_1 %>%
  dplyr::filter(USUBJID=="ABC102-01-004")

## -----------------------------------------------------------------------------
df_impute_2 <- apmx::pk_build(ex, pc_impute,
                              time.rnd = 3, impute = 2)

imputed_events_2 <- df_impute_2 %>%
  dplyr::filter(IMPDV==1 | IMPEX==1)

## -----------------------------------------------------------------------------
times_check_2 <- df_impute_2 %>%
  dplyr::filter(USUBJID=="ABC102-01-004")

## -----------------------------------------------------------------------------
ex_impute <- ex
ex_impute$EXSTDTC[2] <- NA

df_impute <- apmx::pk_build(ex_impute, pc, #no imputation method
                            time.rnd = 3)

## -----------------------------------------------------------------------------
df_impute_1 <- apmx::pk_build(ex_impute, pc, #imputation method 1
                              time.rnd = 3, impute = 1)

imputed_events_1 <- df_impute_1 %>% #imputed records
  dplyr::filter(IMPDV==1 | IMPEX==1)

## -----------------------------------------------------------------------------
df_impute_2 <- apmx::pk_build(ex_impute, pc,
                              time.rnd = 3, impute = 2)

imputed_events_2 <- df_impute_2 %>%
  dplyr::filter(IMPDV==1 | IMPEX==1)

## -----------------------------------------------------------------------------
ex_impute <- ex
ex_impute$EXSTDTC[1] <- NA

df_impute <- apmx::pk_build(ex_impute, pc, # No imputation method, expect a warning
                            time.rnd = 3)

## -----------------------------------------------------------------------------
df_impute_2 <- apmx::pk_build(ex_impute, pc, #imputation method 2
                              time.rnd = 3, impute = 2)

imputed_events_2 <- df_impute_2 %>% #imputed events
  dplyr::filter(IMPDV==1 | IMPEX==1 | IMPFEX==1)

## -----------------------------------------------------------------------------
ex_impute <- ex
ex_impute$EXSTDTC[1:2] <- NA

df_impute <- apmx::pk_build(ex = ex_impute, pc = pc_impute, #no impuation method
                            time.rnd = 3)

## -----------------------------------------------------------------------------
df_impute_2 <- apmx::pk_build(ex = ex_impute, pc = pc_impute, #imputation method 2
                              time.rnd = 3, impute = 2)

## -----------------------------------------------------------------------------
df_full2 <- df_full %>%
  dplyr::filter(DOMAIN!="PD") %>% #remove glucose observations
  dplyr::filter(ID<19) %>% #remove subject 19
  dplyr::group_by(ID) %>%
  dplyr::mutate(NSTUDYC = "ABC103", #update study ID
                USUBJID = gsub("ABC102", "ABC103", USUBJID),
                BAGE = round(rnorm(1, 45, 10)), #re-create all continuous covariates
                BALB = round(rnorm(1, 4, 0.5), 1),
                BALT = round(rnorm(1, 30, 5)),
                BAST = round(rnorm(1, 33, 5)),
                BBILI = round(rnorm(1, 0.7, 0.2), 3),
                BCREAT = round(rnorm(1, 0.85, 0.2), 3),
                TAST = ifelse(NTFD==0, BAST, round(rnorm(1, 33, 5))),
                TALT = ifelse(NTFD==0, BALT, round(rnorm(1, 30, 5)))) %>%
  dplyr::ungroup()

## -----------------------------------------------------------------------------
df_combine <- apmx::pk_combine(df_full, df_full2)

## -----------------------------------------------------------------------------
unique(df_full$DVID)
unique(df_full2$DVID)

## -----------------------------------------------------------------------------
name <- "PK_ABC101_V01.csv"
apmx::pk_write(df_combine, file.path(tempdir(), name))

## -----------------------------------------------------------------------------
vl <- apmx::variable_list_create(variable = c("SEX", "RACE", "ETHNIC", "AGE",
                                              "ALB", "ALT", "AST", "BILI", "CREAT"),
                           categorization = rep("Covariate", 9),
                           description = c("sex", "race", "ethnicity", "age",
                                           "albumin", "alanine aminotransferase",
                                           "aspartate aminotransferase",
                                           "total bilirubin", "serum creatinine"))

## -----------------------------------------------------------------------------
define <- apmx::pk_define(df = df_combine,
                          variable.list=vl)

## -----------------------------------------------------------------------------
define <- apmx::pk_define(df = df_combine,
                          file = file.path(tempdir(), "definition_file.docx"),
                          variable.list=vl,
                          project = "Sponsor Name",
                          data = "Dataset Name")

## -----------------------------------------------------------------------------
vrlg <- apmx::version_log(df = df_combine,
                          name = name,
                          file = file.path(tempdir(), "version_log.docx"),
                          src_data = "original test data")

## -----------------------------------------------------------------------------
sum1 <- apmx::pk_summarize(df = df_combine)

## -----------------------------------------------------------------------------
sum2 <- apmx::pk_summarize(df = df_combine,
                           strat.by = c("NSTUDYC", "NSEXC"),
                           ignore.request = "NRACE == 2")

