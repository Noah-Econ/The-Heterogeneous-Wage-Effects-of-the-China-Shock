library(haven)
library(data.table)
library(ipumsr)
library(sf)
library(dplyr)
library(stringr)
library(xtable)


################## Load data
cw_puma1990_czone <- read_dta("data/cw_puma1990_czone.dta")
cw_puma2000_czone <- read_dta("data/cw_puma2000_czone.dta")
cw_puma2010_czone <- read_dta("data/cw_puma2010_czone.dta")
shocks <- read_dta("data/czone_exposure_by_period_v5_gh.dta")
setDT(cw_puma1990_czone)
setDT(cw_puma2000_czone)
setDT(cw_puma2010_czone)
cz_geo <- st_read("data/cz_geo.shp")
################## Load Census Data
ddi <- read_ipums_ddi("data/usa_00006.xml")
#data <- read_ipums_micro(ddi)
data1 <- read_ipums_micro(ddi, vars = c(
  "YEAR", "AGE", "INCWAGE_CPIU_2010", "UHRSWORK", "SEX", "WKSWORK2", "COUNTYFIP",
  "STATEFIP", "PUMA", "PERWT", "SAMPLE", "INDNAICS", "EDUC", "OCC", "RACE", "MARST", "CITIZEN",
  "METRO"
))

setDT(data1)
data <- data1
data <- data[YEAR >= 2005 & YEAR <= 2019,]
nrow(data)

######### subset age and hours worked variable -2mil obs
data <- subset(
  data,
  AGE >= 18 & AGE <= 64 
)

nrow(data)
############ include commuting zones - 
data_05_11 <- data[YEAR >= 2005 & YEAR <= 2011,]
data_12_19 <- data[YEAR >= 2012 & YEAR <= 2019,]

data_05_11[, PUMAFIP := paste0(
  sprintf("%02d", STATEFIP),
  sprintf("%04d", PUMA)
)]
data_05_11[, PUMAFIP:= as.numeric(PUMAFIP)]

data_12_19[, PUMAFIP := paste0(
  sprintf("%02d", STATEFIP),
  sprintf("%05d", PUMA)
)]
data_12_19[, PUMAFIP:= as.numeric(PUMAFIP)]

data_05_11 <- cw_puma2000_czone[data_05_11,on = .(puma2000 = PUMAFIP), allow.cartesian = T]
data_12_19 <- cw_puma2010_czone[data_12_19,on = .(puma2010 = PUMAFIP), allow.cartesian = T]

data <- rbindlist(list(data_05_11, data_12_19),use.names = TRUE,fill = TRUE)
data[, pweight := PERWT * afactor]


############################ wage construction 
data[, weeks_worked := fifelse(WKSWORK2 == 1, 7,
                               fifelse(WKSWORK2 == 2, 20,
                                       fifelse(WKSWORK2 == 3, 33,
                                               fifelse(WKSWORK2 == 4, 43.5,
                                                       fifelse(WKSWORK2 == 5, 48.5,
                                                               fifelse(WKSWORK2 == 6, 51, NA_real_))))))]

data[, wage := INCWAGE_CPIU_2010 / (weeks_worked * UHRSWORK)]
#### construct lnwage
data <- data %>%
  mutate(lnwage = log(wage))

## education variable ##################################################################################################################
data <- data %>%
  mutate(educ_years = case_when(
    EDUC == 0  ~ 0,
    EDUC == 1  ~ 2.5,
    EDUC == 2  ~ 6.5,
    EDUC == 3  ~ 9,
    EDUC == 4  ~ 10,
    EDUC == 5  ~ 11,
    EDUC == 6  ~ 12,
    EDUC == 7  ~ 13,
    EDUC == 8  ~ 14,
    EDUC == 9  ~ 15,
    EDUC == 10 ~ 16,
    EDUC == 11 ~ 18,
    EDUC == 99 ~ 99,
    TRUE ~ NA_real_
  ))

##experience variable ###########################################################################################
data <- data %>%
  mutate(
    exp = pmax(AGE - educ_years - 6, 0),
    exp_2 = exp^2
  )


### aggregate industry levels on a higher scale
data$INDNAICS <- as.factor(data$INDNAICS)
data <- data %>%
  mutate(
    naics_digits = str_extract(INDNAICS, "^[0-9]+"),   # keep leading digits
    naics2 = substr(naics_digits, 1, 2),
    naics3 = substr(naics_digits, 1, 3)
  )
data$INDNAICS[data$INDNAICS == 3]

data <- data %>%
  mutate(
    naics2 = as.character(naics2),
    naics3 = as.character(naics3)
  )

data <- data %>%
  mutate(
    naics_sector = case_when(
      naics2 == "11" ~ "Agriculture",
      naics2 == "21" ~ "Mining",
      naics2 == "22" ~ "Utilities",
      naics2 == "23" ~ "Construction",
      naics2 %in% c("31","32","33") ~ "Manufacturing",
      naics2 == "42" ~ "Wholesale",
      naics2 %in% c("44","45") ~ "Retail",
      naics2 %in% c("48","49") ~ "Transportation",
      naics2 == "51" ~ "Information",
      naics2 == "52" ~ "Finance",
      naics2 == "53" ~ "RealEstate",
      naics2 == "54" ~ "Professional",
      naics2 == "55" ~ "Management",
      naics2 == "56" ~ "AdminSupport",
      naics2 == "61" ~ "Education",
      naics2 == "62" ~ "HealthCare",
      naics2 == "71" ~ "ArtsRecreation",
      naics2 == "72" ~ "AccommodationFood",
      naics2 == "81" ~ "OtherServices",
      naics2 == "92" ~ "PublicAdmin",
      naics2 %in% c("3","4", "97","98","99") ~ "Unclassified",      
      naics2 %in% c("0",NA) ~ "Unemployed",  
      TRUE ~ NA_character_
    )
  )

data <- data %>%
  mutate(
    naics_sector = as.factor(naics_sector),
  )


##race variable #################
data <- data %>%
  mutate(
    RACE = as.factor(RACE),
  )
summary(data$RACE)

# Citizen status ####
data <- data %>%
  mutate(
    citizen_group = factor(
      ifelse(CITIZEN == 3, "NonCitizen", "Citizen"),
      levels = c("Citizen", "NonCitizen")
    )
  )

summary(data$citizen_group)

# Marital Status #########
data <- data %>%
  mutate(
    married = factor(
      ifelse(MARST %in% c(1,2), "Married", "NotMarried")
    )
  )
summary(data$married)

##data cleaning#################################
n0 <- nrow(data)

valid_cz <- unique(shocks$czone)
data1 <- data[!is.na(czone)]
data1 <- data1 %>%
  filter(czone %in% valid_cz)
n1 <- nrow(data1)

data3 <- data1[!is.na(wage)]
n3 <- nrow(data3)
rm(data, data1)
rm(data_05_11, data_12_19)

wage_cutoffs <- quantile(data3$wage, probs = c(0.01, 0.99), na.rm = TRUE)

data4 <- data3 %>%
  filter(wage > wage_cutoffs[1]) %>%
  filter(wage < wage_cutoffs[2]) %>%
  filter(!is.na(pweight))
n4 <- nrow(data4)


saveRDS(data4, "data/data_robust.rds")
