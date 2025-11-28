
library(readxl)
library(tidyverse)
library(janitor)
library(cowplot)
library(grid)
library(gridExtra)
library(lemon)

##########################################################################################
#SECTION 1: Data set-up and categorising histotypes and calculating age-standardised rates
##########################################################################################
#Before beginning download and have in working directory:
#Cancer datacube 10g: "CDiA-2024-Book-10g-Cancer-incidence-by-histology-vulva-vagina-uterus.xlsx" from AIHW: https://www.aihw.gov.au/reports-data/health-conditions-disability-deaths/cancer/data
#National population data: Population - Australia Population at 30 June, by sex and single year of age, Aust., from 1971 onwards from:https://www.abs.gov.au/statistics/people/population/national-state-and-territory-population/latest-release#data-downloads
#And 'hysterectomyfactions' excel file

dat <- read_excel("CDiA-2024-Book-10g-Cancer-incidence-by-histology-vulva-vagina-uterus.xlsx", 
                  sheet = "Table S10g.3", skip = 4)
dat <- dat %>% janitor::clean_names() %>%
  mutate(
    Subtype= case_when(
      i_d %in% c("U1.01.02.01.01", "U1.01.02.01.02", "U1.01.02.01.03", "U1.01.02.01.04", # Endometrioid adenocarcinomas
                 "U1.01.02.05.01", "U1.01.02.05.02", "U1.01.02.05.03", "U1.01.02.05.04", "U1.01.02.05.05", #mucinous carcinomas
                 "U1.01.02.06.01", "U1.01.02.06.02", "U1.01.02.06.03", "U1.01.02.06.04", "U1.01.03",  "U1.01.02.02.02") ~ "Endometrioid",
      i_d %in% c("U1.01.02.02.01", "U1.01.02.02.03", "U1.01.02.04", "U1.01.06.10") ~ "Serous",
      i_d =="U1.01.02.03" ~ "Clear cell",
      i_d %in% c("U1.01.04.01", "U1.01.04.02", "U1.01.04.03", "U1.01.04.04") ~ "Carcinosarcoma",
      i_d %in% c("U1.02.01.01", "U1.02.02.01", "U1.02.02.02", "U1.02.03.01","U1.02.04.01", "U1.02.04.02", "U1.02.04.03", "U1.02.04.04", "U1.02.05.01", "U1.02.05.02", "U1.02.05.03", "U1.02.05.04", "U1.02.05.05", "U1.02.06.01", "U1.02.06.02", "U1.02.06.03", "U1.02.07", "U1.02.08.01", "U1.02.08.02", "U1.02.08.03", 
                 "U1.02.08.05.01", "U1.02.08.05.02", "U1.02.08.05.03", "U1.02.08.05.04", "U1.02.09.01", "U1.02.10", "U1.03") ~ "Sarcoma",
      i_d %in% c("U1.01.01.01", "U1.01.01.02", "U1.01.01.03", "U1.01.01.04", "U1.01.01.05", "U1.01.01.06", "U1.01.02.06.05", "U1.01.02.06.06", "U1.01.02.06.07", "U1.01.05.01", "U1.01.05.02.01", "U1.01.05.02.02", "U1.01.05.02.03", "U1.01.05.03", "U1.01.06.01", "U1.01.06.02", 
                 "U1.01.06.03", "U1.01.06.04", "U1.01.06.05", "U1.01.06.06", "U1.01.06.07", "U1.01.06.08", "U1.01.06.09", "U1.01.06.11", "U1.02.01.02", "U1.02.08.04", "U1.02.08.05.05", "U1.04.01", "U1.04.02", "U1.04.03",
                 "U1.05", "U1.06") ~ "Other",
      is.na(i_d) ~ "Total",
      .default=" ")) %>%
  rename(Count = count_for_cancer_type,
         Age = age_group_years, 
         ALLUterineCancers = total_count_for_uterine_cancer,
         Year=year) %>%
  filter(sex=="Females" & age_grouping=="10-year age groups" & !Subtype==" "
         & !Age %in% c("All ages (crude rate)", "All ages (age-standardised rate - 2024 Australian population)", 'All ages (age-standardised rate - 2001 Australian Standard Population)')) %>%
  select(c(Subtype, cancer_type, Year, Age, Count, ALLUterineCancers, i_d))

### All uterine cancer counts by histological subtype from 2001-2020
histologicalsubtype <- dat %>%
  filter(!Age=="0–14" & !Age=="15–24") %>%
  group_by(Subtype, cancer_type, i_d) %>% 
  summarise("Number of cases" = sum(Count)) 


### Import population data and calculate 5-year age groups
pop <- read_excel("3101059.xlsx", sheet="Data1") %>%
  rename(Year = '...1') %>%
  select(Year | contains(c("Female"))) %>%
  filter(!Year==c("Unit", "Series Type", "Data Type", "Frequency", "Collection Month", "Series Start", "Series End", "No. Obs", "Series ID")) %>%
  mutate_if(is.character, as.numeric) 

pop$Year <- excel_numeric_to_date(pop$Year)
pop$Year <- lubridate::year(pop$Year)
names(pop) = gsub(pattern="Estimated Resident Population ;  Female  *", replacement="Age", x=names(pop))
names(pop) = gsub(pattern="* ;", replacement="", x=names(pop))
names(pop) = gsub(pattern="*; *", replacement="", x=names(pop))

pop <- pop %>%
  rename("Age100+" = "Age100 and over") %>%
  mutate("0–14"= Age0+Age1+Age2+Age3+Age4+Age5+Age6+Age7+Age8+Age9+Age10+Age11+Age12+Age13+Age14,
         "15–24" = Age15+Age16+Age17+Age18+Age19+Age20+Age21+Age22+Age23+Age24,
         "25–29" = Age25+Age26+Age27+Age28+Age29,
         "30–34" = Age30+Age31+Age32+Age33+Age34,
         "35–39" = Age35+Age36+Age37+Age38+Age39,
         "40–44" = Age40+Age41+Age42+Age43+Age44,
         "45–49" = Age45+Age46+Age47+Age48+Age49,
         "50–54" = Age50+Age51+Age52+Age53+Age54,
         "55–59" = Age55+Age56+Age57+Age58+Age59,
         "60–64" = Age60+Age61+Age62+Age63+Age64,
         "65–69" = Age65+Age66+Age67+Age68+Age69,
         "70–74" = Age70+Age71+Age72+Age73+Age74,
         "75–79" = Age75+Age76+Age77+Age78+Age79,
         "80–84" = Age80+Age81+Age82+Age83+Age84,
         "85+"   = Age85+Age86+Age87+Age88+Age89+Age90+Age91+Age92+Age93+Age94+Age95+Age96+Age97+Age98+Age99+`Age100+`) %>%
  select(-contains("Age")) %>%
  pivot_longer(cols = contains(c("–","+")),
               names_to = "Age",
               values_to = "pop") %>%
  filter(Year %in% seq(from=2001, to=2020))

#calculate population at risk using AIHW hysterectomy fractions then group into 10-year age groups to match with Cancer data 10-year groups


Year <- c(seq(2001, 2020))
fractions <- read_excel("hysterectomy fractions.xlsx") %>% 
  merge(Year) %>%
  rename(Year = y) %>%
  mutate(fraction = case_when(
    Year %in% c(seq(2001,2004)) ~`Applied to the period 2001-2004`,
    Year %in% c(seq(2005,2015)) ~`Applied to the period 2005-2015`,
    Year %in% c(seq(2016,2020)) ~`Applied to the period 2016-2022`),
    Age = case_when(
      Age =="90+" ~ "85+",
      !Age =="90+" ~ Age)) %>%
  select(Year, Age, fraction)
pop <- pop %>%
  merge(fractions, by=c("Year", "Age")) %>%
  mutate(withouthyst = pop*fraction) %>%
  select(c(Year, Age, pop, withouthyst)) %>%
  pivot_wider(names_from=Age,
              values_from = c(pop, withouthyst)) %>%
  mutate("pop_25–34"= `pop_25–29` + `pop_30–34`,
         "pop_35–44"= `pop_35–39` + `pop_40–44`,
         "pop_45–54"= `pop_45–49` + `pop_50–54`,
         "pop_55–64"= `pop_55–59` + `pop_60–64`,
         "pop_65–74"= `pop_65–69` + `pop_70–74`, 
         "pop_75–84"= `pop_75–79` + `pop_80–84`,
         "withouthyst_25–34"= `withouthyst_25–29`+`withouthyst_30–34`,
         "withouthyst_35–44"= `withouthyst_35–39`+`withouthyst_40–44`,
         "withouthyst_45–54"= `withouthyst_45–49`+`withouthyst_50–54`,
         "withouthyst_55–64"= `withouthyst_55–59`+`withouthyst_60–64`,
         "withouthyst_65–74"= `withouthyst_65–69`+`withouthyst_70–74`,
         "withouthyst_75–84"= `withouthyst_75–79`+`withouthyst_80–84`) %>%
  pivot_longer(cols = contains(c("hyst", "pop")),
               names_to = c(".value", "Age"),
               names_pattern = "(.*)_(.*)") %>%
  filter(Age=="0–14" | Age=="15–24" | Age=="25–34"| Age=="35–44"| Age=="45–54"| Age=="55–64"| Age=="65–74" | Age=="75–84" | Age=="85+") %>%
  group_by(Year) %>%
  arrange(Age, .by_group=T) %>%
  ungroup()


### Separate 2001 population to use as the standard population
standardpop <- pop %>%
  filter(Year==2001) %>%
  select(!Year) %>%
  rename(withouthyststandardpop = withouthyst) %>%
  rename(standardpop= pop)

### generate Dataset with event counts for each subtype and age group 
data <- dat %>% 
  group_by(Subtype, Year, Age) %>% 
  summarise(Count = sum(Count), .groups = "keep") %>%
  filter(!Age %in% c("0–14", "15–24")) %>%
  mutate(analysis="primary")


### Create different counts of Endometrioid and Serous after Re-classifying adenocarcinoma NOS to Endometrioid and Serous using a range of proportions for use in sensitivity analysis
###percent of cancers that are Endometrioid vs Serous in 2020


adeno <- dat %>% 
  filter(cancer_type %in% c("Adenocarcinoma, NOS")) %>%
  group_by(Year, Age) %>%
  rename(`100%` = Count) %>%
  select(Year, Age, contains("%"))

dat %>%
  filter(Subtype %in% c("Endometrioid", "Serous") & Year==2020) %>%
  group_by(Subtype) %>%
  summarise(Count = sum(Count)) %>%
  mutate(proportion = Count/(sum(Count)))
### In 2020: endometrioid was 90%, total serous was 10%
for (i in c(90, 75, 50, 25, 10, 0)) {
  adeno <- adeno %>%
    mutate(!!paste0(i, "%") := `100%` *i/100) 
}


### Subtract adenocarcenoma NOS counts from Endometrioid then add the counts corresponding to the different proportions
sensitivity <- data %>%
  filter(Subtype %in% c("Endometrioid", "Serous")) %>%
  merge(adeno, by=c("Year", "Age")) %>%
  mutate(across(where(is.numeric), ~round(., 0)),
         Count = case_when(
           Subtype=="Endometrioid"  ~ Count-`100%`, # subtract adenocarcinoma counts previously added to endometrioid
           .default=Count)) %>%
  pivot_longer(cols=contains("%"), names_to="percent", values_to="adeno") %>%
  mutate(Count = Count + adeno,
         analysis="sensitivity") %>%
  filter(!(percent=="90.0%" & Subtype=="Serous") & !(percent=="10.0%" & Subtype=="Endometrioid"))

data <- rbind(data, sensitivity)%>%
  merge(pop,by=c("Year","Age")) %>%
  merge(standardpop,by="Age") %>%
  arrange(analysis, percent, Subtype, Year, Age) %>%
  mutate(
    atrisk= withouthyst+Count,
    atriskstandardpop= withouthyststandardpop+Count,
    agespecificrate=Count/atrisk,
    agespecificrate_nc=Count/pop) %>%
  select(analysis, percent, Subtype, Year, Age, Count, atrisk, pop, atriskstandardpop, standardpop, agespecificrate, agespecificrate_nc) 


###calculating age-standardised rate 
# First need the 2001 female total population and at risk population for each age group
totalpop <- data %>%
  filter(Year==2001 & Subtype== "Total") %>%
  summarise(totalpop = sum(pop), .groups = "keep")
totalatrisk <- data %>%
  filter(Year==2001 & Subtype== "Total") %>%
  summarise(totalatrisk = sum(atrisk), .groups = "keep")

agestandardised <- data %>% 
  mutate(
    'corrected for hyst' =agespecificrate*atriskstandardpop,
    'not corrected for hyst'= agespecificrate_nc*standardpop) %>%
  pivot_longer(c('corrected for hyst', 'not corrected for hyst'), names_to="corrected", values_to="rate") %>%
  group_by(analysis, percent, Subtype, Year, corrected) %>% 
  summarise(ageadjustedrate = sum(rate), .groups = "keep") %>%
  mutate(ageadjustedrate= case_when(
    corrected=='corrected for hyst' ~ (ageadjustedrate/totalatrisk$totalatrisk)*100000,
    corrected=='not corrected for hyst' ~ (ageadjustedrate/totalpop$totalpop)*100000)) %>%
  ungroup()

data <- data %>%
  filter(analysis=="primary")
### export data with event counts to be imported into Joinpoint software
histotype <- agestandardised %>% 
  filter(analysis=="primary") %>%
  select(corrected, Subtype, Year, ageadjustedrate) %>%
  arrange(corrected, Subtype) 
write.csv(histotype, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/histotype.csv", row.names = T)
sensitivity_year <- histotype %>%
  filter(Year>2005)
write.csv(sensitivity_year, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/sensitivity_year.csv", row.names = T)




total <- histotype %>%
  filter(corrected == "corrected for hyst" & Subtype=="Total") %>%
  rename(Age = Subtype,
         agespecificrate = ageadjustedrate) %>%
  mutate(Count =0,
         atrisk=0) %>%
  select(Age, Year, Count, atrisk, agespecificrate)

Age <- data %>%
  filter(Subtype == "Total") %>%
  arrange(Subtype, Age) %>%
  mutate(agespecificrate = agespecificrate*100000) %>%
  select(Age, Year, Count, atrisk, agespecificrate)
Age <- rbind(Age, total) %>%
  mutate(Age = case_when(
    Age =="Total" ~ "Age–standardised",
    !Age == "Total" ~ Age))

write.csv(Age, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/Age.csv", row.names = T)

adenocarcinoma <- agestandardised %>%
  filter(analysis=="sensitivity" & corrected=="corrected for hyst") %>%
  select(Subtype, percent, Year, ageadjustedrate) %>%
  arrange(Subtype, percent)
write.csv(adenocarcinoma, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/adenocarcinoma.csv")

PAFs <- data %>%
  filter(Subtype %in% c("Endometrioid", "Serous", "Clear cell") & analysis=="primary" & Year %in% c(2001, 2010, 2020)) %>%
  select(Subtype, Age, Year, Count)
write.csv(PAFs, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/For PAF calculation.csv")

checkjoinpoint <- histotype %>%
  filter(Subtype=="Endometrioid" & Year %in% c(2001, 2010, 2020), corrected=="corrected for hyst")
write.csv(checkjoinpoint, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/checktotalrate.csv")


##########################################################
#SECTION 2: Create serous dataset for Sensitivity analyses
##########################################################
###Reassign Carcinosarcoma into a new Serous group, while keeping serous total and keeping serous subgroups (high-grade, low-grade, mixed cell carcinoma and NOS)

###keep original rates of serous uterine cancers without Carcinosarcoma
serous <- agestandardised %>%
  filter(Subtype=="Serous" & analysis=="primary") %>%
  select(Subtype, Year, corrected, ageadjustedrate) %>%
  mutate(Subtype="All Serous exc. Carcinosarcoma")

###create new rates of serous uterine cancers after including Carcinosarcoma as serous
seroustotal <- dat %>%
  filter(Subtype %in% c("Serous", "Carcinosarcoma")) %>%
  mutate(Subtype = "All Serous") %>% 
  group_by(Subtype, Year, Age) %>% 
  summarise(Count = sum(Count), .groups = "keep") %>%
  filter(!Age %in% c("0–14", "15–24")) %>%
  merge(pop,by=c("Year","Age")) %>%
  merge(standardpop,by="Age") %>%
  arrange(Year, Age) %>%
  mutate(
    atrisk= withouthyst+Count,
    atriskstandardpop= withouthyststandardpop+Count,
    agespecificrate=Count/atrisk,
    agespecificrate_nc=Count/pop) %>%
  select(Subtype, Year, Age, Count, atrisk, pop, atriskstandardpop, standardpop, agespecificrate, agespecificrate_nc)

###create  rates of serous uterine cancer subgroups including Carcinosarcoma
seroussubgroups <- dat %>%
  filter(Subtype %in% c("Serous", "Carcinosarcoma")) %>%
  mutate(Subtype = cancer_type,
         Subtype = case_when(
           Subtype=="Serous carcinoma, NOS" ~ "NOS",
           !Subtype=="Serous carcinoma, NOS" ~ Subtype)) %>% 
  group_by(Subtype, Year, Age) %>% 
  summarise(Count = sum(Count), .groups = "keep") %>%
  filter(!Age %in% c("0–14", "15–24") & !Subtype=="Papillary carcinoma, NOS") %>%
  merge(pop,by=c("Year","Age")) %>%
  merge(standardpop,by="Age") %>%
  arrange(Year, Age) %>%
  mutate(
    atrisk= withouthyst+Count,
    atriskstandardpop= withouthyststandardpop+Count,
    agespecificrate=Count/atrisk,
    agespecificrate_nc=Count/pop) %>%
  select(Subtype, Year, Age, Count, atrisk, pop, atriskstandardpop, standardpop, agespecificrate, agespecificrate_nc) 

###Combine new total and subgroup datasets
combined <- rbind(seroustotal, seroussubgroups) %>%
  arrange(Subtype)

### and calculate age-standardised rates for new total and serous subgroups
newserous <- combined %>% 
  mutate(
    'corrected for hyst' =agespecificrate*atriskstandardpop,
    'not corrected for hyst'= agespecificrate_nc*standardpop) %>%
  pivot_longer(c('corrected for hyst', 'not corrected for hyst'), names_to="corrected", values_to="rate") %>%
  group_by(Subtype, Year, corrected) %>% 
  summarise(ageadjustedrate = sum(rate), .groups = "keep") %>%
  mutate(ageadjustedrate= case_when(
    corrected=='corrected for hyst' ~ (ageadjustedrate/totalatrisk$totalatrisk)*100000,
    corrected=='not corrected for hyst' ~ (ageadjustedrate/totalpop$totalpop)*100000))

### combine with original rates of serous (without carcinosarcoma)
serous <- rbind(serous, newserous) %>%
  filter(corrected=="corrected for hyst") %>%
  select(!corrected)

### export data with event counts to be imported into Joinpoint software
write.csv(serous, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/serous.csv", row.names = T)


##############################################################################
# SECTION 3: Preliminary Data Analysis and Plots of cancer counts and rates
##############################################################################

### All uterine cancer counts by histotype from 2001-2020
dat %>%
  filter(!Age=="0–14" & !Age=="15–24") %>%
  group_by(Subtype) %>% 
  summarise("Number of cases" = sum(Count))
### All uterine cancer counts by histotype from 2006-2020
dat %>%
  filter(!Age=="0–14" & !Age=="15–24" &Year>2005) %>%
  group_by(Subtype) %>% 
  summarise("Number of cases" = sum(Count))



### Counts of cases for the sensitivity analysis allocating adenocarcinoma NOS to different groups
sensitivity %>%
  group_by(Subtype, percent) %>% 
  summarise("Number of cases" = sum(Count))

sensitivity %>%
  group_by(Subtype, percent) %>% 
  summarise("Adenocarcinoma NOS" = sum(adeno))

### Cancer counts by age 
dat %>% 
  filter(!Subtype=="Total") %>%
  group_by(Age) %>% 
  summarise("Number of cases" = sum(Count))
### exploring young age groups
dat %>%
  filter(Age=="0–14" & !Subtype=="Total") %>%
  group_by(Year, Age, cancer_type) %>% 
  summarise(Count = sum(Count), .groups = "keep") %>%
  filter(!Count==0)
dat %>%
  filter(Age=="15–24"& !Subtype=="Total") %>%
  group_by(Year, Age, cancer_type) %>% 
  summarise(Count = sum(Count), .groups = "keep")%>%
  filter(!Count==0) %>%
  print(n=30)

### by histotype and year
plot <- histotype %>%
  filter(corrected == "not corrected for hyst")

p <- data %>%
  group_by(Subtype, Year) %>% 
  summarise(Count = sum(Count), .groups = "keep") %>%
  ggplot(aes(x=Year, y=Count, group=Subtype, color=Subtype)) +
  geom_line() +
  facet_wrap(~Subtype, scales="free")
p + theme_minimal()


### plotting age specific rates per 100, 000 women
p <- data %>%
  mutate(rate=agespecificrate*100000,
         rate_nc= agespecificrate_nc*100000) %>%
  pivot_longer(c(rate, rate_nc), names_to="corrected", values_to="rate") %>%
  ggplot(aes(x=Year, y=rate, group=Subtype, color=Subtype)) +
  geom_line() + geom_point() +
  facet_wrap(~Age + corrected) +
  labs(y="age specific rate")
p + theme_minimal() 



###plot counts of Adenocarcinoma NOS over time
adeno %>%
  group_by(Year) %>%
  summarise(Count = sum(`100%`)) %>%
  ggplot(aes(x=Year, y=Count)) + 
  geom_bar(stat="identity", fill="cornflowerblue") +
  labs(y= "      Counts of
       Adenocarcinoma 
       NOS") +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) 



###############################
# Section 4: Mortality analysis
###############################

library(readxl)
library(tidyverse)
library(janitor)
library(cowplot)

dat<- read_excel("CDiA-2024-Book-2a-Cancer-mortality-and-age-standardised-rates-by-age-5-year-groups.xlsx", sheet = "Table S2a.1", skip = 5) %>%
  janitor::clean_names() %>%
  filter(cancer_group_site=="Uterine cancer" & sex=="Females" & data_type=="Actual" & year>2000 & !age_group_years %in% c("00-04", "05-09", "10-14", "15-19", "20-24", "All ages combined")) %>%
  rename(Year= year,
         Age=age_group_years,
         Count = count) %>%
  select(Year, Age, Count)  
dat$Age <- gsub("-", "–", dat$Age)
dat <- dat %>%
  pivot_wider(names_from=Age,
              values_from = Count,
              names_prefix = "Age") %>%
  mutate(age25 = `Age25–29` + `Age30–34`,
         age35 = `Age35–39` + `Age40–44`,
         age45 = `Age45–49` + `Age50–54`,
         age55 = `Age55–59` + `Age60–64`,
         age65 = `Age65–69` + `Age70–74`,
         age75 = `Age75–79` + `Age80–84`,
         age85 = `Age85–89` + `Age90+`) %>%
  select(Year, !contains(c("–", "+"))) %>%
  pivot_longer(cols= contains("age"),
               names_to = "Age",
               values_to = "Count") %>%
  mutate(Age = case_when(
    Age == "age25" ~ "25–34",
    Age == "age35" ~ "35–44",
    Age == "age45" ~ "45–54",
    Age == "age55" ~ "55–64",
    Age == "age65" ~ "65–74",
    Age == "age75" ~ "75–84",
    Age == "age85" ~ "85+"))
write.csv(dat, "mortalitycounts.csv")




### Import population data and calculate 5-year age groups
pop <- read_excel("3101059.xlsx", sheet="Data1") %>%
  rename(Year = '...1') %>%
  select(Year | contains(c("Female"))) %>%
  filter(!Year==c("Unit", "Series Type", "Data Type", "Frequency", "Collection Month", "Series Start", "Series End", "No. Obs", "Series ID")) %>%
  mutate_if(is.character, as.numeric) 

pop$Year <- excel_numeric_to_date(pop$Year)
pop$Year <- lubridate::year(pop$Year)
names(pop) = gsub(pattern="Estimated Resident Population ;  Female  *", replacement="Age", x=names(pop))
names(pop) = gsub(pattern="* ;", replacement="", x=names(pop))
names(pop) = gsub(pattern="*; *", replacement="", x=names(pop))

pop <- pop %>%
  rename("Age100+" = "Age100 and over") %>%
  mutate("0–14"= Age0+Age1+Age2+Age3+Age4+Age5+Age6+Age7+Age8+Age9+Age10+Age11+Age12+Age13+Age14,
         "15–24" = Age15+Age16+Age17+Age18+Age19+Age20+Age21+Age22+Age23+Age24,
         "25–29" = Age25+Age26+Age27+Age28+Age29,
         "30–34" = Age30+Age31+Age32+Age33+Age34,
         "35–39" = Age35+Age36+Age37+Age38+Age39,
         "40–44" = Age40+Age41+Age42+Age43+Age44,
         "45–49" = Age45+Age46+Age47+Age48+Age49,
         "50–54" = Age50+Age51+Age52+Age53+Age54,
         "55–59" = Age55+Age56+Age57+Age58+Age59,
         "60–64" = Age60+Age61+Age62+Age63+Age64,
         "65–69" = Age65+Age66+Age67+Age68+Age69,
         "70–74" = Age70+Age71+Age72+Age73+Age74,
         "75–79" = Age75+Age76+Age77+Age78+Age79,
         "80–84" = Age80+Age81+Age82+Age83+Age84,
         "85+"   = Age85+Age86+Age87+Age88+Age89+Age90+Age91+Age92+Age93+Age94+Age95+Age96+Age97+Age98+Age99+`Age100+`) %>%
  select(-contains("Age")) %>%
  pivot_longer(cols = contains(c("–","+")),
               names_to = "Age",
               values_to = "pop") %>%
  filter(Year %in% seq(from=2001, to=2022))

#calculate population at risk using AIHW hysterectomy fractions then group into 10-year age groups to match with Cancer data 10-year groups


Year <- c(seq(2001, 2022))
fractions <- read_excel("hysterectomy fractions.xlsx") %>% 
  merge(Year) %>%
  rename(Year = y) %>%
  mutate(fraction = case_when(
    Year %in% c(seq(2001,2004)) ~`Applied to the period 2001-2004`,
    Year %in% c(seq(2005,2015)) ~`Applied to the period 2005-2015`,
    Year %in% c(seq(2016,2022)) ~`Applied to the period 2016-2022`),
    Age = case_when(
      Age =="90+" ~ "85+",
      !Age =="90+" ~ Age)) %>%
  select(Year, Age, fraction)
pop <- pop %>%
  merge(fractions, by=c("Year", "Age")) %>%
  mutate(withouthyst = pop*fraction) %>%
  select(c(Year, Age, pop, withouthyst)) %>%
  pivot_wider(names_from=Age,
              values_from = c(pop, withouthyst)) %>%
  mutate("pop_25–34"= `pop_25–29` + `pop_30–34`,
         "pop_35–44"= `pop_35–39` + `pop_40–44`,
         "pop_45–54"= `pop_45–49` + `pop_50–54`,
         "pop_55–64"= `pop_55–59` + `pop_60–64`,
         "pop_65–74"= `pop_65–69` + `pop_70–74`, 
         "pop_75–84"= `pop_75–79` + `pop_80–84`,
         "withouthyst_25–34"= `withouthyst_25–29`+`withouthyst_30–34`,
         "withouthyst_35–44"= `withouthyst_35–39`+`withouthyst_40–44`,
         "withouthyst_45–54"= `withouthyst_45–49`+`withouthyst_50–54`,
         "withouthyst_55–64"= `withouthyst_55–59`+`withouthyst_60–64`,
         "withouthyst_65–74"= `withouthyst_65–69`+`withouthyst_70–74`,
         "withouthyst_75–84"= `withouthyst_75–79`+`withouthyst_80–84`) %>%
  pivot_longer(cols = contains(c("hyst", "pop")),
               names_to = c(".value", "Age"),
               names_pattern = "(.*)_(.*)") %>%
  filter(Age=="0–14" | Age=="15–24" | Age=="25–34"| Age=="35–44"| Age=="45–54"| Age=="55–64"| Age=="65–74" | Age=="75–84" | Age=="85+") %>%
  group_by(Year) %>%
  arrange(Age, .by_group=T) %>%
  ungroup()

### Separate 2001 population to use as the standard population
standardpop <- pop %>%
  filter(Year==2001) %>%
  select(!Year) %>%
  rename(withouthyststandardpop = withouthyst) %>%
  rename(standardpop= pop)
Age <- dat %>%
  merge(pop,by=c("Year","Age")) %>%
  merge(standardpop,by="Age") %>%
  mutate(
    atrisk= withouthyst+Count,
    atriskstandardpop= withouthyststandardpop+Count,
    agespecificrate=Count/atrisk,
    agespecificrate_nc=Count/pop) 


###calculating age-standardised rate 
# First need the 2001 female total population and at risk population for each age group
totalpop <- Age %>%
  filter(Year==2001) %>%
  summarise(totalpop = sum(pop), .groups = "keep")
totalatrisk <- Age %>%
  filter(Year==2001) %>%
  summarise(totalatrisk = sum(atrisk), .groups = "keep")

agestandardised <- Age %>%
  mutate(
    'rate' =agespecificrate*atriskstandardpop,
    'rate_nc'= agespecificrate_nc*standardpop) %>%
  pivot_longer(c('rate', 'rate_nc'), names_to="corrected", values_to="rate") %>%
  group_by(Year, corrected) %>% 
  summarise(ageadjustedrate = sum(rate), .groups = "keep") %>%
  mutate(ageadjustedrate= case_when(
    corrected=='rate' ~ (ageadjustedrate/totalatrisk$totalatrisk)*100000,
    corrected=='rate_nc' ~ (ageadjustedrate/totalpop$totalpop)*100000)) %>%
  ungroup() %>%
  select(corrected, Year, ageadjustedrate) %>%
  arrange(corrected, Year)

write.csv(agestandardised, "mortality.csv")


agestandardised <- agestandardised %>%
  mutate(Age = "Age–standardised") %>%
  rename(rate = ageadjustedrate) %>%
  select(Age, corrected, Year, rate)
Age <- Age %>%
  mutate(rate=agespecificrate*100000,
         rate_nc= agespecificrate_nc*100000) %>%
  pivot_longer(c(rate, rate_nc), names_to="corrected", values_to="rate") %>% 
  select(Age, corrected, Year, rate)

Age <- rbind(Age, agestandardised) %>%
  arrange(Age, corrected, Year)
write.csv(Age, "mortalityAge.csv")

###########################################################################################

dat %>%
  summarise("Number of deaths" = sum(Count))

dat %>%
  group_by(Age) %>%
  summarise("Number of deaths" = sum(Count))
Age %>%
  filter(Year == 2001 & corrected=="rate")
Age %>%
  filter(Year == 2022 & corrected=="rate")

################################################################################
#####CODE WRITTEN BY REBECCA HARRIS#############################################
#####UNIVERSITY OF QUEENSLAND 2025##############################################
################################################################################
