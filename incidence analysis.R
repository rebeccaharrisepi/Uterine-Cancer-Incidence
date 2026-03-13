
library(readxl)
library(tidyverse)
library(janitor)
library(cowplot)
library(grid)
library(gridExtra)
library(lemon)
library(ggsci)
library(ggh4x)
library(patchwork)
library(cowplot)

##########################################################################################
#SECTION 1: Data set-up and categorising histotypes and calculating age-standardised rates
##########################################################################################
#Before beginning download and have in working directory:
#Cancer datacube 10g: "CDiA-2024-Book-10g-Cancer-incidence-by-histology-vulva-vagina-uterus.xlsx" from the main branch of this repository
#National population data: Population - Australia Population at 30 June, by sex and single year of age, Aust., from 1971 onwards from from the main branch of this repository
#And 'hysterectomyfactions' excel file

dat <- read_excel("CDiA-2024-Book-10g-Cancer-incidence-by-histology-vulva-vagina-uterus.xlsx", 
                  sheet = "Table S10g.3", skip = 4)
dat <- dat %>% janitor::clean_names() %>%
  mutate(
    Subtype= case_when(
      i_d %in% c("U1.01.02.01.01", "U1.01.02.01.02", "U1.01.02.01.03", "U1.01.02.01.04", # Endometrioid adenocarcinomas
                 "U1.01.02.05.01", "U1.01.02.05.02", "U1.01.02.05.03", "U1.01.02.05.04", "U1.01.02.05.05", #mucinous carcinomas
                 "U1.01.02.06.02", "U1.01.02.06.03", "U1.01.02.06.04", "U1.01.03",  "U1.01.02.02.02") ~ "Endometrioid",
      i_d %in% c("U1.01.02.02.01", "U1.01.02.02.03", "U1.01.02.04", "U1.01.06.10") ~ "Serous",
      i_d =="U1.01.02.03" ~ "Clear cell",
      i_d %in% c("U1.01.04.01", "U1.01.04.02", "U1.01.04.03", "U1.01.04.04") ~ "Carcinosarcoma",
      i_d %in% c("U1.02.01.01", "U1.02.02.01", "U1.02.02.02", "U1.02.03.01","U1.02.04.01", "U1.02.04.02", "U1.02.04.03", "U1.02.04.04", "U1.02.05.01", "U1.02.05.02", "U1.02.05.03", "U1.02.05.04", "U1.02.05.05", "U1.02.06.01", "U1.02.06.02", "U1.02.06.03", "U1.02.07", "U1.02.08.01", "U1.02.08.02", "U1.02.08.03", 
                 "U1.02.08.05.01", "U1.02.08.05.02", "U1.02.08.05.03", "U1.02.08.05.04", "U1.02.09.01", "U1.02.10.01", "U1.03") ~ "Sarcoma",
      i_d %in% c("U1.01.01.01", "U1.01.01.02", "U1.01.01.03", "U1.01.01.04", "U1.01.01.05", "U1.01.01.06", "U1.01.02.06.05", "U1.01.02.06.06", "U1.01.02.06.07", "U1.01.05.01", "U1.01.05.02.01", "U1.01.05.02.02", "U1.01.05.02.03", "U1.01.05.03", "U1.01.06.01", "U1.01.06.02", 
                 "U1.01.06.03", "U1.01.06.04", "U1.01.06.05", "U1.01.06.06", "U1.01.06.07", "U1.01.06.08", "U1.01.06.09", "U1.01.06.11", "U1.02.01.02", "U1.02.08.04", "U1.02.08.05.05", "U1.04.01", "U1.04.02", "U1.04.03",
                 "U1.05", "U1.06") ~ "Other",
      is.na(i_d) ~ "Total",
      .default=" ")) %>%
  rename(Count = count_for_cancer_type,
         Age = age_group_years, 
         ALLUterineCancers = total_count_for_uterine_cancer,
         Year=year) %>%
  filter(sex=="Females" & age_grouping=="10-year age groups"
         & !Age %in% c("All ages (crude rate)", "All ages (age-standardised rate - 2024 Australian population)", 'All ages (age-standardised rate - 2001 Australian Standard Population)')) %>%
  select(c(Subtype, cancer_type, Year, Age, Count, ALLUterineCancers, i_d))


### Create different counts of Endometrioid, Serous and clear cell after Re-classifying adenocarcinoma NOS to Endometrioid, Serous and clear cell at percent of cancers that are Endometrioid, Serous or clear cell in each year

proportions <- dat %>%
  filter(Subtype %in% c("Endometrioid", "Serous", "Clear cell")) %>%
  group_by(Year, Subtype) %>%
  summarise(Count = sum(Count)) %>%
  mutate(proportion = Count/(sum(Count))) %>%
  select(-Count)

# Reallocate to endometrioid and serous at specific rates for all years
adeno <- dat %>%
  filter(cancer_type == "Adenocarcinoma, NOS") %>%
  select(cancer_type, Year, Age, n = Count) %>%
  merge(proportions, by = c("Year")) %>%
  mutate(raw = n * proportion,
         floor_val = floor(raw),
         decimal = raw - floor_val) %>%
  group_by(Year, Age) %>%
  mutate(
    remainder = n[1] - sum(floor_val),
    add_one = rank(-decimal, ties.method = "first") <= remainder,
    Count = floor_val + as.integer(add_one)
  ) %>%
  ungroup() %>%
  select(cancer_type, Year, Age, Subtype, Count)


# Merge adenocarcinoma NOS to data and add counts to Endo and Ser subtypes as required
dat <- bind_rows(dat, adeno) %>% # add redistributed counts
  filter(!Subtype==" ") # remove uneeded groupings of subtypes (these are higher level groups e.g. all carcinomas)
rm(adeno)


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
    Year %in% c(seq(2005,2010)) ~`Applied to the period 2005-2010`,
    Year %in% c(seq(2011,2015)) ~`Applied to the period 2011-2015`,
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
  filter(!Age %in% c("0–14", "15–24") & !Subtype=="")



data <- data %>%
  merge(pop,by=c("Year","Age")) %>%
  merge(standardpop,by="Age") %>%
  arrange(Subtype, Year, Age) %>%
  mutate(
    atrisk= withouthyst+Count,
    atriskstandardpop= withouthyststandardpop+Count,
    agespecificrate=Count/atrisk,
    agespecificrate_nc=Count/pop) %>%
  select(Subtype, Year, Age, Count, atrisk, pop, atriskstandardpop, standardpop, agespecificrate, agespecificrate_nc) 

###calculating age-standardised rate 
# First need the 2001 female total population and at risk population for each age group
totalpop <- data %>%
  filter(Year==2001 & Subtype== "Total") %>%
  summarise(totalpop = sum(pop), .groups = "keep")
totalatrisk <- data %>%
  filter(Year==2001 & Subtype== "Total") %>%
  summarise(totalatrisk = sum(atrisk), .groups = "keep")

histotype <- data %>% 
  mutate(
    'corrected for hyst' =agespecificrate*atriskstandardpop,    #age specific rates weighted by standard population
    'not corrected for hyst'= agespecificrate_nc*standardpop) %>%
  pivot_longer(c('corrected for hyst', 'not corrected for hyst'), names_to="corrected", values_to="rate") %>%
  group_by(Subtype, Year, corrected) %>% 
  summarise(
    Count = sum(Count), 
    atrisk = sum(atrisk),
    pop = sum(pop), 
    atriskstandardpop = sum(atriskstandardpop),
    standardpop = sum(standardpop),
    ageadjustedrate = sum(rate)) %>%
  mutate(ageadjustedrate= case_when(
    corrected=='corrected for hyst' ~ (ageadjustedrate/totalatrisk$totalatrisk)*100000,
    corrected=='not corrected for hyst' ~ (ageadjustedrate/totalpop$totalpop)*100000)) %>%
  ungroup()


### export data with event counts to be imported into Joinpoint software
histotype <- histotype %>% 
  select(corrected, Subtype, Year, Count, atrisk, ageadjustedrate) %>%
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
rm(total)
write.csv(Age, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/Age.csv", row.names = T)


PAFs <- data %>%
  filter(Subtype %in% c("Endometrioid", "Serous", "Clear cell", "Total") & Year %in% c(2001, 2010, 2020)) %>%
  select(Subtype, Age, Year, Count, atrisk)
write.csv(PAFs, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/For PAF calculation.csv")
rm(PAFs)

Forattributablerates <- histotype %>%
  filter(Subtype%in% c("Endometrioid", "Clear cell", "Serous", "Total") & Year %in% c(2001, 2010, 2020), corrected=="corrected for hyst")
write.csv(Forattributablerates, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/For attributable rates.csv")
rm(Forattributablerates)

##########################################################
#SECTION 2: Create serous dataset for Sensitivity analyses
##########################################################
###Reassign Carcinosarcoma into a new Serous group

###create new rates of serous uterine cancers after including Carcinosarcoma as serous
seroustotal <- dat %>%
  filter(Subtype %in% c("Serous", "Carcinosarcoma")) %>%
  group_by(Year, Age) %>% 
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
  select(Year, Age, Count, atrisk, pop, atriskstandardpop, standardpop, agespecificrate, agespecificrate_nc)


### and calculate age-standardised rates for new total and serous subgroups
seroustotal <- seroustotal %>% 
  mutate(
    'corrected for hyst' =agespecificrate*atriskstandardpop,
    'not corrected for hyst'= agespecificrate_nc*standardpop) %>%
  pivot_longer(c('corrected for hyst', 'not corrected for hyst'), names_to="corrected", values_to="rate") %>%
  group_by(Year, corrected) %>% 
  summarise(
    Count = sum(Count),
    atrisk = sum(atrisk),
    atriskstandardpop = sum(atriskstandardpop),
    standardpop = sum(standardpop),
    ageadjustedrate = sum(rate)
    ) %>%
  mutate(ageadjustedrate= case_when(
    corrected=='corrected for hyst' ~ (ageadjustedrate/totalatrisk$totalatrisk)*100000,
    corrected=='not corrected for hyst' ~ (ageadjustedrate/totalpop$totalpop)*100000)) %>%
  ungroup() %>%
  select(corrected, Year, ageadjustedrate) %>%
  arrange(corrected, Year)


### export data with event counts to be imported into Joinpoint software
write.csv(seroustotal, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/serous.csv", row.names = T)
rm(seroustotal)

##############################################################################
# SECTION 3: Preliminary Data Analysis and Plots of cancer counts and rates
##############################################################################

### All uterine cancer counts by histological subtype from 2001-2020
histologicalsubtype <- dat %>%
  filter(!Age=="0–14" & !Age=="15–24") %>%
  group_by(cancer_type, i_d) %>% 
  summarise("Number of cases" = sum(Count)) 
rm(histologicalsubtype)


### All uterine cancer counts by Subtype from 2001-2020
dat %>%
  filter(!Age=="0–14" & !Age=="15–24") %>%
  group_by(Subtype) %>% 
  summarise("Number of cases" = sum(Count))
### All uterine cancer counts by histotype from 2006-2020
dat %>%
  filter(!Age=="0–14" & !Age=="15–24" &Year>2005) %>%
  group_by(Subtype) %>% 
  summarise("Number of cases" = sum(Count))


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
dat %>%
  filter(!Age=="0–14" & !Age=="15–24" & Year==2020) %>%
  group_by(Subtype) %>% 
  summarise("Number of cases" = sum(Count))

data %>%
  group_by(Subtype, Year) %>% 
  summarise(Count = sum(Count), .groups = "keep") %>%
  ggplot(aes(x=Year, y=Count, group=Subtype, color=Subtype)) +
  geom_line() +
  facet_wrap(~Subtype, scales="free") + 
  theme_minimal()


### plotting age specific rates per 100, 000 women
data %>%
  mutate(rate=agespecificrate*100000,
         rate_nc= agespecificrate_nc*100000) %>%
  pivot_longer(c(rate, rate_nc), names_to="corrected", values_to="rate") %>%
  ggplot(aes(x=Year, y=rate, group=Subtype, color=Subtype)) +
  geom_line() + geom_point() +
  facet_wrap(~Age + corrected) +
  labs(y="age specific rate") + 
  theme_minimal() 



###plot counts of Adenocarcinoma NOS over time

adeno <- dat %>%
  group_by(Year, Age) %>%
  fill(ALLUterineCancers) %>%
  filter(cancer_type =="Adenocarcinoma, NOS") %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    ALLUterineCancers = first(ALLUterineCancers),
    .groups = "drop"  ) %>%
  group_by(Year) %>%
  summarise(
    Count = sum(Count),
    ALLUterineCancers = sum(ALLUterineCancers)) %>%
  mutate(prop = (Count/ALLUterineCancers)*100) 

lancet_col <- pal_lancet()(1)
p <- adeno %>%
  ggplot(aes(x=Year, y=prop)) + 
  geom_line(color = pal_lancet()(2)[1], linewidth = 2) +
  labs(y= "Percent of total") +
  scale_y_continuous(
    breaks = seq(0,25, by=5)) +
  theme_cowplot(12) +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16, margin = margin(r = 30)),       # bigger y-axis label
    axis.title.x = element_text(size = 16, margin = margin(t=20)),        # bigger x-axis label
    axis.text.x = element_text(size = 14),                              # bigger x-axis tick text
    axis.text.y = element_text(size = 14))                               # bigger y-axis tick text
p
ggsave("Adenoproportion.jpeg", plot=p, width=12, height=9, dpi=300)
rm(adeno)


p <- proportions %>%
  mutate(
    proportion = proportion *100,
    Subtype = factor(Subtype, levels = c("Endometrioid", "Serous", "Clear cell"))  ) %>%
  ggplot(aes(x = Year, y = proportion, group = Subtype, color = Subtype)) +
  geom_line(linewidth = 1.5) +
  ylab("Percent of total") +
  scale_color_lancet() +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  theme_cowplot(12) +
  theme(
    axis.title.y = element_text(
      angle = 0,           # horizontal
      vjust = 0.5,         # vertical alignment
      hjust = 1,           # horizontal alignment
      margin = margin(r = 30)), # add space to the right of text, moves it left
    axis.title.x = element_text(
      margin = margin(t=20)))  +
  labs(color = "")
p
ggsave("Histotype proportion.jpeg", plot=p, width=9, height=6, dpi=300)
rm(proportions)
####################################################################################################################################################
# SECTION 4: Plotting Age-standardised rate per 100, 000 women corrected and not corrected for hysterectomy alongside trends from joinpoint analysis
####################################################################################################################################################
############
# By SUBTYPE
############
joinpoint <- read_excel("histotyperesults.xlsx") %>%
  filter(Year %in% c(2001, 2020) | !is.na(`Joinpoint Location`)) %>%
  mutate(Subtype = word(Cohort, 1),
         Subtype = case_when(
           Subtype=="Clear" ~ "Clear cell",
           .default =Subtype),
         corrected = case_when(
           str_detect(Cohort, "/ corrected") ~ "Modeled hysterectomy-corrected trend",
           str_detect(Cohort, "/ not") ~ "Modeled hysterectomy-uncorrected trend"),
         observed="Modeled trend",
         group = paste(Subtype, corrected)) %>%
  select(!Cohort)

joins <- joinpoint %>%
  filter(!is.na(`Joinpoint Location`)) %>%
  mutate(corrected ="Joinpoints")

APC <- read_excel("histotypeAAPC.xlsx") %>%
  mutate(
    corrected = case_when(
      str_detect(Cohort, "/ corrected") ~ "Modeled hysterectomy-corrected trend",
      str_detect(Cohort, "/ not") ~ "Modeled hysterectomy-uncorrected trend"),
    Subtype = str_remove(Cohort, " /.*"),
    signif = ifelse(str_detect(AAPC, "\\*"), "*", "")) %>%
  mutate(
    AAPC = str_remove(AAPC, "\\*$")) %>%
  mutate(
    AAPC =as.numeric(AAPC))
APC$AAPC <- format(round(APC$AAPC, 2), nsmall = 2, trim=TRUE)
APC <- APC %>%
  mutate(label = paste0("APC=",AAPC,signif)) %>%
  mutate(
    label = label %>%
      if_else(str_detect(., "\\*$"), ., paste0(., " ")) %>% # Add space if it doesn't end with "*"
      str_replace("=", function(x) if_else(str_detect(., "-"), x, "= "))) %>% # Add space after "=" if it doesn't contain "-"
  select(Subtype, corrected, label) 

joinpoint <- joinpoint %>%
  merge(APC, by=c("Subtype","corrected"))
labels <- joinpoint %>%
  filter(Year==2020) %>%
  mutate(
    `Modeled ageadjustedrate` = case_when(
      corrected=="Modeled hysterectomy-corrected trend" & Subtype %in% c("Endometrioid") ~ `Modeled ageadjustedrate` +3,
      corrected=="Modeled hysterectomy-corrected trend" & Subtype %in% c("Total") ~ `Modeled ageadjustedrate` +1,
      corrected=="Modeled hysterectomy-corrected trend" & !Subtype %in% c("Total", "Endometrioid")~ `Modeled ageadjustedrate` +0.3,
      corrected=="Modeled hysterectomy-uncorrected trend" & Subtype %in% c("Total", "Endometrioid") ~ `Modeled ageadjustedrate` -1,
      corrected=="Modeled hysterectomy-uncorrected trend" & !Subtype %in% c("Total", "Endometrioid")~ `Modeled ageadjustedrate` -0.3),
    yaxis = case_when(
      corrected=="Modeled hysterectomy-corrected trend" & Subtype %in% c("Endometrioid", "Total") ~ `Modeled ageadjustedrate`,
      corrected=="Modeled hysterectomy-corrected trend" & Subtype %in% c("Serous") ~ `Modeled ageadjustedrate` - 0.3,
      corrected=="Modeled hysterectomy-corrected trend" & Subtype %in% c("Clear cell", "Sarcoma", "Carcinosarcoma", "Other") ~ `Modeled ageadjustedrate`)) %>%
  group_by(Subtype) %>%
  mutate(
    yaxis = case_when(
      corrected == "Modeled hysterectomy-corrected trend" ~ yaxis,
      corrected == "Modeled hysterectomy-uncorrected trend" & Subtype %in% c("Total", "Endometrioid") ~ yaxis[!is.na(yaxis)][1] -8, # take the corrected value for this subtype
      corrected == "Modeled hysterectomy-uncorrected trend" & !Subtype %in% c("Total", "Endometrioid") ~ yaxis[!is.na(yaxis)][1] -0.6,
      TRUE ~ `Modeled ageadjustedrate`
    )
  ) %>%
  ungroup()

histotype <- histotype %>%
  mutate(corrected = case_when(
    corrected == "corrected for hyst" ~ "Observed hysterectomy-corrected rate",
    corrected == "not corrected for hyst" ~ "Observed hysterectomy-uncorrected rate"))




lancet_cols <- pal_lancet()(2)
big_subtypes <- c("Total", "Endometrioid")
small_subtypes <- c("Serous", "Carcinosarcoma", "Sarcoma", "Clear cell", "Other")

p1 <- ggplot() +
  geom_point(data=histotype %>% filter(Subtype %in% big_subtypes), aes(x=Year, y=ageadjustedrate, shape=corrected, color=corrected)) +  
  ylab("") + 
  facet_wrap(~factor(Subtype, levels=c("Total", "Endometrioid")), scales="fixed", ncol=2, axes="all") +
  theme_cowplot(12) +
  theme(
    legend.title = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_blank(),
    panel.spacing = unit(0.3, "cm"),
    panel.spacing.x = unit(0.5, "cm"),
    strip.text = element_text(margin = margin(b=-13), hjust=0.5),
    strip.clip = "off",
    plot.margin = margin(t = 5, r = 5, b = -2, l = 5, unit = "pt")) +       
  geom_line(data=joinpoint %>% filter(Subtype %in% big_subtypes), aes(x=Year, y=`Modeled ageadjustedrate`, color=corrected)) +
  geom_point(data=joins %>% filter(Subtype %in% big_subtypes), aes(x=Year, y=`Modeled ageadjustedrate`, shape=corrected, color=corrected), size=2) +
  geom_label(data=labels %>% filter(Subtype %in% big_subtypes), aes(x=Year, y=yaxis, label=label, color=corrected), fontface="bold", size=4, nudge_x=1, nudge_y=0.02, hjust = 0, show.legend=F) +
  scale_x_continuous(breaks =c(2000, 2005, 2010, 2015, 2020), limits=c(2000, 2025)) +
  scale_y_continuous(limits = c(0, 40), expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values = c(
    "Observed hysterectomy-corrected rate"   = lancet_cols[1],
    "Observed hysterectomy-uncorrected rate" = lancet_cols[2],
    "Modeled hysterectomy-corrected trend"   = lancet_cols[1],
    "Modeled hysterectomy-uncorrected trend" = lancet_cols[2],
    "labels corrected" = "red3",
    "labels uncorrected" = "royalblue4",
    "Joinpoints" = "black")) +
  scale_shape_manual(values = c(
    "Observed hysterectomy-corrected rate"   = 16,
    "Observed hysterectomy-uncorrected rate" = 17,
    "Modeled hysterectomy-corrected trend"   = NA,
    "Modeled hysterectomy-uncorrected trend" = NA,
    "Joinpoints" = 15
  )) + 
  guides(
    color = guide_legend(
      override.aes = list(
        shape = c(15, NA, NA, 16, 17),  # shapes for points/lines in legend
        linetype = c(NA, 1, 1, NA, NA)  # lines for trends
      )
    ),
    shape = "none"  # hides the separate shape legend
  ) +
  labs(x = NULL) + theme(legend.position ="none")

p1


p2 <- ggplot() +
  geom_point(data=histotype %>% filter(Subtype %in% small_subtypes), aes(x=Year, y=ageadjustedrate, shape=corrected, color=corrected)) +  
  labs(y="Rate (per 100, 000 women)") +
  facet_wrap(~factor(Subtype, levels=c("Serous", "Clear cell", "Sarcoma", "Carcinosarcoma", "Other")), scales="fixed", ncol=2, axes="all") +
  theme_cowplot(12) +
  theme(
    legend.title = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    panel.spacing.x = unit(0.5, "cm"),
    strip.placement = "inside",
    strip.text = element_text(margin = margin(b=-13), hjust=0.5),
    strip.clip = "off",
    axis.title.y = element_text(
      hjust = 0.75,
      vjust = 1.5,
      margin = margin(r=15)),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5, unit = "cm")) +
  geom_line(data=joinpoint %>% filter(Subtype %in% small_subtypes), aes(x=Year, y=`Modeled ageadjustedrate`, color=corrected)) +
  geom_point(data=joins %>% filter(Subtype %in% small_subtypes), aes(x=Year, y=`Modeled ageadjustedrate`, shape=corrected, color=corrected), size=2) +
  geom_label(data=labels %>% filter(Subtype %in% small_subtypes), aes(x=Year, y=yaxis, label=label, color=corrected), fontface="bold", size=4, nudge_x=1, nudge_y=0.02, hjust = 0, show.legend=F) +
  scale_x_continuous(breaks =c(2000, 2005, 2010, 2015, 2020), limits=c(2000, 2025)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values = c(
    "Observed hysterectomy-corrected rate"   = lancet_cols[1],
    "Observed hysterectomy-uncorrected rate" = lancet_cols[2],
    "Modeled hysterectomy-corrected trend"   = lancet_cols[1],
    "Modeled hysterectomy-uncorrected trend" = lancet_cols[2],
    "labels corrected" = "red3",
    "labels uncorrected" = "royalblue4",
    "Joinpoints" = "black"
  )) +
  scale_shape_manual(values = c(
    "Observed hysterectomy-corrected rate"   = 16,
    "Observed hysterectomy-uncorrected rate" = 17,
    "Modeled hysterectomy-corrected trend"   = NA,
    "Modeled hysterectomy-uncorrected trend" = NA,
    "Joinpoints" = 15
  )) + 
  guides(
    color = guide_legend(
      override.aes = list(
        shape = c(15, NA, NA, 16, 17),  # shapes for points/lines in legend
        linetype = c(NA, 1, 1, NA, NA)  # lines for trends
      )
    ),
    shape = "none"  # hides the separate shape legend
  ) 

legend <- get_legend(p2)
p2 <- p2 + theme(legend.position ="none")

legend_grob <- wrap_elements(legend)

final_plot <- p1 / p2  +
  plot_layout(heights = c(1, 4)) + # adjust if needed
  inset_element(
    legend_grob,
    left = 0.65,   # move horizontally
    bottom = 0.00, # move lower vertically
    right = 0.95,
    top = 0.30
  )
final_plot




ggsave("Histotype.jpeg", plot=final_plot, width=15, height=9, dpi=300)



###############
# BY AGE GROUP
###############

age_levels <- c("Age–standardised", "25–34", "35–44", "45–54", "55–64", "65–74", "75–84", "85+")

# Get the Lancet palette with the total number of levels
lancet_colors <- pal_lancet()(length(age_levels))

# Create named vector: first color for "Age–standardised", then others for rest
custom_colors <- c("Age–standardised" = lancet_colors[1], setNames(lancet_colors[-1], age_levels[-1]))


### By Age
joinpoint <- read_excel("Ageresults.xlsx") %>%
  rename(Joinpoint = `Joinpoint Location`) %>%
  filter(Year %in% c(2001, 2020) | !is.na(Joinpoint)) %>%
  mutate(Age = str_remove(Cohort, " -.*"),
         observed="Modeled Trend")  %>%
  select(!Cohort) %>%
  mutate(Age = str_replace_all(Age, "â€“", "–"))

joins <- joinpoint %>%
  filter(!is.na(Joinpoint)) %>%
  mutate(observed="Joinpoint")

APC <- read_excel("AgeAAPC.xlsx") %>%
  mutate(Age = str_remove(Cohort, " -.*"),
         signif = ifelse(str_detect(AAPC, "\\*"), "*", "")) %>%
  mutate(Age = str_replace_all(Age, "â€“", "–")) %>% 
  mutate(
    AAPC = str_remove(AAPC, "\\*$")) %>%
  mutate(
    AAPC =as.numeric(AAPC))
APC$AAPC <- format(round(APC$AAPC, 2), nsmall = 2)
APC <- APC %>%
  mutate(label = paste(AAPC,signif)) 


joinpoint <- joinpoint %>%
  merge(APC, by=c("Age"))
labels <- joinpoint %>%
  filter(Year==2020) %>%
  mutate(Year = Year -0.5)
#%>%
#mutate(`Modeled ageadjustedrate` = case_when(
#corrected=="Modeled hysterectomy-uncorrected trend" & Subtype%in% c("Clear cell", "Other") ~ `Modeled agespecificrate` -0.05,
#.default=`Modeled ageadjustedrate`))


Age <- Age %>%
  mutate(observed="Observed Rate")

p <- ggplot() +
  geom_point(data=Age, aes(x=Year, y=agespecificrate, pch=observed, color=Age)) +  
  theme_cowplot(12) +
  ylab("Rate (per 100, 000 women)") + 
  facet_wrap(~factor(Age, levels=c("Age–standardised", "25–34", "35–44", "45–54", "55–64", "65–74", "75–84", "85+")), scales="free") +
  theme_cowplot(12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_blank(),
    panel.spacing = unit(1, "lines")) +
  theme(legend.title = element_blank()) + 
  geom_line(data=joinpoint, aes(x=Year, y=`Modeled agespecificrate`, group=observed, color=Age, linetype=observed)) +
  geom_point(data=joins, aes(x=Year, y=`Modeled agespecificrate`, pch=observed), color="black", size=2) +
  geom_label(data=labels, aes(x=Year, y=`Modeled agespecificrate`, label=label, color=Age), fontface="bold", cex=4, nudge_x=1, nudge_y=0.02, hjust = 0) +
  scale_x_continuous(breaks =c(2000, 2005, 2010, 2015, 2020), limits = c(2000, 2023)) +
  guides(color = "none") +
  scale_shape_manual(values=c("Joinpoint"=15, "Observed Rate"=16, "Modeled Trend"=19)) + 
  scale_color_manual(values = custom_colors)

p <- reposition_legend(p, 'center', panel='panel-3-3')

p
ggsave("Age.jpeg", plot=p, width=15, height=9, dpi=300)



####################################################################################################################################################
# SECTION 5: Calculating age standardised rates for <50 and >50 years
####################################################################################################################################################


### generate Dataset with event counts for each subtype and age group 
earlyvslate <- dat %>% 
  filter(Subtype %in% c("Endometrioid", "Serous", "Clear cell", "Total") & !Age %in% c("0–14", "15–24")) %>% 
  group_by(Subtype, Year, Age) %>% 
  summarise(Count = sum(Count), .groups = "keep")

earlyvslate <- earlyvslate %>%
  merge(pop,by=c("Year","Age"), all.x=T, all.y=F) %>%
  merge(standardpop,by="Age", all.x=T, all.y=F) %>%
  arrange(Subtype, Year, Age) %>%
  group_by(Subtype, Year, Age) %>% 
  mutate(
    atrisk= withouthyst+Count,
    atriskstandardpop= withouthyststandardpop+Count,
    agespecificrate=Count/atrisk,
    agespecificrate_nc=Count/pop) %>%
  select(Subtype, Year, Age, Count, atrisk, pop, atriskstandardpop, standardpop, agespecificrate, agespecificrate_nc) %>%
  mutate(age_group = case_when(
  Age %in% c("25–34", "35–44", "45–54") ~ "<55",
  .default="55+")) 




###calculating age-standardised rate 
# First need the 2001 female total population and at risk population for each age group
totalpop <- earlyvslate %>%
  filter(Year==2001 & Subtype== "Total") %>%
  ungroup() %>%
  summarise(totalpop = sum(pop), .groups = "keep")
totalatrisk <- earlyvslate %>%
  filter(Year==2001 & Subtype== "Total") %>%
  ungroup() %>%
  summarise(totalatrisk = sum(atrisk), .groups = "keep")


earlyvslate <- earlyvslate %>% 
  mutate(
    'corrected for hyst' =agespecificrate*atriskstandardpop,
    'not corrected for hyst'= agespecificrate_nc*standardpop) %>%
  pivot_longer(c('corrected for hyst', 'not corrected for hyst'), names_to="corrected", values_to="rate") %>%
  group_by(Subtype, Year, age_group, corrected) %>% 
  summarise(
    Count=sum(Count),
    atrisk = sum(atrisk),
    pop = sum(pop),
    atriskstandardpop = sum(atriskstandardpop),
    standardpop=sum(standardpop),
    ageadjustedrate = sum(rate)) %>%
  mutate(ageadjustedrate= case_when(
    corrected=='corrected for hyst' ~ ageadjustedrate/(totalatrisk$totalatrisk)*100000,
    corrected=='not corrected for hyst' ~ ageadjustedrate/(totalpop$totalpop)*100000
)) %>%
  ungroup()


### export data with event counts to be imported into Joinpoint software
earlyvslate <- earlyvslate %>%
  filter(corrected=="corrected for hyst") %>% 
  select(corrected, Subtype, age_group, Year, Count, atrisk, ageadjustedrate) %>%
  arrange(corrected, Subtype, age_group) 
write.csv(earlyvslate, "C:/Users/uqrhar23/OneDrive - The University of Queensland/Documents/Uterine cancer/Joinpoint/earlyonset.csv", row.names = T)




################################################################################
#####CODE WRITTEN BY REBECCA HARRIS#############################################
#####UNIVERSITY OF QUEENSLAND 2025##############################################
################################################################################
