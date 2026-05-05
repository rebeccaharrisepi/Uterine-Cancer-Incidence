###############
# Load packages
###############
library(readxl)
library(tidyverse)
library(janitor)
library(cowplot)


###############################
# Section 1: Mortality analysis
###############################

data<- read_excel("CDiA-2024-Book-2a-Cancer-mortality-and-age-standardised-rates-by-age-5-year-groups.xlsx", sheet = "Table S2a.1", skip = 5) %>%
  janitor::clean_names() %>%
  filter(cancer_group_site=="Uterine cancer" & sex=="Females" & data_type=="Actual" & year>2000 & !age_group_years %in% c("00-04", "05-09", "10-14", "15-19", "20-24", "All ages combined")) %>%
  rename(Year= year,
         Age=age_group_years,
         Count = count) %>%
  select(Year, Age, Count)  
data$Age <- gsub("-", "–", data$Age)
data <- data %>%
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


### Import population data and calculate 10-year age groups
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
    Year %in% c(seq(2005,2010)) ~`Applied to the period 2005-2010`,
    Year %in% c(seq(2011,2015)) ~`Applied to the period 2011-2015`,
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
Age <- data %>%
  merge(pop,by=c("Year","Age")) %>%
  merge(standardpop,by="Age") %>%
  mutate(
    atrisk= withouthyst+Count,
    atriskstandardpop= withouthyststandardpop+Count,
    agespecificrate=Count/atrisk,
    agespecificrate_nc=Count/pop) %>%
  arrange(Year, Age)


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

rm(totalatrisk, totalpop, standardpop, pop, fractions, Year)
##################################
# Section 2: Summary of uterine cancer deaths
##################################

data %>%
  summarise("Number of deaths" = sum(Count))

data %>%
  group_by(Age) %>%
  summarise("Number of deaths" = sum(Count))


Age %>%
  filter(Year == 2001 & corrected=="rate")

Age %>%
  filter(Year == 2022 & corrected=="rate")


##################################################################################################
# Section 3: Plotting hysterectomy-corrected mortality rates for each age group and overall (Age–standardised)
##################################################################################################

options(OutDec = "·")

joinpoint <- read_excel("Mortalityresults.xlsx") %>%
  rename(Joinpoint = `Joinpoint Location`) %>%
  filter(Year %in% c(2001, 2022) | !is.na(Joinpoint)) %>%
  mutate(Age = str_remove(Cohort, " /.*"),
         observed="Modeled Trend",
         Joinpoint = case_when(
           !is.na(Joinpoint) ~ "Joinpoint "))  %>%
  select(!Cohort)
joinpoint$Age <- gsub("â€“", "–", joinpoint$Age)
joins <- joinpoint %>%
  filter(!is.na(Joinpoint)) %>%
  mutate(observed="Joinpoint")

APC <- read_excel("MortalityAAPC.xlsx") %>%
  mutate(Age = str_remove(Cohort, " /.*"),
         signif = ifelse(str_detect(AAPC, "\\*"), "*", " "),
         AAPC = str_remove(AAPC, "\\*$"),
         AAPC = as.numeric(AAPC))
APC$AAPC <- round(APC$AAPC, digits = 2) 
APC <- APC %>%
  mutate(label = paste(AAPC,signif)) %>%
  select(Age, label) 

APC$Age <- gsub("â€“", "–", APC$Age) 

joinpoint <- joinpoint %>%
  merge(APC, by="Age")
labels <- joinpoint %>%
  filter(Year==2022)

Age <- Age %>%
  mutate(observed="Observed Rate") %>%
  filter(corrected=="rate")

age_levels <- c("Age–standardised", "25–34", "35–44", "45–54", "55–64", "65–74", "75–84", "85+")

# Get the Lancet palette with the total number of levels
lancet_colors <- pal_lancet()(length(age_levels))

# Create named vector: first color for "Age–standardised", then others for rest
custom_colors <- c("Age–standardised" = lancet_colors[1], setNames(lancet_colors[-1], age_levels[-1]))



p <- ggplot() +
  geom_point(data=Age, aes(x=Year, y=rate, color=Age, pch=observed)) +
  ylab("Rate (per 100, 000 women)") +
  theme_cowplot(12) + 
  facet_wrap(
    ~factor(Age, 
            levels = c("Age–standardised", "25–34", "35–44", 
                       "45–54", "55–64", "65–74", 
                       "75–84", "85+")),
    scales="free",
    nrow = 4,
    labeller = as_labeller(
      c(
        "Age–standardised" = "'Age–standardised'",
        "25–34" = "'25–34'^a",
        "35–44" = "'35–44'",
        "45–54" = "'45–54'",
        "55–64" = "'55–64'",
        "65–74" = "'65–74'",
        "75–84" = "'75–84'",
        "85+" = "'85+'"
      ),
      label_parsed))+
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks =c(2000, 2005, 2010, 2015, 2020), limits = c(2000, 2025)) +
  geom_line(data=joinpoint, aes(x=Year, y=`Modeled rate`, group=observed, linetype=observed, color=Age), linewidth = 1) +
  geom_point(data=joins, aes(x=Year, y=`Modeled rate`, pch="Joinpoints"), color="black", size=2) +
  geom_label(data=labels, aes(x=Year, y=`Modeled rate`, label=label), fontface="bold", cex=4, nudge_x=1, nudge_y=0.02, hjust = 0) +
  scale_shape_manual(values=c("Joinpoints"=16, "Observed Rate"=17, "Modeled Trend"=19)) +
  guides(color="none") + 
  scale_color_manual(values = custom_colors)

p

ggsave("Mortality.jpeg", plot=p, width=15, height=9, dpi=300)



###########################################
# SECTION 4: PAF calculations for mortality
###########################################

library(deltapif)

# PAFs assuming 10-year latent period

dat <- data %>%
  select(Age, Year)

prev <- read_excel("OW OB prevalence.xlsx")
RRs <- read_excel("Relative risks_Mortality.xlsx")

dat <- dat %>%
  merge(prev, by = c("Year", "Age")) %>% 
  merge(RRs, by="Exposure") %>%
  arrange(Year, Age) %>%
  mutate(P = Prevalence / 100) %>%
  group_by(Year, Age) %>%
  summarise(
    P_vec = list(P),
    RR_vec = list(RR),
    RRLL_vec = list(RRLL),
    RRUL_vec = list(RRUL),
    .groups = "drop"
  )
results_df <- dat %>%
  pmap_dfr(function(Year, Age, P_vec, RR_vec, RRLL_vec, RRUL_vec) {
    
    # Convert lists to vectors
    P_vec <- unlist(P_vec)
    RR_vec <- unlist(RR_vec)
    RRLL_vec <- unlist(RRLL_vec)
    RRUL_vec <- unlist(RRUL_vec)
    
    # Compute log-RR and variance
    beta <- log(RR_vec)
    var_beta <- ((log(RRUL_vec) - log(RRLL_vec)) / (2*1.96))^2
    var_beta[1] <- 0
    var_beta <- diag(var_beta)
    
    # Compute PAF
    res <- paf(p = P_vec, beta = beta, var_beta = var_beta)
    
    # Extract point estimate and CI
    PAF_val <- coef(res)
    CI <- confint(res)
    
    data.frame(
      Age = Age,
      Year = Year,
      PAF = round(PAF_val, 2),
      PAF_LL = round(CI[1], 2),
      PAF_UL = round(CI[2], 2)
    )
  })
results_df


Nobs <- data %>%
  select(Age, Year, Count) %>%
  rename(Nobs = Count)

results_df <- results_df %>%
  merge(Nobs, by=c("Year", "Age")) %>%
  mutate(
    Nexc = round((Nobs*PAF), 0),
    Nexc_LL = round((Nobs*PAF_LL),0),
    Nexc_UL = round((Nobs*PAF_UL),0))

totals <- results_df %>%
  group_by(Year) %>%
  summarise(
    Age = "Total",
    Nobs = sum(Nobs, na.rm = TRUE),
    Nexc = sum(Nexc, na.rm = TRUE),
    Nexc_LL = sum(Nexc_LL, na.rm = TRUE),
    Nexc_UL = sum(Nexc_UL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(PAF = Nexc/Nobs,
         PAF_LL = Nexc_LL/Nobs,
         PAF_UL = Nexc_UL/Nobs)
results_df <- bind_rows(results_df, totals)
results_df


results_wide <- results_df %>%
  mutate(
    PAF_pct = round(PAF * 100, 0),
    LL_pct = round(PAF_LL * 100, 0),
    UL_pct = round(PAF_UL * 100, 0),
    PAF_CI = paste0(PAF_pct, " (", LL_pct, ", ", UL_pct, ")"),
    Nexc_CI = paste0("(", Nexc_LL, ", ", Nexc_UL, ")")
  ) %>%
  select(Age, Year, PAF_CI, Nobs, Nexc, Nexc_CI) %>%
  pivot_wider(
    names_from = Year,
    values_from = c(PAF_CI, Nobs, Nexc, Nexc_CI),
    names_glue = "{Year}_{.value}"
  )

# Reorder columns by year 
years <- c(2001, 2010, 2020)

results_wide <- results_wide %>%
  select(
    Age,
    unlist(lapply(years, function(y) {
      c(
        paste0(y, "_PAF_CI"),
        paste0(y, "_Nobs"),
        paste0(y, "_Nexc"),
        paste0(y, "_Nexc_CI")
      )
    }))
  ) %>%
  arrange(Age)
results_wide

write.csv(results_wide, "PAF mortality.csv")


# Sensitivity analysis assuming 5-year latent period

dat <- data %>%
  select(Age, Year)

prev <- read_excel("OW OB prevalence.xlsx")
RRs <- read_excel("Relative risks_Mortality.xlsx")

dat <- dat %>%
  merge(prev, by = c("Year", "Age")) %>% 
  merge(RRs, by="Exposure") %>%
  arrange(Year, Age) %>%
  mutate(P = Prevalence_5year / 100) %>%
  group_by(Year, Age) %>%
  summarise(
    P_vec = list(P),
    RR_vec = list(RR),
    RRLL_vec = list(RRLL),
    RRUL_vec = list(RRUL),
    .groups = "drop"
  )
results_df <- dat %>%
  pmap_dfr(function(Year, Age, P_vec, RR_vec, RRLL_vec, RRUL_vec) {
    
    # Convert lists to vectors
    P_vec <- unlist(P_vec)
    RR_vec <- unlist(RR_vec)
    RRLL_vec <- unlist(RRLL_vec)
    RRUL_vec <- unlist(RRUL_vec)
    
    # Compute log-RR and variance
    beta <- log(RR_vec)
    var_beta <- ((log(RRUL_vec) - log(RRLL_vec)) / (2*1.96))^2
    var_beta[1] <- 0
    var_beta <- diag(var_beta)
    
    # Compute PAF
    res <- paf(p = P_vec, beta = beta, var_beta = var_beta)
    
    # Extract point estimate and CI
    PAF_val <- coef(res)
    CI <- confint(res)
    
    data.frame(
      Age = Age,
      Year = Year,
      PAF = round(PAF_val, 2),
      PAF_LL = round(CI[1], 2),
      PAF_UL = round(CI[2], 2)
    )
  })
results_df


Nobs <- data %>%
  select(Age, Year, Count) %>%
  rename(Nobs = Count)

results_df <- results_df %>%
  merge(Nobs, by=c("Year", "Age")) %>%
  mutate(
    Nexc = round((Nobs*PAF), 0),
    Nexc_LL = round((Nobs*PAF_LL),0),
    Nexc_UL = round((Nobs*PAF_UL),0))

totals <- results_df %>%
  group_by(Year) %>%
  summarise(
    Age = "Total",
    Nobs = sum(Nobs, na.rm = TRUE),
    Nexc = sum(Nexc, na.rm = TRUE),
    Nexc_LL = sum(Nexc_LL, na.rm = TRUE),
    Nexc_UL = sum(Nexc_UL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(PAF = Nexc/Nobs,
         PAF_LL = Nexc_LL/Nobs,
         PAF_UL = Nexc_UL/Nobs)
results_df <- bind_rows(results_df, totals)
results_df


results_wide <- results_df %>%
  mutate(
    PAF_pct = round(PAF * 100, 0),
    LL_pct = round(PAF_LL * 100, 0),
    UL_pct = round(PAF_UL * 100, 0),
    PAF_CI = paste0(PAF_pct, " (", LL_pct, ", ", UL_pct, ")"),
    Nexc_CI = paste0("(", Nexc_LL, ", ", Nexc_UL, ")")
  ) %>%
  select(Age, Year, PAF_CI, Nobs, Nexc, Nexc_CI) %>%
  pivot_wider(
    names_from = Year,
    values_from = c(PAF_CI, Nobs, Nexc, Nexc_CI),
    names_glue = "{Year}_{.value}"
  )

# Reorder columns by year 
years <- c(2001, 2010, 2020)

results_wide <- results_wide %>%
  select(
    Age,
    unlist(lapply(years, function(y) {
      c(
        paste0(y, "_PAF_CI"),
        paste0(y, "_Nobs"),
        paste0(y, "_Nexc"),
        paste0(y, "_Nexc_CI")
      )
    }))
  ) %>%
  arrange(Age)
results_wide


write.csv(results_wide, "PAF mortality 5 year lag.csv")




################################################################################
#####CODE WRITTEN BY REBECCA HARRIS#############################################
#####UNIVERSITY OF QUEENSLAND 2025##############################################
################################################################################

