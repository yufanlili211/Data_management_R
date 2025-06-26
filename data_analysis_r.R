# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Examination project in R 
# Author: Yufan Yao
# Email: yufanyao426@gmail.com
# Submission date: 2025/6/12
# Version: 3
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  

## install and library packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)

# 1 Data Management ---------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: data import, variable assignment, dataset reorganisation (merge + long format),
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
setwd("/Users/yaoyufan")
## 1.1 data import ----
pk_data <- read.csv(file = "Desktop/semester1/R_exam/BPI889_PK_44.csv", header = TRUE) # import pk data assigned to the name of pk_data
cov_data <- read.table(file = "Desktop/semester1/R_exam/BPI889_cov_44.txt", header = TRUE, sep = ",") # import cov data seperated by comma

## 1.2 variables assignment and data classification ----
### PK ----
pk_data <- pk_data %>% rename(ID = pat) # change the name of column from pat to ID
# Pivot to long format
long_pk_data <- pivot_longer(pk_data, 
                             cols = starts_with("X"),
                             names_to = "Time(h)", 
                             values_to = "C(mg/L)") # concentration
long_pk_data$`Time(h)` <- sub("^X", "", long_pk_data$`Time(h)`) # remove X in the column of time
long_pk_data$`Time(h)` <- sub("h", "", long_pk_data$`Time(h)`) # remove h
long_pk_data <- long_pk_data %>%
  mutate(across(-ID, as.numeric)) # convert concentration and time to numeric data
str(long_pk_data) # check data classification
long_pk_data$`C(mg/L)`[long_pk_data$`C(mg/L)` == "."] <- NA # replace missing value shown by "." to NA
long_pk_data <- long_pk_data[!is.na(long_pk_data$`C(mg/L)`), ] # remove NA

### cov ----
# change the names of columns
cov_data <- cov_data %>% rename(ID = patid, 
                                `sex(M/F)` = Sex, 
                                `age(yrs)` = Age..yrs.,
                                `height(cm)` = Height..cm.,
                                `weight(kg)` = Weight..kg.,
                                `BMI(kg/m2)` = BMI..kg.m2.,
                                `CLcr(mL/min)` = CLcr..mL.h.) 
# data classification
cov_data <- cov_data %>%
  mutate(across(-c(`sex(M/F)`, X2D6, X3A4, X2C9, X2C19), as.numeric)) %>% # classify numeric variables - except for sex and enzymes
  mutate(across(c(`sex(M/F)`, X2D6, X3A4, X2C9, X2C19), as.factor)) %>%  # classify categorical variables - sex and enzymes
  rename(
    CYP2D6 = `X2D6`,
    CYP3A4 = `X3A4`,
    CYP2C9 = `X2C9`,
    CYP2C19 = `X2C19`
  )
str(cov_data) # check data classification

## 1.3 merge pk and cov dataframes----
merged_data <- long_pk_data %>%
  left_join(cov_data, by = "ID") # merged by ID with full data of concentration - time 
head(merged_data) # check the merged dataframe

# 2 Variable calculations ---------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: calculation of body size measurement, categorization of body size measurement, 
# PK variable calculation
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
## 2.1 calculation and categorization of BFP----
merged_data <- merged_data %>%
  mutate(`height(cm)` = `height(cm)` / 100) %>%  # Convert height to meters
  mutate(`BFP(%)` = ifelse(`sex(M/F)` == "M",
                           1.20 * (`weight(kg)` / `height(cm)`^2) + 0.23 * `age(yrs)` - 16.2,
                           1.20 * (`weight(kg)` / `height(cm)`^2) + 0.23 * `age(yrs)` - 5.4)) %>% # Create BFP column based on sex
  mutate(BFP_group = ifelse(`BFP(%)` > 25, "Above 25%", "25% or Below")) # categorize into two groups

## 2.2 PK calculation -----
### Cmax, tmax ----
cmax_tmax <- merged_data %>%
  group_by(ID) %>%
  filter(`C(mg/L)` == max(`C(mg/L)`)) %>%
  slice(1) %>%  # In case multiple time points have the same Cmax
  select(ID, `Time(h)`, `C(mg/L)`) %>%
  rename(Tmax = `Time(h)`, Cmax = `C(mg/L)`)
merged_data <- merged_data %>%
  left_join(cmax_tmax, by = "ID") # Join back to merged_data

### Vd -----
dose <- 200
vd <- merged_data %>%
  filter(`Time(h)` >= Tmax) %>% # Keep only time-points after tmax
  group_by(ID) %>%
  summarise(
    ln_C0 = coef(lm(log(`C(mg/L)`) ~ `Time(h)`))["(Intercept)"],  # ln(C0)
    k = -coef(lm(log(`C(mg/L)`) ~ `Time(h)`))["`Time(h)`"],       # k
    C0 = exp(ln_C0),                                              # C0
    Vd = dose / C0                                                # Vd
  )

### CL(AUC) -----
# calculate half life
vd <- vd %>%
  mutate(t_half = log(2) / k) 
# Calculate trapezoidal AUC for each ID
auc_trap <- merged_data %>% 
  group_by(ID) %>% 
  arrange(ID, `Time(h)`) %>% 
  summarise(
    AUC_trap = sum(
      diff(`Time(h)`) * (head(`C(mg/L)`, -1) + tail(`C(mg/L)`, -1)) / 2
    ),
    C_last = last(`C(mg/L)`)
  ) 
# Join with vd to get k for each ID
auc_full <- auc_trap %>%
  left_join(vd %>% select(ID, k), by = "ID") %>%
  mutate(
    Residual = ifelse(!is.na(k), C_last / k, NA_real_),  # residual area Cn/k
    AUC_full = AUC_trap + Residual
  )
# CL
cl_data <- auc_full %>%
  mutate(CL = dose / AUC_full) %>%   # Calculate clearance = Dose / AUC
  select(ID, CL, everything()) 

# 3 Data Exploration --------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: numerical summary of PK variables, graphical assessment of 1) Concentration versus time,
# 2) PK variable correlations, 3) PK variable by phenotype of enzyme, 
# 4) PK variable-body size measurement correlation with linear regression
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
## 3.1 PK numerical summary ---------
pk_results <- cmax_tmax %>% 
  left_join(vd %>% select(ID, Vd), by = "ID") %>% 
  left_join(cl_data %>% select(ID, CL), by = "ID") %>% 
  ungroup() # merge tmax, cmax, cl and vd

str(pk_results)

pk_summary <- pk_results %>%   
  ungroup() %>% 
  summarise(
    Cmax_mean = mean(Cmax, na.rm = TRUE),
    Cmax_median = median(Cmax, na.rm = TRUE),
    Cmax_sd = sd(Cmax, na.rm = TRUE),
    Cmax_range = paste0(min(Cmax, na.rm = TRUE), " - ", max(Cmax, na.rm = TRUE)),
    
    Tmax_mean = mean(Tmax, na.rm = TRUE),
    Tmax_median = median(Tmax, na.rm = TRUE),
    Tmax_sd = sd(Tmax, na.rm = TRUE),
    Tmax_range = paste0(min(Tmax, na.rm = TRUE), " - ", max(Tmax, na.rm = TRUE)),
    
    Vd_mean = mean(Vd, na.rm = TRUE),
    Vd_median = median(Vd, na.rm = TRUE),
    Vd_sd = sd(Vd, na.rm = TRUE),
    Vd_range = paste0(min(Vd, na.rm = TRUE), " - ", max(Vd, na.rm = TRUE)),
    
    CL_mean = mean(CL, na.rm = TRUE),
    CL_median = median(CL, na.rm = TRUE),
    CL_sd = sd(CL, na.rm = TRUE),
    CL_range = paste0(min(CL, na.rm = TRUE), " - ", max(CL, na.rm = TRUE))
  )

## 3.2 Spaghetti plot --------
ggplot(merged_data, aes(x = `Time(h)`, y = `C(mg/L)`, group = ID)) + 
  geom_line(linewidth = 0.3, alpha = 0.3) +  
  labs(
    title = "Spaghetti Plot of BPI889 Concentrations vs. Time",
    x = "Time (hours)",
    y = "BPI889 Concentration (mg/L)"
  ) +
  theme_minimal() + 
  theme(legend.position = "none")  

## 3.3 Correlation plot ------
# between Cmax, tmax, Vd, and CL
corr_data <- pk_results %>% 
  select(Tmax, Cmax, Vd, CL)
ggpairs(corr_data,
        title = "Pairwise Correlation Plot of Cmax, Tmax, Vd, and CL")  # scatter plots

## 3.4  box-and-whisker plots -----
merged_data <- merged_data %>%
  left_join(pk_results %>% select(ID, CL, Vd), by = "ID") %>% 
  group_by(ID) %>% 
  slice(1) 

long_data <- merged_data %>%
  pivot_longer(cols = c(CYP2D6, CYP3A4, CYP2C9, CYP2C19),
               names_to = "Enzyme",
               values_to = "Phenotype")

# Vd plot
ggplot(long_data, aes(x = Phenotype, y = Vd, fill = Phenotype)) +
  geom_boxplot() +
  facet_wrap(~ Enzyme, scales = "free_x") +
  scale_fill_manual(
    name = "Phenotype",
    values = c("0" = "skyblue", "1" = "lightgreen", "2" = "orange"),
    labels = c("0 = Poor metabolizer", "1 = Normal metabolizer", "2 = Extensive metabolizer")
  ) +
  labs(title = "Vd vs Enzyme Phenotypes",
       x = "Phenotype", y = "Volume of distribution") +
  theme_minimal()
# CL plot
ggplot(long_data, aes(x = Phenotype, y = CL, fill = Phenotype)) +
  geom_boxplot() +
  facet_wrap(~ Enzyme, scales = "free_x") +
  scale_fill_manual(
    name = "Phenotype",
    values = c("0" = "skyblue", "1" = "lightgreen", "2" = "orange"),
    labels = c("0 = Poor metabolizer", "1 = Normal metabolizer", "2 = Extensive metabolizer")
  ) +
  labs(title = "CL vs Enzyme Phenotypes",
       x = "Phenotype", y = "Clearance (CL)") +
  theme_minimal()

## 3.5 Scatter plot -------
ggplot(merged_data, aes(x = `BFP(%)`, y = Vd)) +
  geom_point(color = "blue", alpha = 0.6, size = 2) +         
  geom_smooth(method = "lm", color = "red", se = FALSE) +      
  labs(
    title = "Correlation between Body Fat Percentage (BFP) and Vd",
    x = "Body Fat Percentage (BFP)",
    y = "Volume of Distribution (Vd)"
  ) +
  theme_minimal()

# 4 Statistical testing -----------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: ANOVA of PK variables for SNPs, linear regression of PK variable for body size measurement 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
## 4.1 ANOVA on Cmax and CL ------
str(long_data) # check the phenotype is categorical
enzymes <- c("CYP2D6", "CYP3A4", "CYP2C9", "CYP2C19")

run_enzyme_anova <- function(data, enzymes, response_vars = c("Cmax", "CL")) {
  for (enzyme in enzymes) {
    for (response in response_vars) {
      formula <- as.formula(paste(response, "~", enzyme))
      cat("--------------------------------------------------\n")
      cat("ANOVA result for", response, "by", enzyme, "\n")
      print(summary(aov(formula, data = data)))
    }
  }
}
run_enzyme_anova(data = merged_data, enzymes = enzymes)

## 4.2 Linear regression on Vd ~ BFP -----
lm_result <- lm(Vd ~ `BFP(%)`, data = merged_data)
summary_lm <- summary(lm_result)
coef_BFP <- summary_lm$coefficients["`BFP(%)`", "Estimate"]
p_value <- summary_lm$coefficients["`BFP(%)`", "Pr(>|t|)"]
r_squared <- summary_lm$r.squared
# formatting p value
formatted_p <- ifelse(p_value < 0.0001, "< 0.0001", round(p_value, 4))
# report
cat(paste0(
  "A linear regression was conducted to evaluate the association between body fat percentage (BFP) and volume of distribution (Vd). ",
  "The analysis showed that BFP is a statistically significant predictor of Vd (β = ",
  round(coef_BFP, 2), ", p = ", formatted_p, "). ",
  "Specifically, for every 1% increase in BFP, the Vd increases by approximately ",
  round(coef_BFP, 2), " units. ",
  "The model explains about ", round(r_squared * 100, 1), "% of the variance in Vd (R² = ",
  round(r_squared, 3), "), indicating a moderate but meaningful linear relationship."
))

