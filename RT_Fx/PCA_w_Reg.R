setwd("~/Documents/Moffitt/Project/MATLAB/Scripts/RT_Fx")


# Loading in Data ---------------------------------------------------------
# Datatable
ti_all <- readxl::read_excel("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_Anonymized.xlsx")
ti_outcome <- readxl::read_excel("/Users/saketpandit/Documents/Moffitt/Project/Data/TumorAndImmuneCounts_LungOutcomeSubset_Anonymized.xlsx")
ti.sub <- dplyr::left_join(ti_outcome[1:(nrow(ti_outcome)-9),], ti_all)

# Variables of Interest
lrf <- ti.sub$`Local Failure` # 0 = Tumor resolution, 1 = Tumor escape
SF2 <- ti.sub$`SF (2)`
rsi <- ti.sub$RSI
t <- ti.sub$`Cancer Cells`
at <- ti.sub$`Total Anti-Tumor Case 3`
rc <- ti.sub$`Total Pro-Tumor Case 3`

# PCA Plot ----------------------------------------------------------------
# Data Cleaning
t <- t/max(t)
at <- at/max(at)
rc <- rc/max(rc)
dat <- data.frame(t, at, rc)
dat <- dat * 10^6

# Calculating Principle Components
pcr1 <- prcomp(dat)
biplot(pcr1)
summary(pcr1)
pcr1
