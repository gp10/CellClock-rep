rm(list = ls()) # clear the console

### REQUIRED PACKAGES:
library("readxl")
library(dplyr)      #data_frame, %>%, filter, summarise, group_by
library(emmeans)    #emmeans, contrast
library(phia)       #testInteractions
library(tidyr)      #spread
library(ARTool)     #art, artlm
library(ggplot2)    #ggplot, stat_..., geom_..., etc
library("rstatix")
library("dplyr")
library(vtree)
library(reshape2)

################################################################################


### LOADING THE DATA (FREQUENCY TABLE):
mydata <- read_excel("Summary_EdU_BrdU_data.xlsx", sheet = "Sheet1")
mydata$FoV = as.factor(mydata$FoV)
mydata$Animal = as.factor(mydata$Animal)
mydata$Genotype = as.factor(mydata$Genotype)
mydata["bEdU"] <- mydata$EdU_blank + mydata$EdU_BrdU + mydata$EdU_pH3G2 + mydata$EdU_pH3M
mydata <- data.frame(mydata)
mydata$tEdU <- mydata$bEdU + mydata$sbEdU
mydata$bBrdU <- mydata$EdU_BrdU + mydata$Blank_BrdU
mydata$tEdU_bBrdU <- mydata$tEdU / mydata$bBrdU
mydata$sbEdU_tEdU <- mydata$sbEdU / mydata$tEdU
mydata$bEdU_bBrdU <- mydata$bEdU / mydata$bBrdU
head(mydata) # dataframe



### PLOTTING FULL DATASET TO SHOW SOURCES OF VARIABILITY:
source("./Plot_Rscripts/Plot_variability_EdU_dotplot.R")
Plot_variability_EdU_dotplot()
source("./Plot_Rscripts/Plot_variability_All_dotplot.R")
Plot_variability_All_dotplot()



### RECONSTRUCTING INDIVIDUAL CELL DATA (RAW COUNT TABLE):
# rbind(melt(mydata[,c(2,3,1,4:7)]), melt(mydata[,c(2,3,1,8:11)]), melt(mydata[,c(2,3,1,13)]))
mydataval <- cbind(mydata[rep(seq_len(nrow(mydata)), times = 9), c(3,2,1)],
                   Loc = c(rep("Basal",nrow(mydata)*8), rep("SB",nrow(mydata))),
                   Lineage = c(rep("EdU",nrow(mydata)*4), rep("Blank",nrow(mydata)*4), rep("EdU",nrow(mydata))),
                   Phase = c(rep("Blank",nrow(mydata)), rep("BrdU",nrow(mydata)), rep("pH3G2",nrow(mydata)), rep("pH3M",nrow(mydata)),
                             rep("Blank",nrow(mydata)), rep("BrdU",nrow(mydata)), rep("pH3G2",nrow(mydata)), rep("pH3M",nrow(mydata)),
                             rep("Blank",nrow(mydata))),
                   Freq = c(mydata$EdU_blank, mydata$EdU_BrdU, mydata$EdU_pH3G2, mydata$EdU_pH3M,
                            mydata$Blank_blank, mydata$Blank_BrdU, mydata$Blank_pH3G2, mydata$Blank_pH3M,
                            mydata$sbEdU)
                   )
mydataval <- mydataval[rep(seq_len(nrow(mydataval)), times = mydataval$Freq), -ncol(mydataval)]
mydataval$Loc <- as.factor(mydataval$Loc)
mydataval$Lineage <- as.factor(mydataval$Lineage)
mydataval$Phase <- as.factor(mydataval$Phase)
mydataval$Genotype <- relevel(mydataval$Genotype, "WT")

table(mydataval)
nrow(mydataval)
head(mydataval)

# Overview of all covariates:
vtree(mydataval, c("Genotype","Animal","FoV","Loc","Lineage"), sortfill = TRUE, horiz = FALSE)




### SUBSET INDIVIDUAL BASAL CELLS FOR CELL-CYCLE ANALYSIS:
mydataval_b <- subset(mydataval, Loc == "Basal")
nrow(mydataval_b)
sum(mydata$sbEdU) == nrow(mydataval) - nrow(mydataval_b) # check

# Overview of covariates within basal data:
vtree(mydataval_b, c("Genotype","Animal","FoV","Lineage"), sortfill = TRUE, horiz = FALSE)

# SUBSET INDIVIDUAL {BASAL, EDU+} or {BASAL, Blank} CELL DATA:
mydataval_bEdU <- subset(mydataval_b, Lineage == "EdU", select = names(mydataval))
mydataval_bBlank <- subset(mydataval_b, Lineage == "Blank", select = names(mydataval))
nrow(mydataval_bEdU) + nrow(mydataval_bBlank) == nrow(mydataval_b) # check

# Overview of covariates within {basal, EdU+} or {basal, Blank} data:
vtree(mydataval_bEdU, c("Genotype","Animal","FoV"), sortfill = TRUE, horiz = FALSE)





### MULTINOMIAL LOGISTIC REGRESSION MODELS (GLM): FACTORS (PREDICTORS) INFLUENCING CELL-CYCLE PHASE DISTRIBUTION:
# GLM on basal cells:
glm_b <- glm(Phase ~ Genotype + Animal + Lineage + FoV, family=binomial(link='logit'), data=mydataval_b)
summary(glm_b)
#contrasts(mydataval_b$Animal)
alias(glm_b)$Complete # effects that cannot be estimated independently of terms (levels) which occur earlier in the model... they sum 7 co-linear levels
# (note that linkage of Animal and FoV with Genotype leaves 1 & 6 dependent levels, respectively)
# check degrees of freedom of test: independent levels or rank - 1:
glm_b_model <- mydataval_b %>% dplyr::count(Genotype, Animal, FoV, Lineage, Phase)
glm_b_model <- model.matrix(~ Genotype + Animal + FoV + Lineage, glm_b_model[,-5])
library(Matrix)
rankMatrix(glm_b_model[c(1:nrow(glm_b_model)),]) - 1 # Df = (2-1) + (7-1)-1 + (21-1)-6 + (2-1) = 21 levels
# difference between the null deviance and the model's deviance is distributed as a chi-squared
# under the hypothesis that the model as a whole is no better than the null model, 
# with degrees of freedom equal to the null df minus the model's df:
# p-val = 1 - pchisq( Null_Deviance - Residual_Deviance, df= (Null_df - Residual_df)).
anova(glm_b, test="Chisq") # most variability explained by:  Lineage > Animal > Genotype > FoV

# GLM on {basal, EdU+} cells:
glm_bEdU <- glm(Phase ~ Genotype + Animal + FoV, family=binomial(link='logit'), data=mydataval_bEdU)
#glm_bEdU <- glm(as.factor(mydataval_bEdU$Phase=="Blank") ~ Genotype + Animal + FoV, family=binomial(link='logit'), data=mydataval_bEdU)
summary(glm_bEdU)
anova(glm_bEdU, test="Chisq") # most variability explained by: Animal > Genotype > FoV

# GLM on {basal, Blank} cells:
glm_bBlank <- glm(Phase ~ Genotype + Animal + FoV, family=binomial(link='logit'), data=mydataval_bBlank)
summary(glm_bBlank)
anova(glm_bBlank, test="Chisq") # most variability explained by: Animal > Genotype > FoV





### SUBSET AND OBTAIN {BASAL, EDU+} & {BASAL, ALL} FREQUENCY TABLE:
freq_bAll <- mydataval_b %>% dplyr::count(Genotype, Animal, FoV, Phase)
freq_bEdU <- mydataval_bEdU %>% dplyr::count(Genotype, Animal, FoV, Phase)
nrow(freq_bAll) == 3*4*4 + 3*3*4  #  Check it corresponds to:  FoV*AnimalWT*Phase + FoV*AnimalMut*Phase
nrow(freq_bEdU) == 3*4*4 + 3*3*4  #  Check it corresponds to:  FoV*AnimalWT*Phase + FoV*AnimalMut*Phase

# Calculate sum {Basal, EdU+} and {Basal, All} cells per FoV:
for (i in 1:nrow(freq_bAll))
{
  freq_bAll$bCells[i] = sum(freq_bAll$n[ freq_bAll$Genotype==freq_bAll$Genotype[i] & freq_bAll$Animal==freq_bAll$Animal[i] & freq_bAll$FoV==freq_bAll$FoV[i] ])
  freq_bEdU$bEdU[i] = sum(freq_bEdU$n[ freq_bEdU$Genotype==freq_bEdU$Genotype[i] & freq_bEdU$Animal==freq_bEdU$Animal[i] & freq_bEdU$FoV==freq_bEdU$FoV[i] ])
}
freq_bAll$Fraction <- freq_bAll$n / freq_bAll$bCells * 100
freq_bEdU$Fraction <- freq_bEdU$n / freq_bEdU$bEdU * 100





### AVERAGING OVER FoV:
freq_bAll_avg <- aggregate(Fraction ~ Genotype + Animal + Phase, data=freq_bAll, FUN=mean)
freq_bAll_avg <- merge(freq_bAll_avg, aggregate(n ~ Genotype + Animal + Phase, data=freq_bAll, FUN=sum))
freq_bAll_avg <- merge(freq_bAll_avg, aggregate(bCells ~ Genotype + Animal + Phase, data=freq_bAll, FUN=sum))
freq_bAll_avg$FracEst <- freq_bAll_avg$n / freq_bAll_avg$bCells * 100
(freq_bAll_avg$FracEst - freq_bAll_avg$Fraction) / freq_bAll_avg$FracEst # max deviation with both approximations
freq_bAll_avg$n_FoV <- round(freq_bAll_avg$n/3)
freq_bEdU_avg <- aggregate(Fraction ~ Genotype + Animal + Phase, data=freq_bEdU, FUN=mean)
freq_bEdU_avg <- merge(freq_bEdU_avg, aggregate(n ~ Genotype + Animal + Phase, data=freq_bEdU, FUN=sum))
freq_bEdU_avg <- merge(freq_bEdU_avg, aggregate(bEdU ~ Genotype + Animal + Phase, data=freq_bEdU, FUN=sum))
freq_bEdU_avg$FracEst <- freq_bEdU_avg$n / freq_bEdU_avg$bEdU * 100
# Calculate max deviation with both approximations:
(freq_bEdU_avg$FracEst - freq_bEdU_avg$Fraction) / freq_bEdU_avg$FracEst
freq_bEdU_avg$n_FoV <- round(freq_bEdU_avg$n/3)





### PLOT Basal EdU or total cell fractions in different cell-cycle phases:
source("./Plot_Rscripts/Plot_All_cellcycle_circular_barplot.R")
Plot_All_cellcycle_circular_barplot()
source("./Plot_Rscripts/Plot_EdU_cellcycle_circular_barplot.R")
Plot_EdU_cellcycle_circular_barplot()





### AVERAGING OVER MICE:
freq_bAll_avgmice <- aggregate(Fraction ~ Genotype + Phase, data=freq_bAll_avg, FUN=mean)
freq_bAll_avgmice <- merge(freq_bAll_avgmice, aggregate(n_FoV ~ Genotype + Phase, data=freq_bAll_avg, FUN=mean))
freq_bAll_avgmice$n_FoV <- round(freq_bAll_avgmice$n_FoV)
toappend <- aggregate(Fraction ~ Genotype + Phase, data=freq_bAll_avg, FUN=sd)
colnames(toappend)[ncol(toappend)] <- "Frac_ci"
freq_bAll_avgmice <- merge(freq_bAll_avgmice, toappend)
freq_bAll_avgmice$up <- freq_bAll_avgmice$Fraction + freq_bAll_avgmice$Frac_ci
freq_bAll_avgmice$dn <- freq_bAll_avgmice$Fraction - freq_bAll_avgmice$Frac_ci

freq_bEdU_avgmice <- aggregate(Fraction ~ Genotype + Phase, data=freq_bEdU_avg, FUN=mean)
freq_bEdU_avgmice <- merge(freq_bEdU_avgmice, aggregate(n_FoV ~ Genotype + Phase, data=freq_bEdU_avg, FUN=mean))
freq_bEdU_avgmice$n_FoV <- round(freq_bEdU_avgmice$n_FoV)
toappend <- aggregate(Fraction ~ Genotype + Phase, data=freq_bEdU_avg, FUN=sd)
colnames(toappend)[ncol(toappend)] <- "Frac_ci"
freq_bEdU_avgmice <- merge(freq_bEdU_avgmice, toappend)
freq_bEdU_avgmice$up <- freq_bEdU_avgmice$Fraction + freq_bEdU_avgmice$Frac_ci
freq_bEdU_avgmice$dn <- freq_bEdU_avgmice$Fraction - freq_bEdU_avgmice$Frac_ci

# Merge EdU and All elements in a single dataframe:
freq_bAllandbEdU_avgmice <- rbind(freq_bAll_avgmice, freq_bEdU_avgmice)
freq_bAllandbEdU_avgmice$Lineage <- as.factor(c(rep("All", 8), rep("EdU", 8)))
# write.csv(freq_bAllandbEdU_avgmice, file="freq_CellCyclePhases_R.xlsx")

# Compare EdU vs Blank in a single dataframe:
freq_bEdUvsBlank_avgmice <- freq_bEdU_avgmice[,c(1:2,4)]
names(freq_bEdUvsBlank_avgmice)[3] <- "n_FoV_EdU"
freq_bEdUvsBlank_avgmice$n_FoV_Blank <- freq_bAll_avgmice$n_FoV - freq_bEdUvsBlank_avgmice$n_FoV_EdU

freq_bEdUvsBlank_avg <- freq_bEdU_avg[,c(1:3,8)]
names(freq_bEdUvsBlank_avg)[4] <- "n_FoV_EdU"
freq_bEdUvsBlank_avg$n_FoV_Blank <- freq_bAll_avg$n_FoV - freq_bEdUvsBlank_avg$n_FoV_EdU



### PLOT Cell-cycle Errorbars:
source("./Plot_Rscripts/Plot_cellcycle_errorbars.R")
Plot_cellcycle_errorbars()


### PLOT Donut Charts:
source("./Plot_Rscripts/Plot_cellcycle_Donut.R")
Plot_cellcycle_Donut()


### PLOT Piechart:
source("./Plot_Rscripts/Plot_cellcycle_Piechart.R")
Plot_cellcycle_Piechart()





### CONTINGENCY TESTS: GENOTYPE ASSOCIATION WITH CELL-CYCLE PHASE DISTRIBUTION IN {Basal, EdU} POPULATION:
# 2x4  Contingency table: Genotype x Phase in {Basal, EdU}:
contin_EdU <- melt(freq_bEdU_avgmice, c("Genotype","Phase"), "n_FoV")
contin_EdU <- dcast(contin_EdU, Genotype ~ Phase)
rownames(contin_EdU) <- contin_EdU$Genotype
contin_EdU <- contin_EdU[,c(2:ncol(contin_EdU))]
contin_EdU
mosaicplot(contin_EdU, color=TRUE)
chisq.test(contin_EdU)$expected

# Chi2 test on 2x4 contingency table:
chisq.test(contin_EdU)

# Fisher's exact test on 2x4 contingency table:
test <- fisher.test(contin_EdU)
test$p.value # no sig. association between genotype and cell cycle distribution

# 2x2 Binomial contingency tables:
contin_EdU_Blank <- cbind(Blank = contin_EdU[,1], Others = rowSums(contin_EdU[,c(2:4)]))
contin_EdU_BrdU <- cbind(BrdU = contin_EdU[,2], Others = rowSums(contin_EdU[,c(1,3,4)]))
contin_EdU_pH3G2 <- cbind(pH3G2 = contin_EdU[,3], Others = rowSums(contin_EdU[,c(1,2,4)]))
contin_EdU_pH3M <- cbind(pH3M = contin_EdU[,4], Others = rowSums(contin_EdU[,c(1:3)]))

# Post-hoc tests (pairwise Fisher's exact tests with FDR correction):
library(rstatix)
pairwise_fisher_test(as.matrix(contin_EdU), p.adjust.method = "fdr")

# Fisher's exact tests on phase-specific comparisons:
pval_EdU_Blank <- fisher.test(contin_EdU_Blank)$p.value
pval_EdU_BrdU <- fisher.test(contin_EdU_BrdU)$p.value
pval_EdU_pH3G2 <- fisher.test(contin_EdU_pH3G2)$p.value
pval_EdU_pH3M <- fisher.test(contin_EdU_pH3M)$p.value
padj_groups <- p.adjust(c(pval_EdU_Blank,pval_EdU_BrdU,pval_EdU_pH3G2,pval_EdU_pH3M), method = "BH")





### CONTINGENCY TESTS: LINEAGE ASSOCIATION WITH CELL-CYCLE PHASE DISTRIBUTION IN {WT, Basal} POPULATION:
# 2x4  Contingency table: Lineage x Phase in {WT, Basal}:
freq_bEdUandBlank_avgmice <- cbind(freq_bEdUvsBlank_avgmice[,c(1:2)],
                                    n_FoV = freq_bEdUvsBlank_avgmice$n_FoV_EdU,
                                    Lineage = rep("EdU", nrow(freq_bEdUvsBlank_avgmice)))
freq_bEdUandBlank_avgmice <- rbind(freq_bEdUandBlank_avgmice, cbind(freq_bEdUvsBlank_avgmice[,c(1:2)],
                                                                      n_FoV = freq_bEdUvsBlank_avgmice$n_FoV_Blank,
                                                                      Lineage = rep("Blank", nrow(freq_bEdUandBlank_avgmice))))
contin_WT <- melt(freq_bEdUandBlank_avgmice[freq_bEdUandBlank_avgmice$Genotype=="WT",], c("Lineage", "Phase"), "n_FoV")
contin_WT <- dcast(contin_WT, Lineage ~ Phase)
rownames(contin_WT) <- contin_WT$Lineage
contin_WT <- contin_WT[,c(2:ncol(contin_WT))]
contin_WT
mosaicplot(contin_WT, color=TRUE)
chisq.test(contin_WT)$expected

# Chi2 test on 2x4 contingency table:
chisq.test(contin_WT)

# Fisher's exact test on 2x4 contingency table:
test <- fisher.test(contin_WT)
test$p.value # there IS strong association between lineage and cell cycle distribution

# 2x2 Binomial contingency tables:
contin_WT_Blank <- cbind(Blank = contin_WT[,1], Others = rowSums(contin_WT[,c(2:4)]))
contin_WT_BrdU <- cbind(BrdU = contin_WT[,2], Others = rowSums(contin_WT[,c(1,3,4)]))
contin_WT_pH3G2 <- cbind(pH3G2 = contin_WT[,3], Others = rowSums(contin_WT[,c(1,2,4)]))
contin_WT_pH3M <- cbind(pH3M = contin_WT[,4], Others = rowSums(contin_WT[,c(1:3)]))

# Fisher's exact tests on phase-specific comparisons:
pval_WT_Blank <- fisher.test(contin_WT_Blank)$p.value
pval_WT_BrdU <- fisher.test(contin_WT_BrdU)$p.value
pval_WT_pH3G2 <- fisher.test(contin_WT_pH3G2)$p.value
pval_WT_pH3M <- fisher.test(contin_WT_pH3M)$p.value
padj_groups <- p.adjust(c(pval_WT_Blank,pval_WT_BrdU,pval_WT_pH3G2,pval_WT_pH3M), method = "BH")

# odds ratios from sig. tests:
1 / fisher.test(contin_WT_Blank)$estimate
1 / fisher.test(contin_WT_BrdU)$estimate

# 2x4x4  Contingency table: Lineage x Phase x Animal in {WT, Basal}:
names(freq_bEdUvsBlank_avg)[4:5] <- c("EdU","Blank")
freq_bEdUandBlank_avg <- melt(freq_bEdUvsBlank_avg, c(1:3), c(4:5))
names(freq_bEdUandBlank_avg)[4:5] <- c("Lineage","n_FoV")
contin_WT_mice <- melt(freq_bEdUandBlank_avg[freq_bEdUandBlank_avg$Genotype=="WT",], c("Lineage","Phase","Animal"), "n_FoV")
contin_WT_mice <- xtabs(value ~ Lineage + Phase + Animal, data = contin_WT_mice)[,,1:4]
ftable(contin_WT_mice)
library(vcd)
oddsratio(contin_WT_mice, log =FALSE)

# Check homogeneous association across animals (odds ratios across animal strata being identical?):
# Woolf test:
library(vcd)
woolf_test(contin_WT_mice) # homogeneity checked -> Cochran-Mantel-Haenszel test is appropriate

# Cochran-Mantel-Haenszel (CMH) test (given homogeneous association through different strata):
mantelhaen.test(contin_WT_mice) # cell cycle distribution depends on lineage across mice

# CMH post-hoc tests on individual mice: (Fisher)
library(rcompanion)
groupwiseCMH(contin_WT_mice)





### CONTINGENCY TESTS: LINEAGE ASSOCIATION WITH CELL-CYCLE PHASE DISTRIBUTION IN {Mut, Basal} POPULATION:
# 2x4  Contingency table: Lineage x Phase in {Mut, Basal}:
contin_Mut <- melt(freq_bEdUandBlank_avgmice[freq_bEdUandBlank_avgmice$Genotype=="Mut",], c("Lineage", "Phase"), "n_FoV")
contin_Mut <- dcast(contin_Mut, Lineage ~ Phase)
rownames(contin_Mut) <- contin_Mut$Lineage
contin_Mut <- contin_Mut[,c(2:ncol(contin_Mut))]
contin_Mut
mosaicplot(contin_Mut, color=TRUE)
chisq.test(contin_Mut)$expected

# Chi2 test on 2x4 contingency table:
chisq.test(contin_Mut)

# Fisher's exact test on 2x4 contingency table:
test <- fisher.test(contin_Mut)
test$p.value # there IS strong association between lineage and cell cycle distribution

# 2x2 Binomial contingency tables:
contin_Mut_Blank <- cbind(Blank = contin_Mut[,1], Others = rowSums(contin_Mut[,c(2:4)]))
contin_Mut_BrdU <- cbind(BrdU = contin_Mut[,2], Others = rowSums(contin_Mut[,c(1,3,4)]))
contin_Mut_pH3G2 <- cbind(pH3G2 = contin_Mut[,3], Others = rowSums(contin_Mut[,c(1,2,4)]))
contin_Mut_pH3M <- cbind(pH3M = contin_Mut[,4], Others = rowSums(contin_Mut[,c(1:3)]))

# Fisher's exact tests on phase-specific comparisons:
pval_Mut_Blank <- fisher.test(contin_Mut_Blank)$p.value
pval_Mut_BrdU <- fisher.test(contin_Mut_BrdU)$p.value
pval_Mut_pH3G2 <- fisher.test(contin_Mut_pH3G2)$p.value
pval_Mut_pH3M <- fisher.test(contin_Mut_pH3M)$p.value
padj_groups <- p.adjust(c(pval_Mut_Blank,pval_Mut_BrdU,pval_Mut_pH3G2,pval_Mut_pH3M), method = "BH")

# odds ratios from sig. tests:
1 / fisher.test(contin_Mut_Blank)$estimate
1 / fisher.test(contin_Mut_BrdU)$estimate

# 2x4x4  Contingency table: Lineage x Phase x Animal in {WT, Basal}:
contin_Mut_mice <- melt(freq_bEdUandBlank_avg[freq_bEdUandBlank_avg$Genotype=="Mut",], c("Lineage","Phase","Animal"), "n_FoV")
contin_Mut_mice <- xtabs(value ~ Lineage + Phase + Animal, data = contin_Mut_mice)[,,5:7]
ftable(contin_Mut_mice)
oddsratio(contin_Mut_mice, log =FALSE)

# Check homogeneous association across animals (odds ratios across animal strata being identical?):
# Woolf test:
library(vcd)
woolf_test(contin_Mut_mice) # homogeneity checked -> Cochran-Mantel-Haenszel test is appropriate

# Cochran-Mantel-Haenszel (CMH) test (given homogeneous association through different strata):
mantelhaen.test(contin_Mut_mice) # cell cycle distribution depends on lineage across mice

# CMH post-hoc tests on individual mice: (Fisher)
library(rcompanion)
groupwiseCMH(contin_Mut_mice)






### AVERAGING OVER FoV:
mydata_avg <- aggregate(tEdU_bBrdU ~ Genotype + Animal, data=mydata, FUN=mean) # tEdU_bBrdU
mydata_avg <- merge(mydata_avg, aggregate(tEdU ~ Genotype + Animal, data=mydata, FUN=sum))
mydata_avg <- merge(mydata_avg, aggregate(bBrdU ~ Genotype + Animal, data=mydata, FUN=sum))
mydata_avg <- merge(mydata_avg, aggregate(sbEdU_tEdU ~ Genotype + Animal, data=mydata, FUN=mean)) # sbEdU_tEdU
mydata_avg <- merge(mydata_avg, aggregate(bEdU_bBrdU ~ Genotype + Animal, data=mydata, FUN=mean)) # bEdU_bBrdU
head(mydata_avg) # dataframe




### ANALYZE STRATIFICATION RATE:  {SB EdU / TOTAL EdU}

# equality of variances (F-test):
var.test(sbEdU_tEdU ~ Genotype, data = mydata_avg) # checked (n.s.)
#var.test(sbEdU_tEdU ~ Genotype, data = mydata) # checked

# Compare genotypes: two-sample t-test (two-tailed)
ttest_strat_gt <- t.test(sbEdU_tEdU ~ Genotype, data = mydata_avg, var.equal = TRUE)
#wilcox.test(sbEdU_tEdU ~ Genotype, data = mydata_avg)
cat("Stratification rate Mut vs WT (t-test): p-val =", ttest_strat_gt$p.value, "\n")

# Plot:
source("./Plot_Rscripts/Plot_sbEdU_tEdU_dotplot.R")
Plot_sbEdU_tEdU_dotplot()




### ANALYZE RENEWAL RATE:  {TOTAL EdU / BASAL BrdU}

# equality of variances (F-test):
var.test(tEdU_bBrdU ~ Genotype, data = mydata_avg) # checked (n.s.)

# Compare genotypes: two-sample t-test (two-tailed)
ttest_renew_gt <- t.test(tEdU_bBrdU ~ Genotype, data = mydata_avg, var.equal = TRUE)
#wilcox.test(tEdU_bBrdU ~ Genotype, data = mydata_avg)
cat("Renewal rate Mut vs WT (t-test): p-val =", ttest_renew_gt$p.value, "\n")

# Compare to value=1: one-sample t-test (one-tailed)
ttest_renew_WTvs1 <- t.test(mydata_avg[mydata_avg$Genotype=="WT",]$tEdU_bBrdU,
                            mu = 1, alternative = "greater")
ttest_renew_Mutvs1 <- t.test(mydata_avg[mydata_avg$Genotype=="Mut",]$tEdU_bBrdU,
                             mu = 1, alternative = "greater")
cat("Renewal rate WT > 1 (t-test): p-val =", ttest_renew_WTvs1$p.value, "\n")
cat("Renewal rate Mut > 1 (t-test): p-val =", ttest_renew_Mutvs1$p.value, "\n")

# Compare to BrdU ratios across rnd pairs of samples (permutation):
library(gtools)
perm <- permutations(length(mydata$bBrdU),2, mydata$bBrdU, repeats.allowed = T)
perm_ratio <- perm[,1] / perm[,2]
# two-sample t-test (two-tailed) against permuted sample of bBrdU values:
ttest_renew_WTvsPerm <- t.test(mydata_avg[mydata_avg$Genotype=="WT",]$tEdU_bBrdU,
                               perm_ratio, alternative = "greater", var.equal = TRUE)
ttest_renew_MutvsPerm <- t.test(mydata_avg[mydata_avg$Genotype=="Mut",]$tEdU_bBrdU,
                                perm_ratio, alternative = "greater", var.equal = TRUE)

# Plot:
source("./Plot_Rscripts/Plot_tEdU_bBrdU_dotplot.R")
Plot_tEdU_bBrdU_dotplot()

# Plot BASAL EdU vs. BASAL BrdU:
source("./Plot_Rscripts/Plot_bEdU_bBrdU_dotplot.R")
Plot_bEdU_bBrdU_dotplot()


sessionInfo()

