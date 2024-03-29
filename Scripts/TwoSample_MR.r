###Two-Sample MR of platelet phenotypes and platelet proteins on risk of arterial and venous thrombotic events###
##Platelet phenotypes: Platelet count and mean platelet volume##
##Platelet proteins: P-selectin and P-selectin glycoprotein ligand 1##
##Disease outcomes: Acute MI, Mortality after MI, DVT, PE 

#Install relevant packages 
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages("ggplot2")

#Load packages
library(TwoSampleMR)
library(MRPRESSO)

#Extract platelet exposures 
ao<-available_outcomes()
id.platelet.exp <- c("ukb-d-30080_irnt", "ukb-d-30100_irnt", "prot-a-2667", "prot-a-2668")
exposure <- extract_instruments(id.platelet.exp)

#OR 

#Read in exposure data (e.g. RALA and RALB)
data <- read.table("RALB_eqtl_instruments_FINAL", header = T)
exposure <- format_data(data, type ="exposure", header = TRUE, snp_col = "SNP", gene_col = "Gene", chr_col = "Chr", pos_col = "Pos", effect_allele_col = "EA", other_allele_col = "OA", pval_col = "Pvalue", samplesize_col = "N", eaf_col = "EAF", beta_col = "Beta", se_col = "SE")

#OR 

#Read in platelet GWAS data 
PLT_data <- read.table("ieu-a-1008_PLT_SD", header = T)
exposure <- format_data(PLT_data, type ="exposure", header = TRUE, snp_col = "SNP", chr_col = "chr", pos_col = "pos", beta_col = "beta_SD", se_col = "se_SD", pval_col = "pval", samplesize_col = "samplesize", effect_allele_col = "effect_allele", other_allele_col = "other_allele")

MPV_data <- read.table("ieu-a-1006_MPV_SD", header = T)
exposure <- format_data(MPV_data, type ="exposure", header = TRUE, snp_col = "SNP", chr_col = "chr", pos_col = "pos", beta_col = "beta_SD", se_col = "se_SD", pval_col = "pval", samplesize_col = "samplesize", effect_allele_col = "effect_allele", other_allele_col = "other_allele")

#LD clumping 
try(exposure <- clump_data(exposure)) 

#Extract CVD outcomes 
ao<-available_outcomes()
id.cvd.out <- c("finn-b-I9_MI_EXNONE", "finn-b-I9_PHLETHROMBDVTLOW_EXNONE", "finn-b-I9_PULMEMB_EXNONE") #FinnGen MI, DVT and PE
outcome_dat <- extract_outcome_data(snps = exposure$SNP, outcomes = id.cvd.out)

#OR 

#Read in unadjusted and adjusted outcome data

#Unadjusted 
outcome_dat <- read_outcome_data(snps = exposure$SNP, filename = "Mortality_EA_jointgender_25Feb2021_1_processed.csv", sep = " ", snp_col = "SNP", effect_allele_col = "Allele1", other_allele_col = "Allele2", eaf_col = "Freq1", se_col = "StdErr", pval_col = "PVALUE", beta_col = "Effect")

#Dudbridge et al approach adjusted 
outcome_dat <- read_outcome_data(snps = exposure$SNP, filename = "Mortality_EA_jointgender_25Feb2021_Dudbridgecorrected.txt", sep = " ", snp_col = "SNP", beta_col = "beta.corrected", se_col = "se.corrected", pval_col = "pval.corrected", effect_allele_col = "EA.prognosis", other_allele_col = "OA.prognosis", eaf_col = "EAF.prognosis")

#SH 0.08 adjustment factor 
outcome_dat <- read_outcome_data(snps = exposure$SNP, filename = "Mortality_EA_jointgender_25Feb2021_SHcorrected_0.08adjustmentfac_xp_thresh_0.001.txt", snp_col = "SNP", beta_col = "beta.corrected", se_col = "se.corrected", pval_col = "pval.corrected", effect_allele_col = "EA.prognosis", other_allele_col = "OA.prognosis", eaf_col = "EAF.prognosis")

#SH 0.59 adjustment factor
outcome_dat <- read_outcome_data(snps = exposure$SNP, filename = "Mortality_EA_jointgender_25Feb2021_SHcorrected_0.59adjustmentfac_xp_thresh_0.1.txt", snp_col = "SNP", beta_col = "beta.corrected", se_col = "se.corrected", pval_col = "pval.corrected", effect_allele_col = "EA.prognosis", other_allele_col = "OA.prognosis", eaf_col = "EAF.prognosis")

#Harmonise exposure and outcome data 
dat <- harmonise_data(exposure_dat = exposure, outcome_dat = outcome_dat)

#MR analysis 
Res <- mr(dat) #IVW (default), Egger, WM, MODE, Wald ratio 
Res_pleio <- mr_pleiotropy_test(dat) # MR-Egger intercept test 
Res_hetero <- mr_heterogeneity(dat) #
Res_single <- mr_singlesnp(dat) #single SNP analysis
or_results_res <- generate_odds_ratios(Res) #Odds ratio 
Res_rucker_jackknife <- mr_rucker_jackknife(dat, parameters = default_parameters())
Res_MR_PRESSO <- run_mr_presso(dat) #Run if there is evidence of heterogeneity in Res_hetero findings 
Res_mr_leaveoneout <- mr_leaveoneout(dat, parameters = default_parameters(), method = mr_ivw)

#MR-PRESSO
Res_MR_PRESSO <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)


#Plot MR-leave one out 
png("MR_leaveoneout_exposure_outcome.png", width = 4, height = 4, units = 'in', res = 300)
mr_leaveoneout_plot(Res_mr_leaveoneout)
dev.off()

#Write results 
write.table(Res,"MR_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_pleio,"MR_pleio_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_hetero,"MR_hetero_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_single,"MR_single_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(or_results_res,"MR_OR_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_rucker_jackknife, "MR_rucker_jackknife_exposure_outcome", sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_MR_PRESSO, "MR_PRESSO_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F) 
write.table(Res_mr_leaveoneout, "Res_mr_leaveoneout_exposure_outcome", sep="\t",col.names=T,row.names=F,quote=F) 

##Check instrument strength 

#Rename required columns
exposure$BetaXG<- exposure$beta.exposure
exposure$seBetaXG<- exposure$se.exposure
BetaXG   = exposure$BetaXG
seBetaXG = exposure$seBetaXG 
seBetaYG<-dat$se.outcome

BXG             = abs(BetaXG)         # gene-exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted = Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

mF  # Want mF to be high for good IVW performance
unIsq # Want Isq to be close to 1 for good MR-Egger performance. This can then be used to indicate whether SIMEX correction is needed for MR-Egger. If I2 is <90%, then this is typically recommended. The code for running SIMEX is below: 

##SIMEX correction 
#Run the simex correction for each phenotype
install.packages("simex")
#load package
library(simex)

#create empty dataframe to store output
simexegger<-c()

#run simex
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)        

# MR-Egger regression (weighted)
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE)

# Simulation extrapolation
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE")
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE")
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

write.table(mod1, "SIMEX_correction_weighted_exposure_outcome", sep ="\t",col.names=T,row.names=F,quote=F) 
write.table(mod2, "SIMEX_correction_unweighted_exposure_outcome", sep ="\t",col.names=T,row.names=F,quote=F) 
