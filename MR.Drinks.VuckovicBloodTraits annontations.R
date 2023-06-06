#Liu et al DrinksPerWeek on VUCKOVIC 2020 CELL blood traits

###########################################
##  Below is the iterative run for all blood traits, requires pre-filtered and clumped data for exposure
###############################


###########   Filter


library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
#library(MVMR)
library(TwoSampleMR)
#library(MRPRESSO)

setwd("/home/thomc/thomlab_MR/thomc/Behaviors_Blood/MR.Drinks.VuckovicBloodTraits")     ##DONT FORGET change

### Load data and filter. Have arbitrarily chosen MPV as the blood trait of choice here to initially filter (other blood traits done in 2nd part)

##########assuming this works, but an if/else loop at top so it'd happen automatically if no filtered data sets exist

if( file.exists("Drinks.filtered.txt") & file.exists("MPV.filtered.txt")) {                                    
  print("Working from pre-tabulated files")
  e1_comm <- fread("Drinks.filtered.txt", header=T) ## exposure 1
  #out_comm <- fread("MPV.filtered.txt", header=T)  ## Outcome common
  
## MPV means that it is mean platelet volume
  
} else {
  print("Loading summ stats and making filtered files")
  
  e1<- fread("/mnt/isilon/thom_lab/thomlab_datasets/GWAS/52_Liu_Smoking/DrinksPerWeek.txt.gz", header=T)  ## looking at the file in DrinksPerWeek, adding a column with Drinks for the Phenotype
  e1$Phenotype<- "Drinks"
  e1$SNP<- e1$RSID                 ##the column from e1 that contains the column RSID is assigned to the SNP column
  e1$ChrPosRefAlt<- paste0(e1$CHROM, ":", e1$POS, "_", e1$REF, "_", e1$ALT)   ##creating two separate columns where they have the chromosome:Position_ ref allele_ alternate allele (then switched) so that each time they get 50% of the combinations
  e1$ChrPosAltRef<- paste0(e1$CHROM, ":", e1$POS, "_", e1$ALT, "_", e1$REF)   ##paste0 means no spaces
  
  ## $ means it is accessing a column

  #e2<- fread("/project/voight_datasets/GWAS/04_giant/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T)        ## e1 was drinks, e2 is BMI, and e3 is MPV
  #e2$Phenotype<- "BMI"

  e3<- fread("/mnt/isilon/thom_lab/thomlab_datasets/GWAS/59_2020BloodTraits/BCX2_MPV_EA_GWAMA.out.gz", header=T)  #loading data from the MPV GWAS
  e3$Phenotype<- "MPV"
  
  #e1 <- fread("/project/voight_datasets/GWAS/61_dbGAP_MVP/79613/MVP.Lipids/MVP.EUR.TC.gwas.dbGAP.txt.gz", header=T)   ##loading data from the TC GWAS
  #e1$Phenotype<- "TC"
  #e1$ChrPosRefAlt<- paste0(e1$Chromosome, ":", e1$Position, "_", e1$Allele1, "_", e1$Allele2)
  #e1$ChrPosAltRef<- paste0(e1$Chromosome, ":", e1$Position, "_", e1$Allele2, "_", e1$Allele1)
  
  print("made it to line 47 loaded data sets\n")
  
  
  #e1_comm<- e1[e1$SNP %in% e2$SNP,]                   
  #e1_comm<- e1_comm[e1_comm$SNP %in% out$SNP_ID,]     ## making e1_comm equal to the matched SNPs from exposure 1 and exposure 2, then matching the common SNPs to the SNP_ID in the outcome
  
  e1_a<- e1[e1$ChrPosRefAlt %in% e3$rs_number,]        ## matching the drinks ChrPosRefAlt to the MPV rs_number
  e1_b<- e1[e1$ChrPosAltRef %in% e3$rs_number,]
  e1_comm<- bind_rows(e1_a,e1_b)                       ## bind both the combinations of alleles and bind rows so that you can find the common SNPs
  
  #e2_comm<- e2[e2$SNP %in% e1_comm$SNP_ID,]
  #out_comm<- out[out$SNP_ID %in% e1_comm$SNP,]
  
  #if using GRCh17_to_rsid files instead of primary Vuckovic files
  #e1_comm<- e1[e1$RSID %in% e3$name,]
  #e3_comm<- e3[e3$name %in% e1_comm$RSID,]
  
  e3_joina<- left_join(e1_a, e3, by= c("ChrPosRefAlt" = "rs_number"))   ## LEFT join  merge  two data frames where the merge returns all of the rows from the left side and any matching rows from the second table,  a merge operation between two data frames where the merge returns all of the rows from one table (the left side) and any matching rows from the second table(all rows from the matching SNPs between original and filtered - drinks per week number with HGB GWAS SNPs, and the matching rows from the Outcome
  
  e3_joinb<- left_join(e1_b, e3, by= c("ChrPosAltRef" = "rs_number"))
  e3_comm<- bind_rows(e3_joina, e3_joinb)                               ##combines rows from a and b on the bottom
  
  write.table(e1_comm, "Drinks.filtered.txt", sep="\t", quote= FALSE, row.names= FALSE, col.names= TRUE)  ## Write table prints the first argument from the second argument which is the file
  #write.table(e2_comm, "BMI.filtered.txt", sep="\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
  write.table(e3_comm, "MPV.filtered.txt", sep="\t", quote= FALSE, row.names= FALSE, col.names= TRUE) # this should match what is later output as 'MPV.filtered'
  #write.table(out_comm, "TC.filtered.txt", sep="\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
}  


##quit()
  
######   ClumpingFilteredSNPsets.SigOnly.R
######   Change outcome phenotype to alter working directory!

#library(TwoSampleMR)
###library(MRPRESSO)
#library(dplyr)
#library(stringr)
#library(ggplot2)
#library(data.table)

####Change WD!!
#setwd("/Users/thomc/Downloads/MVMR.Blood.sig/")


if( file.exists("Drinks.Clumped")) {
  print("Working from pre-tabulated clump files")
  e1_clump <- fread("Drinks.Clumped", header=T)
  e1_sig <- subset(e1_comm, PVALUE <= 5E-8)                ##selects only values from e1 that are statistically significant
  #out_clump <- fread("HGB.Clumped", header=T)
 
} else {
  
  e1_sig <- subset(e1_comm, PVALUE <= 5E-8)
  e1_data <- format_data(e1_sig, type="exposure", snps=e1_sig$RSID, header=T, phenotype_col = "Phenotype",
                         snp_col="RSID",
                         beta_col="BETA", se_col="SE",
                         eaf_col="AF", effect_allele_col="ALT",          ##Pval.exposure
                         other_allele_col="REF", pval_col="PVALUE")     ##formatting the data, telling the program what each column stands for, transferring the values to what 2 sample MR needs
  e1_clump <- clump_data(e1_data, clump_kb = 500, clump_r2 = 0.01)
  write.table(e1_clump, "Drinks.Clumped", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  print("made it to line 97 clumped exposure\n")
}

################## MR

# load Lifetime Smoking data and an example Astle trait to filter and clump/prepare the instrument. This is already formatted
####if doing on a new trait, re-do the stuff above in order to format and clump!
#Drinks <- fread("Drinks.filtered.txt", header=T)
#Drinks_data <- fread("Drinks.Clumped", header=T)
#Drinks$Phenotype <- "Drinks"

#Drinks_data <- format_data(Drinks, type="exposure", snps=Drinks$SNP, header=T,
#                       snp_col="SNP", phenotype_col = "Phenotype",
#                       beta_col="BETA", se_col="SE",
#                       eaf_col="AF", effect_allele_col="ALT",
#                       other_allele_col="REF", pval_col="PVALUE")

#Load and run MR for Drinks->Each blood trait
files <- list.files(path="/mnt/isilon/thom_lab/thomlab_datasets/GWAS/59_2020BloodTraits/", pattern="*EA_GWAMA.out.gz", full.names=FALSE, recursive=FALSE)

print("made it to line 118 got file list\n")      ## finds the pattern and creates a list of files for each blood trait. * wild card,  and removes the path

lapply(files, function(x) {                       ##  data frame as input and gives output in the form of a list object, everything after bracket is "x" is done for each blood trait

  t <- fread(paste0("/mnt/isilon/thom_lab/thomlab_datasets/GWAS/59_2020BloodTraits/", x), header=TRUE)
  t1 <- inner_join(e1_sig, t, by=c("ChrPosRefAlt" = "rs_number")) ## finds matching values for both combination of alleles
  t2 <- inner_join(e1_sig, t, by=c("ChrPosAltRef" = "rs_number"))
  t <- bind_rows(t1, t2)	                                        ##lists all combinations under t
  t_clumped <- t[t$RSID %in% e1_sig$RSID,]                        ## filters so only significant RS_IDs are listed, all the t data (rows) where the value present in t rsid is also found in e1 sig

  #splt at _ and use 1st sapply for real ; use   . and 2 for quick plt    ## ##  data frame as input and gives output in the form of a vector or matrix

  spl <- strsplit(x, "_" ,fixed=TRUE)                                     ##splits the input string vector into sub-strings
  trait <- sapply(spl, "[", 2) #first element is BCX2, 2nd is the trait   ## []for every data point, looking for second variable
  t_clumped$Phenotype <- trait
    
  # format blood trait as outcome
  t_data <- format_data(t_clumped, type="outcome", snps=t$SNP, header=T,
                       snp_col="SNP", phenotype_col = "Phenotype",
                       beta_col="beta", se_col="se",   # when already formatted might be 'beta'
                       eaf_col="eaf", effect_allele_col="reference_allele", #or EAF, A1  ... or effect_allele, eaf
                       other_allele_col="other_allele", pval_col="p-value")  # A2, P .... or other_allele and pval
    
  # clump and harmonize Drinks_data to blood trait
  #Drinks_data <- Drinks_data %>% filter(SNP %in% t_data$SNP)
  #Drinks_data_clump <- clump_data(Drinks_data)
  dat1 <- harmonise_data(exposure_dat = e1_clump, outcome_dat = t_data)     ##harmonise data so it is the same alleles 

  #print out instrument stats (for supp table) and F statistic (based on Shim et al PMID 25898129). 
  #print("this is R2 for the F-statistic based on Shim et al - plug in to mRnd") 
  # sample size for Astle is ~ 131564, Wootton et al LfSmk = 462690, Liu et al Drinks = 1232091
  dat1$r2 <- 2*(dat1$beta.exposure)^2 * dat1$eaf.exposure * (1-dat1$eaf.exposure) /
    (2*(dat1$beta.exposure)^2 * dat1$eaf.exposure * (1-dat1$eaf.exposure) +
       (dat1$se.exposure)^2*2*1232091*dat1$eaf.exposure * (1-dat1$eaf.exposure))
  r2 <- sum(dat1$r2)
  
  write.table(dat1, paste0("Drinks.", trait, ".IV.txt"), quote=F,col.names=T,row.names=F,sep="\t")  ##instrumental variable that has all the harmonized data that is used in results
  #run MR
  res <- mr(dat1) 
  
  #generate HTML file
  #mr_report(dat1)
  
  #print output
  write.table(res, paste0("Drinks.", trait, ".MR.txt"), quote=F,col.names=T,row.names=F,sep="\t")
  
  #heterogeneity test
  het <- mr_heterogeneity(dat1)  ##how variable each of the SNPs are
  write.table(het, paste0("Drinks.", trait, ".het.txt"), quote=F,col.names=T,row.names=F,sep="\t")
  #horiz pleiotropy test
  hp <- mr_pleiotropy_test(dat1)  ## checks to see if SNPs are independently affecting the exposure and outcome or if they are SNPs - exposure- outcome (check pvalue)
  write.table(hp, paste0("Drinks.", trait, ".hp.txt"), quote=F,col.names=T,row.names=F,sep="\t")
  



  #Plots

  #Scatterplot
  pdf(paste0("Drinks.", trait, ".Scatterplot.pdf"))  ##create a pdf first?
  p1 <- mr_scatter_plot(res, dat1)
  print(p1[[1]])                            ## purpose of the [1]?
  dev.off()                                 ##DEV.OFF function?
  
  #Forest plot
  pdf(paste0("Drinks.", trait, ".Forestplot.pdf"))  
  res_single <- mr_singlesnp(dat1)
  p2 <- mr_forest_plot(res_single)
  print(p2[[1]])
  dev.off()
  
  pdf(paste0("Drinks.", trait, ".OtherForestplot.pdf"))
  res_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_two_sample_ml"))
  p3 <- mr_forest_plot(res_single)
  print(p3[[1]])
  dev.off()
  
  #Leave One Out Plot
  pdf(paste0("Drinks.", trait, ".LeaveOneOutplot.pdf"))
  res_loo <- mr_leaveoneout(dat1)
  p4 <- mr_leaveoneout_plot(res_loo)
  print(p4[[1]])
  dev.off()
  
  #Good forest plot -- taking this out for now, as seems to fail. This worked well with dplyr/ggplot2 in past though
  #pdf(paste0("Drinks.", trait, ".GoodForest.pdf"))
  #res_single <- mr_singlesnp(dat1, all_method = c("mr_ivw",
  #                                                "mr_egger_regression",
  #                                                "mr_weighted_median"))
  #singlesnp_results <- res_single
  #exponentiate <- FALSE
  # 
  #requireNamespace("ggplot2", quietly = TRUE)
  #requireNamespace("plyr", quietly = TRUE)
  #res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d) {
  #  d <- plyr::mutate(d)
  #  
  #  if (sum(!grepl("All", d$SNP)) < 2) {
  #    return(blank_plot("Insufficient number of SNPs"))
  #    }
  #  
  #  levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "Inverse variance weighted"
  #  levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "MR Egger"
  #  levels(d$SNP)[levels(d$SNP) == "All - Weighted median"] <- "Weighted median"
  #  #d[d$SNP == "All - Inverse variance weighted", "SNP"] <- "Inverse variance weighted"
  #  #d[d$SNP == "All - MR Egger", "SNP"] <- "MR Egger"
  #  #d[d$SNP == "All - Weighted median", "SNP"] <- "Inverse variance weighted"
  #  am <- grep("All", d$SNP, value = TRUE)
  #  d$up <- d$b + 1.96 * d$se #or 
  #  d$lo <- d$b - 1.96 * d$se
  #                   
  #  # change unit
  #  # continuous gets exp(beta), binary is exp(ln(2)*beta)
  #  d$b <- exp(d$b) # d$b * log(2) for binary exposure
  #  d$up <- exp(d$up) # d$b * log(2) for binary exposure
  #                   d$lo <- exp(d$lo) # d$b * log(2) for binary exposure
  #                    
  #                   d$tot <- 0.01
  #                   d$tot[d$SNP %in% am] <- 1
  #                   d$SNP <- as.character(d$SNP)
  #                   nom <- d$SNP[!d$SNP %in% am]
  #                   nom <- nom[order(d$b)]
  #                   d <- rbind(d, d[nrow(d), ])
  #                   #d$SNP[nrow(d) - 1] <- ""
  #                   #d$b[nrow(d) - 1] <- NA
  #                   #d$up[nrow(d) - 1] <- NA
  #                   #d$lo[nrow(d) - 1] <- NA
  #                   d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
  #                   xint <- 0
  #                   if (exponentiate) {
  #                     d$b <- exp(d$b)
  #                     d$up <- exp(d$up)
  #                     d$lo <- exp(d$lo)
  #                     xint <- 1
  #                   }
  #                   #print(tail(d, 4))
  #                   d <- tail(d, 4)
  #                   d <- head(d, 3)
  #                   d[d$SNP == "Inverse variance weighted", "samplesize"] <- 3
  #                   d[d$SNP == "Weighted median", "samplesize"] <- 2
  #                   d[d$SNP == "MR Egger", "samplesize"] <- 1
  #                   d$SNP2 <- reorder(d$SNP, d$samplesize)
  #                   print(d)
  #                   ggplot2::ggplot(d, aes(y = SNP2, x = b)) + 
  #                     ggplot2::geom_vline(xintercept = 1, linetype = "dotted") + 
  #                     ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = up, size=as.factor(tot),colour = as.factor(tot)), size=4,
  #                              height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot)), size=16)  + 
  #                     ggplot2::scale_colour_manual(values = c("black", "red")) + 
  #                     ggplot2::scale_size_manual(values = c(0.3, 1)) + 
  #                     ggplot2::theme(legend.position = "none", 
  #                      axis.text.y = ggplot2::element_text(size = 14), 
  #                      axis.ticks.y = ggplot2::element_line(size = 0), 
  #                      axis.title.x = ggplot2::element_text(size = 14)) + 
  #                     ggplot2::labs(y = "", x = paste0("MR odds ratio for\n'", 
  #                                                      d$exposure[1], "' on '", d$outcome[1], "'"))
  #                 })
  #  print(res) #this prints if within Rstudio
    
})