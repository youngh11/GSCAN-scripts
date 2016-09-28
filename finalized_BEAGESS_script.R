beagess_data <- read.table("/work/KellerLab/GSCAN/dbGaP/Barretts/PhenoGenotypeFiles/RootStudyConsentSet_phs000869.BEAGESS.v1.p1.c1.GRU-MDS/PhenotypeFiles/phs000869.v1.pht004610.v1.p1.c1.BEAGESS_Subject_Phenotypes.GRU-MDS.txt.gz", header=TRUE, sep="\t", stringsAsFactors=F)
geno_data <- read.table("/work/KellerLab/GSCAN/dbGaP/Barretts/PhenoGenotypeFiles/RootStudyConsentSet_phs000869.BEAGESS.v1.p1.c1.GRU-MDS/GenotypeFiles/phg000580.v1.NCI_BEAGESS.genotype-calls-matrixfmt.c1.GRU-MDS.update/BEAGESS_dbGaP_29Jan2015.fam", header=FALSE, stringsAsFactors=F)
geno_data <- geno_data[-c(6)]
geno_data$V3[geno_data$V3 =="0"] = "x"
geno_data$V4[geno_data$V4 =="0"] = "x"

### IMPORTANT: All columns in original phenotype file are shifted one column to the the right, so the label of the column does not match the data in that column!!!
si <- beagess_data$BMI
sc <- beagess_data$BMI
age <- beagess_data$sex

beagess_data <- cbind(beagess_data, beagess_data[,6])
colnames(beagess_data)[13] <- "sc"


### SMOKING INITIATION 
###
### BEAGESS variable name is "cig_smk_status" and is listed under the column "BMI".
###    "Cigarette smoking status."
###    Response options are "-99", "-9", "0", "1" and "2".
###       -99 = not consented
###       -9 = Missing
###        0 = Never 
###        1 = Former
###        2 = Current
###    
### Descriptives: 
###
### > table(beagess_data$BMI)
### 
### -99   -9    0    1    2 
### 494 2051 1515 2288  575 
###
### > summary(beagess_data$BMI)
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
### -99.000  -9.000   0.000  -9.234   1.000   2.000 

beagess_data$BMI[beagess_data$BMI == -99] = "x"
beagess_data$BMI[beagess_data$BMI == -9] = "x"
beagess_data$BMI[beagess_data$BMI == 1] = 2
beagess_data$BMI[beagess_data$BMI == 0] = 1

### SMOKING Cessation 
###
### BEAGESS variable name is "cig_smk_status" and is listed under the column "BMI".
###    "Cigarette smoking status."
###    Response options are "-99", "-9", "0", "1" and "2".
###       -99 = not consented
###       -9 = Missing
###        0 = Never 
###        1 = Former
###        2 = Current
### 
### Descriptives: 
###
### > table(beagess_data$BMI)
### 
### -99   -9    0    1    2 
### 494 2051 1515 2288  575 
###
### > summary(beagess_data$BMI)
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
### -99.000  -9.000   0.000  -9.234   1.000   2.000 

beagess_data$sc[beagess_data$sc == -99] = "x"
beagess_data$sc[beagess_data$sc == -9] = "x"
beagess_data$sc[beagess_data$sc == 0] = "x"

### AGE 
###
### BEAGESS variable name is "agegpcat" and is listed under the column "sex".
###   "Age in 5 year categories"
###   Response options are integers 1 through 14 or -9. 
###       -9 = Missing
###       1 = 15-29 years of age
###       2 = 30-34 years of age
###       3 = 35-39 years of age
###       4 = 40-44 years of age
###       5 = 45-49 years of age
###       6 = 50-54 years of age
###       7 = 55-59 years of age
###       8 = 60-64 years of age
###       9 = 65-69 years of age
###       10 = 70-74 years of age
###       11 = 75-79 years of age
###       12 = 80-84 years of age
###       13 = 85-89 years of age
###       14 = 90-100 years of age
###
### Descriptives: 
###
### > table(beagess_data$sex)
### 
### -9    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
### 27   19   48  112  190  378  657  958 1134 1238 1018  794  246   87   17 
###
### > summary(beagess_data$BMI)
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
### -9.000   7.000   8.000   8.227  10.000  14.000 

beagess_data$sex[beagess_data$sex== 1] = 22
beagess_data$sex[beagess_data$sex== 2] = 32
beagess_data$sex[beagess_data$sex== 3] = 37
beagess_data$sex[beagess_data$sex== 4] = 42
beagess_data$sex[beagess_data$sex== 5] = 47
beagess_data$sex[beagess_data$sex== 6] = 52
beagess_data$sex[beagess_data$sex== 7] = 57
beagess_data$sex[beagess_data$sex== 8] = 62
beagess_data$sex[beagess_data$sex== 9] = 67
beagess_data$sex[beagess_data$sex== 10] = 72
beagess_data$sex[beagess_data$sex== 11] = 77
beagess_data$sex[beagess_data$sex== 12] = 82
beagess_data$sex[beagess_data$sex== 13] = 87
beagess_data$sex[beagess_data$sex== 14] = 95
beagess_data$sex[beagess_data$sex== -9] = "x"
beagess_data$sex <- as.numeric(beagess_data$sex)

### rename SUBJID column to dbGaP_Subject_ID in genotype file so it matches phenotype file where SUBJID is misslabelled as dbGaP_Subject_ID
colnames(geno_data) <- c("famid", "dbGaP_Subject_ID", "patid", "matid", "sex")

### merge geno and pheno files
phen <- merge(geno_data,beagess_data, by="dbGaP_Subject_ID", all =TRUE)
phen <- phen[,c(1:5,7,10,17)]

phen <- phen[,c(2,1,3,4,5,6,7,8)]
colnames(phen) <- c("famid", "dbGaP_Subject_ID", "patid", "matid", "sex", "age","si","sc")


phenotypes <- data.frame(famid=phen$famid, dbGaP_Subject_ID=phen$dbGaP_Subject_ID, patid=phen$patid, matid=phen$matid, sex=phen$sex, si=phen$si, sc=phen$sc, currentformersmoker=phen$sc, age=phen$age, age2=phen$age^2) 
colnames(phenotypes) <- c("famid", "dbGaP_Subject_ID", "patid", "matid", "sex", "si", "sc", "currentformersmoker", "age", "age2")

pcs <- read.table("/rc_scratch/hayo0753/BEAGESS/BEAGESS_pcs_and_ancestrys.txt",head=TRUE, stringsAsFactors=F, sep=" ")
colnames(pcs) [1] <- "dbGaP_Subject_ID"

gscan.phenotypes <- merge(phenotypes, pcs, by="dbGaP_Subject_ID", all=TRUE)
gscan.phenotypes[is.na(gscan.phenotypes)] <- "x" 

### EUROPEANS - entire sample is EUR
phenotypes.EUR.ped <- subset(gscan.phenotypes, ancestry == "EUR", select=c("famid","dbGaP_Subject_ID","patid","matid", "sex", "si", "sc"))
write.table(phenotypes.EUR.ped, file="BEAGESS.EUR.phenotypes.ped", quote=F, col.names=T, row.names=F, sep="\t")

covariates.EUR.ped <- subset(gscan.phenotypes, ancestry == "EUR", select=c("famid","dbGaP_Subject_ID","patid","matid", "sex", "age", "age2",  "currentformersmoker", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
write.table(covariates.EUR.ped, file="BEAGESS.EUR.covariates.ped", quote=F, col.names=T, row.names=F, sep="\t")

### TODO
### make age^2 values and currentformersmoker
