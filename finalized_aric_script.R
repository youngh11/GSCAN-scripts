phenotypes <- read.table("/work/KellerLab/GSCAN/dbGaP/ARIC/PhenoGenotypeFiles/ChildStudyConsentSet_phs000090.ARIC_RootStudy.v3.p1.c1.HMB-IRB/PhenotypeFiles/phs000090.v3.pht000114.v2.p1.c1.GENEVA_ARIC_Subject_Phenotypes.HMB-IRB.txt.gz",header=T,sep="\t",stringsAsFactors=F)

phenotypes <- subset(phenotypes, select=c("geneva_id", "racegrp", "gender","v1age01", "anta01","anta04", "drnkr01", "hom29", 'hom35', "hom32", 'cigt01','evrsmk01', 'dtia90','dtia96', 'dtia97','dtia98', 'cursmk01','forsmk01'))

### rename phenotypes to be readable
names(phenotypes)[c(1,2,3,4,5,6)] <- c("geneva_id", "race", "sex", "age", "height", "weight")

### To connect sample ids to geneva ids, take SAMPID (the ID used in the genotype fam file, and SUBJID (aka geneva_id) )
id_map <- read.table(gzfile("/work/KellerLab/GSCAN/dbGaP/ARIC/PhenoGenotypeFiles/ChildStudyConsentSet_phs000090.ARIC_RootStudy.v3.p1.c1.HMB-IRB/GenotypeFiles/phg000035.v1.ARIC_GEI.genotype-qc.MULTI/geno-qc/samp-subj-mapping.csv.gz"),header=T,sep=",",stringsAsFactors=F)[,c(2,1)]
names(id_map) <- c("SAMPID", "geneva_id")
phenotypes <- merge(phenotypes, id_map, by="geneva_id", all.x=TRUE)


### import genotype data to get family info
fam_data <- read.table("/work/KellerLab/GSCAN/dbGaP/ARIC/PhenoGenotypeFiles/ChildStudyConsentSet_phs000090.ARIC_RootStudy.v3.p1.c1.HMB-IRB/GenotypeFiles/phg000035.v1.ARIC_GEI.genotype-calls-matrixfmt.c1.GRU.update1/Genotypes_with_flagged_chromosomal_abnormalities_zeroed_out/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out.fam", col.names = c("fam_id", "SAMPID", "patid", "matid", "sex", "dummy"))

### Replace 0's with "x" for rvTest preferred formatting
fam_data$patid <- fam_data$matid <- fam_data$fam_id[fam_data$fam_id == 0] <- "x"


########################################
###---- Derive GSCAN phenotypes -----###
########################################

### DRINKER VERSUS NON-DRINKER
### ARIC variable name is "drnkr01".
### Combination of "Do you presently drink alcoholic beverages?" and "Have you ever consumed alcoholic beverages?" 
###   Response option for both questions are "yes" or "no", which are turned into the options below.  
###   Response options:
###           1 = Current Drinker
###           2 = Former Drinker
###           3 = Never Drinker
###           4 = Unknown
###
### Descriptives: 
### table(phenotypes$drnkr01)
###    1    2    3    4 
### 7257 2309 3153    6 
###
###  To obtain GSCAN "DND" collapse across Former and Never Drinkers
###  and make "Non-Drinkers". Current Drinkers will be made "Drinkers"

dnd <- phenotypes$drnkr01
dnd[dnd == 1] <- "Current Drinker"
dnd[dnd == 2 | dnd == 3] <- 1
dnd[dnd == "Current Drinker"] <- 2
dnd[dnd == 4 | is.na(dnd)] <- "x"

### AGE OF INITIATION OF SMOKING
###
### ARIC variable name is "hom29".
###    "How old were you when you first started regular cigarette smoking?"
###    Response option is an integer value.
###
### Descriptives: 
###
### > table(phenotypes$hom29)
###     0    1    4    5    6    7    8    9   10   11   12   13   14   15   16   17
###    19    1    2   10   11   15   22   32   65   44  154  187  302  659  941  715
###    18   19   20   21   22   23   24   25   26   27   28   29   30   31   32   33
###  1219  567  703  447  275  129   88  247   56   45   59   22  100    8   34    8
###    34   35   36   37   38   39   40   41   42   43   44   45   46   47   48   49
###    17   51    9    9   10    8   35    3   11    5    5   16    3    4    2    3
###    50   51   52   55   57   59   60   62
###     7    1    2    1    1    2    2    1
###
### > summary(phenotypes$hom29)
###   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
###   0.00   16.00   18.00   18.77   20.00   62.00    5377 
###
ai <- phenotypes$hom29
### remove ages older than 35 and younger than 10
ai[ai > 35 | ai < 10 | is.na(ai)] <- "x"


### CIGARETTES PER DAY
### ARIC variable name is "hom35"
###    "On the average of the entire time you smoked, how many cigarettes did you usually smoke per day?"
###    Response option is integer, or "0" for <1 cigarette per day
###
### > table(phenotypes$hom35)
###     0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
###    46   54  118  187  133  254  146   84   89   17  990   23  106   30   12  520
###    16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31
###    31   30   62    5 2559    1    5   10    5  193    9    3    7    6  771    3
###    32   33   34   35   36   37   38   40   42   43   45   50   51   54   55   58
###     1    1    1   44    3    1    1  572    1    3   14   72    1    1    3    1
###    60   65   70   75   80   86   90   99
###   100    1    4    2   10    1    1    3
### > summary(phenotypes$hom35)
###    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
###    0.00   10.00   20.00   19.67   24.00   99.00    5420
### Responses are binned in accordance with the GSCAN Analysis Plan. 
cpd <- phenotypes$hom35
cpd[cpd <=  5 & cpd >=  1] <- 1
cpd[cpd <= 15 & cpd >=  6] <- 2
cpd[cpd <= 25 & cpd >= 16] <- 3
cpd[cpd <= 35 & cpd >= 26] <- 4
cpd[cpd >= 36 & cpd <= 60] <- 5
cpd[cpd > 60 | is.na(cpd)] <- "x"


### DRINKS PER WEEK
### ARIC variable names are "dtia96", "dtia97", and "dtia98"
###     "dtia96" - "How many glasses of wine do you usualy have per week? (4oz. glasses; round down)." 
###     "dtia97" - "How many bottles of cans or beer do you usualy have per week? (12oz. bottles or cans; round down)." 
###     "dtia98" - "How many drinks of hard liquor do you usualy have per week? (4oz. glasses; round down)."
###      Response option for all three is integer. 
###
###  Descriptives:
###
### >table(phenotypes$dtia96)
###  0    1    2    3    4    5    6    7    8    9   10   11   12   14   15   16 
###  5226  844  461  255  147   90   75   50   27    3   34    1   15   28    9    2 
###  17   18   20   21   25   28   30   32   33   35   40 
###   1    3    7    5    1    2    3    1    1    1    1 
###
### >summary(phenotypes$dtia96)
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
###  0.000   0.000   0.000   0.868   1.000  40.000    5478 
###
### >table(phenotypes$dtia97)
### 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
### 4356  674  528  312  214  107  297   93   76   10   97    3  186    4   40   28 
###  16   18   19   20   21   22   23   24   25   28   30   32   33   35   36   40 
###   5   28    1   36   19    1    1   89    8    8   12    2    1   10    8    6 
###  42   45   48   49   50   56   60   63   70   72   80   92 
###  13    2    6    1    3    2    4    1    1    2    1    1 
###
### >summary(phenotypes$dtia97)
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
###  0.000   0.000   0.000   2.609   2.000  92.000    5474 
###
### >table(phenotypes$dtia98)
###  0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
### 4387  735  551  295  239  190  148  138   70   11  151   14   57    3  103   33 
###  16   17   18   20   21   24   25   26   27   28   30   32   33   34   35   36 
###   9    9    7   36   30    1   10    1    1    9    6    2    1    2    4    1 
###  39   40   44   45   47   48   50   51   52   54   55   56   63   64   75   77 
###   1    7    1    1    1    2    5    1    1    1    1    3    1    2    1    2 
###  90   99 
###   1    2 
###
### >summary(phenotypes$dtia98)
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
###  0.000   0.000   0.000   2.227   2.000  99.000    5483 

wine <- phenotypes$dtia96 # 1 drink = 4oz
beer <- phenotypes$dtia97 # 1 drink = 12oz
spirits <- phenotypes$dtia98 # 1 drink = 1.5oz

### muliply wine by 4 and divide wine by 5 to normalize to standard drink of 5 oz for wine
### Combine all alcohol types, left-anchor at 1, and log
dpw <- log((wine*4/5 + beer + spirits) + 1) 
dpw[is.na(dpw)] <- "x"


### SMOKING INITIATION
### ARIC variable name is "evrsmk01"
### The variable checks answers to "Have you ever smoked cigarettes?" and "Do you now smoke cigarettes?".
###    Response options are "yes" or "no". 
###
### Descriptives:
###
### >table(phenotypes$evrsmk01)
### 0    1 
### 5328 7434 
###
### >summary(phenotypes$evrsmk01)
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
### 0.0000  0.0000  1.0000  0.5825  1.0000  1.0000       9 

si <- phenotypes$evrsmk01
si[si == 1] <- 2
si[si == 0] <- 1
si[is.na(si)] <- "x"


### SMOKING CESSATION
### ARIC variable names are "cursmk01" and "forsmk01"
### Both varaiables take into account the questions: "Have you ever smoked cigarettes?" and "Do you now smoke cigarettes?"
###    Response options are "yes" or "no". 
###    Smoking inititation (si) is coded as "2" for "Smoker" if "yes" to "Have you ever smoked cigarettes?"
###            If a subsequent "yes" to "Do you now smoke cigarettes?", smoking cessation (sc) is coded as "2" for "Current Smoker".
###            If a subsequent "no" to "Do you now smoke cigarettes?",  smoking cessation (sc) is coded as "1" for "Former Smoker".       
###    Smoking inititation (si) is coded as "1" for "Non Smoker" if "no" to "Have you ever smoked cigarettes?" 

current.smoker <- subset(phenotypes, select=c("cursmk01"))
former.smoker <- subset(phenotypes, select=c("forsmk01"))
N <- nrow(phenotypes)
sc <- rep(NA, N)
for(i in 1:N){
  if(is.na(current.smoker[i,1]) | is.na(former.smoker[i,1])){
    sc[i] <- NA
  }
  else if (current.smoker[i,1]  == 0 & former.smoker[i,1] == 0){
    sc[i] <- NA
  }
  else if (current.smoker[i,1]  == 0 & former.smoker[i,1] == 1){
    sc[i] <- 1 ### former smokers are coded as 1
  }
  else if (current.smoker[i,1]  == 1 & former.smoker[i,1] == 0){
    sc[i] <- 2 ### current smokers are coded as 2
  }
}
sc[is.na(sc)] <- "x"



### Create dataframe with our new GSCAN variables
N <- nrow(phenotypes)
NAs <- rep("x", N)
gscan.phenotypes <- data.frame(SAMPID = phenotypes$SAMPID, famid = NAs,geneva_id = phenotypes$geneva_id,patid = NAs,matid = NAs,sex = ifelse(phenotypes$sex == "M", 1, 2),cpd = cpd,ai = ai,si = si,sc = sc,dnd = dnd,dpw = dpw,age = phenotypes$age,age2 = phenotypes$age^2,height = phenotypes$height,weight = phenotypes$weight,currentformersmoker = sc)
                               
### Merge in the SAMPID, which is used in the genotype files
gscan.phenotypes <- merge(gscan.phenotypes, id_map, by="geneva_id", all.x=TRUE) 
                               
### Reorder phenotype file to make pedigree file consistent with genotype IDs
gscan.phenotypes <- gscan.phenotypes[c(17,2,16,3:15)]
colnames(gscan.phenotypes) [2] <- "SAMPID" 
                               
### Read in PCs and add to pedigree file, then write out to a phenotype and covariate file
### [ here read in PCs and merge into phenotype file (probably by the SAMPID) ]
pcs <- read.table("/rc_scratch/hayo0753/aric/aric_ancestry_and_pcs", head=TRUE, stringsAsFactors=F)
colnames(pcs) [1] <- "SAMPID"
pcs <- merge(pcs, gscan.phenotypes, by="SAMPID", all.x=TRUE)
                               
                               
############# PRELIMINARY #####################
                               
### Write to file [NOTE TO HANNAH: will have to be changed once we
### have PCs and ancestry groups identified. PCs will have to be read
### in like with read.table() and we'll have to subset the dataset
### into European and African ancestry, and then write out one
### phenotype and covariate file per ancestry group.
                               

### EUROPEANS
phenotypes.EUR.ped <- subset(pcs, ancestry == "EUR",select=c("famid","SAMPID","patid","matid", "sex","cpd", "ai","si", "sc", "dnd","dpw"))
write.table(phenotypes.EUR.ped, file="ARIC.EUR.phenotypes1.ped", quote=F, col.names=T, row.names=F, sep="\t")
                               
covariates.EUR.ped <- subset(pcs, ancestry == "EUR", select=c("famid","SAMPID","patid","matid", "sex","age", "age2", "height", "weight", "currentformersmoker","PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8","PC9", "PC10"))
write.table(covariates.EUR.ped, file="ARIC.EUR.covariates1.ped", quote=F, col.names=T, row.names=F, sep="\t")

phenotypes.AFR.ped <- subset(pcs, ancestry == "AFR",select=c("famid","SAMPID","patid","matid", "sex","cpd", "ai","si", "sc", "dnd","dpw"))
write.table(phenotypes.AFR.ped, file="ARIC.AFR.phenotypes1.ped", quote=F, col.names=T, row.names=F, sep="\t")
 
covariates.AFR.ped <- subset(pcs, ancestry == "AFR", select=c("famid","SAMPID","patid","matid", "sex","age", "age2", "height", "weight", "currentformersmoker", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8","PC9", "PC10"))
write.table(covariates.AFR.ped, file="ARIC.AFR.covariates1.ped", quote=F, col.names=T, row.names=F, sep="\t")
                               

### Must remove duplicates from all files, use UNIX
### sort ARIC.EUR.phenotypes1.ped | uniq > ARIC.EUR.phenotypes.ped
### sort ARIC.EUR.covariates.ped | uniq > ARIC.EUR.covariates.ped
### sort ARIC.AFR.phenotypes1.ped | uniq > ARIC.AFR.phenotypes.ped
### sort ARIC.AFR.covariates.ped | uniq > ARIC.AFR.covariates.ped




