



############### STEP 1 #########################################################

#### extracting repeat length from VCF files using python and reading in 
#### the phenotypes of interest from the RAP

import os
import gzip
import dxpy
import multiprocessing
import numpy as np
import pandas as pd
from tqdm import tqdm
from glob import glob


pip install numpy==1.25.2 --no-cache-dir
pip install --no-cache-dir --force-reinstall scipy
pip install --no-cache-dir --force-reinstall dxdata
pip install tables --no-cache-dir

finaltypes = {'chr': str,
  'refstart': int,
  'refend': int,
  'refcount': int,
  'reflength': int,
  'repeatunit': str,
  'allele1': float,
  'allele2':float,
  'eid': int}


def find_headerline(path):
  with gzip.open(path, 'r') as file:
  for idx, line in enumerate(file, start=0):
  if line.startswith(b'#CHROM'):
  return idx

def parse_vcfinfo(info: pd.Series):
  varinfo = info.str.split(';', expand=True)[[0,1,2,3,4]]
varinfo = varinfo.apply(lambda x: x.str.split('=', expand=True)[1], axis=0)
return varinfo.rename(columns={0:'refend',1:'refcount',2:'reflength',3:'repeatunit',4:'varid'})

def generate_alleles(vcfsample: pd.Series):
  genotype = vcfsample.str.split(':', expand=True)[2]
alleles = genotype.str.split('/', expand=True)
return alleles.rename(columns={0: 'allele1', 1: 'allele2'})

def generate_allele_models(alleles: pd.Series, varinfo: pd.Series):
  out = pd.DataFrame()
out['totalcount'] = alleles['allele1'] + alleles['allele2']
out['maxcount'] = alleles[['allele1','allele2']].max(axis=1)
out['totalratio'] = out['totalcount'] / (2 * varinfo['refcount'])



def parse_vcf_to_tabular(eid, path):
  data = pd.read_csv(path, sep='\t', header=find_headerline(path))
os.remove(path)
data.rename(columns={'#CHROM':'chr', 'POS':'refstart'}, inplace=True)

varinfo = parse_vcfinfo(data['INFO'])
alleles = generate_alleles(data[data.columns[-1]])

out = data[['chr','refstart']].join(varinfo).join(alleles)
out.replace('.', np.nan, inplace=True)
out['eid'] = eid
autosomal = out[out['chr']!='chrX'].dropna(how='any')
xchrom = out[out['chr']=='chrX']########################
return pd.concat([autosomal, xchrom]).astype(finaltypes).set_index('varid', drop=True)



def get_vcflist(vcflist):
  with open(vcflist, 'r') as file:
  return [x.strip('\n') for x in file.readlines()]

def initiate_vcffile(dxid):
  eid = dxpy.DXFile(dxid).name.split('_')[0]
localname = f"{eid}_calls.vcf.gz"
# Download to instance, overwrite last copy and verify that change
dxpy.download_dxfile(dxid, localname)
return eid, localname


def initiate_vcf_and_parse(dxid):
  eid, localname = initiate_vcffile(dxid)
return parse_vcf_to_tabular(eid, localname)

def process_vcflist(path: str):
  vcflist = get_vcflist(path)
proc = multiprocessing.cpu_count()
print(f'Available CPUs for pool: {proc}')
with multiprocessing.Pool(proc) as pool:
  return list(tqdm(pool.imap(initiate_vcf_and_parse, vcflist), total=len(vcflist)))


def store_to_hdf(results, eidblock):
  hdfpath = f'strcalls_clinsubset_dragen_male_eid{eidblock}.hdf'

if os.path.exists(hdfpath):
  os.remove(hdfpath)

with pd.HDFStore(hdfpath, 'w') as store:
  for x in tqdm(results):
  store.put(f'eids_{eidblock}_prefix', x, 
            format='table', append=True, 
            data_columns=['eid'], 
            complib='blosc', 
            complevel=9)
return hdfpath


for i in range(50,61):
  results = process_vcflist(f'/mnt/project/resources/strcalls/dragen_exphunter_vcflists/male/strcalls_dragen_vcflist_male_eid{i}.txt')
localhdf = store_to_hdf(results, i)
#dxpy.upload_local_file(localhdf)



def store_to_hdf(results, eidblock):
  # Define the path for the HDF file (you can add a local directory like /tmp/ if needed)
  hdfpath = f'strcalls_clinsubset_dragen_female_eid{eidblock}.hdf'

# If the file exists, remove it
if os.path.exists(hdfpath):
  os.remove(hdfpath)

# Open the HDFStore and write the data
with pd.HDFStore(hdfpath, 'w') as store:
  for x in tqdm(results):
  store.put(f'eids_{eidblock}_prefix', x, 
            format='table', append=True, 
            data_columns=['eid'], 
            complib='blosc', 
            complevel=5)

return hdfpath


for i in range(11,12):
  results = process_vcflist(f'/mnt/project/resources/strcalls/dragen_exphunter_vcflists/female/strcalls_dragen_vcflist_female_eid{i}.txt')
localhdf = store_to_hdf(results, i)
#dxpy.upload_local_file(localhdf)


import dxpy
import os
from glob import glob
import pandas as pd
import numpy as np
import tables



# Load up FMR1 STR genotypes from processed HDFs (must already be dx downloaded to location shown)
# fmr1 dataframe produced shows a1 and a2 repeat counts, NOT base counts (if it were base counts, they'd need to be a multiple of 3)
frames = list()
for x in glob(f'/mnt/project/resources/strcalls/dragen_exphunter_hdfs_male/*.hdf'):
  with pd.HDFStore(x, 'r') as store:
  i = os.path.splitext(os.path.basename(x))[0].split('_')[4].lstrip('eid')
subframe = store.select(f"eids_{i}_prefix", 
                        where='index=FMR1', 
                        columns=['eid','allele1','allele2','refcount'])
frames.append(subframe)
fmr1_male = pd.concat(frames).set_index('eid', drop=True)


data=fmr1_male


data.to_csv('fmr1_male.txt', sep='\t', index=True)



#### extracting the relevant phenotypes 


# Pull in ovarian ageing phenotypes from group projects and own covariates
thisproj = 'project-GgB82BQJK6Y3469x45fkP727'
phensandbox = 'project-G7G4YgjJYVkPQ93jPpFY2pxg'
dxpy.download_folder(thisproj, 'dragen_exphunter_hdfs', '/resources/strcalls/dragen_exphunter_hdfs')
dxpy.download_dxfile("file-G82b4xjJYVk24bfJ47QP4F87", 'anm.txt', project = phensandbox)
dxpy.download_dxfile("file-G7p0kZ0JYVk89zK628525Vkp", 'aam.txt', project = phensandbox)
dxpy.download_dxfile("file-G82b58QJYVkJb1J03Bq2Xpj3", 'poi.txt', project = phensandbox)
dxpy.download_dxfile("file-Gjgp1X0JX6GPxG7Jvkb762k7", 'births.txt')
dxpy.download_dxfile("file-GkZ6yGQJ7vyX57JpJjY7Fj8J", 'covars.txt')



# Read in group phenotypes and index on EID
anm = pd.read_csv('anm.txt', delim_whitespace=True, usecols=['IID','nat_mage_min34']).set_index('IID', drop=True)
aam = pd.read_csv('aam.txt', delim_whitespace=True, usecols=['IID','AAM_all']).set_index('IID', drop=True)
poi = pd.read_csv('poi.txt', delim_whitespace=True, usecols=['IID','nat_mage_POI']).set_index('IID', drop=True)

# Derive early and late menarche categorical
aam['relative_aam'] = pd.cut(aam['AAM_all'], [0,11,15,27], right=False, labels=['early','normal','late'])


# Read in and process covars, index on EID
covars = pd.read_csv('covars.txt', sep='\t').set_index('eid', drop=True)
covars.rename(columns={'p22001': 'sex', 
  'p21022': 'age', 
  'p22006': 'eur', 
  'p22021': 'kinship',
  'p32060': 'coverage',
  'p32051': 'seqprovider',
  'p22000': 'arraybatch'}, inplace=True)



# Read in birth phenotypes and get derived phenotypes, index on EID
births = pd.read_csv('births.txt', sep='\t').set_index('eid', drop=True)
births['live_births'] = births[['p2734_i0','p2734_i1','p2734_i2','p2734_i3']].max(axis=1)
births['no_children'] = pd.cut(births['live_births'], [0,1,23], right=False, labels=['No children','At least one child'])
births['miscarriages'] = births[['p3839_i0','p3839_i1','p3839_i2','p3839_i3']].max(axis=1)
births['stillbirths'] = births[['p3829_i0','p3829_i1','p3829_i2','p3829_i3']].max(axis=1)
births['lost_pregnancy'] = np.where(((births['miscarriages']>0) | (births['stillbirths']>0)),1, np.NaN)
births['lost_pregnancy'] = np.where(((births['miscarriages']==0) & (births['stillbirths']==0)), 0, births['lost_pregnancy'])
births = births[['live_births','no_children','lost_pregnancy','miscarriages','stillbirths']]



# Load up FMR1 STR genotypes from processed HDFs (must already be dx downloaded to location shown)
# fmr1 dataframe produced shows a1 and a2 repeat counts, NOT base counts (if it were base counts, they'd need to be a multiple of 3)
frames = list()
for x in glob('dragen_exphunter_hdfs/*.hdf'):
  with pd.HDFStore(x, 'r') as store:
  i = os.path.splitext(os.path.basename(x))[0].split('_')[4].lstrip('eid')
subframe = store.select(f"eids_{i}_prefix", 
                        where='index=FMR1', 
                        columns=['eid','allele1','allele2','refcount'])
frames.append(subframe)
fmr1 = pd.concat(frames).set_index('eid', drop=True)


# There is one person with a full mutation - omit them 
fmr1 = fmr1[((fmr1['allele1']<200) & (fmr1['allele2']<200))]


# Create derived allele exposures
fmr1['allele_sum'] = fmr1['allele1'] + fmr1['allele2']
fmr1['max_allele'] = fmr1[['allele1','allele2']].max(axis=1)
premut1 = np.where(fmr1['allele1'].between(55,200), 1, 0)
premut2 = np.where(fmr1['allele2'].between(55,200), 1, 0)
fmr1['has_premut'] = np.where(((premut1==1) | (premut2==1)), 1, 0)
highrisk1 = np.where(fmr1['allele1'].between(70,100), 1, 0)
highrisk2 = np.where(fmr1['allele2'].between(70,100), 1, 0)
fmr1['has_highrisk'] = np.where(((highrisk1==1) | (highrisk2==1)), 1, 0)



data = fmr1.join([covars, anm, aam, poi, births], how='inner', validate='1:1')
data.index.name = 'IID'
assert all(data['sex']==0)


data.to_csv('fmr1_main_analytic.txt', sep='\t', index=True)


dxpy.upload_local_file('fmr1_main_analytic.txt', folder='/fmr1')












####### STEP 2 #################################################################


### phenotypic analysis with extracted repeat length in R

library(stringr)  
library(dplyr)   
library(tidyr)    
library(readr) 
library(ggplot2)
library(broom)
library(data.table)
library(ggExtra)
library(GGally)
library(ggeffects)
library(splines)
library(ggbreak)



# read in main dataframe
main_analytic <- fread("fmr1_main_analytic.txt")

# filter to keep only unrelated individuals and also remove those who have pulled out
include <- fread("INCLUDEFOR_EUR_Unrelated.txt")
include <- include %>% 
  rename(IID=V1)
main_analytic <- main_analytic %>% 
  filter(IID %in% include$IID)


# looking at the distribution of menopause age
ggplot(main_analytic, aes(x=nat_mage_min34)) +
  geom_histogram()


# making sure sequence provider is a factor variable not numeric
main_analytic$seqprovider <- as.factor(main_analytic$seqprovider)
main_analytic$seqprovider <- as.character(main_analytic$seqprovider)


# looking at some of the outcome variables a bit more closely
table(main_analytic$eur, useNA = "ifany")
table(main_analytic$kinship, useNA = "ifany")
table(main_analytic$sex, useNA = "ifany")
table(main_analytic$arraybatch, useNA = "ifany")
table(main_analytic$nat_mage_POI, useNA = "ifany")
table(main_analytic$max_allele, useNA = "ifany")
table(main_analytic$nat_mage_min34, useNA = "ifany")
table(main_analytic$p21000_i0, useNA = "ifany")
table(main_analytic$p22006, useNA = "ifany")
table(main_analytic$max_allele, useNA = "ifany")

table(main_analytic$has_premut, main_analytic$nat_mage_min34, useNA = "ifany")
sum(main_analytic$has_premut)


had_menopause <- fread("had_menopause_participant.tsv")

table(had_menopause$p2724_i0)

had_menopause <- had_menopause %>% 
  rename(IID=eid)

main_analytic_meno_info <- main_analytic %>% 
  select(IID, has_premut, nat_mage_min34)


meno_check <- left_join(main_analytic_meno_info, had_menopause, by="IID")

meno_check_premut <- meno_check %>% 
  filter(has_premut==1)

table(meno_check_premut$nat_mage_min34, meno_check_premut$p2724_i0, useNA = "ifany")

nrow(meno_check_premut)



filtered_data <- meno_check_premut %>%
  filter(p2724_i0 == "Yes" & is.na(nat_mage_min34))


### checking this against non premutation carriers


meno_check_nopremut <- meno_check %>% 
  filter(has_premut==0)

table(meno_check_nopremut$nat_mage_min34, meno_check_nopremut$p2724_i0, useNA = "ifany")

nrow(meno_check_nopremut)



filtered_data <- meno_check_nopremut %>%
  filter(p2724_i0 == "Yes" & is.na(nat_mage_min34))


### creating the clinical cut offs 

main_analytic <- main_analytic %>% 
  mutate(clinical_category = case_when(
    max_allele<44 ~ "Normal",
    max_allele>43&max_allele<55 ~ "Intermediate",
    max_allele>54&max_allele<200 ~ "Premutation",
    max_allele>199 ~ "Full mutation"
  ))


### checking for homozygosity

(sum(main_analytic$allele1==main_analytic$allele2))/nrow(main_analytic)



### reading in data on other outcome variables and confounders 

# menstrual cycles 

cycles <- fread("regular_cycles_20-50_days_inc_contraception.tsv")

contraception_covars <- fread("BOLT_500k_Covariates_200k_batch_contraception.tsv")

contraception_covars <- contraception_covars %>% 
  select(FID, IID, contraception_binary)

cycles <- merge(cycles, contraception_covars, by="FID", all.x = TRUE)

cyles <- cycles %>% 
  select(FID, Cycle_length, contraception_binary)

cycles <- cycles %>% 
  rename(eid=FID)

cycles <- cycles %>% 
  rename(IID=eid)

main_analytic <- merge(main_analytic, cycles, by="IID", all.x=TRUE)



### BMI

bmi <- fread("height_bmi_pheno_participant.tsv")

bmi <- bmi %>% rename(IID=eid,
                      BMI=p21001_i0)

main_analytic <- left_join(main_analytic, bmi, by="IID")


# SMOKING


smoking <- fread("smoking_participant.tsv")

smoking <- smoking %>% 
  rename(IID=eid,
         smoking_status=p20116_i0,
         smoking_pack_years=p20161_i0)

main_analytic <- left_join(main_analytic, smoking, by="IID")

main_analytic <- main_analytic %>%
  mutate(ever_smoked = ifelse(smoking_status=="Current"| smoking_status=="Previous", 1,0)) 




main_analytic <- main_analytic %>%
  mutate(clinical_category = factor(clinical_category, levels = c("Normal", "Intermediate", "Premutation")))

descriptive_info <- main_analytic %>%
  group_by(clinical_category) %>%
  summarise(
    n_women = n(),
    
    # Age
    mean_age = mean(age, na.rm = TRUE),
    age_sd = sd(age, na.rm = TRUE),
    age_se = age_sd / sqrt(n_women),
    age_ci_lower = mean_age - 1.96 * age_se,
    age_ci_upper = mean_age + 1.96 * age_se,
    age_median = median(age, na.rm = TRUE),
    age_iqr_low = quantile(age, 0.25, na.rm = TRUE),
    age_iqr_high = quantile(age, 0.75, na.rm = TRUE),
    
    # POI
    n_POI = sum(nat_mage_POI == 1, na.rm = TRUE),
    n_POI_total = sum(nat_mage_POI %in% c(0,1), na.rm = TRUE),
    prop_POI = n_POI / n_POI_total,
    poi_ci_lower = prop.test(n_POI, n_POI_total)$conf.int[1],
    poi_ci_upper = prop.test(n_POI, n_POI_total)$conf.int[2],
    
    # Repeat length
    mean_repeat_length = mean(max_allele, na.rm = TRUE),
    repeat_sd = sd(max_allele, na.rm = TRUE),
    repeat_se = repeat_sd / sqrt(n_women),
    repeat_ci_lower = mean_repeat_length - 1.96 * repeat_se,
    repeat_ci_upper = mean_repeat_length + 1.96 * repeat_se,
    repeat_median = median(max_allele, na.rm = TRUE),
    repeat_iqr_low = quantile(max_allele, 0.25, na.rm = TRUE),
    repeat_iqr_high = quantile(max_allele, 0.75, na.rm = TRUE),
    
    # Sum of alleles
    mean_sum_alleles = mean(allele_sum, na.rm = TRUE),
    sum_sd = sd(allele_sum, na.rm = TRUE),
    sum_se = sum_sd / sqrt(n_women),
    sum_ci_lower = mean_sum_alleles - 1.96 * sum_se,
    sum_ci_upper = mean_sum_alleles + 1.96 * sum_se,
    sum_median = median(allele_sum, na.rm = TRUE),
    sum_iqr_low = quantile(allele_sum, 0.25, na.rm = TRUE),
    sum_iqr_high = quantile(allele_sum, 0.75, na.rm = TRUE),
    
    # Age at menopause
    mean_meno_age = mean(nat_mage_min34, na.rm = TRUE),
    meno_sd = sd(nat_mage_min34, na.rm = TRUE),
    meno_se = meno_sd / sqrt(sum(!is.na(nat_mage_min34))),
    meno_ci_lower = mean_meno_age - 1.96 * meno_se,
    meno_ci_upper = mean_meno_age + 1.96 * meno_se,
    meno_median = median(nat_mage_min34, na.rm = TRUE),
    meno_iqr_low = quantile(nat_mage_min34, 0.25, na.rm = TRUE),
    meno_iqr_high = quantile(nat_mage_min34, 0.75, na.rm = TRUE),
    
    # Age at menarche
    mean_age_menarche = mean(AAM_all, na.rm = TRUE),
    menarche_sd = sd(AAM_all, na.rm = TRUE),
    menarche_se = menarche_sd / sqrt(sum(!is.na(AAM_all))),
    menarche_ci_lower = mean_age_menarche - 1.96 * menarche_se,
    menarche_ci_upper = mean_age_menarche + 1.96 * menarche_se,
    menarche_median = median(AAM_all, na.rm = TRUE),
    menarche_iqr_low = quantile(AAM_all, 0.25, na.rm = TRUE),
    menarche_iqr_high = quantile(AAM_all, 0.75, na.rm = TRUE),
    
    # Number of children
    mean_n_children = mean(live_births, na.rm = TRUE),
    children_sd = sd(live_births, na.rm = TRUE),
    children_se = children_sd / sqrt(sum(!is.na(live_births))),
    children_ci_lower = mean_n_children - 1.96 * children_se,
    children_ci_upper = mean_n_children + 1.96 * children_se,
    
    
    # Cycle length
    mean_cycle_length = mean(Cycle_length, na.rm = TRUE),
    cycle_sd = sd(Cycle_length, na.rm = TRUE),
    cycle_se = cycle_sd / sqrt(sum(!is.na(Cycle_length))),
    cycle_ci_lower = mean_cycle_length - 1.96 * cycle_se,
    cycle_ci_upper = mean_cycle_length + 1.96 * cycle_se,
    
    # Lost pregnancies
    mean_lost_pregnancy = mean(lost_pregnancy, na.rm = TRUE),
    lost_sd = sd(lost_pregnancy, na.rm = TRUE),
    lost_se = lost_sd / sqrt(sum(!is.na(lost_pregnancy))),
    lost_ci_lower = mean_lost_pregnancy - 1.96 * lost_se,
    lost_ci_upper = mean_lost_pregnancy + 1.96 * lost_se,
    lost_median = median(lost_pregnancy, na.rm = TRUE),
    lost_iqr_low = quantile(lost_pregnancy, 0.25, na.rm = TRUE),
    lost_iqr_high = quantile(lost_pregnancy, 0.75, na.rm = TRUE),
    
    # BMI
    mean_bmi = mean(BMI, na.rm = TRUE),
    bmi_sd = sd(BMI, na.rm = TRUE),
    bmi_se = bmi_sd / sqrt(sum(!is.na(bmi))),
    bmi_ci_lower = mean_bmi - 1.96 * bmi_se,
    bmi_ci_upper = mean_bmi + 1.96 * bmi_se,
    
    # Smoking
    n_smokers = sum(ever_smoked == 1, na.rm = TRUE),
    n_smoke_total = sum(ever_smoked %in% c(0, 1), na.rm = TRUE),
    prop_smoked = n_smokers / n_smoke_total,
    smoke_ci_lower = prop.test(n_smokers, n_smoke_total)$conf.int[1],
    smoke_ci_upper = prop.test(n_smokers, n_smoke_total)$conf.int[2]
  )


library(tidyr)

descriptive_info_wide <- descriptive_info %>%
  pivot_longer(-clinical_category, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = clinical_category, values_from = value) %>% 
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE, digits = 3)))


write_tsv(descriptive_info_wide, "updated_descriptive_info.tsv")




### adding pairwise comparisons with significance levels 

main_analytic <- main_analytic %>%
  filter(clinical_category %in% c("Normal", "Intermediate", "Premutation")) %>%
  mutate(clinical_category = factor(clinical_category, levels = c("Normal", "Intermediate", "Premutation")))

continuous_vars <- c("age", "BMI", "Cycle_length", "AAM_all", "live_births", "allele_sum", "nat_mage_min34")
binary_vars <- c("nat_mage_POI", "ever_smoked")

library(dplyr)

run_pairwise_tests <- function(var, data, type = c("continuous", "binary")) {
  type <- match.arg(type)
  
  res <- data.frame(
    variable = var,
    comparison = c("Intermediate vs Normal", "Premutation vs Normal"),
    p_value = NA_real_
  )
  
  for (i in 2:3) {  # compare levels 2 and 3 vs. level 1 (Normal)
    group1 <- data %>% filter(clinical_category == "Normal") %>% pull(var)
    group2 <- data %>% filter(clinical_category == levels(data$clinical_category)[i]) %>% pull(var)
    
    p <- tryCatch({
      if (type == "continuous") {
        t.test(group1, group2)$p.value
      } else {
        tbl <- table(c(rep("Normal", length(group1)), rep("Group", length(group2))),
                     c(group1, group2))
        if (any(tbl < 5)) fisher.test(tbl)$p.value else chisq.test(tbl)$p.value
      }
    }, error = function(e) NA)
    
    res$p_value[i - 1] <- p
  }
  
  return(res)
}

results_list <- list()

for (v in continuous_vars) {
  results_list[[v]] <- run_pairwise_tests(v, main_analytic, type = "continuous")
}

for (v in binary_vars) {
  results_list[[v]] <- run_pairwise_tests(v, main_analytic, type = "binary")
}

pairwise_results <- bind_rows(results_list) %>%
  mutate(p_value = signif(p_value, 3),
         significance = cut(p_value,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                            labels = c("***", "**", "*", "ns"),
                            right = FALSE))





### Looking at the distribution of the allele lengths 


ggplot(main_analytic, aes(x=allele1, y=allele2)) +
  geom_point()


### creating a new "min allele" variable to deal with the genotyping order issue
### (described in email to isobel)

main_analytic <- main_analytic %>% 
  mutate(min_allele=pmin(allele1, allele2))

main_analytic %>% 
  filter(nat_mage_POI==1) %>% 
  count(min_allele<45 & max_allele<45)

main_analytic %>% 
  filter(nat_mage_POI==1) %>% 
  count(max_allele<45)

ggplot(main_analytic, aes(x=min_allele, y=max_allele)) +
  geom_point() +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60) ) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200))



main_analytic <- main_analytic %>% 
  mutate(concurrence=ifelse(allele2==max_allele | allele1==allele2, 1, 0))

table(main_analytic$concurrence)

concurrence <- main_analytic %>% 
  filter(concurrence==0)


sum(main_analytic$max_allele==30)

ggplot(main_analytic, aes(x=allele_sum)) +
  geom_density()

table(main_analytic$allele_sum)

poi <- main_analytic %>% 
  filter(nat_mage_POI==1)

prop_below_45 <- poi %>%
  filter(max_allele < 45) %>%
  summarise(proportion = n() / nrow(poi)) %>%
  pull(proportion)

prop_below_45

prop_both_below_45 <- poi %>%
  filter(max_allele < 45, min_allele < 45) %>%
  summarise(proportion = n() / nrow(poi)) %>%
  pull(proportion)

prop_both_below_45

prop_max_above_55 <- poi %>%
  filter(max_allele >= 55) %>%
  summarise(proportion = n() / nrow(poi)) %>%
  pull(proportion)



poi_max_gt_45_min_lt_45 <- poi %>%
  filter(max_allele > 45, min_allele < 45)
nrow(poi_max_gt_45_min_lt_45)




library(dplyr)

poi_summary <- poi %>%
  mutate(
    premutation_min = min_allele >= 45,
    premutation_max = max_allele >= 45,
    premutation_category = case_when(
      premutation_min & premutation_max ~ "Two alleles in premutation",
      (premutation_min | premutation_max) & !(premutation_min & premutation_max) ~ "One allele in premutation",
      TRUE ~ "No alleles in premutation"
    )
  ) %>%
  group_by(premutation_category) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

poi_summary




# Subset of women without POI
non_poi <- main_analytic %>%
  filter(nat_mage_POI != 1)

# Summarize premutation status for non-POI group
non_poi_summary <- non_poi %>%
  mutate(
    premutation_min = min_allele >= 55,
    premutation_max = max_allele >= 55,
    premutation_category = case_when(
      premutation_min & premutation_max ~ "Two alleles in premutation",
      (premutation_min | premutation_max) & !(premutation_min & premutation_max) ~ "One allele in premutation",
      TRUE ~ "No alleles in premutation"
    )
  ) %>%
  group_by(premutation_category) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

non_poi_summary

library(dplyr)

# Proportion of POI women with at least one allele between 45 and 54 (inclusive)
poi_with_intermediate <- poi %>%
  filter((min_allele >= 45 & min_allele < 55) | (max_allele >= 45 & max_allele < 55))

# Calculate proportion
prop_intermediate <- nrow(poi_with_intermediate) / nrow(poi)

prop_intermediate




### MAX ALLELE
mean_menopause_age <- main_analytic %>%
  group_by(max_allele) %>%
  summarise(mean_age = mean(nat_mage_min34, na.rm = TRUE),
            n_women = n())

# Plot the results
ggplot(mean_menopause_age, aes(x = max_allele, y = mean_age)) +
  geom_smooth(method = "loess", color = "blue") +   # Smooth line with loess
  geom_point(color = "red") +                        # Points on the line
  theme_minimal() +
  labs(title = "Smooth Line: Mean Age at Menopause by Max Allele",
       x = "Repeat Length", 
       y = "Mean Age at Menopause") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = seq(0, max(mean_menopause_age$max_allele), by = 10))




#### SAME FOR MENARCHE


### MAX ALLELE
mean_menarche_age <- main_analytic %>%
  group_by(max_allele) %>%
  summarise(mean_age = mean(AAM_all, na.rm = TRUE),
            n_women = n())

# Plot the results
ggplot(mean_menarche_age, aes(x = max_allele, y = mean_age)) +
  geom_smooth(method = "loess", color = "blue") +   # Smooth line with loess
  geom_point(color = "red") +                        # Points on the line
  theme_minimal() +
  labs(title = "Smooth Line: Mean Age at Menopause by Max Allele",
       x = "Repeat Length", 
       y = "Mean Age at Menarche") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = seq(0, max(mean_menopause_age$max_allele), by = 10))





### plotting the distribution of allele length

table(main_analytic$max_allele, useNA = "ifany")


# Calculate mean menopause age for each max_allele group
mean_menopause_age <- main_analytic %>%
  group_by(max_allele) %>%
  summarise(mean_mage = mean(nat_mage_min34, na.rm = TRUE))


# Plot the density of max_allele
ggplot() +
  # Plot density of max_allele
  geom_density(data = main_analytic, aes(x = max_allele), adjust = 10, color = "dodgerblue3", linewidth = 0.7, alpha = 0.1) +
  theme_classic() +
  labs(
    title = "Density of Maximum FMR1 Repeat Length in Women",
    x = "Max Allele Length",
    y = "Density (Max Allele)"
  ) +
  scale_x_continuous(breaks = seq(10, 200, by=10))




ggplot() +
  # Plot density of max_allele with a legend
  geom_density(data = main_analytic, aes(x = max_allele, color = "Max Allele"), adjust = 10, linewidth = 0.7, alpha = 0.1) +
  # Plot density of min_allele with a legend
  geom_density(data = main_analytic, aes(x = min_allele, color = "Min Allele"), adjust = 10, linewidth = 0.7, alpha = 0.1) +
  theme_classic() +
  labs(
    x = "Allele Length",
    y = "Density"
  ) +
  scale_x_continuous(breaks = seq(10, 200, by = 10)) +
  scale_color_manual(values = c("Max Allele" = "dodgerblue3", "Min Allele" = "seagreen")) +
  theme(legend.title = element_blank())  # Removes the title for the legend


summary(main_analytic$min_allele)
summary(main_analytic$max_allele)




POI_cases <- main_analytic %>% 
  filter(nat_mage_POI==1)

POI_controls <- main_analytic %>% 
  filter(nat_mage_POI==0)


ggplot() +
  # Plot density of max_allele with a legend
  geom_density(data = POI_cases, aes(x = max_allele, color="POI"), adjust = 10, linewidth = 0.7, alpha = 0.1) +
  # Plot density of min_allele with a legend
  geom_density(data = POI_controls, aes(x = max_allele, color="Normal"), adjust = 10, linewidth = 0.7, alpha = 0.1) +
  theme_classic() +
  labs(
    x = "Allele Length",
    y = "Density"
  ) +
  scale_x_continuous(breaks = seq(10, 200, by = 10)) +
  scale_color_manual(values = c("Normal" = "dodgerblue3", "POI" = "darkred")) +
  theme(legend.title = element_blank())  # Removes the title for the legend



### splitting the x axis 



n_cases <- nrow(POI_cases)
n_controls <- nrow(POI_controls)
total_individuals <- n_cases + n_controls

cat("Number of POI cases:", n_cases, "\n")
cat("Number of controls:", n_controls, "\n")
cat("Total individuals represented:", total_individuals, "\n")


### adding this on to the plot as a label 

# Count individuals
n_cases <- nrow(POI_cases)
n_controls <- nrow(POI_controls)
total_individuals <- n_cases + n_controls
label_text <- paste("N =", total_individuals)




ggplot() +
  geom_density(data = POI_cases, aes(x = max_allele, color = "POI"), adjust = 10, linewidth = 0.7, alpha = 0.1) +
  geom_density(data = POI_controls, aes(x = max_allele, color = "Normal"), adjust = 10, linewidth = 0.7, alpha = 0.1) +
  theme_classic() +
  labs(x = "Allele Length", y = "Density") +
  scale_color_manual(values = c("Normal" = "dodgerblue3", "POI" = "darkred")) +
  scale_x_break(c(70, 160), scales = 0.3) +  # squashes instead of removing
  annotate("text", x = 170, y = 0.16, label = label_text, size = 4) +
  theme(legend.title = element_blank())


### making the plot a histogram instead

ggplot() +
  geom_histogram(data = POI_cases, aes(x = max_allele, fill = "POI"), 
                 bins = 50, position = "identity", alpha = 0.5, color = "darkred") +
  geom_histogram(data = POI_controls, aes(x = max_allele, fill = "Normal"), 
                 bins = 50, position = "identity", alpha = 0.5, color = "dodgerblue3") +
  theme_classic() +
  labs(x = "Allele Length", y = "Count") +
  scale_fill_manual(values = c("Normal" = "dodgerblue3", "POI" = "darkred")) +
  theme(legend.title = element_blank()) +
  scale_x_break(c(70, 160)) +
  annotate("text", x = 170, y = 15, label = label_text, size = 4)


# Combine both datasets and add group labels
POI_cases$Group <- "POI"
POI_controls$Group <- "Normal"
hist_data <- rbind(POI_cases, POI_controls)

ggplot(hist_data, aes(x = max_allele, fill = Group)) +
  geom_histogram(bins = 50, color = "black", alpha = 0.7) +
  facet_wrap(~ Group, scales = "free_y") +
  theme_classic() +
  labs(x = "Allele Length", y = "Count") +
  scale_fill_manual(values = c("Normal" = "dodgerblue3", "POI" = "darkred")) +
  theme(legend.position = "none") 


library(ggbreak)

ggplot(hist_data, aes(x = max_allele, fill = Group)) +
  geom_histogram(bins = 50, color = "black", alpha = 0.7) +
  facet_wrap(~ Group, scales = "free_y") +
  theme_classic() +
  labs(x = "Allele Length", y = "Count") +
  scale_fill_manual(values = c("Normal" = "dodgerblue3", "POI" = "darkred")) +
  theme(legend.position = "none")   # <-- Break between 70 and 120


library(ggplot2)

# Subset data by Group
data_POI <- subset(hist_data, Group == "POI")
data_Normal <- subset(hist_data, Group == "Normal")

# Plot for POI
p_POI <- ggplot(data_POI, aes(x = max_allele, fill = Group)) +
  geom_histogram(bins = 50, color = "black", alpha = 0.7) +
  theme_classic() +
  labs(title = "POI Group", x = "Allele Length", y = "Count") +
  scale_fill_manual(values = c("POI" = "darkred")) +
  theme(legend.position = "none")  +  scale_x_break(c(70, 120)) 

# Plot for Normal
p_Normal <- ggplot(data_Normal, aes(x = max_allele, fill = Group)) +
  geom_histogram(bins = 50, color = "black", alpha = 0.7) +
  theme_classic() +
  labs(title = "Normal Group", x = "Allele Length", y = "Count") +
  scale_fill_manual(values = c("Normal" = "dodgerblue3")) +
  theme(legend.position = "none") + scale_x_break(c(70, 160)) 

# Combine plots side by side
library(patchwork)
combined <- p_Normal + p_POI + plot_layout(ncol = 2)
print(combined)





table(main_analytic$AAM_all, useNA="ifany")

library(ggplot2)

ggplot() +
  # Plot density of menopause age
  geom_bar(data = main_analytic, aes(x = AAM_all), adjust = 10, color = "dodgerblue3", linewidth = 0.7, alpha = 0.1) +
  theme_classic() +
  labs(
    title = "Histogram of age at Menarche in Women",
    x = "Age at menarche",
    y = "Frequency"
  ) 





library(data.table)
irregular_periods <- fread("irregular_cycles_self_report_inc_contraception.tsv")

ggplot(irregular_periods, aes(x=factor(irregular_cycles))) +
  geom_bar()

main_analytic <- merge(main_analytic, irregular_periods, by="IID", all.x=TRUE)




#### REGRESSIONS WITH PROTEIN LEVELS



library(tidyverse)
protein <- fread("fmr1_protein_olink_instance_0.tsv")
protein <- protein %>% 
  rename(IID=`Participant ID(olink_instance_0 - eid)`)
protein <- protein %>% 
  rename(fmr1_protein=`FMR1;Synaptic functional regulator FMR1`)

main_analytic <- left_join(main_analytic, protein, by="IID")

sum(!is.na(main_analytic$fmr1_protein))


ggplot(main_analytic, aes(x=max_allele, y=fmr1_protein)) +
  geom_point() +
  theme_classic()


ggplot(main_analytic, aes(y=fmr1_protein)) +
  geom_density() +
  theme_classic()



has_mutation <- main_analytic %>% 
  filter(has_premut==1)

no_mutation <- main_analytic %>% 
  filter(has_premut==0)


ggplot() +
  geom_density(data = has_mutation, aes(x = fmr1_protein, color = "has_mutation"), adjust = 10, linewidth = 0.6, alpha = 0.1) +
  geom_density(data = no_mutation, aes(x = fmr1_protein, color = "no_mutation"), adjust = 10, linewidth = 0.6, alpha = 0.1) +
  scale_color_manual(name = "Legend", values = c("has_mutation" = "red", "no_mutation" = "dodgerblue")) +
  theme_minimal()





poi_cases <- main_analytic %>% 
  filter(nat_mage_POI==1)

poi_controls <- main_analytic %>% 
  filter(nat_mage_POI==0)




ggplot() +
  geom_density(data = poi_cases, aes(x = fmr1_protein, color = "poi_cases"), adjust = 10, linewidth = 0.6, alpha = 0.1) +
  geom_density(data = poi_controls, aes(x = fmr1_protein, color = "poi_controls"), adjust = 10, linewidth = 0.6, alpha = 0.1) +
  scale_color_manual(name = "Legend", values = c("poi_cases" = "red", "poi_controls" = "dodgerblue")) +
  labs(
    x = "FMR1 Protein Level",
    y = "Density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )


ggplot(main_analytic, aes(x = max_allele, y = fmr1_protein)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "#0072B2", fill = "#0072B2", alpha = 0.2) +
  labs(
    x = "Max allele (CGG repeats)",
    y = "FMR1 protein level"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )






#### SIMPLE LOGISTIC REGRESSION OF PREMUTATION WITH POI



#filter anyone that doesn't have menopause age

POI <- main_analytic  %>% 
  filter(nat_mage_min34!="NA")


### HAS PREMUTATION BINARY EXPOSURE

model_premut_POI <- glm(nat_mage_POI ~ has_premut + age + seqprovider + coverage +
                          p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = POI, family="binomial")



### MAX ALLELE CONTINUOUS EXPOSURE

model_max_POI <- glm(nat_mage_POI ~ max_allele + age + seqprovider + coverage + 
                       p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = POI, family="binomial")



### SUM OF ALLELES CONTINUOUS EXPOSURE 


model_sum_POI <- glm(nat_mage_POI ~ allele_sum + age + seqprovider + coverage +
                       p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = POI, family="binomial")






### BASIC LINEAR REGRESSION FOR AGE AT MENOPAUSE

### HAS PREMUT
model_premut_ANM <- lm(nat_mage_min34 ~ has_premut + age + seqprovider + coverage + 
                         p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)




### MAX ALLELE

model_max_ANM <- lm(nat_mage_min34 ~ max_allele + age + coverage + seqprovider +
                      p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)



### SUM ALLELES

model_sum_ANM <- lm(nat_mage_min34 ~ allele_sum + age + coverage + seqprovider +
                      p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)
summary(model_sum_ANM)








### ISOLATING THE AGE AT MENOPAUSE SAMPLE TO THOSE WITH MENOPAUSE ABOVE 40

POI_ANM_40 <- POI %>% 
  filter(nat_mage_min34>=40)

model_premut_ANM40 <- lm(nat_mage_min34 ~ has_premut + age + coverage + seqprovider +
                           p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = POI_ANM_40)




model_max_ANM40 <- lm(nat_mage_min34 ~ max_allele + age + coverage + seqprovider +
                        p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = POI_ANM_40)



### REPEATING FOR SUM OF ALLELES


model_sum_ANM40 <- lm(nat_mage_min34 ~ allele_sum + age + coverage + seqprovider +
                        p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = POI_ANM_40)





### BASIC LINEAR REGRESSION FOR AGE AT MENARCHE

### HAS PREMUT
model_premut_AAM <- lm(AAM_all ~ has_premut + age + seqprovider + coverage + 
                         p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)



### MAX ALLELE

model_max_AAM <- lm(AAM_all ~ max_allele + age + seqprovider + coverage + 
                      p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)




### SUM ALLELES

model_sum_AAM <- lm(AAM_all ~ allele_sum + age + seqprovider + coverage + 
                      p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)




### BASIC LINEAR REGRESSION FOR CYCLE LENGTH

### HAS PREMUT
model_premut_CYCLE <- lm(Cycle_length ~ has_premut + age + seqprovider + coverage + 
                           p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)




### MAX ALLELE

model_max_CYCLE <- lm(Cycle_length ~ max_allele + age + coverage + seqprovider +
                        p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)



### SUM ALLELES

model_sum_CYCLE <- lm(Cycle_length ~ allele_sum + age + coverage + seqprovider +
                        p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)




### number of children 

### HAS PREMUTATION BINARY EXPOSURE

main_analytic$live_births[main_analytic$live_births == -3] <- NA


model_premut_BIRTHS <- glm(live_births ~ has_premut + age + seqprovider + coverage +
                             p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic, family="poisson")


summary(model_premut_BIRTHS)

coef(summary(model_premut_BIRTHS))

### MAX ALLELE CONTINUOUS EXPOSURE

model_max_BIRTHS <- glm(live_births ~ max_allele + age + seqprovider + coverage + 
                          p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic, family="poisson")



### SUM OF ALLELES CONTINUOUS EXPOSURE 


model_sum_BIRTHS <- glm(live_births ~ allele_sum + age + seqprovider + coverage +
                          p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic, family="poisson")





### irregular periods 

### HAS PREMUTATION BINARY EXPOSURE


model_premut_IRREG <- glm(irregular_cycles ~ has_premut + age + seqprovider + coverage +
                            p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic, family="binomial")
summary(model_premut_IRREG)




### MAX ALLELE CONTINUOUS EXPOSURE

model_max_IRREG <- glm(irregular_cycles ~ max_allele + age + seqprovider + coverage + 
                         p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic, family="binomial")



### SUM OF ALLELES CONTINUOUS EXPOSURE 


model_sum_IRREG <- glm(irregular_cycles ~ allele_sum + age + seqprovider + coverage +
                         p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic, family="binomial")





### miscarriage


#### lost pregnancy 

model_premut_MISCARRY <- glm(lost_pregnancy ~ has_premut + age + seqprovider + coverage + 
                               p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic, family="binomial")
summary(model_premut_MISCARRY)


model_max_MISCARRY <- glm(lost_pregnancy ~ max_allele + age + seqprovider + coverage + 
                            p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic, family="binomial")



model_sum_MISCARRY <- glm(lost_pregnancy ~ allele_sum + age + seqprovider + coverage + 
                            p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic, family="binomial")





### HAS PREMUT
model_premut_PROT <- lm(fmr1_protein ~ has_premut + age + seqprovider + coverage + 
                          p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)



### MAX ALLELE

model_max_PROT <- lm(fmr1_protein ~ max_allele + age + seqprovider + coverage + 
                       p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)
summary(model_max_PROT)



### SUM ALLELES

model_sum_PROT <- lm(fmr1_protein ~ allele_sum + age + seqprovider + coverage + 
                       p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)





### extracting the results to get them into a table

# Load libraries
library(dplyr)

# Helper functions
get_or_ci_p <- function(model, var) {
  coefs <- summary(model)$coefficients
  beta <- coefs[var, "Estimate"]
  se <- coefs[var, "Std. Error"]
  or <- exp(beta)
  ci_low <- exp(beta - 1.96 * se)
  ci_high <- exp(beta + 1.96 * se)
  p <- coefs[var, "Pr(>|z|)"]
  return(c(sprintf("%.4f (%.4f–%.4f)", or, ci_low, ci_high), signif(p, 3)))
}

get_beta_ci_p <- function(model, var) {
  coefs <- summary(model)$coefficients
  beta <- coefs[var, "Estimate"]
  se <- coefs[var, "Std. Error"]
  ci_low <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  p <- coefs[var, "Pr(>|t|)"]
  return(c(sprintf("%.4f (%.4f–%.4f)", beta, ci_low, ci_high), signif(p, 3)))
}

get_beta_pois_ci_p <- function(model, var) {
  coefs <- summary(model)$coefficients
  beta <- coefs[var, "Estimate"]
  se <- coefs[var, "Std. Error"]
  ci_low <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  p <- coefs[var, "Pr(>|z|)"]
  return(c(sprintf("%.4f (%.4f–%.4f)", beta, ci_low, ci_high), signif(p, 3)))
}


# Initialize table
results <- data.frame(
  Exposure = c("Premutation", "Max allele", "Sum alleles"),
  POI_OR_CI = NA, POI_P = NA,
  ANM_Beta_CI = NA, ANM_P = NA,
  ANM40_Beta_CI = NA, ANM40_P = NA,
  AAM_Beta_CI = NA, AAM_P = NA,
  Cycle_Beta_CI = NA, Cycle_P = NA,
  LiveBirths_Beta_CI = NA, LiveBirths_P = NA,
  IRREG_OR_CI = NA, IRREG_P = NA,
  MC_OR_CI = NA, MC_P = NA,
  PROT_Beta_CI = NA, PROT_P = NA
)

### === Fill in results for each exposure ===

## Row 1: Premutation
results[1, c("POI_OR_CI", "POI_P")] <- get_or_ci_p(model_premut_POI, "has_premut")
results[1, c("ANM_Beta_CI", "ANM_P")] <- get_beta_ci_p(model_premut_ANM, "has_premut")
results[1, c("ANM40_Beta_CI", "ANM40_P")] <- get_beta_ci_p(model_premut_ANM40, "has_premut")
results[1, c("AAM_Beta_CI", "AAM_P")] <- get_beta_ci_p(model_premut_AAM, "has_premut")
results[1, c("Cycle_Beta_CI", "Cycle_P")] <- get_beta_ci_p(model_premut_CYCLE, "has_premut")
results[1, c("LiveBirths_Beta_CI", "LiveBirths_P")] <- get_beta_pois_ci_p(model_premut_BIRTHS, "has_premut")
results[1, c("IRREG_OR_CI", "IRREG_P")] <- get_or_ci_p(model_premut_IRREG, "has_premut")
results[1, c("MC_OR_CI", "MC_P")] <- get_or_ci_p(model_premut_MISCARRY, "has_premut")
results[1, c("PROT_Beta_CI", "PROT_P")] <- get_beta_ci_p(model_premut_PROT, "has_premut")

## Row 2: Max allele
results[2, c("POI_OR_CI", "POI_P")] <- get_or_ci_p(model_max_POI, "max_allele")
results[2, c("ANM_Beta_CI", "ANM_P")] <- get_beta_ci_p(model_max_ANM, "max_allele")
results[2, c("ANM40_Beta_CI", "ANM40_P")] <- get_beta_ci_p(model_max_ANM40, "max_allele")
results[2, c("AAM_Beta_CI", "AAM_P")] <- get_beta_ci_p(model_max_AAM, "max_allele")
results[2, c("Cycle_Beta_CI", "Cycle_P")] <- get_beta_ci_p(model_max_CYCLE, "max_allele")
results[2, c("LiveBirths_Beta_CI", "LiveBirths_P")] <- get_beta_pois_ci_p(model_max_BIRTHS, "max_allele")
results[2, c("IRREG_OR_CI", "IRREG_P")] <- get_or_ci_p(model_max_IRREG, "max_allele")
results[2, c("MC_OR_CI", "MC_P")] <- get_or_ci_p(model_max_MISCARRY, "max_allele")
results[2, c("PROT_Beta_CI", "PROT_P")] <- get_beta_ci_p(model_max_PROT, "max_allele")

## Row 3: Sum of alleles
results[3, c("POI_OR_CI", "POI_P")] <- get_or_ci_p(model_sum_POI, "allele_sum")
results[3, c("ANM_Beta_CI", "ANM_P")] <- get_beta_ci_p(model_sum_ANM, "allele_sum")
results[3, c("ANM40_Beta_CI", "ANM40_P")] <- get_beta_ci_p(model_sum_ANM40, "allele_sum")
results[3, c("AAM_Beta_CI", "AAM_P")] <- get_beta_ci_p(model_sum_AAM, "allele_sum")
results[3, c("Cycle_Beta_CI", "Cycle_P")] <- get_beta_ci_p(model_sum_CYCLE, "allele_sum")
results[3, c("LiveBirths_Beta_CI", "LiveBirths_P")] <- get_beta_pois_ci_p(model_sum_BIRTHS, "allele_sum")
results[3, c("IRREG_OR_CI", "IRREG_P")] <- get_or_ci_p(model_sum_IRREG, "allele_sum")
results[3, c("MC_OR_CI", "MC_P")] <- get_or_ci_p(model_sum_MISCARRY, "allele_sum")
results[3, c("PROT_Beta_CI", "PROT_P")] <- get_beta_ci_p(model_sum_PROT, "allele_sum")


## print and export
print(results)
write_tsv(results, "summary_results_table_regressions.tsv")




### gathering numbers of observations 

nobs(model_premut_POI)
nobs(model_max_POI)
nobs(model_sum_POI)

nobs(model_premut_ANM)
nobs(model_max_ANM)
nobs(model_sum_ANM)

nobs(model_premut_ANM40)
nobs(model_max_ANM40)
nobs(model_sum_ANM40)

nobs(model_premut_AAM)
nobs(model_max_AAM)
nobs(model_sum_AAM)

nobs(model_premut_CYCLE)
nobs(model_max_CYCLE)
nobs(model_sum_CYCLE)

nobs(model_premut_BIRTHS)
nobs(model_max_BIRTHS)
nobs(model_sum_BIRTHS)

nobs(model_premut_IRREG)
nobs(model_max_IRREG)
nobs(model_sum_IRREG)

nobs(model_premut_MISCARRY)
nobs(model_max_MISCARRY)
nobs(model_sum_MISCARRY)

nobs(model_premut_PROT)
nobs(model_max_PROT)
nobs(model_sum_PROT)






### plotting the mean allele length or the proportion of POI cases across the 
### menopause ages 


summary_df <- main_analytic %>%
  filter(!is.na(nat_mage_min34)) %>% 
  mutate(
    age_floor = floor(nat_mage_min34),
    age_group = if_else(age_floor >= 60, 60L, age_floor)  # Group 60+ as 60
  ) %>%
  group_by(age_group) %>%
  summarise(
    n = n(),
    carriers = sum(has_premut),
    proportion = carriers / n
  ) %>%
  mutate(
    binom_res = purrr::pmap(list(carriers, n), ~ binom::binom.confint(..1, ..2, method = "wilson")),
    lower = purrr::map_dbl(binom_res, ~ .x$lower),
    upper = purrr::map_dbl(binom_res, ~ .x$upper)
  )

summary_df <- summary_df %>%
  mutate(age_label = if_else(age_group == 60, "60+", as.character(age_group)))

summary_df$age_label <- as.numeric(summary_df$age_label)


summary_df$color_group <- cut(
  summary_df$age_label,
  breaks = c(-Inf, 40, 45, Inf),
  labels = c("POI", "Early Menopause", "Normal Menopause")
)

ggplot(summary_df, aes(x = age_label, y = proportion)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4) +
  labs(
    x = "Menopause Age",
    y = "Proportion with Premutation (95% CI)",
    title = "Proportion of Premutation Carriers by Menopause Age"
  ) +
  theme_classic()

ggplot(summary_df, aes(x = age_label, y = proportion, color = color_group)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4) +
  labs(
    x = "Menopause Age",
    y = "Proportion with Premutation (95% CI)",
    title = "Proportion of Premutation Carriers by Menopause Age",
    color = "Age Group"
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("POI" = "darkred", "Early Menopause" = "orange", "Normal Menopause" = "deepskyblue2"),
    na.translate = FALSE
  )


ggplot(summary_df, aes(x = age_label, y = proportion, color = color_group)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4) +
  geom_hline(yintercept = 0.006, linetype = "dashed", color = "black") +
  annotate("text", x = 60, y = 0.005, label = "Sample\npremutation\nprevalence=0.005", 
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "black") +
  labs(
    x = "Menopause Age",
    y = "Proportion with Premutation (95% CI)",
    title = "Proportion of Premutation Carriers by Menopause Age",
    color = "Age Group"
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("POI" = "darkred", "Early Menopause" = "orange", "Normal Menopause" = "deepskyblue2"),
    na.translate = FALSE
  )



library(ggplot2)
library(cowplot)

p <- ggplot(summary_df, aes(x = age_label, y = proportion, color = color_group)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4) +
  geom_hline(yintercept = 0.005, linetype = "dashed", color = "black") +
  labs(
    x = "Menopause Age",
    y = "Proportion with Premutation (95% CI)",
    color = "Age Group"
  ) +
  theme_classic() +
  scale_color_manual(
    values = c("POI" = "darkred", "Early Menopause" = "orange", "Normal Menopause" = "deepskyblue2"),
    na.translate = FALSE
  )

# Combine with label below
final_plot <- ggdraw() +
  draw_plot(p, x = 0, y = 0.1, width = 1, height = 0.9) +
  draw_label("Population\npremutation\nprevalence=0.005", 
             x = 0.96, y = 0.30, hjust = 1, size = 10)

final_plot






#### spline regressions at 36 and 55 repeats 

### manually setting the knots
library(splines)
library(broom)
spline_adj <- lm(nat_mage_min34 ~ ns(max_allele, knots=c(36, 55))
                 + age + seqprovider + coverage + 
                   p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 +
                   p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 +
                   p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 +
                   p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = main_analytic)
summary(spline_adj)


summary(spline_adj)

sum_mod <- summary(spline_adj)

# Extract coefficients table with p-values
coef_table <- sum_mod$coefficients

# Extract p-values (numeric)
p_values <- coef_table[, "Pr(>|t|)"]

# Display p-values
print(p_values)

model_results <- tidy(spline_adj, conf.int = TRUE)

# View the first few rows of the results
head(model_results)

# Filter for spline terms only (those that contain 'ns(max_allele)')
spline_results <- model_results %>%
  filter(grepl("ns\\(max_allele", term)) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(
    `β (95% CI)` = paste0(round(estimate, 3), " (", round(conf.low, 3), ", ", round(conf.high, 3), ")"),
    p = round(p.value, 3)
  ) %>%
  select(spline = term, `β (95% CI)`, p)

# View the filtered result table
print(spline_results)


ggplot(main_analytic, aes(x = max_allele, y = nat_mage_min34)) +
  geom_smooth(method = "lm", formula = y ~ ns(x, knots = c(36, 55)) , se = TRUE, color = "blue") +
  labs(title = "a)",
       x = "Maximum Allele Length",
       y = "Age at Menopause") +
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 180, by = 10))



library(splines)

# Fit model
fit <- lm(nat_mage_min34 ~ ns(max_allele, knots = c(36, 55)), data = main_analytic)

# Create a new dataset over the range of x
newdat <- data.frame(max_allele = seq(0, 180, by = 1))
newdat$pred <- predict(fit, newdata = newdat)

# Plot predicted curve
ggplot(main_analytic, aes(x = max_allele, y = nat_mage_min34)) +
  geom_line(data = newdat, aes(x = max_allele, y = pred), color = "blue", size = 1) +
  labs(
    title = "a)",
    x = "Maximum Allele Length",
    y = "Predicted Age at Menopause"
  ) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 180, by = 10))







spline_adj <- lm(
  AAM_all ~ ns(max_allele, knots=c(36,55)) + age + seqprovider + coverage + 
    p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 +
    p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 +
    p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 +
    p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20,
  data = main_analytic
)

summary(spline_adj)

model_results <- tidy(spline_adj, conf.int = TRUE)

plot(spline_adj)

# View the first few rows of the results
head(model_results)

print(model_results)

# Filter for spline terms only (those that contain 'ns(max_allele)')
spline_results <- model_results %>%
  filter(grepl("ns\\(max_allele", term)) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(
    `β (95% CI)` = paste0(round(estimate, 3), " (", round(conf.low, 3), ", ", round(conf.high, 3), ")"),
    p = round(p.value, 3)
  ) %>%
  select(spline = term, `β (95% CI)`, p)

# View the filtered result table
print(spline_results)



ggplot(main_analytic, aes(x = max_allele, y = AAM_all)) +
  geom_smooth(method = "lm", formula = y ~ ns(x, knots=c(36,55)), se = TRUE, color = "blue") +
  labs(title="b)",
       x = "Maximum Allele Length",
       y = "Age at Menarche") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 180, by = 10))





#### looking at the association between the levels of the protein and age at menopause


protein <- merge(POI, protein, by="IID", all.x=TRUE)
protein <- protein %>% 
  select(-fmr1_protein.y) %>% 
  rename(fmr1_protein=fmr1_protein.x)

lm <- lm(nat_mage_min34 ~ fmr1_protein + age, data = protein)

summary(lm)


ggplot(protein, aes(x=fmr1_protein, y=nat_mage_min34)) +
  geom_point() +
  geom_smooth(method = "lm")







#### looking at the protein levels in women with POI versus those without

POI_cases <- protein %>% 
  filter(nat_mage_POI==1)

controls <- protein %>% 
  filter(nat_mage_POI==0)


summary(POI_cases$fmr1_protein)
summary(controls$fmr1_protein)



ggplot() +
  # Plot density of max_allele with a legend
  geom_density(data = POI_cases, aes(x = fmr1_protein, color="POI"), adjust = 10, linewidth = 0.7, alpha = 0.1) +
  # Plot density of min_allele with a legend
  geom_density(data = controls, aes(x = fmr1_protein, color="Normal"), adjust = 10, linewidth = 0.7, alpha = 0.1) +
  theme_classic() +
  labs(
    x = "FMR1 protein",
    y = "Density"
  )  +
  scale_color_manual(values = c("Normal" = "dodgerblue3", "POI" = "darkred")) +
  theme(legend.title = element_blank())  # Removes the title for the legend




ggplot(protein, aes(x=factor(nat_mage_POI), y=fmr1_protein)) +
  geom_boxplot() +
  theme_classic()



protein_normal <- protein %>% 
  filter(max_allele<55)

ggplot(protein_normal, aes(x=max_allele, y=fmr1_protein)) +
  geom_point()





#### adjusting for the PRS when looking at the association between fmr1 and the 
#### level of the protein 
anm_score <- fread("anm_score_fmr1_info.txt")
anm_score <- anm_score %>% 
  rename(IID=V1)


anm_score <- anm_score %>% 
  select(IID, anmscore)

protein <- merge(protein, anm_score, by="IID", all.x=TRUE)

protein <- protein %>% 
  rename(anm_increasing_score=anmscore)


#### calculating menopause age decreasing score 

protein <- protein %>% 
  mutate(anm_decreasing_score = 580-anm_increasing_score)


lm <- lm(nat_mage_min34 ~ fmr1_protein + anm_decreasing_score + age + seqprovider + coverage + 
           p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data = protein)

summary(lm)









#### creating a dummy variable for premutation status based on different cut off
#### ppints and plotting the odds ratio of POI based on these varying 
#### definitions of cases and controls 


library(cowplot)

POI <- main_analytic %>%  filter(!is.na(nat_mage_min34))

pcs <- paste("p22009_a", 1:20, sep="")
POI$age <- as.numeric(POI$age)

# Add columns for the new grouped categories
POI <- POI %>%
  mutate(
    has_premut_75_plus = ifelse(max_allele >= 75, 1, 0)
  )

# Loop to create 'has_premut_*' columns for 25 to 84
for (i in 25:74) {
  POI <- POI %>%
    mutate(!!paste0("has_premut_", i) := ifelse(max_allele >= i, 1, 0))
}

has_premut_vars <- c(paste0("has_premut_", 25:74), "has_premut_75_plus")

# Initialize an empty list to store the results
results_list <- list()

# Run logistic regression for each variable
for (var in has_premut_vars) {
  
  formula_string <- paste("nat_mage_POI ~", var, "+ age + seqprovider + coverage + ", paste(pcs, collapse=" + "))
  
  logistic_model <- glm(as.formula(formula_string), data = POI, family = binomial())
  
  tidy_model <- tidy(logistic_model)
  tidy_model$variable <- var
  
  results_list[[var]] <- tidy_model
}

# Combine all results into a single dataframe
results_df <- bind_rows(results_list)

# Filter results to include only coefficients for 'has_premut_*'
results_df_betas <- results_df %>% 
  filter(str_starts(term, "has_premut_"))

# Extract numeric cutoff values, handling the custom groups
results_df_betas <- results_df_betas %>%
  mutate(pre_mutation_cut_off = case_when(
    variable == "has_premut_75_plus" ~ 75,  # Representing '75+' as 75
    TRUE ~ as.numeric(gsub("has_premut_", "", variable))
  ))

### **Adding case/control counts**
results_list <- list()

# Loop through each cutoff (including the grouped categories)
for (cutoff in c(25:74, 75)) {
  
  cutoff_var <- case_when(
    cutoff == 75 ~ "has_premut_75_plus",
    TRUE ~ paste0("has_premut_", cutoff)
  )
  
  contingency_table <- table(POI[[cutoff_var]], POI$nat_mage_POI)
  
  premut_POI_count <- ifelse("1" %in% rownames(contingency_table) && "1" %in% colnames(contingency_table), contingency_table["1", "1"], 0)
  premut_noPOI_count <- ifelse("1" %in% rownames(contingency_table) && "0" %in% colnames(contingency_table), contingency_table["1", "0"], 0)
  nopremut_POI_count <- ifelse("0" %in% rownames(contingency_table) && "1" %in% colnames(contingency_table), contingency_table["0", "1"], 0)
  nopremut_noPOI_count <- ifelse("0" %in% rownames(contingency_table) && "0" %in% colnames(contingency_table), contingency_table["0", "0"], 0)
  
  result_df <- data.frame(
    premut_POI = premut_POI_count,
    premut_noPOI = premut_noPOI_count,
    nopremut_POI = nopremut_POI_count,
    nopremut_noPOI = nopremut_noPOI_count,
    pre_mutation_cut_off = cutoff
  )
  
  results_list[[paste0("cutoff_", cutoff)]] <- result_df
}

# Combine all results into a single dataframe
final_results_df <- do.call(rbind, results_list)

### **Final results table**
results <- merge(results_df_betas, final_results_df, by = "pre_mutation_cut_off")

results <- results %>% 
  mutate(OR = exp(estimate),
         UCI = exp(estimate + (2 * std.error)),
         LCI = exp(estimate - (2 * std.error)))

write_tsv(results, "logistic_regression_POI_premut_25-75.tsv")



results$OR <- format(results$OR, scientific = FALSE)
results$UCI <- format(results$UCI, scientific = FALSE)
results$LCI <- format(results$LCI, scientific = FALSE)

results$OR <- as.numeric(results$OR)
results$UCI <- as.numeric(results$UCI)
results$LCI <- as.numeric(results$LCI)

## adding the scientifically defined POI classification to the table

results <- results %>% 
  mutate(clinical_category = case_when(
    pre_mutation_cut_off<44 ~ "Normal",
    pre_mutation_cut_off>43&pre_mutation_cut_off<55 ~ "Intermediate",
    pre_mutation_cut_off>54&pre_mutation_cut_off<200 ~ "Premutation",
    pre_mutation_cut_off>199 ~ "Full mutation"
  ))

results$clinical_category <- 
  factor(results$clinical_category,
         levels = c("Normal", "Intermediate", "Premutation", "Full mutation"))


ggplot(results, aes(x = pre_mutation_cut_off, y = OR)) +
  geom_point(aes(color = clinical_category)) + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI, width = 0.1, color = clinical_category)) + 
  theme_minimal() +
  labs(
    x = "Pre-mutation Cut-off",
    y = "POI OR",
    color = "Clinical category\nof repeat length",
    title = "Odds ratio of POI by Pre-mutation Cut-off"
  ) +
  scale_x_continuous(
    breaks = seq(25, 90, by = 5), 
    labels = as.character(seq(25, 90, by = 5))  
  ) + 
  scale_y_continuous(
    breaks = seq(0, 15, by = 1),
    minor_breaks = seq(0, 15, 0.5),
    labels = as.character(seq(0, 15, by = 1))
  ) +
  scale_color_manual(values=c("dodgerblue",
                              "darkorange",
                              "red4")) +
  theme(panel.background = element_rect(fill = NA, colour="gray2", linewidth = 1))

# Find the first significant point
sig_point <- results[(results$LCI > 1 | results$UCI < 1), ][1, ]

p1 <- ggplot(results, aes(x = pre_mutation_cut_off, y = OR)) +
  geom_point(aes(color = clinical_category)) + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI, width = 0.1, color = clinical_category)) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  geom_text(
    data = sig_point,
    aes(label = "*"),
    vjust = -1,
    color = "black",
    size = 3.5
  ) +
  theme_minimal() +
  labs(
    x = "Pre-mutation Cut-off",
    y = "POI OR",
    color = "Clinical category\nof repeat length",
    title = "Odds ratio of POI by Pre-mutation Cut-off"
  ) +
  scale_x_continuous(
    breaks = seq(25, 90, by = 5), 
    labels = as.character(seq(25, 90, by = 5))  
  ) + 
  scale_y_continuous(
    breaks = seq(0, 15, by = 1),
    minor_breaks = seq(0, 15, 0.5),
    labels = as.character(seq(0, 15, by = 1))
  ) +
  scale_color_manual(values = c("dodgerblue", "darkorange", "red4")) +
  theme(panel.background = element_rect(fill = NA, colour = "gray2", linewidth = 1))


### adding a zoomed in version of the area we are interested in 

library(ggplot2)
library(cowplot)

# Create main plot
main_plot <- ggplot(results, aes(x = pre_mutation_cut_off, y = OR)) +
  geom_point(aes(color = clinical_category)) + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI, width = 0.1, color = clinical_category)) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  theme_minimal() +
  labs(
    x = "Pre-mutation Cut-off",
    y = "OR of POI",
    color = "Clinical category\nof repeat length",
  ) +
  scale_x_continuous(
    breaks = seq(25, 90, by = 5), 
    labels = as.character(seq(25, 90, by = 5))  
  ) + 
  scale_y_continuous(
    breaks = seq(0, 15, by = 1),
    minor_breaks = seq(0, 15, 0.5),
    labels = as.character(seq(0, 15, by = 1))
  ) +
  scale_color_manual(values = c("dodgerblue", "darkorange", "red4")) +
  theme(panel.background = element_rect(fill = NA, colour = "gray2", linewidth = 1))

# Create zoomed-in plot
zoom_plot <- ggplot(results, aes(x = pre_mutation_cut_off, y = OR)) +
  geom_point(aes(color = clinical_category)) + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI, width = 0.1, color = clinical_category)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  coord_cartesian(xlim = c(32, 40), ylim = c(0.8, 1.5)) +
  scale_x_continuous(breaks = 32:40) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    plot.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
  ) +
  scale_color_manual(values = c("dodgerblue", "darkorange", "red4"))

# Combine plots with zoom inset
final_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(zoom_plot, x = 0.075, y = 0.45, width = 0.4, height = 0.4) +
  draw_line(x = c(0.18, 0.075), y = c(0.16, 0.45), color = "gray", linewidth = 0.5) +  # x ≈ 32
  draw_line(x = c(0.3, 0.475), y = c(0.16, 0.45), color = "gray", linewidth = 0.5) 



print(final_plot)










#### investigating the AUC curves for each of the logistic regression models
# Load necessary library for ROC plot
library(pROC)
library(dplyr)

complete_POI <- POI %>%
  dplyr::select(
    matches("^has_premut_(24|25|26|27|28|29|30|31|32|33|34|35|36|37|38|39|40|41|42|43|44|45|46|47|48|49|50|51|52|53|54|55|56|57|58|59|60|61|62|63|64|65|66|67|68|69|70|71|72|73|74|75)$"),
    nat_mage_POI, age, coverage, seqprovider, 
    p22009_a1:p22009_a20
  )



complete_POI <- na.omit(complete_POI)


# Create the formula string dynamically
formula_string <- paste("nat_mage_POI ~ has_premut_45 + age + coverage +", paste(pcs, collapse = " + "))

# Fit the logistic regression model
mylogit <- glm(as.formula(formula_string), data = complete_POI, family = binomial())

# Get the predicted probabilities
prob <- predict(mylogit,type = c("response"))

test_roc <- roc(complete_POI$nat_mage_POI ~ prob, plot=TRUE, print.auc=TRUE)






#### PUTTING THIS INTO A LOOP 



# Define the range of pre-mutation variables
premut_vars <- 45:67

# Create an empty list to store the ROC objects
roc_list <- list()

# Create an empty vector to store the AUC values
auc_values <- numeric(length(premut_vars))

# Loop through the pre-mutation variables from 45 to 65
for (i in premut_vars) {
  # Create the formula string dynamically for the current pre-mutation variable
  formula_string <- paste("nat_mage_POI ~ has_premut_", i, " + age + seqprovider + coverage + ", paste(pcs, collapse=" + "), sep = "")
  
  # Fit the logistic regression model
  mylogit <- glm(as.formula(formula_string), data = complete_POI, family = binomial())
  
  # Get the predicted probabilities
  prob <- predict(mylogit, type = "response")
  complete_POI$prob <- prob  # Store the probabilities in the data
  
  # Create the ROC curve for the current pre-mutation variable and store it in the list
  test_roc <- roc(complete_POI$nat_mage_POI ~ prob, plot = FALSE)  # plot = FALSE to avoid plotting
  roc_list[[paste("test_roc_", i, sep = "")]] <- test_roc
  
  # Calculate and store the AUC for the current model
  auc_values[i - 44] <- auc(test_roc)
}

# Print the AUC values for each model
print(auc_values)

# Plot the ROC curves
plot(roc_list[["test_roc_45"]], col = 1, lty = 2, main = "ROC", xlim = c(1, 0), ylim = c(0, 1))

# Add the other ROC curves for test_roc_46, test_roc_47, ..., test_roc_65
for (i in 46:67) {
  # Plot each subsequent ROC curve with different color and line type
  plot(roc_list[[paste("test_roc_", i, sep = "")]], col = i - 45, lty = 1 + (i %% 3), add = TRUE)
}

# Optionally, add a legend to label each ROC curve
legend("bottomleft", legend = paste(45:67, "AUC =", round(auc_values, 3)), col = 1:21, lty = 1 + (1:21 %% 3), cex = 0.65)






### extracting the sensitivity and specificity for all the models

# Create the formula string 
formula_string <- paste("nat_mage_POI ~ has_premut_45 + age + coverage +", paste(pcs, collapse = " + "))

# Fit the logistic regression model
mylogit <- glm(as.formula(formula_string), data = complete_POI, family = binomial())

# Get the predicted probabilities
prob <- predict(mylogit,type = c("response"))

test_roc <- roc(complete_POI$nat_mage_POI ~ prob, plot=TRUE, print.auc=TRUE)

coords(test_roc, x = 0.54, input = "threshold", ret = c("sensitivity", "specificity"))

summary(test_roc)

best_coords <- coords(test_roc, "best", ret = c("threshold", "sensitivity", "specificity"))
print(best_coords)


formula_string <- paste("nat_mage_POI ~ has_premut_55 + age + coverage +", paste(pcs, collapse = " + "))

# Fit the logistic regression model
mylogit <- glm(as.formula(formula_string), data = complete_POI, family = binomial())

# Get the predicted probabilities
prob <- predict(mylogit,type = c("response"))

test_roc <- roc(complete_POI$nat_mage_POI ~ prob, plot=TRUE, print.auc=TRUE)

coords(test_roc, x = 0.54, input = "threshold", ret = c("sensitivity", "specificity"))

summary(test_roc)

best_coords <- coords(test_roc, "best", ret = c("threshold", "sensitivity", "specificity"))
print(best_coords)



#### putting it in a loop

library(dplyr)
library(pROC)

# Define the premutation columns
premut_columns <- paste0("has_premut_", 25:74)

# Loop through each premutation column
for (premut in premut_columns) {
  
  # Create the formula string dynamically
  formula_string <- paste("nat_mage_POI ~", premut, "+ age + coverage +", paste(pcs, collapse = " + "))
  
  # Fit the logistic regression model
  mylogit <- glm(as.formula(formula_string), data = complete_POI, family = binomial())
  
  # Get the predicted probabilities
  prob <- predict(mylogit, type = "response")
  
  # Create the ROC curve
  test_roc <- roc(complete_POI$nat_mage_POI ~ prob, plot = TRUE, print.auc = TRUE)
  
  # Get the best threshold based on Youden's index
  best_coords <- coords(test_roc, "best", ret = c("threshold", "sensitivity", "specificity"))
  print(best_coords)
  
  # Optionally, print the summary of the ROC object for each loop iteration
  summary(test_roc)
}




library(dplyr)
library(pROC)

# premutation columns
premut_columns <- paste0("has_premut_", 25:74)

results <- data.frame(
  premut_column = character(),
  threshold = numeric(),
  sensitivity = numeric(),
  specificity = numeric(),
  auc = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each premutation column
for (premut in premut_columns) {
  
  # Create the formula string dynamically
  formula_string <- paste("nat_mage_POI ~", premut, "+ age + coverage +", paste(pcs, collapse = " + "))
  
  # Fit the logistic regression model
  mylogit <- glm(as.formula(formula_string), data = complete_POI, family = binomial())
  
  # Get the predicted probabilities
  prob <- predict(mylogit, type = "response")
  
  # Create the ROC curve
  test_roc <- roc(complete_POI$nat_mage_POI ~ prob, plot = FALSE)
  
  
  # Get the best threshold based on Youden's index
  best_coords <- coords(test_roc, "best", ret = c("threshold", "sensitivity", "specificity"))
  
  # Get the AUC value
  auc_value <- auc(test_roc)
  
  # Append results to the table
  results <- rbind(results, data.frame(
    premut_column = premut,
    threshold = best_coords$threshold,
    sensitivity = best_coords$sensitivity,
    specificity = best_coords$specificity,
    auc = auc_value
  ))
}

# View the results table
print(results)

write_tsv(results, "sensitivty_specificity.tsv")










#### repeating this but changing the POI threshold rather than the premutation threshold


# Set POI thresholds
POI_thresholds <- 34:60

# Create POI indicator variables
for (threshold in POI_thresholds) {
  col_name <- paste0("POI_less_than_", threshold)
  POI[[col_name]] <- ifelse(POI$nat_mage_min34 <= threshold, 1,
                            ifelse(!is.na(POI$nat_mage_min34), 0, NA))
}


library(broom)

# Covariates
pcs <- paste0("p22009_a", 1:20)

# Initialize list for storing results
regression_results <- list()

# Loop through each threshold-specific POI column
for (threshold in POI_thresholds) {
  
  outcome_var <- paste0("POI_less_than_", threshold)
  
  formula_str <- paste0(outcome_var, " ~ has_premut + age + seqprovider + coverage + ", 
                        paste(pcs, collapse = " + "))
  
  model <- try(glm(as.formula(formula_str), data = POI, family = binomial()), silent = TRUE)
  
  if (inherits(model, "try-error") || !model$converged) {
    warning(paste("Model did not converge or failed for threshold", threshold))
    next
  }
  
  tidy_model <- tidy(model)
  tidy_model$threshold <- threshold
  regression_results[[as.character(threshold)]] <- tidy_model
}

# Combine and filter for has_premut coefficient
results_df <- bind_rows(regression_results) %>%
  filter(term == "has_premut") %>%
  mutate(
    OR = exp(estimate),
    LCI = exp(estimate - 2 * std.error),
    UCI = exp(estimate + 2 * std.error)
  )


write_tsv(results_df, "POI_threshold_sweep.tsv")

ggplot(results_df, aes(x = threshold, y = OR)) +
  geom_point() +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  scale_x_continuous(breaks = 34:60) +
  theme_classic() +
  labs(
    x = "POI defined as menopause before age...",
    y = "Odds Ratio (Premutation vs. Normal)",
    title = "Effect of Premutation on POI across Age Thresholds"
  )












##### now doing what we did with the premutation but changing the reference group to 
##### everyone with 30 repeats

# Updated logistic regression loop comparing to repeat length 30
results_list <- list()

for (i in c(30, 30:74, 75)) {
  
  var_name <- paste0("compare_to_30_", i)
  
  POI <- POI %>%
    mutate(!!sym(var_name) := case_when(
      max_allele == 30 ~ 0,
      max_allele >= i ~ 1,
      TRUE ~ NA_real_
    ))
  
  # Drop NAs so we only include relevant groups (exactly 30 vs. ≥ i)
  subset_data <- POI %>% filter(!is.na(!!sym(var_name)))
  
  formula_string <- paste("nat_mage_POI ~", var_name, "+ age + seqprovider + coverage +", paste(pcs, collapse = " + "))
  
  logistic_model <- glm(as.formula(formula_string), data = subset_data, family = binomial())
  
  tidy_model <- tidy(logistic_model)
  tidy_model$variable <- var_name
  tidy_model$cutoff <- i
  
  results_list[[var_name]] <- tidy_model
}

# Combine all regression results
results_df <- bind_rows(results_list)

# Filter to keep only the rows for the cutoff variable
results_df <- results_df %>%
  filter(term == variable) %>%
  mutate(
    OR = exp(estimate),
    LCI = exp(estimate - 1.96 * std.error),
    UCI = exp(estimate + 1.96 * std.error)
  )


library(ggplot2)

ggplot(results_df, aes(x = cutoff, y = OR)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  theme_minimal() +
  labs(
    x = "CGG Repeat Cutoff (≥ vs. 30)",
    y = "Odds Ratio for POI",
    title = "Odds of POI by CGG Repeat Threshold (Compared to 30 Repeats)"
  ) +
  scale_y_continuous(trans = "log10",  # optional: log scale if large OR spread
                     breaks = c(0.1, 0.5, 1, 2, 5, 10),
                     labels = c("0.1", "0.5", "1", "2", "5", "10")) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

results_df <- results_df %>%
  mutate(clinical_category = case_when(
    cutoff < 44 ~ "Normal",
    cutoff >= 44 & cutoff <= 54 ~ "Intermediate",
    cutoff > 54 ~ "Premutation"
  ))

results_df$clinical_category <- factor(results_df$clinical_category,
                                       levels = c("Normal", "Intermediate", "Premutation"))


ggplot(results_df, aes(x = cutoff, y = OR, color = clinical_category)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  theme_minimal() +
  labs(
    x = "CGG Repeat Cutoff (≥ vs. 30)",
    y = "Odds Ratio for POI",
    title = "Odds of POI by CGG Repeat Threshold (Compared to 30 Repeats)",
    color = "Clinical Category"
  ) +
  scale_color_manual(values = c("Normal" = "dodgerblue",
                                "Intermediate" = "darkorange",
                                "Premutation" = "red4")) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

sum(main_analytic$max_allele==30)

#### based on log odds ratio and with all the relevant information 

# Main plot
main_plot <- ggplot(results_df, aes(x = cutoff, y = log(OR), color = clinical_category)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = log(LCI), ymax = log(UCI)), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal() +
  labs(
    x = "CGG Repeat Cutoff (≥ vs. 30)",
    y = "Log Odds Ratio for POI",
    color = "Clinical Category"
  ) +
  scale_color_manual(values = c(
    "Normal" = "dodgerblue",
    "Intermediate" = "darkorange",
    "Premutation" = "red4"
  )) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5)) +
  scale_x_continuous(breaks = seq(30, 75, by = 5)) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

main_plot <- ggplot(results_df, aes(x = cutoff, y = log(OR), color = clinical_category)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = log(LCI), ymax = log(UCI)), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal() +
  labs(
    x = "CGG Repeat Cutoff (≥ vs. 30)",
    y = "Odds Ratio for POI",
    color = "Clinical Category"
  ) +
  scale_color_manual(values = c(
    "Normal" = "dodgerblue",
    "Intermediate" = "darkorange",
    "Premutation" = "red4"
  )) +
  scale_y_continuous(
    breaks = log(c(0.25, 0.5, 1, 2, 4, 8)),
    labels = c(0.25, 0.5, 1, 2, 4, 8)
  ) +
  scale_x_continuous(breaks = seq(30, 75, by = 5)) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )



asterisk_y <- results_df %>%
  filter(cutoff == 36) %>%
  pull(UCI) %>%
  log() + 0.02  # small offset above the CI

# Zoomed-in plot
# Add the asterisk annotation
zoom_plot <- ggplot(results_df, aes(x = cutoff, y = log(OR))) +
  geom_point(aes(color = clinical_category)) + 
  geom_errorbar(aes(ymin = log(LCI), ymax = log(UCI), color = clinical_category), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_text(data = results_df %>% filter(cutoff == 36),
            aes(x = 36, y = asterisk_y), label = "* p<0.05", size = 2) +
  coord_cartesian(xlim = c(32, 40), ylim = c(0, 0.4)) +
  scale_x_continuous(breaks = 32:40) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    plot.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
  ) +
  scale_color_manual(values = c(
    "Normal" = "dodgerblue",
    "Intermediate" = "darkorange",
    "Premutation" = "red4"
  ))


final_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(zoom_plot, x = 0.08, y = 0.6, width = 0.35, height = 0.35) +
  draw_line(x = c(0.12, 0.078), y = c(0.16, 0.60), color = "gray30", linewidth = 0.2, linetype=2) +
  draw_line(x = c(0.25, 0.432), y = c(0.16, 0.60), color = "gray30", linewidth = 0.2, linetype=2)

# Print the final combined plot
print(final_plot)


write_tsv(results_df, "newcomparing_with_30_results.tsv")













#### calculating the roc curves with 30 as the reference group


library(pROC)

roc_list <- list()
auc_list <- data.frame(cutoff = numeric(), auc = numeric())

for (i in c(30, 30:74, 75)) {
  
  var_name <- paste0("compare_to_30_", i)
  
  POI <- POI %>%
    mutate(!!sym(var_name) := case_when(
      max_allele == 30 ~ 0,
      max_allele >= i ~ 1,
      TRUE ~ NA_real_
    ))
  
  subset_data <- POI %>% filter(!is.na(!!sym(var_name)))
  
  formula_string <- paste("nat_mage_POI ~", var_name, "+ age + seqprovider + coverage +", paste(pcs, collapse = " + "))
  
  logistic_model <- glm(as.formula(formula_string), data = subset_data, family = binomial())
  
  tidy_model <- broom::tidy(logistic_model)
  tidy_model$variable <- var_name
  tidy_model$cutoff <- i
  results_list[[var_name]] <- tidy_model
  
  # Predict probabilities
  subset_data$prob <- predict(logistic_model, type = "response")
  
  # ROC and AUC
  roc_obj <- roc(subset_data$nat_mage_POI, subset_data$prob)
  roc_list[[var_name]] <- roc_obj
  auc_list <- rbind(auc_list, data.frame(cutoff = i, auc = auc(roc_obj)))
}


ggplot(auc_list, aes(x = cutoff, y = auc)) +
  geom_line(color = "black") +
  geom_point(color = "blue", size = 2) +
  theme_minimal() +
  labs(
    x = "CGG Repeat Cutoff (≥ vs. 30)",
    y = "AUC (Discrimination of POI)",
    title = "AUC of Models Comparing Repeat ≥ Cutoff vs 30"
  ) +
  ylim(0.5, 1) +
  theme(
    panel.border = element_rect(color = "gray30", fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )


# Extract AUCs
auc_45 <- round(auc(roc_list[["compare_to_30_45"]]), 3)
auc_55 <- round(auc(roc_list[["compare_to_30_55"]]), 3)
auc_65 <- round(auc(roc_list[["compare_to_30_65"]]), 3)

# Plot ROC curves
plot(roc_list[["compare_to_30_45"]], col = "red", lty = 1)
plot(roc_list[["compare_to_30_55"]], add = TRUE, col = "blue", lty = 2)
plot(roc_list[["compare_to_30_65"]], add = TRUE, col = "darkgreen", lty = 3)

# Add AUCs to legend
legend("bottomright", 
       legend = c(paste0("≥45 (AUC = ", auc_45, ")"),
                  paste0("≥55 (AUC = ", auc_55, ")"),
                  paste0("≥65 (AUC = ", auc_65, ")")),
       col = c("red", "blue", "darkgreen"), lty = 1:3)



# Compute AUC
auc_55 <- round(auc(roc_list[["compare_to_30_55"]]), 3)

# Set up empty plot with axis labels
plot(
  roc_list[["compare_to_30_55"]],
  col = "blue",
  lwd = 2.5,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  legacy.axes = TRUE,
  print.auc = FALSE
)

# Add diagonal line for reference
segments(x0 = 1, y0 = 0, x1 = 0, y1 = 1, lty = 2, col = "gray50")

# Add legend with AUC
legend("bottomright",
       legend = paste0("CGG ≥55 vs. 30 (AUC = ", auc_55, ")"),
       col = "blue",
       lwd = 4,
       bty = "n")




write_tsv(auc_list, "auc_by_repeat_threshold.tsv")








### repeating this for the sum of alleles with the reference group as 60 repeats

results_list_allele_sum <- list()

# Loop starting at 60 instead of 30
for (i in c(60, 60:125, 125)) {
  
  var_name <- paste0("compare_to_60_", i)
  
  POI <- POI %>%
    mutate(!!sym(var_name) := case_when(
      allele_sum == 60 ~ 0,
      allele_sum >= i ~ 1,
      TRUE ~ NA_real_
    ))
  
  # Filter to relevant groups (exactly 60 vs. ≥ i)
  subset_data <- POI %>% filter(!is.na(!!sym(var_name)))
  
  formula_string <- paste("nat_mage_POI ~", var_name, "+ age + seqprovider + coverage +", paste(pcs, collapse = " + "))
  
  logistic_model <- glm(as.formula(formula_string), data = subset_data, family = binomial())
  
  tidy_model <- tidy(logistic_model)
  tidy_model$variable <- var_name
  tidy_model$cutoff <- i
  
  results_list_allele_sum[[var_name]] <- tidy_model
}

# Combine all regression results
results_df_allele_sum <- bind_rows(results_list_allele_sum)

# Filter rows for the cutoff variable
results_df_allele_sum <- results_df_allele_sum %>%
  filter(term == variable) %>%
  mutate(
    OR = exp(estimate),
    LCI = exp(estimate - 1.96 * std.error),
    UCI = exp(estimate + 1.96 * std.error)
  )

# Define clinical categories based on cutoff for allele_sum (adjust if needed)
results_df_allele_sum <- results_df_allele_sum %>%
  mutate(clinical_category = case_when(
    cutoff < 90 ~ "Normal",        # Adjust these ranges as per allele_sum biology or analysis
    cutoff >= 90 & cutoff <= 110 ~ "Intermediate",
    cutoff > 110 ~ "Premutation"
  ))

results_df_allele_sum$clinical_category <- factor(results_df_allele_sum$clinical_category,
                                                  levels = c("Normal", "Intermediate", "Premutation"))

# Then you can plot with the same plotting code, just replace results_df with results_df_allele_sum and adjust labels accordingly

ggplot(results_df_allele_sum, aes(x = cutoff, y = OR, color = clinical_category)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  theme_minimal() +
  labs(
    x = "Allele Sum Cutoff (≥ vs. 60)",
    y = "Odds Ratio for POI",
    title = "Odds of POI by Allele Sum Threshold (Compared to 60)",
    color = "Clinical Category"
  ) +
  scale_color_manual(values = c("Normal" = "dodgerblue",
                                "Intermediate" = "darkorange",
                                "Premutation" = "red4")) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

fwrite(results_df_allele_sum, "allele_sum_vs_60.tsv", sep = "\t")









#### repeating this process but with age at menopause as the outcome

ANM <- main_analytic 

## filtering for women who have not had POI (menopause before 40)



# Loop to create 'has_premut_*' columns
for (i in 25:100) {
  ANM <- ANM %>%
    mutate(!!paste0("has_premut_", i) := ifelse(max_allele >= i, 1, 0))
}


has_premut_vars <- paste0("has_premut_", 25:100)


# Initialize an empty list to store the results
results_list <- list()

# Loop over each variable and run the logistic regression
for (var in has_premut_vars) {
  
  formula_string <- paste("nat_mage_min34~", var, "+ age + seqprovider + coverage +", paste(pcs, collapse=" + "))
  
  # Fit the linear regression model
  linear_model <- lm(as.formula(formula_string), data = ANM)
  
  # Tidy the model output and add the variable name
  tidy_model <- tidy(linear_model)
  tidy_model$variable <- var
  
  # Append the results to the list
  results_list[[var]] <- tidy_model
}

# Combine all results into a single data frame
results_df <- bind_rows(results_list)

# Filter the results to only include coefficients for nat_mage_POI
results_df_betas <- results_df %>% 
  filter(str_starts(term, "has_premut_"))

# Create a new column for the pre-mutation cutoff values
results_df_betas <- results_df_betas %>% 
  mutate(pre_mutation_cut_off = as.numeric(gsub("has_premut_", "", variable)))

# Display the results
print(results_df_betas)


write_tsv(results_df_betas, "linear_regression_ANM_premut.tsv")

```
```{r}
ANM <- main_analytic 

# Create 'has_premut_<25' group
ANM <- ANM %>%
  mutate(has_premut_25 = ifelse(max_allele < 25, 1, 0))

# Create binary indicators from 25 to 75
for (i in 26:75) {
  ANM <- ANM %>%
    mutate(!!paste0("has_premut_", i) := ifelse(max_allele >= i, 1, 0))
}

# Create 'has_premut_>75' group
ANM <- ANM %>%
  mutate(has_premut_76 = ifelse(max_allele > 75, 1, 0))

# Combine variable names
has_premut_vars <- c("has_premut_25", paste0("has_premut_", 26:75), "has_premut_76")

# Run regression loop
results_list <- list()

for (var in has_premut_vars) {
  formula_string <- paste("nat_mage_min34 ~", var, "+ age + seqprovider + coverage +", paste(pcs, collapse = " + "))
  
  linear_model <- lm(as.formula(formula_string), data = ANM)
  tidy_model <- tidy(linear_model)
  tidy_model$variable <- var
  results_list[[var]] <- tidy_model
}

# Combine and filter results
results_df <- bind_rows(results_list)

results_df_betas <- results_df %>% 
  filter(str_starts(term, "has_premut_")) %>%
  mutate(pre_mutation_cut_off = case_when(
    variable == "has_premut_25" ~ "25",
    variable == "has_premut_76" ~ "76",
    TRUE ~ gsub("has_premut_", "", variable)
  ))

# Save
write_tsv(results_df_betas, "linear_regression_ANM_premut_grouped.tsv")


results_df_betas$pre_mutation_cut_off <- as.numeric(results_df_betas$pre_mutation_cut_off)

results_df_betas <- results_df_betas %>% 
  mutate(clinical_category = case_when(
    pre_mutation_cut_off < 44 ~ "Normal",
    (pre_mutation_cut_off >= 44 & pre_mutation_cut_off < 55) ~ "Intermediate",
    (pre_mutation_cut_off >= 55 & pre_mutation_cut_off < 200) ~ "Premutation",
    pre_mutation_cut_off == ">75" ~ "Premutation"
  ))


results_df_betas$clinical_category <- 
  factor(results_df_betas$clinical_category,
         levels = c("Normal", "Intermediate", "Premutation"))


ggplot(results_df_betas, aes(x = pre_mutation_cut_off, y = estimate, colour = clinical_category)) +
  geom_point() + 
  geom_errorbar(aes(ymin = estimate - 2*std.error, ymax = estimate + 2*std.error), width = 0.1) + 
  theme_minimal() +
  labs(
    x = "Pre-mutation Cut-off",
    y = "Change in ANM (years)",
    color = "Clinical category\nof repeat length",
    title = "ANM by Pre-mutation Cut-off"
  ) +
  ylim(-4, 1) +
  scale_color_manual(values=c("dodgerblue",
                              "darkorange",
                              "red4")) +
  theme(panel.background = element_rect(fill = NA, colour="gray2", linewidth = 1),
        panel.grid.minor.x = element_blank())








##### looking at the same thing but only in women who do not have POI


ANM <- main_analytic 

ANM <- ANM %>% filter(nat_mage_min34>39)

# Create 'has_premut_<25' group
ANM <- ANM %>%
  mutate(has_premut_25 = ifelse(max_allele < 25, 1, 0))

# Create binary indicators from 25 to 75
for (i in 26:75) {
  ANM <- ANM %>%
    mutate(!!paste0("has_premut_", i) := ifelse(max_allele >= i, 1, 0))
}

# Create 'has_premut_>75' group
ANM <- ANM %>%
  mutate(has_premut_76 = ifelse(max_allele > 75, 1, 0))

# Combine variable names
has_premut_vars <- c("has_premut_25", paste0("has_premut_", 26:75), "has_premut_76")

# Run regression loop
results_list <- list()

for (var in has_premut_vars) {
  formula_string <- paste("nat_mage_min34 ~", var, "+ age + seqprovider + coverage +", paste(pcs, collapse = " + "))
  
  linear_model <- lm(as.formula(formula_string), data = ANM)
  tidy_model <- tidy(linear_model)
  tidy_model$variable <- var
  results_list[[var]] <- tidy_model
}

# Combine and filter results
results_df <- bind_rows(results_list)

results_df_betas <- results_df %>% 
  filter(str_starts(term, "has_premut_")) %>%
  mutate(pre_mutation_cut_off = case_when(
    variable == "has_premut_25" ~ "25",
    variable == "has_premut_76" ~ "76",
    TRUE ~ gsub("has_premut_", "", variable)
  ))

# Save
write_tsv(results_df_betas, "linear_regression_ANM_premut_grouped_non_POI.tsv")

results_df_betas$pre_mutation_cut_off <- as.numeric(results_df_betas$pre_mutation_cut_off)

results_df_betas <- results_df_betas %>% 
  mutate(clinical_category = case_when(
    pre_mutation_cut_off < 44 ~ "Normal",
    (pre_mutation_cut_off >= 44 & pre_mutation_cut_off < 55) ~ "Intermediate",
    (pre_mutation_cut_off >= 55 & pre_mutation_cut_off < 200) ~ "Premutation",
    pre_mutation_cut_off == ">75" ~ "Premutation"
  ))


results_df_betas$clinical_category <- 
  factor(results_df_betas$clinical_category,
         levels = c("Normal", "Intermediate", "Premutation"))


ggplot(results_df_betas, aes(x = pre_mutation_cut_off, y = estimate, colour = clinical_category)) +
  geom_point() + 
  geom_errorbar(aes(ymin = estimate - 2*std.error, ymax = estimate + 2*std.error), width = 0.1) + 
  theme_minimal() +
  labs(
    x = "Pre-mutation Cut-off",
    y = "Change in ANM (years)",
    color = "Clinical category\nof repeat length",
    title = "ANM by Pre-mutation Cut-off"
  ) +
  ylim(-4, 1) +
  scale_color_manual(values=c("dodgerblue",
                              "darkorange",
                              "red4")) +
  theme(panel.background = element_rect(fill = NA, colour="gray2", linewidth = 1),
        panel.grid.minor.x = element_blank())






#### doing the same thing but with age at menarche as the outcome



AAM <- main_analytic 



# Loop to create 'has_premut_*' columns
for (i in 30:90) {
  AAM <- AAM %>%
    mutate(!!paste0("has_premut_", i) := ifelse(max_allele >= i, 1, 0))
}


has_premut_vars <- paste0("has_premut_", 30:90)


# Initialize an empty list to store the results
results_list <- list()

# Loop over each variable and run the logistic regression
for (var in has_premut_vars) {
  # Construct the formula
  formula_string <- paste("AAM_all ~", var, "+ age + seqprovider + coverage +", paste(pcs, collapse=" + "))
  
  # Fit the linear regression model
  linear_model <- lm(as.formula(formula_string), data = AAM)
  
  
  # Tidy the model output and add the variable name
  tidy_model <- tidy(linear_model)
  tidy_model$variable <- var
  
  # Append the results to the list
  results_list[[var]] <- tidy_model
}

# Combine all results into a single data frame
results_df <- bind_rows(results_list)

# Filter the results to only include coefficients for nat_mage_POI
results_df_betas <- results_df %>% 
  filter(str_starts(term, "has_premut_"))

# Create a new column for the pre-mutation cutoff values
results_df_betas <- results_df_betas %>% 
  mutate(pre_mutation_cut_off = as.numeric(gsub("has_premut_", "", variable)))

# Display the results
print(results_df_betas)


write_tsv(results_df_betas, "linear_regression_AAM_premut.tsv")

```


### Plot the results 

```{r}

results_df_betas <- results_df_betas %>% 
  mutate(clinical_category = case_when(
    pre_mutation_cut_off<44 ~ "Normal",
    pre_mutation_cut_off>43&pre_mutation_cut_off<55 ~ "Intermediate",
    pre_mutation_cut_off>54&pre_mutation_cut_off<200 ~ "Premutation",
    pre_mutation_cut_off>199 ~ "Full mutation"
  ))

results_df_betas$clinical_category <- 
  factor(results_df_betas$clinical_category,
         levels = c("Normal", "Intermediate", "Premutation", "Full mutation"))


ggplot(results_df_betas, aes(x = pre_mutation_cut_off, y = estimate, colour = clinical_category)) +
  geom_point() + 
  geom_errorbar(aes(ymin = estimate - 2*std.error, ymax = estimate + 2*std.error), width = 0.1) + 
  theme_minimal() +
  labs(
    x = "Pre-mutation Cut-off",
    y = "Change in AAM (years)",
    color = "Clinical category\nof repeat length",
    title = "AAM by Pre-mutation Cut-off"
  ) +
  scale_x_continuous(
    breaks = seq(25, 100, by = 5),  # Set the breaks at every integer between 45 and 69
    labels = as.character(seq(25, 100, by = 5))  # Ensure the labels are characters
  ) +
  ylim(-2,1) +
  scale_color_manual(values=c("dodgerblue",
                              "darkorange",
                              "red4")) +
  theme(panel.background = element_rect(fill = NA, colour="gray2", linewidth = 1),
        panel.grid.minor.x = element_blank())













#### investigation of the interaction with the PRS
## We want to model the interaction between the PRS and the effect of the FMR1 mutation
## We can also see if any specific significant SNPs in the PRS have a mediating effect 



## reading in the ANM scores 

anm_score <- fread("anm_score_fmr1_info.txt")
anm_score <- anm_score %>% 
  rename(IID=V1)


anm_score <- anm_score %>% 
  dplyr::select(IID, anmscore)

main_analytic <- merge(main_analytic, anm_score, by="IID", all.x=TRUE)

main_analytic <- main_analytic %>% 
  rename(anm_increasing_score=anmscore)

### 

ggplot(anm_score, aes(x=anmscore)) +
  geom_density()



### plotting the effect of the score on ANM

main_analytic$age_bin <- floor(main_analytic$nat_mage_min34)
main_analytic$genetic_score_scaled <- scale(main_analytic$anm_increasing_score)
main_analytic$genetic_score_scaled <- as.numeric(scale(main_analytic$anm_increasing_score))
main_analytic$age_bin <- floor(pmin(main_analytic$nat_mage_min34, 60))



mean_scores <- main_analytic %>%
  group_by(age_bin) %>%
  summarise(
    mean_genetic_score = mean(genetic_score_scaled, na.rm = TRUE),
    sd_genetic_score = sd(genetic_score_scaled, na.rm = TRUE),
    n = n()
  )


library(ggplot2)
ggplot(mean_scores, aes(x = age_bin, y = mean_genetic_score)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 2) +
  labs(
    x = "Age at Menopause (binned)",
    y = "Mean Genetic Score",
    title = "Mean Genetic Score by Age at Menopause"
  ) +
  theme_minimal()





#### having a look at this by premutation status 

#### doing this by premutation status


# Cap age and create bins
main_analytic$age_bin <- floor(pmin(main_analytic$nat_mage_min34, 60))

# Standardize genetic score
main_analytic$genetic_score_scaled <- as.numeric(scale(main_analytic$anm_increasing_score))

# Group by both age bin and premutation status
mean_scores <- main_analytic %>%
  group_by(age_bin, has_premut) %>%
  summarise(
    mean_genetic_score = mean(genetic_score_scaled, na.rm = TRUE),
    sd_genetic_score = sd(genetic_score_scaled, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Plot
library(ggplot2)

ggplot(mean_scores, aes(x = age_bin, y = mean_genetic_score, color = factor(has_premut))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Age at Menopause (binned, max 60)",
    y = "Standardized Genetic Score (mean = 0, SD = 1)",
    color = "Premutation Status",
    title = "Mean Genetic Score by Age at Menopause, Split by Premutation Status"
  ) +
  theme_classic()









#### calculating menopause age decreasing score 

main_analytic <- main_analytic %>% 
  mutate(anm_decreasing_score = 580-anm_increasing_score)


ggplot(main_analytic, aes(x=anm_decreasing_score)) +
  geom_density()

ggplot(main_analytic, aes(x=anm_decreasing_score, y=nat_mage_min34)) +
  geom_point()


##### looking at it descriptively 

# Create the plot with different line colors for each group
ggplot(main_analytic, aes(x = anm_decreasing_score, color = factor(has_premut))) +
  geom_density(size = 1) +  # Adjust the line thickness with size argument
  scale_color_manual(values = c("blue", "red"), labels = c("No Mutation", "Has Mutation")) +
  labs(title = "Density of ANM Decreasing Score by Mutation Status",
       x = "ANM Decreasing Score",
       y = "Density",
       color = "Mutation Status") +
  theme_minimal()


main_analytic_mutation <- main_analytic %>% 
  filter(max_allele>54)

main_analytic_no_mutation <- main_analytic %>% 
  filter(max_allele<55)

ggplot(main_analytic_mutation, aes(x=anm_decreasing_score, y=nat_mage_min34)) +
  geom_point()

ggplot(main_analytic_no_mutation, aes(x=anm_decreasing_score, y=nat_mage_min34)) +
  geom_point()





### running logistic regression of PRS on POI in premutation and non-premutation carriers


# Fit models
glm_mutation <-glm(nat_mage_POI ~ anm_decreasing_score + age + seqprovider + coverage + arraybatch +
                     p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + 
                     p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + 
                     p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + 
                     p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, 
                   data = main_analytic_mutation, family = "binomial")


summary(glm_mutation)

nobs(glm_mutation)



glm_no_mutation <- glm(nat_mage_POI ~ anm_decreasing_score + age + seqprovider + coverage + arraybatch +
                         p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + 
                         p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + 
                         p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + 
                         p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, 
                       data = main_analytic_no_mutation, family = "binomial")


summary(glm_no_mutation)

nobs(glm_no_mutation)

p_values <- summary(glm_no_mutation)$coefficients[, 4]  # Extract the fourth column (Pr(>|z|))
print(p_values)




### now doing this for ANM
lm_no_mutation <- lm(nat_mage_min34 ~ anm_decreasing_score + age + seqprovider + coverage + arraybatch +
                       p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + 
                       p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + 
                       p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + 
                       p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, 
                     data = main_analytic_no_mutation)



summary(lm_no_mutation)


lm_mutation <- lm(nat_mage_min34 ~ anm_decreasing_score + age + seqprovider + coverage + arraybatch +
                    p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + 
                    p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + 
                    p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + 
                    p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, 
                  data = main_analytic_mutation)



summary(lm_mutation)







### adding PRS as an interaction term

glm_interaction <- glm(nat_mage_POI ~ has_premut*anm_decreasing_score + age + seqprovider + coverage + arraybatch + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data=main_analytic, family="binomial")

summary(glm_interaction)





# Extract results and filter for "anm_decreasing_score"
results_mutation <- tidy(glm_mutation, conf.int = TRUE) %>%
  filter(term == "anm_decreasing_score") %>%
  mutate(Group = "Mutation")

results_no_mutation <- tidy(glm_no_mutation, conf.int = TRUE) %>%
  filter(term == "anm_decreasing_score") %>%
  mutate(Group = "No Mutation")

# Combine results
results <- bind_rows(results_mutation, results_no_mutation)





ggplot(results, aes(x = Group, y = estimate, ymin = conf.low, ymax = conf.high, color = Group)) +
  geom_pointrange(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "ANM Decreasing Score Effect on POI in mutation\nversus non mutation carriers\nusing 290 genome-wide significant snps",
       x = "Group",
       y = "Log Odds Estimate",
       color = "Group") +
  theme_minimal()




#### DOING THE SAME BUT FOR ANM


### just adjusting for PRS

lm <- lm(nat_mage_min34 ~ has_premut + anm_decreasing_score + age + seqprovider + coverage + arraybatch +
           p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20 , data=main_analytic)
summary(lm)






### plotting the potential interaction between premutation and PRS on ANM


ggplot(main_analytic_mutation, aes(x=anm_decreasing_score, y=nat_mage_min34)) +
  geom_point()


lm <- lm(nat_mage_min34 ~ anm_decreasing_score + age + seqprovider + coverage + arraybatch + p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data=main_analytic_mutation)

summary(lm)




# Fit the model
lm <- lm(nat_mage_min34 ~ anm_decreasing_score + age + arraybatch + coverage + 
           p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + 
           p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + 
           p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + 
           p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, 
         data = main_analytic_no_mutation)

summary(lm)

# Given values
t_val <- -95.040
df <- 91840

# Calculate two-sided p-value
p_val <- 2 * pt(abs(t_val), df=df, lower.tail=FALSE)
p_val

# Fit the model
lm <- lm(nat_mage_min34 ~ anm_decreasing_score + age + arraybatch + coverage + 
           p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + 
           p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + 
           p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + 
           p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, 
         data = main_analytic_mutation)

# View the model summary
summary(lm)

### adding PRS as an interaction term

lm_interaction <- lm(nat_mage_min34 ~ has_premut*anm_decreasing_score + age + seqprovider + coverage  + arraybatch + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data=main_analytic)

summary(lm_interaction)



### without the dots 

# Create the scatter plot with regression lines for both groups (mutation and no mutation)
ggplot() +
  
  # Regression line for the no mutation group (line in red)
  geom_smooth(data = main_analytic_no_mutation, aes(x = anm_decreasing_score, y = nat_mage_min34, color = "No Mutation"), 
              method = "lm", se = FALSE) +
  
  # Regression line for the mutation group (line in green)
  geom_smooth(data = main_analytic_mutation, aes(x = anm_decreasing_score, y = nat_mage_min34, color = "Mutation"), 
              method = "lm", se = FALSE) + 
  
  # Add titles and labels
  labs(title = "Effect of ANM PRS on age at natural menopause based on\npremutation status",
       x = "ANM Decreasing Score (290 significant SNPs)",
       y = "Age at natural menopause", 
       color = "Group") +  # Label for the legend
  
  # Minimal theme for a cleaner look
  theme_minimal()





ggplot() +
  
  # Regression line for the no mutation group with CI (red)
  geom_smooth(data = main_analytic_no_mutation, aes(x = anm_decreasing_score, y = nat_mage_min34, color = "No Premutation"), 
              method = "lm", se = TRUE, color = "skyblue", fill = "skyblue", alpha = 0.2) +
  
  # Regression line for the mutation group with CI (green)
  geom_smooth(data = main_analytic_mutation, aes(x = anm_decreasing_score, y = nat_mage_min34, color = "Premutation"), 
              method = "lm", se = TRUE, color = "red", fill = "red", alpha = 0.2) + 
  
  # Titles and labels
  labs(x = "ANM Decreasing Score (290 significant SNPs)",
       y = "Age at natural menopause", 
       color = "Group") +  
  
  # Minimal theme
  theme_minimal() 






ggplot() +
  
  geom_smooth(data = main_analytic_no_mutation,
              aes(x = anm_decreasing_score, y = nat_mage_min34, color = "No Premutation", fill = "No Premutation"),
              method = "lm", se = TRUE, alpha = 0.2) +
  
  geom_smooth(data = main_analytic_mutation,
              aes(x = anm_decreasing_score, y = nat_mage_min34, color = "Premutation", fill = "Premutation"),
              method = "lm", se = TRUE, alpha = 0.2) +
  
  labs(x = "ANM Decreasing Score (290 significant SNPs)",
       y = "Age at natural menopause",
       color = "Group",
       fill = "Group") +
  
  scale_color_manual(values = c("No Premutation" = "skyblue", "Premutation" = "red")) +
  scale_fill_manual(values = c("No Premutation" = "skyblue", "Premutation" = "red")) +
  
  theme_minimal()








##### looping through all the SNPs individually 
anm_score <- fread("anm_ruth2021_variants_ukb.txt")

snp_data <- anm_score


### flipping the alleles to the menopause decreasing alignment 

library(dplyr)

# Assuming the first column is "IID" and the rest are SNPs
snp_data <- snp_data %>%
  mutate(across(-IID, ~ 2 - .))




main_analytic <- left_join(main_analytic, snp_data, by="IID")




# Identify SNP columns (those starting with "rs")
snp_cols <- grep("^(rs|chr)", names(main_analytic), value = TRUE)

# Create an empty list to store results
results_list <- list()

# Loop through each SNP and run the logistic regression
for (snp in snp_cols) {
  formula <- as.formula(paste("nat_mage_POI ~ has_premut *", snp, 
                              "+ age + seqprovider + coverage + arraybatch +", 
                              paste(paste0("p22009_a", 2:20), collapse = " + ")))
  
  model <- glm(formula, data = main_analytic, family = "binomial")
  
  # Extract the interaction term result
  interaction_term <- paste("has_premut:", snp, sep = "")
  results_list[[snp]] <- tidy(model) %>% filter(term == interaction_term)
}

# Combine results into a single dataframe
final_results_poi <- bind_rows(results_list, .id = "SNP")

# View results
head(final_results)

write_tsv(final_results_poi, "snp:premut_interactions.tsv")

final_results_sorted <- final_results %>%
  arrange(p.value)









#### looking at them with ANM

# Identify SNP columns (those starting with "rs" or "chr")
snp_cols <- grep("^(rs|chr)", names(main_analytic), value = TRUE)

# Create an empty list to store results
results_list <- list()

# Loop through each SNP and run the logistic regression
for (snp in snp_cols) {
  formula <- as.formula(paste("nat_mage_min34 ~ has_premut *", snp, 
                              "+ age + seqprovider + coverage + arraybatch +", 
                              paste(paste0("p22009_a", 2:20), collapse = " + ")))
  
  model <- lm(formula, data = main_analytic)
  
  # Extract the interaction term result
  interaction_term <- paste("has_premut:", snp, sep = "")
  results_list[[snp]] <- tidy(model) %>% filter(term == interaction_term)
}

# Combine results into a single dataframe
final_results_anm <- bind_rows(results_list, .id = "SNP")

# View results
head(final_results)

final_results_sorted <- final_results %>% 
  arrange(p.value)

write_tsv(final_results_anm, "snp:premut_interactions_anm.tsv")





### merging all those results into one table

# Already created:
final_results_poi   # From POI logistic models
final_results_anm   # From ANM linear models

final_results_poi <- final_results_poi %>%
  select(SNP, estimate, std.error, p.value) %>%
  rename(
    estimate_POI = estimate,
    std.error_POI = std.error,
    p.value_POI = p.value
  )

final_results_anm <- final_results_anm %>%
  select(SNP, estimate, std.error, p.value) %>%
  rename(
    estimate_ANM = estimate,
    std.error_ANM = std.error,
    p.value_ANM = p.value
  )

combined_results <- left_join(final_results_poi, final_results_anm, by = "SNP")

combined_results <- combined_results %>%
  arrange(p.value_ANM)

combined_results <- combined_results %>% 
  select(SNP, estimate_ANM, std.error_ANM, p.value_ANM, estimate_POI, std.error_POI, p.value_POI)

combined_results <- combined_results %>% 
  mutate(conf_int_ANM=paste0("(",round(estimate_ANM-1.96*std.error_ANM, 2),"-",round(estimate_ANM+1.96*std.error_ANM, 2),")"),
         conf_int_POI=paste0("(",round(estimate_POI-1.96*std.error_POI, 2),"-",round(estimate_POI+1.96*std.error_POI,2),")"))


combined_results <- combined_results %>% 
  select(SNP, estimate_ANM, conf_int_ANM, p.value_ANM, estimate_POI, conf_int_POI, p.value_POI)

combined_results <- combined_results %>% 
  mutate(estimate_ANM_confint=paste0(round(estimate_ANM, 2), " ", conf_int_ANM),
         estimate_POI_confint=paste0(round(estimate_POI, 2), " ", conf_int_POI))

combined_results <- combined_results %>% 
  select(SNP, estimate_ANM_confint, p.value_ANM, estimate_POI_confint, p.value_POI)

combined_results <- combined_results %>%
  mutate(
    p.value_ANM = round(p.value_ANM, 6),
    p.value_POI = round(p.value_POI, 6)
  )


write_tsv(combined_results, "combined_snp_premut_interaction_results.tsv")



combined_results_filtered <- combined_results %>% 
  filter(p.value_ANM<0.07 & p.value_POI<0.07)

write_tsv(combined_results_filtered, "filtered_combined_snp_premut_interaction_results.tsv")






#### checking if there is an interaction with the full ANM score released by UKBB - with all the variants in the 
#### genome and not just the genome wide significant ones


## reading in the ANM scores 

anm_score <- fread("ANM_PRS_participant.tsv")
anm_score <- anm_score %>% 
  rename(IID=`Participant ID`)


main_analytic <- fread("fmr1_main_analytic.txt")


main_analytic_meno <- main_analytic %>% 
  filter(nat_mage_min34!="NA")

main_analytic <- left_join(main_analytic, ethnicity, by="IID")


### filter to keep only unrelated individuals and also remove those who have pulled out

include <- fread("INCLUDEFOR_EUR_Unrelated.txt")

include <- include %>% 
  rename(IID=V1)


main_analytic_meno <- merge(include, main_analytic_meno, by="IID")


main_analytic <- main_analytic_meno


main_analytic$seqprovider <- as.factor(main_analytic$seqprovider)
main_analytic$seqprovider <- as.character(main_analytic$seqprovider)



main_analytic <- merge(main_analytic, anm_score, by="IID", all.x=TRUE)


ggplot(main_analytic, aes(x=`Standard PRS for age at menopause (AAM)`, y=nat_mage_min34)) +
  geom_point()

main_analytic <- main_analytic %>% 
  rename(anm_increasing_score=`Standard PRS for age at menopause (AAM)`)

### 

ggplot(main_analytic, aes(x=anm_increasing_score, y=nat_mage_min34)) +
  geom_point()


#### calculating menopause age decreasing score 

main_analytic <- main_analytic %>% 
  mutate(anm_decreasing_score = -(anm_increasing_score))


ggplot(main_analytic, aes(x=anm_decreasing_score)) +
  geom_density()

ggplot(main_analytic, aes(x=anm_decreasing_score, y=nat_mage_min34)) +
  geom_point()


##### looking at it descriptively 

# Create the plot with different line colors for each group
ggplot(main_analytic, aes(x = anm_decreasing_score, color = factor(has_premut))) +
  geom_density(size = 1) +  # Adjust the line thickness with size argument
  scale_color_manual(values = c("blue", "red"), labels = c("No Mutation", "Has Mutation")) +
  labs(title = "Density of ANM Decreasing Score by Mutation Status",
       x = "ANM Decreasing Score",
       y = "Density",
       color = "Mutation Status") +
  theme_minimal()


main_analytic_mutation <- main_analytic %>% 
  filter(max_allele>54)

main_analytic_no_mutation <- main_analytic %>% 
  filter(max_allele<55)

ggplot(main_analytic_mutation, aes(x=anm_decreasing_score, y=nat_mage_min34)) +
  geom_point()

ggplot(main_analytic_no_mutation, aes(x=anm_decreasing_score, y=nat_mage_min34)) +
  geom_point()





#### running the PRS alone 

glm <- glm(nat_mage_POI ~ anm_decreasing_score + age + arraybatch + coverage + p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data=main_analytic, family="binomial")

summary(glm)

p_values <- summary(glm)$coefficients[, 4]  # Extract the fourth column (Pr(>|z|))
print(p_values)



lm <- lm(nat_mage_min34 ~ anm_decreasing_score + age + arraybatch + coverage + p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data=main_analytic)

summary(lm)

# Extract p-values from the glm summary
p_values <- summary(lm)$coefficients[2, 4]


# Print formatted p-values
print(p_values)




### running logistic regression of PRS on POI in premutation and non-premutation carriers


glm <- glm(nat_mage_POI ~ anm_decreasing_score + age + seqprovider + coverage +
             p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20 , data=main_analytic_mutation, family="binomial")
summary(glm)



glm <- glm(nat_mage_POI ~ anm_decreasing_score + age + seqprovider + coverage +
             p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20 , data=main_analytic_no_mutation, family="binomial")
summary(glm)





library(ggplot2)
library(dplyr)
library(broom)

# Fit models
glm_mutation <- glm(nat_mage_POI ~ anm_decreasing_score + age + seqprovider + coverage +
                      p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + 
                      p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + 
                      p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + 
                      p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, 
                    data = main_analytic_mutation, family = "binomial")


summary(glm_mutation)

glm_no_mutation <- glm(nat_mage_POI ~ anm_decreasing_score + age + seqprovider + coverage +
                         p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + 
                         p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + 
                         p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + 
                         p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, 
                       data = main_analytic_no_mutation, family = "binomial")


summary(glm_no_mutation)

# Extract results and filter for "anm_decreasing_score"
results_mutation <- tidy(glm_mutation, conf.int = TRUE) %>%
  filter(term == "anm_decreasing_score") %>%
  mutate(Group = "Mutation")

results_no_mutation <- tidy(glm_no_mutation, conf.int = TRUE) %>%
  filter(term == "anm_decreasing_score") %>%
  mutate(Group = "No Mutation")

# Combine results
results <- bind_rows(results_mutation, results_no_mutation)


ggplot(results, aes(x = Group, y = estimate, ymin = conf.low, ymax = conf.high, color = Group)) +
  geom_pointrange(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "ANM Decreasing Score Effect on POI in mutation\nversus non mutation carriers\nusing whole genome score",
       x = "Group",
       y = "Log Odds Estimate",
       color = "Group") +
  theme_minimal()




### just adjusting for PRS

glm <- glm(nat_mage_POI ~ has_premut + anm_decreasing_score + age + seqprovider + coverage +
             p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20 , data=main_analytic, family="binomial")
summary(glm)


### adding PRS as an interaction term

glm_interaction <- glm(nat_mage_POI ~ anm_decreasing_score*has_premut + age + seqprovider + coverage + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data=main_analytic, family="binomial")

summary(glm_interaction)






#### DOING THE SAME BUT FOR ANM


### just adjusting for PRS

lm <- lm(nat_mage_min34 ~ has_premut + anm_decreasing_score + age + seqprovider + coverage +
           p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20 , data=main_analytic)
summary(lm)


### adding PRS as an interaction term

lm_interaction <- lm(nat_mage_min34 ~ anm_decreasing_score*has_premut + age + seqprovider + coverage + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data=main_analytic)

summary(lm_interaction)












### plotting the potential interaction between premutation and PRS on ANM


ggplot(main_analytic_mutation, aes(x=anm_decreasing_score, y=nat_mage_min34)) +
  geom_point()


lm <- lm(nat_mage_min34 ~ anm_decreasing_score + age + seqprovider + coverage + p22009_a1 +  p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, data=main_analytic_mutation)

summary(lm)




# Fit the model
lm <- lm(nat_mage_min34 ~ anm_decreasing_score + age + arraybatch + coverage + 
           p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + 
           p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 + 
           p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 + 
           p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20, 
         data = main_analytic_no_mutation)

# View the model summary
summary(lm)


p_values <- summary(lm)$coefficients[, 4]  # Extract the fourth column (Pr(>|z|))
print(p_values)


# Create the scatter plot with regression lines for both groups (mutation and no mutation)
ggplot() +
  # Scatter plot for the mutation group (points in black)
  geom_point(data = main_analytic_mutation, aes(x = anm_decreasing_score, y = nat_mage_min34), color = "black") + 
  
  
  # Scatter plot for the no mutation group (points in black)
  geom_point(data = main_analytic_no_mutation, aes(x = anm_decreasing_score, y = nat_mage_min34), color = "black") + 
  
  # Regression line for the no mutation group (line in red)
  geom_smooth(data = main_analytic_no_mutation, aes(x = anm_decreasing_score, y = nat_mage_min34, color = "No Mutation"), 
              method = "lm1", se = FALSE) +
  
  # Regression line for the mutation group (line in green)
  geom_smooth(data = main_analytic_mutation, aes(x = anm_decreasing_score, y = nat_mage_min34, color = "Mutation"), 
              method = "lm2", se = FALSE) + 
  
  # Add titles and labels
  labs(title = "Regression of nat_mage_min34 on ANM Decreasing Score\nusing whole genome score",
       x = "ANM Decreasing Score",
       y = "nat_mage_min34", 
       color = "Group") +  # Label for the legend
  
  # Minimal theme for a cleaner look
  theme_minimal()





### without the black dots

# Create the scatter plot with regression lines for both groups (mutation and no mutation)
ggplot() +
  
  
  # Regression line for the no mutation group (line in red)
  geom_smooth(data = main_analytic_no_mutation, aes(x = anm_decreasing_score, y = nat_mage_min34, color = "No Mutation"), 
              method = "lm", se = TRUE) +
  
  # Regression line for the mutation group (line in green)
  geom_smooth(data = main_analytic_mutation, aes(x = anm_decreasing_score, y = nat_mage_min34, color = "Mutation"), 
              method = "lm", se = TRUE) + 
  
  # Add titles and labels
  labs(title = "Regression of age at menopause on ANM Decreasing Score\nusing whole genome score",
       x = "ANM Decreasing Score",
       y = "Age at menopause", 
       color = "Group") +  # Label for the legend
  
  # Minimal theme for a cleaner look
  theme_minimal()







ggplot() +
  # Regression + CI for no mutation
  geom_smooth(data = main_analytic_no_mutation,
              aes(x = anm_decreasing_score, y = nat_mage_min34, color = "No Mutation", fill = "No Mutation"),
              method = "lm", se = TRUE) +
  
  # Regression + CI for mutation
  geom_smooth(data = main_analytic_mutation,
              aes(x = anm_decreasing_score, y = nat_mage_min34, color = "Mutation", fill = "Mutation"),
              method = "lm", se = TRUE) +
  
  # Titles and legend
  labs(
    x = "ANM Decreasing Score",
    y = "Age at Menopause",
    color = "Group",
    fill = "Group"  # Add fill legend label
  ) +
  theme_minimal() +
  
  # Custom colors
  scale_color_manual(values = c("No Mutation" = "skyblue", "Mutation" = "red3")) +
  scale_fill_manual(values = c("No Mutation" = "skyblue", "Mutation" = "palered")) +
  theme(
    legend.title=element_text(size=16),
    legend.text = element_text(size=14),
    legend.key.size = unit(0.7, "cm")
  )












##### investigating environmental exposures



library(ggeffects)
library(ggplot2)

main_analytic$BMI_cat <- cut(main_analytic$BMI,
                             breaks = c(18.5, 25, 30, Inf),
                             labels = c("Normal", "Overweight", "Obese"))

table(main_analytic$BMI_cat, main_analytic$has_premut)

bmi_model <- glm(nat_mage_POI ~ BMI_cat + has_premut + age + seqprovider + coverage + 
                   p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 +
                   p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 +
                   p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 +
                   p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20,
                 data = main_analytic, family = "binomial")

summary(bmi_model)

model_binned <- glm(nat_mage_POI ~ has_premut * BMI_cat + age + seqprovider + coverage + 
                      p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 +
                      p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 +
                      p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 +
                      p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20,
                    data = main_analytic, family = "binomial")




summary(model_binned)


preds_binned <- ggpredict(model_binned, terms = c("BMI_cat", "has_premut"))

ggplot(preds_binned, aes(x = x, y = predicted, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.6), width = 0.2) +
  labs(
    x = "BMI category",
    y = "Predicted Probability of POI",
    fill = "FMR1 Premutation"
  ) +
  theme_minimal() +
  scale_fill_manual(values =c("0" = "skyblue", "1" = "darkred")) +
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "bottom")








smoking <- fread("smoking_participant.tsv")
smoking <- smoking %>% 
  rename(IID=eid,
         smoking_status=p20116_i0,
         smoking_pack_years=p20161_i0)

main_analytic <- left_join(main_analytic, smoking, by="IID")

table(main_analytic$smoking_status, main_analytic$has_premut)

main_analytic <- main_analytic %>% 
  mutate(smoker=ifelse(smoking_status=="Current" | smoking_status=="Previous", 1,0))

main_analytic <- main_analytic %>% 
  filter(smoking_status=="Never"|  smoking_status=="Previous"| smoking_status=="Current")


main_analytic$smoking_status <- factor(main_analytic$smoking_status, levels = c("Never", "Previous", "Current"))
table(main_analytic$smoking_status)



model_binned <- glm(nat_mage_POI ~ smoking_status + has_premut + age + seqprovider + coverage +  p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 +
                      p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 +
                      p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 +
                      p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20,
                    data = main_analytic, family = "binomial")


summary <- summary(model_binned)

summary(model_binned)

p_values <- coef(summary)[, "Pr(>|z|)"]

model_binned <- glm(nat_mage_POI ~ has_premut * smoking_status + age + seqprovider + coverage +  p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 +
                      p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10 +
                      p22009_a11 + p22009_a12 + p22009_a13 + p22009_a14 + p22009_a15 +
                      p22009_a16 + p22009_a17 + p22009_a18 + p22009_a19 + p22009_a20,
                    data = main_analytic, family = "binomial")


summary(model_binned)


model_binned$fitted.values


preds_binned <- ggpredict(model_binned, terms = c("smoking_status", "has_premut"))

ggplot(preds_binned, aes(x = x, y = predicted, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.6), width = 0.2) +
  labs(
    x = "Smoking status",
    y = "Predicted Probability of POI",
    fill = "FMR1 Premutation"
  ) +
  theme_minimal() +
  scale_fill_manual(values =c("0" = "skyblue", "1" = "darkred")) +
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "bottom")







prs_sex <- fread("anm_prs_with_sex_participant.tsv")

prs_sex <- prs_sex %>% 
  rename(PRS=`Standard PRS for age at menopause (AAM)`)



lm <- lm(PRS~Sex, data=prs_sex)
summary(lm)

prs_male <- prs_sex %>%  
  filter(Sex=="Male")

prs_female <- prs_sex %>%  
  filter(Sex=="Female")


ggplot() +
  geom_density(data = prs_male, aes(x = PRS, color = "prs_male"), adjust = 10, linewidth = 0.6, alpha = 0.1) +
  geom_density(data = prs_female, aes(x = PRS, color = "prs_female"), adjust = 10, linewidth = 0.6, alpha = 0.1) +
  scale_color_manual(name = "Legend", values = c("prs_female" = "red", "prs_male" = "dodgerblue")) +
  theme_minimal() +
  labs(
    x = "ANM PRS",
    y = "Density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )


t.test(prs_male$PRS, prs_female$PRS)



prs_sex$Sex <- factor(prs_sex$Sex, levels = c("Male", "Female"))

# Fit the linear regression model
model <- lm(PRS ~ Sex, data = prs_sex)

# View the model summary
summary(model)




# Fit the linear regression model
prs_sex$Sex <- factor(prs_sex$Sex, levels = c("Male", "Female"))
model <- lm(PRS ~ Sex, data = prs_sex)

# Extract the regression results (coefficients)
model_summary <- summary(model)
coefficients <- model_summary$coefficients
intercept <- coefficients["(Intercept)", "Estimate"]
sex_female_coeff <- coefficients["SexFemale", "Estimate"]
p_value <- model_summary$coefficients["SexFemale", "Pr(>|t|)"]

# Create the plot with added text annotation for the model results

ggplot() +
  geom_density(data = prs_male, aes(x = PRS, color = "prs_male"), adjust = 10, linewidth = 0.6, alpha = 0.1) +
  geom_density(data = prs_female, aes(x = PRS, color = "prs_female"), adjust = 10, linewidth = 0.6, alpha = 0.1) +
  scale_color_manual(name = "Legend", values = c("prs_female" = "red", "prs_male" = "dodgerblue")) +
  theme_minimal() +
  labs(
    x = "ANM PRS",
    y = "Density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  # Add the model results as text annotation
  annotate("text", x = 3, y = 0.25, 
           label = paste( "β Female : ", "0.002", "\n",
                          "P-value: ", "0.328"),
           size = 4, color = "black", hjust = 0.5)

































