#########################
### Uro LFMM analysis ###
### per reviewer recs. ##
### Mol Ecol Jul-2024 ###
#########################

#devtools::install_github("bcm-uga/lfmm")
#BiocManager::install("LEA")

library(lfmm)
library(LEA)
library(vcfR)
library(mice)

### Genotype:
# Read VCF file (replace 'your_file.vcf' with your actual VCF file)
vcf <- read.vcfR("Re_analysis_2023/PopGen/Data/atlantic.vcf")

# Extract the genotype matrix
geno <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
geno <- sub("0/0", "0", geno)
geno <- sub("0/1", "1", geno)
geno <- sub("1/1", "2", geno)

# Convert the genotype matrix to a data frame
geno_df <- as.data.frame(geno)

# Sanitize column names to ensure they are valid R identifiers
colnames(geno_df) <- make.names(colnames(geno_df))

# Replace non-numeric entries (e.g., missing values represented by ".")
#geno_df[geno_df == "."] <- NA

# Convert to numeric matrix
geno_matrix <- as.matrix(geno_df)
geno_matrix <- apply(geno_matrix, 2, as.numeric)

# Print the dimensions and a sample of the genotype matrix
print(dim(geno_matrix))
print(head(geno_matrix))

# Perform imputation using mice
imputed_data <- mice(geno_matrix, m = 1, method = 'pmm', maxit = 5)

# Extract the completed data
imputed_geno <- complete(imputed_data, 1)

# Convert back to matrix if needed
imputed_geno_matrix <- as.matrix(imputed_geno)

# Print the imputed genotype matrix
print(head(imputed_geno_matrix))

### Environment:

Env<-read.csv(file="./Re_analysis_2023/PopGen/environment_data/In_situ/environmental_summary_stats.csv",
              header = T)
rownames(Env)<-Env$Site
Env$Site<-NULL
Env <- scale(Env, center=TRUE, scale=TRUE)
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')
Env <- as.data.frame(Env)
Env <- dplyr::select(Env, annual_min, annual_mean, annual_max, annual_range)

geno_matrix2<-t(imputed_geno_matrix)
row.names(geno_matrix2)<-gsub("_.*","",row.names(geno_matrix2))

pop_order<-as.data.frame(gsub("\\..*","",row.names(geno_matrix2)))
colnames(pop_order)<-"V1"

env<-merge(Env, pop_order, by.x="row.names", by.y = "V1")
env$Row.names<-NULL

# Convert environmental data to matrix
env_data_matrix <- as.matrix(env)

# Write the outputs
write.table(file="geno_matrix2_labels.lfmm", quote = F, geno_matrix2)
write.table(file="env_data_matrix_labels.env", quote = F, env_data_matrix)

# Set the number of latent factors (K). This is often determined by cross-validation.
K <- 6
reps <- 5

# Prepare and run the LFMM project
project <- lfmm(input.file = "Re_analysis_2023/PopGen/LFMM/geno_matrix2.lfmm",
                env = "Re_analysis_2023/PopGen/LFMM/env_data_matrix.env",
                     K = K, repetitions = reps, project = "new", CPU = 1)

#Record z-scores from the 5 runs in the zs matrix
loci_list<-list()
for(q in 1:4){
zs = z.scores(project, K = 6, d = q)

#Combine z-scores using the median
zs.median = apply(zs, MARGIN = 1, median)

#Compute the GIF
lambda = median(zs.median^2)/qchisq(0.5, df = 1)
lambda

# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)

#histogram of p-values
hist(adj.p.values, col = "red")

# compute adjusted p-values from the combined z-scores
#adj.p.values = pchisq(zs.median^2/.55, df = 9, lower = FALSE)

#histogram of p-values
hist(adj.p.values, col = "green")

## FDR control: Benjamini-Hochberg at level q
## L = number of loci
L = 7671
#fdr level q
q = 0.001
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]

## FDR control: Storey's q-values 
library(qvalue)
#plot(qvalue(adj.p.values))
candidates.qv = which(qvalue(adj.p.values, fdr = .001)$signif)

loci_list[q]<-candidates.qv

df<-read.table(file="Re_analysis_2023/PopGen/LFMM/geno_matrix2_labels.lfmm", sep = " ")
colnames(df)<-rownames(geno_df)
sig_loci<-df[,candidates.qv]

write.table(file="Re_analysis_2023/PopGen/LFMM/output_LFMM_Atlantic_rangeSST.txt", sig_loci, quote = F)

}

#min,mean,max,range

####################
### Pacific LFMM ###
####################

### Genotype:
# Read VCF file (replace 'your_file.vcf' with your actual VCF file)
vcf <- read.vcfR("Re_analysis_2023/PopGen/Data/pacific.vcf")

# Extract the genotype matrix
geno <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
geno <- sub("0/0", "0", geno)
geno <- sub("0/1", "1", geno)
geno <- sub("1/1", "2", geno)

# Convert the genotype matrix to a data frame
geno_df <- as.data.frame(geno)

# Sanitize column names to ensure they are valid R identifiers
colnames(geno_df) <- make.names(colnames(geno_df))

# Replace non-numeric entries (e.g., missing values represented by ".")
#geno_df[geno_df == "."] <- NA

# Convert to numeric matrix
geno_matrix <- as.matrix(geno_df)
geno_matrix <- apply(geno_matrix, 2, as.numeric)

# Print the dimensions and a sample of the genotype matrix
print(dim(geno_matrix))
print(head(geno_matrix))

# Perform imputation using mice
imputed_data <- mice(geno_matrix, m = 1, method = 'pmm', maxit = 5)

# Extract the completed data
imputed_geno <- complete(imputed_data, 1)

# Convert back to matrix if needed
imputed_geno_matrix <- as.matrix(imputed_geno)

# Print the imputed genotype matrix
print(head(imputed_geno_matrix))

### Environment:

Env<-read.csv(file="./Re_analysis_2023/PopGen/environment_data/In_situ/environmental_summary_stats.csv",
              header = T)
rownames(Env)<-Env$Site
Env$Site<-NULL
Env <- scale(Env, center=TRUE, scale=TRUE)
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')
Env <- as.data.frame(Env)
Env <- dplyr::select(Env, annual_min, annual_mean, annual_max, annual_range)

geno_matrix2<-t(imputed_geno_matrix)
row.names(geno_matrix2)<-gsub("_.*","",row.names(geno_matrix2))

pop_order<-as.data.frame(gsub("\\..*","",row.names(geno_matrix2)))
colnames(pop_order)<-"V1"

env<-merge(Env, pop_order, by.x="row.names", by.y = "V1")
env$Row.names<-NULL

# Convert environmental data to matrix
env_data_matrix <- as.matrix(env)

# Write the outputs
write.table(file="Re_analysis_2023/PopGen/LFMM/geno_matrix_labels_pacific.lfmm", quote = F, geno_matrix2)
write.table(file="Re_analysis_2023/PopGen/LFMM/env_data_matrix_labels_pacific.env", quote = F, env_data_matrix)

# Set the number of latent factors (K). This is often determined by cross-validation.
K <- 3
reps <- 5

# Prepare and run the LFMM project
project <- lfmm(input.file = "Re_analysis_2023/PopGen/LFMM/geno_matrix_pacific.lfmm",
                env = "Re_analysis_2023/PopGen/LFMM/env_data_matrix_pacific.env",
                K = K, repetitions = reps, project = "new", CPU = 1)

#Record z-scores from the 5 runs in the zs matrix
loci_list<-list()
for(x in 1:4){
  zs = z.scores(project, K = 3, d = x)
  
  #Combine z-scores using the median
  zs.median = apply(zs, MARGIN = 1, median)
  
  #Compute the GIF
  lambda = median(zs.median^2)/qchisq(0.5, df = 1)
  lambda
  
  # compute adjusted p-values from the combined z-scores
  adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)
  
  #histogram of p-values
  hist(adj.p.values, col = "red")
  
  # compute adjusted p-values from the combined z-scores
  #adj.p.values = pchisq(zs.median^2/.55, df = 9, lower = FALSE)
  
  #histogram of p-values
  hist(adj.p.values, col = "green")
  
  ## FDR control: Benjamini-Hochberg at level q
  ## L = number of loci
  L = 7671
  #fdr level q
  q = 0.001
  w = which(sort(adj.p.values) < q * (1:L)/L)
  candidates.bh = order(adj.p.values)[w]
  
  ## FDR control: Storey's q-values 
  library(qvalue)
  #plot(qvalue(adj.p.values))
  candidates.qv = which(qvalue(adj.p.values, fdr = .001)$signif)
  
  loci_list[q]<-candidates.qv
  
  df<-read.table(file="Re_analysis_2023/PopGen/LFMM/geno_matrix_labels_pacific.lfmm", sep = " ")
  colnames(df)<-rownames(geno_df)
  sig_loci<-df[,candidates.qv]
  
  write.table(file="Re_analysis_2023/PopGen/LFMM/output_LFMM_Pacific_minSST.txt", sig_loci, quote = F)
  
}
