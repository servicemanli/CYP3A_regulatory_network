#-----------------------------quality control-----------------------------------

#-------------------------------------read data---------------------------------

SampleCon = read.csv("Confounder.csv", sep = ",", header = T)
dim(SampleCon)  
## 912 samples
## 1 column of ID
## 5 columns of samples information
SampleCon[1:3, ]

SamplePheo = read.csv("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.csv", header = T)
dim(SamplePheo) 
## 980 samples
## 1 column of donor ID
## 3 column of phenotype
SamplePheo[1:3, ]

##---------------------------- extract samples with RIN>=6 ---------------------

table(SampleCon$SMTS) 
## 226 samples of liver
## 187 samples of small intestine

table(SampleCon$SMRIN >= 6)
SampleCon = SampleCon[SampleCon$SMRIN >= 6, ] 
## extract samples with RIN>=6

table(SampleCon$SMTS) 
## 193 samples of liver after filtering
## 175 samples of small intestine after filtering

SampleCon[1:3, ]
LiverCon = SampleCon[SampleCon$SMTS %in% "Liver", ]
SmallIntestineCon = SampleCon[SampleCon$SMTS %in% "Small Intestine", ]

dim(LiverCon) 
## 193   6
dim(SmallIntestineCon) 
## 175   6
LiverCon[1:3, ]
SmallIntestineCon[1:3, ]

library(stringr)

## add donors ID for samples of liver
LiverCon$SUBJID = sapply(LiverCon$SAMPID, 
                         function(SAMPID)
                           {
                            paste(unlist(str_split(SAMPID, "-"))[1:2], collapse = "-")
                           }
                         )
LiverCon[1:3, ]

## add donors ID for samples of small intestine
SmallIntestineCon$SUBJID = sapply(SmallIntestineCon$SAMPID, 
                                  function(SAMPID)
                                    {
                                     paste(unlist(str_split(SAMPID, "-"))[1:2], collapse = "-")
                                    }
                                  )
SmallIntestineCon[1:3, ]


## integrate information phenotype
LiverPheo = SamplePheo[match(LiverCon$SUBJID, SamplePheo$SUBJID), ]
dim(LiverPheo)
## 193   4
LiverPheo[1:3, ]

SmallIntestinePheo = SamplePheo[match(SmallIntestineCon$SUBJID, SamplePheo$SUBJID), ]
dim(SmallIntestinePheo) 
##175   4
SmallIntestinePheo[1:3, ]

Liver_conf_pheno = merge(LiverCon, LiverPheo, by = "SUBJID")  
dim(Liver_conf_pheno) 
##193  10
Liver_conf_pheno[1:3, ]

SmallIntestine_conf_pheno = merge(SmallIntestineCon, SmallIntestinePheo, by = "SUBJID") 
SmallIntestine_conf_pheno[1:3, ]

#----------------- density plot of RIN and discrete SMTSISCH and RIN -----------
library(ggplot2)

## liver samples
fivenum(Liver_conf_pheno$SMTSISCH)
fivenum(Liver_conf_pheno$SMRIN) 
## range 6-9.8

x = as.data.frame(Liver_conf_pheno$SMRIN)
ggplot(x, aes(x = Liver_conf_pheno$SMRIN)) + geom_density(color = "black", fill = "gray")


break1 = seq(0, 1500, by = 300)
labels1 = c("1-300", "301-600", "601-900", "901-1200", "1201-1500")
discrete_SMTSISCH = cut(Liver_conf_pheno$SMTSISCH, breaks = break1, 
                        labels = labels1, right = F,             
                        ordered_result = T, include.lowest=T)  
table(discrete_SMTSISCH)

Liver_conf_pheno$discrete_SMTSISCH = discrete_SMTSISCH
Liver_conf_pheno[1:3, ]

break2 = seq(6, 10, by = 1)
labels2 = c("6-7", "7-8", "8-9", "9-10")
discrete_RIN = cut(Liver_conf_pheno$SMRIN, breaks = break2,    
                   labels = labels2, right = F,             
                   ordered_result = T, include.lowest=T)  
table(discrete_RIN)

Liver_conf_pheno$discrete_RIN = discrete_RIN
Liver_conf_pheno[1:3, ]

## small intestine samples
fivenum(SmallIntestine_conf_pheno$SMTSISCH) 
## range 84-1682
fivenum(SmallIntestine_conf_pheno$SMRIN) 
## range 6-8.8

x = as.data.frame(SmallIntestine_conf_pheno$SMTSISCH)
ggplot(x, aes(x = SmallIntestine_conf_pheno$SMTSISCH)) + geom_density(color = "black", fill = "gray")

break1 = seq(0, 1800, by = 300)
labels1 = c("1-300", "301-600", "601-900", "901-1200", "1201-1500", "1501-1800")
discrete_SMTSISCH = cut(SmallIntestine_conf_pheno$SMTSISCH, breaks = break1,    
                        labels = labels1, right = F,             
                        ordered_result = T, include.lowest=T) 
table(discrete_SMTSISCH)

SmallIntestine_conf_pheno$discrete_SMTSISCH = discrete_SMTSISCH
SmallIntestine_conf_pheno[1:3, ]

break2 = seq(6, 9, by = 1)
labels2 = c("6-7", "7-8", "8-9")
discrete_RIN = cut(SmallIntestine_conf_pheno$SMRIN, breaks = break2,    
                   labels = labels2, right = F,             
                   ordered_result = T, include.lowest=T)
table(discrete_RIN)

SmallIntestine_conf_pheno$discrete_RIN = discrete_RIN
SmallIntestine_conf_pheno[1:3, ]

#saveRDS(SmallIntestine_conf_pheno, "SmallIntestineAnn.rds")
#saveRDS(Liver_conf_pheno,  "LiverAnn.rds")
#write.csv(Liver_conf_pheno,  "LiverAnn.csv", row.names = F)
#write.csv(SmallIntestine_conf_pheno, "SmallIntestineAnn.csv", row.names = F)

#--------------------------- read gene expression data -------------------------
tpm = read.table("fourtissue_exprs.txt", header = T, row.names = 1, check.names = F)
dim(tpm)  
##  56200    genes
##    912  samples
tpm[1:3, 1:4]

tpm_sum = apply(tpm, 1, sum)
Liver_tpm = tpm[, match(Liver_conf_pheno$SAMPID, colnames(tpm))]
table(colnames(Liver_tpm) == Liver_conf_pheno$SAMPID)

SmallIntestine_tpm = tpm[, match(SmallIntestine_conf_pheno$SAMPID, colnames(tpm))]
table(colnames(SmallIntestine_tpm) == SmallIntestine_conf_pheno$SAMPID)

dim(Liver_tpm) 
## 56200   genes
##   193 samples

dim(SmallIntestine_tpm) 
## 56200   genes
##   175 samples

Liver_tpm[1:3, 1:3]
SmallIntestine_tpm[1:3, 1:3]

#----------------------- calculate mean expression of genes --------------------
Livermean = apply(Liver_tpm, 1, mean) 
SmallIntestinemean = apply(SmallIntestine_tpm, 1, mean)

quantile(Livermean)
quantile(SmallIntestinemean)

hist(SmallIntestinemean)
hist(Livermean)

#saveRDS(SmallIntestine_tpm, "Original_SmallIntestine_tpm.rds")
#saveRDS(Liver_tpm, "Original_Liver_tpm.rds")


#-------------------------- expression data preprocessing ----------------------

#------------------------------------- genes filter ----------------------------

# genes with TPM > quantile(m.mean)[2] seen in more than 80% of the samples 
filter_5 <- function(dat, top_rate = 2, zero_th = 0.8)
  {
  m.mean = apply(dat, 1, mean) 
  condition = rowSums(dat > quantile(m.mean)[2])  
  dat_filter = dat[condition/dim(dat)[2] > zero_th, ] 
  return(dat_filter)
  }

dat = Liver_tpm
Liver_tpm_filter = filter_5(Liver_tpm)
SmallIntestine_tpm_filter = filter_5(SmallIntestine_tpm)

dim(Liver_tpm_filter)  
## 22304   genes
##   193 samples
Liver_tpm_filter[1:3, 1:3]

dim(SmallIntestine_tpm_filter) 
## 25464   genes
##   175 samples
SmallIntestine_tpm_filter[1:3, 1:3]


#--------------------------------- quantile normalization ----------------------
library(preprocessCore)

## liver samples
dim(Liver_tpm_filter) 
## 22304   193
Liver_tpm_filter[1:3, 1:3] 

Liver_quantile = normalize.quantiles(as.matrix(Liver_tpm_filter))
Liver_quantile[1:3, 1:3]

colnames(Liver_quantile) = colnames(Liver_tpm_filter)
rownames(Liver_quantile) = rownames(Liver_tpm_filter)
Liver_quantile[1:3, 1:3]

## small intestine samples
dim(SmallIntestine_tpm_filter) 
## 25464   175
SmallIntestine_tpm_filter[1:3, 1:3]

SmallIntestine_quantile = normalize.quantiles(as.matrix(SmallIntestine_tpm_filter))
SmallIntestine_quantile[1:3, 1:3]

colnames(SmallIntestine_quantile) = colnames(SmallIntestine_tpm_filter)
rownames(SmallIntestine_quantile) = rownames(SmallIntestine_tpm_filter)
SmallIntestine_quantile[1:3, 1:3]



#--------------------------------- log normalization ---------------------------
Liver_quantile[1:3, 1:3]
Liver_log = log2(Liver_quantile + 1)
Liver_log[1:3, 1:3]

SmallIntestine_quantile[1:3, 1:3]
SmallIntestine_log = log2(SmallIntestine_quantile + 1)
SmallIntestine_log[1:3, 1:3]


#---------------------------- check data before WGCNA --------------------------
## goodSamplesGenes function to check exp matrix before doing WGCNA

library(WGCNA)

Liver_log[1:3, 1:3]
Liver_log_t = as.data.frame(t(Liver_log))
dim(Liver_log_t) 
## 193 22304

Liver_log_t[1:3, 1:3]
gsg = goodSamplesGenes(Liver_log_t, verbose = 3)
gsg$allOK   
## all return T means samples and genes are ready


SmallIntestine_log[1:3, 1:3]
SmallIntestine_log_t = as.data.frame(t(SmallIntestine_log))
dim(SmallIntestine_log_t) 
## 175 25464

SmallIntestine_log_t[1:3, 1:3]
gsg = goodSamplesGenes(SmallIntestine_log_t, verbose = 3)
gsg$allOK  
## all return T means samples and genes are ready



## if there is F returned, could try these:
if (!gsg$allOK)
  {
   if (sum(!gsg$goodGenes)>0)
     printFlush(paste("Removing genes:", paste(names(Liver_log_t)[!gsg$goodGenes], collapse= ", ")));
  
   if (sum(!gsg$goodSamples)>0)
     printFlush(paste("Removing samples:", paste(rownames(Liver_log_t)[!gsg$goodSamples], collapse=", ")));
  
   Liver_log_t= Liver_log_t[gsg$goodSamples, gsg$goodGenes]
  }



#------------------------------ data visualization -----------------------------

dat = as.data.frame(t(Liver_log_t))
dim(dat) 
dat[1:3, 1:3]

batch_df = Liver_conf_pheno
dim(batch_df)
batch_df[1:3,]

identical(batch_df$SAMPID, colnames(dat)) 
## must be TRUE


## PCA plot
if(F)
  {
   library(factoextra)
   pca = prcomp(t(dat), scale=TRUE)
   p1 = fviz_eig(pca, addlabels = TRUE)  
   print(p1)
  
  
  xlab <- paste("PC1","(",round((summary(pca))$importance[2,1]*100,1),"%)",sep="")
  ylab <- paste("PC2","(",round((summary(pca))$importance[2,2]*100,1),"%)",sep="")
  # PCA plot of SMTSISCH
  SMTSISCH = batch_df$discrete_SMTSISCH
  pca1 = ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = SMTSISCH)) + 
    geom_point() + stat_ellipse() + xlab(xlab) + ylab(ylab) 
  print(pca1)
  
  # PCA plot of RIN
  RIN = batch_df$discrete_RIN
  pca2 = ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = RIN)) + 
    geom_point() + stat_ellipse() + xlab(xlab) + ylab(ylab) 
  print(pca2)
  
  # PCA plot of gender
  sex = as.factor(batch_df$SEX)
  pca3 = ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = sex)) + 
    geom_point() + stat_ellipse() + xlab(xlab) + ylab(ylab)
  print(pca3)
  
  # PCA plot of age
  age = as.factor(batch_df$AGE)
  pca4 = ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = age)) + 
    geom_point() + stat_ellipse() + xlab(xlab) + ylab(ylab)
  print(pca4)
  
  # PCA plot of DTHHRDY
  DTHHRDY = as.factor(batch_df$DTHHRDY)
  pca5 = ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, color = DTHHRDY)) + 
    geom_point() + stat_ellipse() + xlab(xlab) + ylab(ylab)
  print(pca5)
  
  }

## boxplot
if(F)
  {
  library(reshape2)
  library(ggplot2)
  dat[1:3, 1:3]
  dim(dat) 
  
  dat_m = data.frame(gene_id = rownames(dat), dat[, 1:ncol(dat)])
  dim(dat_m)
  dat_m[1:3, 1:3]
  dat_L = melt(dat_m, id = "gene_id")
  colnames(dat_L)[2:3] = c("sample", "exprs")
  dim(dat_L)
  dat_L[1:4, ]
  dat_L$sample = gsub(".", "-", dat_L$sample, fixed = TRUE) 
  dat_L[1:4, ]
  
  
  dat_L$SUBJID = batch_df$SUBJID[match(dat_L$sample, batch_df$SAMPID)]  
  dat_L$SMTSISCH = as.factor(batch_df$discrete_SMTSISCH[match(dat_L$sample, batch_df$SAMPID)])  
  dat_L[1:5, ] 
  
  
  boxplot = ggplot(data = dat_L, aes(x = SUBJID, y = exprs), alpha = 0.5) + 
    geom_boxplot(alpha=0.5) + theme(axis.text.x = element_text(angle = 90, size = 5))
  
  }


## cluster plot
if(F)
  {
  sampleTree = hclust(dist(t(dat)), method = "average")
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  }


#----------------------- adjust influence of SMTSISCH --------------------------

dat = as.data.frame(t(SmallIntestine_log_t))
dim(dat)

dat[1:3, 1:3]

batch_df = SmallIntestine_conf_pheno
dim(batch_df)

batch_df[1:3,]

identical(batch_df$SAMPID, colnames(dat)) 
## must be TRUE

gene_adjust =  apply(dat, 1, 
                     function(gene)
                       {
                        fit = lm(gene ~ batch_df$SMTSISCH)
                        xx = fit$residuals
                       }
                     )


gene_adjust[1:3, 1:3]

dat_adjust = t(gene_adjust)
dim(dat_adjust)

dat_adjust[1:3, 1:3]

identical(colnames(dat_adjust), batch_df$SAMPID)

#saveRDS(dat_adjust, "Liver_adjust.rds")
#saveRDS(dat_adjust,  "SmallIntestine_adjust.rds")


#------------------ visualization of data after adjustment --------------------- 

## PCA plot
if(F)
  {
  library(factoextra)
  pca_adjust = prcomp(t(dat_adjust), scale=TRUE)
  p1 = fviz_eig(pca_adjust, addlabels = TRUE)  
  print(p1)
  
  
  xlab <- paste("PC1","(",round((summary(pca_adjust))$importance[2,1]*100,1),"%)",sep="")
  ylab <- paste("PC2","(",round((summary(pca_adjust))$importance[2,2]*100,1),"%)",sep="")
  # PCA plot of SMTSISCH
  SMTSISCH = batch_df$discrete_SMTSISCH
  pca1 = ggplot(as.data.frame(pca_adjust$x), aes(x = PC1, y = PC2, color = SMTSISCH)) + 
    geom_point() + stat_ellipse() + xlab(xlab) + ylab(ylab) 
  print(pca1)
  
  # PCA plot of RIN
  RIN = batch_df$discrete_RIN
  pca2 = ggplot(as.data.frame(pca_adjust$x), aes(x = PC1, y = PC2, color = RIN)) + 
    geom_point() + stat_ellipse() + xlab(xlab) + ylab(ylab) 
  print(pca2)
  
  # PCA plot of gender
  sex = as.factor(batch_df$SEX)
  pca3 = ggplot(as.data.frame(pca_adjust$x), aes(x = PC1, y = PC2, color = sex)) + 
    geom_point() + stat_ellipse() + xlab(xlab) + ylab(ylab)
  print(pca3)
  
  # PCA plot of age
  age = as.factor(batch_df$AGE)
  pca4 = ggplot(as.data.frame(pca_adjust$x), aes(x = PC1, y = PC2, color = age)) + 
    geom_point() + stat_ellipse() + xlab(xlab) + ylab(ylab)
  print(pca4)
  
  #PCA plot of DTHHRDY
  DTHHRDY = as.factor(batch_df$DTHHRDY)
  pca5 = ggplot(as.data.frame(pca_adjust$x), aes(x = PC1, y = PC2, color = DTHHRDY)) + 
    geom_point() + stat_ellipse() + xlab(xlab) + ylab(ylab)
  print(pca5)
  
}

## boxplot
if(F)
  {
  library(reshape2)
  library(ggplot2)
  dat_adjust[1:3, 1:3]
  dim(dat_adjust) 
  
  dat_adjust_m = data.frame(gene_id = rownames(dat_adjust), dat[, 1:ncol(dat_adjust)])
  dim(dat_adjust_m)
  dat_adjust_m[1:3, 1:3]
  dat_adjust_L = melt(dat_adjust_m, id = "gene_id")
  colnames(dat_adjust_L)[2:3] = c("sample", "exprs")
  dim(dat_adjust_L)
  dat_adjust_L[1:4, ]
  dat_adjust_L$sample = gsub(".", "-", dat_adjust_L$sample, fixed = TRUE)  
  dat_adjust_L[1:4, ]
  
  
  dat_adjust_L$SUBJID = batch_df$SUBJID[match(dat_adjust_L$sample, batch_df$SAMPID)]  
  dat_adjust_L$SMTSISCH = as.factor(batch_df$discrete_SMTSISCH[match(dat_adjust_L$sample, batch_df$SAMPID)])  
  dat_adjust_L[1:5, ] 
  
  
  boxplot = ggplot(data = dat_adjust_L, aes(x = SUBJID, y = exprs), alpha = 0.5) + 
    geom_boxplot(alpha=0.5) + theme(axis.text.x = element_text(angle = 90, size = 5))
}


## heatmap
if(F)
  {
  library(pheatmap)
  select_dat = t(scale(t(dat_adjust)))
  dim(select_dat)
  select_dat[1:4,1:4] 
  ac = data.frame
  (
    age = as.factor(batch_df$AGE),
    gender = as.factor(batch_df$SEX),
    dead_type = as.factor(batch_df$DTHHRDY)
  )  
  
  ac[1:3, ]
  dim(ac)
  
  rownames(ac) = batch_df$SAMPID
  ac[1:3, ]
  
  lapply(ac, function(x){ class(x)})
  
  
  pheatmap(select_dat, cluster_col = TRUE, cluster_row = FALSE,
           show_rownames = F, show_colnames = F, cellwidth  = 2,   
           annotation_col = ac[, 3:6], annotation_legend = T, fontsize = 10)
  
  pheatmap(select_dat, cluster_col = TRUE, cluster_row = FALSE,
           show_rownames = F, show_colnames = F, cellwidth  = 2,   
           annotation_col = ac, annotation_legend = T, fontsize = 10)
  }

## cluster plot

if(F)
  {   
  sampleTree = hclust(dist(t(dat_adjust)), method = "average")
  }

sampleTree = readRDS("Liver_sampleTree.rds")

par(cex = 0.4)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)



## Construction of the gene network and identification of modules
#----------------------- preparing packages and functions ----------------------

library(WGCNA)
library(reshape2)
library(stringr)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Allow multi-threading within WGCNA. This helps speed up certain calculations. 
# At present this call is necessary. 
# Any error here may be ignored but you may want to update WGCNA if you see one. 
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()

Liver_exprs <- readRDS("Liver_adjust.rds")
dim(Liver_exprs) 
##  22304   genes
##    193 samples

Liver_exprs[1:3, 1:4]
Liver_exprs = t(Liver_exprs)

Liver_exprs[1:3, 1:3]

LiverAnn <- readRDS(paste0(step2_input_dir, "LiverAnn.rds"))
dim(LiverAnn)  
## 193 samples

head(LiverAnn)


# module extraction
get_CYP3A4.5.7_moduleinfo <- function(net, original_matrix)
  {
  moduleLabels = net$colors
  CYP3A4_ind = which(names(net$colors) == "ENSG00000160868.14")
  CYP3A5_ind = which(names(net$colors) == "ENSG00000106258.13")
  CYP3A7_ind = which(names(net$colors) == "ENSG00000160870.12")
  
  CYP3A4_moduleID = net$colors[CYP3A4_ind]
  CYP3A5_moduleID = net$colors[CYP3A5_ind]
  CYP3A7_moduleID = net$colors[CYP3A7_ind]
  
  CYP3A4_module = labels2colors(CYP3A4_moduleID)
  CYP3A5_module = labels2colors(CYP3A5_moduleID)
  CYP3A7_module = labels2colors(CYP3A7_moduleID)
  
  CYP3A4_moduleGene_num = length(rownames(original_matrix)[moduleLabels == CYP3A4_moduleID])
  CYP3A5_moduleGene_num = length(rownames(original_matrix)[moduleLabels == CYP3A5_moduleID])
  CYP3A7_moduleGene_num = length(rownames(original_matrix)[moduleLabels == CYP3A7_moduleID])
  
  CYP3A4.5.7 = data.frame(gene = c("CYP3A4", 'CYP3A5', 'CYP3A7'), 
                          ind = c(CYP3A4_ind, CYP3A5_ind, CYP3A7_ind), 
                          moduleID = c(CYP3A4_moduleID, CYP3A5_moduleID, CYP3A7_moduleID),
                          moduleColor = c(CYP3A4_module, CYP3A5_module, CYP3A7_module),
                          moduleGene_num = c(CYP3A4_moduleGene_num, CYP3A5_moduleGene_num, CYP3A7_moduleGene_num))
  CYP3A4.5.7$gene = as.character(CYP3A4.5.7$gene)
  return(CYP3A4.5.7)
  }




## extract exp matrix of modules
get_CYP3A4.5.7_module_matrix <- function(original_matrix = Liver_exprs, net = bwnet, moduleinfo = CYP3A4.5.7)
  {
  original_matrix = as.matrix(original_matrix)
  moduleLabels = net$colors
  module_exprs_list = list()
  for(i in 1:dim(moduleinfo)[1])
    {
    module_exprs_list[[i]] = original_matrix[, moduleLabels == moduleinfo$moduleID[i]];
    names(module_exprs_list)[i] = moduleinfo$gene[i];
    }
  
  return(module_exprs_list)
  }


#------------------------- WGCNA analysis of liver samples ---------------------

#------------------------------ pick soft threshold ----------------------------

if(F)
  {
  ## bicor,unsigned
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  corType = "bicor"
  networkType = "unsigned"
  maxPOutliers = 0.1
  sft = pickSoftThreshold(Liver_exprs, 
                          powerVector = powers, 
                          networkType = networkType,
                          corFnc = corType,  
                          corOptions = list(maxPOutliers = 0.1), 
                          verbose = 5, 
                          RsquaredCut = 0.8)
  saveRDS(sft, file = "unsigned_bicor_sft.rds")
  }

sft = readRDS("unsigned_bicor_sft.rds")
sft$powerEstimate	

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, 
     cex = cex1, 
     col = "red")

# this line corresponds to using an R^2 cut-off of h
abline(h = 0.80, col = "red") 
# Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, 
     cex = cex1, 
     col = "red")



#--------------------------- Block-wise network construction -------------------

if(F)
  {
  corType = "bicor"
  networkType = "unsigned"
  maxPOutliers = 0.1
  power = 4
  setwd(step2_output_dir)
  bwnet = blockwiseModules(Liver_exprs, 
                           maxBlockSize = 30000,
                           power = power, 
                           maxPOutliers = maxPOutliers, 
                           networkType = networkType, 
                           corType = corType, 
                           TOMType = networkType,
                           minModuleSize = 20, 
                           mergeCutHeight = 0.1,
                           reassignThreshold = 0, 
                           nThreads = 8, 
                           deepSplit = 4, 
                           numericLabels = TRUE, 
                           saveTOMs = T,
                           saveTOMFileBase = "bicor_unsigned_TOM-blockwise",
                           verbose = 3)
  saveRDS("unsigned_bicor_bwnet.rds")
  }

bwnet = readRDS("unsigned_bicor_bwnet.rds")

moduleLabels = bwnet$colors

table(moduleLabels)

# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels)

# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)

table(bwModuleColors)

# Plot the dendrogram and the module colors underneath

CYP3A4.5.7 = get_CYP3A4.5.7_moduleinfo(net = bwnet, 
                                       original_matrix = Liver_exprs)

CYP3A4.5.7



#------------------------------------ visualization ----------------------------

if(F)
  {
  plotDendroAndColors(bwnet$dendrograms[[1]], 
                      bwModuleColors[bwnet$blockGenes[[1]]],
                      "Module colors", 
                      main = "Gene dendrogram and module colors of Liver tissue dataset",
                      dendroLabels = FALSE, 
                      hang = 0.03,
                      addGuide = TRUE, 
                      guideHang = 0.05)
  
  
  MEs = bwnet$MEs
  
  MEs_col = MEs
  colnames(MEs_col) = paste0("ME", 
                             labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME","")))
                             )
  MEs_col = orderMEs(MEs_col)
  
  # marDendro/marHeatmap 
  plotEigengeneNetworks(MEs_col, 
                        "Eigengene adjacency heatmap of Liver tissue dataset", 
                        marDendro = c(3,3,2,4),
                        marHeatmap = c(3,4,2,2), 
                        plotDendrograms = T, 
                        xLabelsAngle = 90)  
  
  
  CYP3A4_module =  CYP3A4.5.7$moduleColor[1]
  ME = MEs_col[, paste("ME",CYP3A4_module, sep="")]
  
  par(mfrow = c(2,1), mar = c(0, 4.1, 4, 2.05))
  
  plotMat(t(scale(Liver_exprs)),
          nrgcols = 30, 
          rlabels = F, 
          rcols = CYP3A4_module,
          main = CYP3A4_module, 
          cex.main = 2)
  
  par(mar = c(2, 2.3, 0.5, 0.8))
  
  barplot(ME, 
          col = CYP3A4_module,
          main = "", cex.main = 2,
          ylab = "eigengene expression", 
          xlab = "sample")
  
  }



#--------------------------- calculate kME of genes ----------------------------

if(F)
  {
  maxPOutliers = 0.1
  power = 4
  
  connet = abs(bicor(Liver_exprs, maxPOutliers = maxPOutliers , nThreads = 8)^power)
  
  bwModuleColors = labels2colors(bwnet$colors)
  
  Alldegrees = intramodularConnectivity(connet, bwModuleColors)
  
  head(Alldegrees)
  
  AllKME = signedKME(Liver_exprs, MEs_col, corFnc = "bicor")
  
  dim(AllKME) 
  #22304    96
  
  AllKME[1:3, 1:3]
  
  saveRDS(AllKME, "AllKME.rds")
  }



#------------------------ get infomation of target modules ---------------------

#### add Description for target modules

ENSEMBL_Description = read.table("ENSEMBL_Description.txt", header = T)
head(ENSEMBL_Description)

table(bwModuleColors == "turquoise")

ModuleGeneName = colnames(Liver_exprs)[bwModuleColors == "turquoise"]
length(ModuleGeneName)  
## 3103 genes

head(ModuleGeneName)

table(ModuleGeneName %in% ENSEMBL_Description$ENSEMBL)

Description = ENSEMBL_Description[match(ModuleGeneName, ENSEMBL_Description$ENSEMBL), 2]
head(Description)

standard_ENSEMBL = sapply(ModuleGeneName, 
                          function(gene)
                            {unlist(str_split(gene, '[.]'))[1]}
                          ) 

ENSEMBL_df = data.frame(gene_id = ModuleGeneName, standard_ENSEMBL, Description)
head(ENSEMBL_df)

identical(rownames(ENSEMBL_df), as.character(ENSEMBL_df$gene_id)) 
## must be TRUE



## add genes connectivity for target modules

Alldegrees = readRDS("Alldegrees.rds")
dim(Alldegrees) 
## 22304  genes 
##     4 column

head(Alldegrees)

ModuleGene_degrees = Alldegrees[match(ModuleGeneName, rownames(Alldegrees)), "kWithin"]
length(ModuleGene_degrees) 
## 3103 genes



## add kME for target modules

AllKME = readRDS( "AllKME.rds" )
ModuleGene_KME = AllKME[match(ModuleGeneName, rownames(AllKME)), which("kMEturquoise" ==  colnames(AllKME))]

length(ModuleGene_KME) 
## 3103 genes

table(abs(ModuleGene_KME) > 0.8) 
## 265 genes |KME| > 0.8



module_geneinfo = data.frame(gene_id = ENSEMBL_df$gene_id,
                             ENSEMBL = ENSEMBL_df$standard_ENSEMBL,
                             Description = ENSEMBL_df$Description,
                             kME = ModuleGene_KME,
                             KIN = ModuleGene_degrees
                             #ensembl_symbol = ENSEMBL_df$ensembl_symbol,
                             #ensembl_entrezid = ENSEMBL_df$ensembl_entrezid,
                             #KIN = ModuleGene_degree$kWithin,
                             #MMturquoise = ModuleGene_geneModuleMembership,
                             #p.MMturquoise = ModuleGene_MMPvalue
                             #GS.Liver = ModuleGene_GS.Liver,
                             #GSPvalue.Liver = ModuleGene_GSPvalue,
                             )

module_geneinfo = module_geneinfo[order(module_geneinfo$kME, decreasing = T), ]                     
head(module_geneinfo)

dim(module_geneinfo)  
## 3103 genes

fivenum(abs(module_geneinfo$kME)) 
# range 0.1988422-0.9182030

module_hubgene_info = module_geneinfo[abs(module_geneinfo$kME) > 0.8, ] 
## |kME| > 0.8 regard to be hub genes

dim(module_hubgene_info) 
# 265 genes

#saveRDS(module_geneinfo, "turquoise_geneinfo.rds")
#write.csv(module_geneinfo, "turquoise_geneinfo.csv", row.names = F)

#saveRDS(module_hubgene_info,  "turquoise_hubgene_info.rds")
#write.csv(module_hubgene_info,  "turquoise_hubgene_info.csv", row.names = F)



#---------------------- import target modules to cytoscape ---------------------

## extract TOM matrix of target modules

bwLabels = matchLabels(bwnet$colors, moduleLabels) 
bwModuleColors = labels2colors(bwLabels) 
table(bwModuleColors)

if(F)
  {
  TOM = load("/blockwiseModules_output/bicor_OneBlock_mergeCutHeight0.1_deepSplit4_TOM-blockwise-block.1.RData")
  TOM = as.matrix(TOM)
  # Select module
  select_module = "turquoise" 
  select_module_flag = as.vector(bwModuleColors == select_module)
  modTOM = TOM[select_module_flag, select_module_flag]
  }


## get TOM matrix of all genes and interaction weight value > 0.02

modTOM = readRDS("turquoise_TOM.rds", sep = "")
dim(modTOM) 
# 3370 genes
modTOM[1:3, 1:3]

select_module = "turquoise" 
select_module_flag = as.vector(bwModuleColors == select_module)

select_modgene = colnames(Liver_exprs)[select_module_flag]

ENSEMBL = sapply(select_modgene, function(x){unlist(str_split(x, '[.]'))[1]})

ENSEMBL = as.data.frame(ENSEMBL)
head(ENSEMBL)

turquoise_geneinfo = readRDS("turquoise_geneinfo.rds")
head(turquoise_geneinfo)

select_modDescription = turquoise_geneinfo$Description[match(rownames(ENSEMBL), turquoise_geneinfo$gene_id)]

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(select_module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(select_module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = select_modgene,
                               altNodeNames = select_modDescription,
                               nodeAttr = bwModuleColors[select_module_flag])



#--------------------- WGCNA analysis of small intestine samples ---------------

#------------------------------ pick soft threshold ----------------------------

enableWGCNAThreads()

SmallIntestine_exprs <- readRDS("SmallIntestine_adjust.rds")

dim(SmallIntestine_exprs) 
##  22304   genes
##    193 samples

SmallIntestine_exprs[1:3, 1:4]

SmallIntestine_exprs = t(SmallIntestine_exprs)
SmallIntestine_exprs[1:3, 1:3]

SmallIntestineAnn <- readRDS( "SmallIntestineAnn.rds")
dim(SmallIntestineAnn)  
## 193 samples

head(SmallIntestineAnn)


if(F)
  {
  ###### bicor,unsigned
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  corType = "bicor"
  networkType = "unsigned"
  maxPOutliers = 0.1
  sft = pickSoftThreshold(SmallIntestine_exprs, 
                          powerVector = powers, 
                          networkType = networkType,
                          corFnc = corType,  
                          corOptions = list(maxPOutliers = 0.1), 
                          verbose = 5, 
                          RsquaredCut = 0.8)
  
  saveRDS(sft, file = "unsigned_bicor_sft.rds")
  
  }

sft = readRDS("unsigned_bicor_sft.rds")
sft$powerEstimate

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, 
     cex = cex1, 
     col = "red")

# this line corresponds to using an R^2 cut-off of h
abline(h = 0.80, col = "red") 
# Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, 
     cex = cex1, 
     col = "red")




#------------------------Block-wise network construction -----------------------

if(F)
  {
  corType = "bicor"
  networkType = "unsigned"
  maxPOutliers = 0.1
  power = 7

  bwnet = blockwiseModules(SmallIntestine_exprs, 
                           maxBlockSize = 30000,
                           power = power, 
                           maxPOutliers = maxPOutliers, 
                           networkType = networkType, 
                           corType = corType, 
                           TOMType = networkType,
                           minModuleSize = 20,
                           mergeCutHeight = 0.1,
                           reassignThreshold = 0, 
                           nThreads = 8, 
                           deepSplit = 4, 
                           numericLabels = TRUE, 
                           saveTOMs = T,
                           saveTOMFileBase = "bicor_unsigned_TOM-blockwise",
                           verbose = 3)
  }

bwnet = readRDS("maxBlockSize30000_power4/bwnet.rds")

moduleLabels = bwnet$colors

table(moduleLabels) 

# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels)

# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)

table(bwModuleColors)

# Plot the dendrogram and the module colors underneath

CYP3A4.5.7 = get_CYP3A4.5.7_moduleinfo(net = bwnet, original_matrix = SmallIntestine_exprs)

CYP3A4.5.7




#-------------------------------- visualization --------------------------------

if(F)
  {
  plotDendroAndColors(bwnet$dendrograms[[1]], 
                      bwModuleColors[bwnet$blockGenes[[1]]],
                      "Module colors", 
                      main = "Gene dendrogram and module colors of Small Intestine",
                      dendroLabels = FALSE, 
                      hang = 0.03,
                      addGuide = TRUE, 
                      guideHang = 0.05)
  
  
  MEs = bwnet$MEs
  
  MEs_col = MEs
  
  colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  
  MEs_col = orderMEs(MEs_col)
  
  # marDendro/marHeatmap 
  plotEigengeneNetworks(MEs_col, 
                        "Eigengene adjacency heatmap of Small Intestine", 
                        marDendro = c(3,3,2,4),
                        marHeatmap = c(3,4,2,2),
                        plotDendrograms = T, 
                        xLabelsAngle = 90)  
  
  
  CYP3A4_module =  CYP3A4.5.7$moduleColor[1]
  
  ME = MEs_col[, paste("ME",CYP3A4_module, sep="")]
  
  par(mfrow = c(2,1), mar = c(0, 4.1, 4, 2.05))
  
  plotMat(t(scale(SmallIntestine_exprs)),
          nrgcols = 30, 
          rlabels = F, 
          rcols = CYP3A4_module,
          main = CYP3A4_module, 
          cex.main = 2)
  
  par(mar = c(2, 2.3, 0.5, 0.8))
  
  barplot(ME, 
          col = CYP3A4_module, main = "", 
          cex.main = 2,
          ylab = "eigengene expression", 
          xlab = "sample")
  }




#-------------------------- calculate kME of genes -----------------------------

if(F)
  {
  maxPOutliers = 0.1
  power = 7
  
  connet = abs(bicor(SmallIntestine_exprs, maxPOutliers = maxPOutliers , nThreads = 8)^power)  
  
  bwModuleColors = labels2colors(bwnet$colors)
  
  Alldegrees = intramodularConnectivity(connet, bwModuleColors)
  
  head(Alldegrees)
  
  AllKME = signedKME(SmallIntestine_exprs, MEs_col, corFnc = "bicor")
  
  dim(AllKME) 
  ## 22304    96
  
  AllKME[1:3, 1:3]
  
  saveRDS(AllKME, file = "maxBlockSize30000_power4/AllKME.rds")
  }



#---------------------- get information of target modules ----------------------

## add Description for target modules

ENSEMBL_Description = read.table("ENSEMBL_Description.txt", header = T)
head(ENSEMBL_Description)

table(bwModuleColors == "red")

ModuleGeneName = colnames(SmallIntestine_exprs)[bwModuleColors == "red"]
length(ModuleGeneName)  
## 3370 genes

head(ModuleGeneName)

table(ModuleGeneName %in% ENSEMBL_Description$ENSEMBL)

Description = ENSEMBL_Description[match(ModuleGeneName, ENSEMBL_Description$ENSEMBL), 2]
head(Description)

standard_ENSEMBL = sapply(ModuleGeneName, function(gene){unlist(str_split(gene, '[.]'))[1]})

ENSEMBL_df = data.frame(gene_id = ModuleGeneName, standard_ENSEMBL, Description)
head(ENSEMBL_df)

identical(rownames(ENSEMBL_df), as.character(ENSEMBL_df$gene_id)) 
## must be TRUE




## add KME for target modules

AllKME = readRDS("maxBlockSize30000_power4/AllKME.rds")

ModuleGene_KME = AllKME[match(ModuleGeneName, rownames(AllKME)), which("kMEred" ==  colnames(AllKME))] 
length(ModuleGene_KME) 
## 3370 genes

table(abs(ModuleGene_KME) > 0.8) 
## 267 genes |KME| > 0.8



module_geneinfo = data.frame(gene_id = ENSEMBL_df$gene_id,
                             ENSEMBL = ENSEMBL_df$standard_ENSEMBL,
                             Description = ENSEMBL_df$Description,
                             kME = ModuleGene_KME
                             #ensembl_symbol = ENSEMBL_df$ensembl_symbol,
                             #ensembl_entrezid = ENSEMBL_df$ensembl_entrezid,
                             #KIN = ModuleGene_degree$kWithin,
                             #MMred = ModuleGene_geneModuleMembership, 
                             #p.MMred = ModuleGene_MMPvalue
                             #GS.SmallIntestine = ModuleGene_GS.SmallIntestine,
                             #GSPvalue.SmallIntestine = ModuleGene_GSPvalue,
                             )

module_geneinfo = module_geneinfo[order(module_geneinfo$kME, decreasing = T), ]                     
head(module_geneinfo) 

dim(module_geneinfo)  
## 3103 genes

fivenum(abs(module_geneinfo$kME)) 
## range 0.1988422-0.9182030

module_hubgene_info = module_geneinfo[abs(module_geneinfo$kME) > 0.8, ] 
## |kME| > 0.8 regard to be hub genes

dim(module_hubgene_info) 

#saveRDS(module_geneinfo, file = "maxBlockSize30000_power4/red_geneinfo.rds")
#write.csv(module_geneinfo, file = "maxBlockSize30000_power4/red_geneinfo.csv", row.names = F)
#saveRDS(module_hubgene_info, file = "maxBlockSize30000_power4/red_hubgene_info.rds")
#write.csv(module_hubgene_info, file = "maxBlockSize30000_power4/red_hubgene_info.csv", row.names = F)



#-------------------- import target modules to cytoscape -----------------------

## extract TOM matrix of target modules

bwLabels = matchLabels(bwnet$colors, moduleLabels) 
bwModuleColors = labels2colors(bwLabels) 
table(bwModuleColors)

if(F)
  {
  TOM = load("/blockwiseModules_output/bicor_OneBlock_mergeCutHeight0.1_deepSplit4_TOM-blockwise-block.1.RData") 
  TOM = as.matrix(TOM)
  
  # Select module
  select_module = "red" 
  select_module_flag = as.vector(bwModuleColors == select_module)
  modTOM = TOM[select_module_flag, select_module_flag]
  }

## get TOM matrix of all genes and interaction > 0.02

modTOM = readRDS("maxBlockSize30000_power4/red_TOM.rds")
dim(modTOM) 
# 3370 genes

modTOM[1:3, 1:3]

select_module = "red" 
select_module_flag = as.vector(bwModuleColors == select_module)
select_modgene = colnames(SmallIntestine_exprs)[select_module_flag]

ENSEMBL = sapply(select_modgene, function(x) unlist(str_split(x, '[.]'))[1])
ENSEMBL = as.data.frame(ENSEMBL)

head(ENSEMBL)

red_geneinfo = readRDS("maxBlockSize30000_power4/red_geneinfo.rds")
head(red_geneinfo)

select_modDescription = red_geneinfo$Description[match(rownames(ENSEMBL), red_geneinfo$gene_id)]

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(select_module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(select_module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = select_modgene,
                               altNodeNames = select_modDescription,
                               nodeAttr = bwModuleColors[select_module_flag])


#-------------------------- GO and KEGG pathway analysis -----------------------


#---------------------------------- data preparing -----------------------------

turquoise_geneinfo = readRDS("turquoise_geneinfo.rds")
dim(turquoise_geneinfo)  
## 3370 genes 

head(turquoise_geneinfo)
table(!is.na(turquoise_geneinfo$ensembl_symbol)) 
## 2864 genes have SYMBOL in ensemble


#---------------------------- transform ID into ENTREZID -----------------------

library(org.Hs.eg.db)
library(clusterProfiler)
library(topGO)  


if(F)
  {
  keytypes(org.Hs.eg.db)
  
  eg = bitr(turquoise_geneinfo$ENSEMBL, 
            fromType = "ENSEMBL", 
            toType = c("ENTREZID", "SYMBOL"), 
            OrgDb = "org.Hs.eg.db")
  
  dim(eg)  
  # 2553 succeed
  
  head(eg)
  
  saveRDS(eg, file = "clusterProfiler_ENSEMBL2ENTREZID2SYMBOL.rds")
  write.csv(eg, file = "clusterProfiler_ENSEMBL2ENTREZID2SYMBOL.csv", row.names = F)
  
  
  head(turquoise_geneinfo)
  
  turquoise_geneinfo$org.Hs.eg.db_Symbol = eg$SYMBOL[match(turquoise_geneinfo$ENSEMBL, eg$ENSEMBL)]
  turquoise_geneinfo$org.Hs.eg.db_Entrezid = eg$ENTREZID[match(turquoise_geneinfo$ENSEMBL, eg$ENSEMBL)]
  head(turquoise_geneinfo)
  
  turquoise_geneinfo_add = data.frame(gene_id = turquoise_geneinfo$gene_id,
                                      Description = turquoise_geneinfo$Description,
                                      ENSEMBL = turquoise_geneinfo$ENSEMBL,
                                      org.Hs.eg.db_Symbol = turquoise_geneinfo$org.Hs.eg.db_Symbol,
                                      org.Hs.eg.db_Entrezid = turquoise_geneinfo$org.Hs.eg.db_Entrezid,
                                      kME = turquoise_geneinfo$kME)
  head(turquoise_geneinfo_add)
  
  saveRDS(turquoise_geneinfo_add, file =  "turquoise_geneinfo_ID_transform.rds")
  write.csv(turquoise_geneinfo_add, file = "turquoise_geneinfo_ID_transform.csv", row.names = F)
  }



#--------------------------------- GO pathway ----------------------------------

eg = readRDS("clusterProfiler_ENSEMBL2ENTREZID2SYMBOL.rds")
dim(eg) 
# 2760 genes have symbol and entrezid ID


GO_result = enrichGO(eg$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "ALL",   
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,  
                     qvalueCutoff = 0.05,
                     keyType = 'ENTREZID',
                     readable = T)

dim(GO_result) 
## 287 terms
head(GO_result)

#saveRDS(GO_result, file =  "/GO/turquoise_enrich_GO_original.rds")
#write.csv(GO_result, file = "/GO/turquoise_enrich_GO_original.csv", row.names = F)



#--------------------------------- KEGG pathway --------------------------------

KEGG_result <- enrichKEGG(gene = eg$ENTREZID,
                          organism = "hsa",
                          pvalueCutoff = 0.05, 
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 0.05)

dim(KEGG_result) 
## 31 terms

head(KEGG_result)

kk = setReadable(KEGG_result, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(kk)

#saveRDS(kk, file = "/KEGG/turquoise_enrich_KEGG_original.rds")
#write.csv(kk, file = "/KEGG/turquoise_enrich_KEGG_original.csv", row.names = F)



#-------------- get GO and KEGG pathway include target CYP3A genes -------------
library(stringr)

GO_result = readRDS("/GO/turquoise_enrich_GO_original.rds")
turquoise_geneinfo = readRDS("turquoise_geneinfo_ID_transform.rds")

head(turquoise_geneinfo)

dim(turquoise_geneinfo) 
## 3103 6


### CYP3A4 GO pathway

CYP3A4_flag = sapply(GO_result@result$geneID, 
                     function(x)
                       {"CYP3A4" %in% unlist(str_split(x, "/"))}
                     )

table(CYP3A4_flag) 
## 42 GO terms include CYP3A4

CYP3A4_GO = GO_result[CYP3A4_flag, asis = T] 

GO_geneID = unlist(str_split(CYP3A4_GO$geneID, "/")) 
GO_geneID_uniq = GO_geneID[!duplicated(GO_geneID)] 
length(GO_geneID_uniq) 

head(GO_geneID_uniq)
CYP3A4_GOterm_gene = turquoise_geneinfo[match(GO_geneID_uniq, turquoise_geneinfo$org.Hs.eg.db_Symbol), ]

identical(length(GO_geneID_uniq), nrow(CYP3A4_GOterm_gene))
head(CYP3A4_GOterm_gene)


## CYP3A5 GO pathway

CYP3A5_flag = sapply(GO_result@result$geneID, 
                     function(x){"CYP3A5" %in% unlist(str_split(x, "/"))}
                     )

table(CYP3A5_flag) 
## 22 GO terms include CYP3A5

CYP3A5_GO = GO_result[CYP3A5_flag, asis = T] 
GO_geneID = unlist(str_split(CYP3A5_GO$geneID, "/")) 
GO_geneID_uniq = GO_geneID[!duplicated(GO_geneID)]

length(GO_geneID_uniq) 
head(GO_geneID_uniq)

CYP3A5_GOterm_gene = turquoise_geneinfo[match(GO_geneID_uniq, turquoise_geneinfo$org.Hs.eg.db_Symbol), ]

identical(length(GO_geneID_uniq), nrow(CYP3A5_GOterm_gene))
head(CYP3A5_GOterm_gene)


## CYP3A7 GO pathway

CYP3A7_flag = sapply(GO_result@result$geneID, 
                     function(x){"CYP3A7" %in% unlist(str_split(x, "/"))}
                     )

table(CYP3A7_flag) 
## 21 GO terms include CYP3A7

CYP3A7_GO = GO_result[CYP3A7_flag, asis = T] 
GO_geneID = unlist(str_split(CYP3A5_GO$geneID, "/"))  
GO_geneID_uniq = GO_geneID[!duplicated(GO_geneID)] 
length(GO_geneID_uniq)

head(GO_geneID_uniq)
CYP3A7_GOterm_gene = turquoise_geneinfo[match(GO_geneID_uniq, turquoise_geneinfo$org.Hs.eg.db_Symbol), ]

identical(length(GO_geneID_uniq), nrow(CYP3A7_GOterm_gene))
head(CYP3A7_GOterm_gene)



library(stringr)

KEGG_result = readRDS("KEGG/turquoise_enrich_KEGG_original.rds")
turquoise_geneinfo = readRDS("turquoise_geneinfo_ID_transform.rds")

head(turquoise_geneinfo)
dim(turquoise_geneinfo) 
## 3370 8


### CYP3A4 KEGG term
CYP3A4_flag = sapply(KEGG_result@result$geneID, 
                     function(x){"CYP3A4" %in% unlist(str_split(x, "/"))}
                     )

table(CYP3A4_flag) 
## 8 KEGG terms include CYP3A4

CYP3A4_KEGG = KEGG_result[CYP3A4_flag, asis = T] 
which(is.na(CYP3A4_KEGG@result$ID))
CYP3A4_KEGG_new = CYP3A4_KEGG[!is.na(CYP3A4_KEGG@result$ID), asis = T]
nrow(CYP3A4_KEGG_new@result)  
## 6 terms

KEGG_geneID = unlist(str_split(CYP3A4_KEGG_new$geneID, "/"))  
KEGG_geneID_uniq = KEGG_geneID[!duplicated(KEGG_geneID)] 
length(KEGG_geneID_uniq) 

head(KEGG_geneID_uniq)

CYP3A4_KEGGterm_gene = turquoise_geneinfo[match(KEGG_geneID_uniq, turquoise_geneinfo$org.Hs.eg.db_Symbol), ]
identical(length(KEGG_geneID_uniq), nrow(CYP3A4_KEGGterm_gene))
head(CYP3A4_KEGGterm_gene)

## CYP3A5 KEGG terms
CYP3A5_flag = sapply(KEGG_result@result$geneID, 
                     function(x){"CYP3A5" %in% unlist(str_split(x, "/"))}
                     )

table(CYP3A5_flag) 
## 5 KEGG terms include CYP3A5

CYP3A5_KEGG = KEGG_result[CYP3A5_flag, asis = T] 
which(is.na(CYP3A5_KEGG@result$ID)) 
KEGG_geneID = unlist(str_split(CYP3A5_KEGG$geneID, "/")) 
KEGG_geneID_uniq = KEGG_geneID[!duplicated(KEGG_geneID)] 
length(KEGG_geneID_uniq)

head(KEGG_geneID_uniq)
CYP3A5_KEGGterm_gene = turquoise_geneinfo[match(KEGG_geneID_uniq, turquoise_geneinfo$org.Hs.eg.db_Symbol), ]

identical(length(KEGG_geneID_uniq), nrow(CYP3A5_KEGGterm_gene))
head(CYP3A5_KEGGterm_gene)


## CYP3A7 KEGG term
CYP3A7_flag = sapply(KEGG_result@result$geneID, 
                     function(x){"CYP3A7" %in% unlist(str_split(x, "/"))}
                     )

table(CYP3A7_flag) 
## 3 KEGG terms include CYP3A7

CYP3A7_KEGG = KEGG_result[CYP3A7_flag, asis = T] 
which(is.na(CYP3A7_KEGG@result$ID)) 
KEGG_geneID = unlist(str_split(CYP3A5_KEGG$geneID, "/"))  
KEGG_geneID_uniq = KEGG_geneID[!duplicated(KEGG_geneID)] 
length(KEGG_geneID_uniq)

head(KEGG_geneID_uniq)
CYP3A7_KEGGterm_gene = turquoise_geneinfo[match(KEGG_geneID_uniq, turquoise_geneinfo$org.Hs.eg.db_Symbol), ]

identical(length(KEGG_geneID_uniq), nrow(CYP3A7_KEGGterm_gene))
head(CYP3A7_KEGGterm_gene)


#----------------------- GO KEGG pathway visualization -------------------------

## GO terms ¡ª¡ª CYP3A4(42),CYP3A5(22), CYP3A7(21)
## KEGG terms ¡ª¡ª CYP3A4(6),CYP3A5(5), CYP3A7(3)

library(clusterProfiler)
library(topGO)  
library(ggplot2)

## bubble plot

p1 = dotplot(CYP3A4_GO, 
             x = "Count",
             showCategory = 20, 
             color = "qvalue", 
             title="Enrichment GO Top20")
print(p1)


p1 = dotplot(CYP3A5_GO, 
             x = "Count", 
             showCategory = nrow(CYP3A5_GO), 
             color = "qvalue", 
             title="Enrichment GO All")

p1 + scale_y_discrete(labels = function(x) str_wrap(x))


p1 = dotplot(CYP3A7_GO, 
             x = "Count", 
             showCategory = nrow(CYP3A7_GO), 
             color = "qvalue", 
             title="Enrichment GO All")

p1 + scale_y_discrete(labels = function(x) str_wrap(x))



p1 = dotplot(CYP3A4_KEGG, 
             x = "Count", 
             showCategory = nrow(CYP3A4_KEGG), 
             color = "qvalue", 
             title="Enrichment KEGG All")
print(p1)


p1 = dotplot(CYP3A5_KEGG, 
             x = "Count", 
             showCategory = nrow(CYP3A5_KEGG), 
             color = "qvalue", 
             title="Enrichment KEGG All")

p1 + scale_y_discrete(labels = function(x) str_wrap(x))


p1 = dotplot(CYP3A7_KEGG, 
             x = "Count", 
             showCategory = nrow(CYP3A7_KEGG), 
             color = "qvalue", 
             title="Enrichment KEGG All")

p1 + scale_y_discrete(labels = function(x) str_wrap(x))




## barplot

p1 = barplot(CYP3A4_GO[order(CYP3A4_GO@result$qvalue), asis = T], 
             showCategory = 20, 
             color = "qvalue", 
             title="Enrichment GO Top20") 

p1 + scale_x_discrete(labels = function(x) str_wrap(x))


p1 = barplot(CYP3A5_GO[order(CYP3A5_GO@result$qvalue), asis = T], 
             showCategory = nrow(CYP3A5_GO), 
             color = "qvalue", 
             title="Enrichment GO All")

p1 + scale_x_discrete(labels = function(x) str_wrap(x))



p1 = barplot(CYP3A7_GO[order(CYP3A7_GO@result$qvalue), asis = T], 
             showCategory = nrow(CYP3A7_GO), 
             color = "qvalue", 
             title="Enrichment GO All")

p1 + scale_x_discrete(labels = function(x) str_wrap(x))



p1 = barplot(CYP3A4_KEGG[order(CYP3A4_KEGG@result$qvalue), asis = T], 
             showCategory = nrow(CYP3A4_KEGG), 
             color = "qvalue", 
             title="Enrichment KEGG All") 

p1 + scale_x_discrete(labels = function(x) str_wrap(x))


p1 = barplot(CYP3A5_KEGG[order(CYP3A5_GO@result$qvalue), asis = T], 
             showCategory = nrow(CYP3A5_KEGG), 
             color = "qvalue", 
             title="Enrichment KEGG All")

p1 + scale_x_discrete(labels = function(x) str_wrap(x))



p1 = barplot(CYP3A7_KEGG[order(CYP3A7_KEGG@result$qvalue), asis = T], 
             showCategory = nrow(CYP3A7_KEGG), 
             color = "qvalue", 
             title="Enrichment KEGG All")

p1 + scale_x_discrete(labels = function(x) str_wrap(x))



dim(GO_result) 
## 287 GO terms

dim(KEGG_result) 
## 31 KEGG terms


turquoise_geneinfo = readRDS("turquoise_geneinfo_ID_transform.rds")
head(turquoise_geneinfo)

dim(turquoise_geneinfo) 
## 3103 6


####### bubble plot
p1 = dotplot(GO_result[order(GO_result@result$Count, decreasing = T), asis = T],
             x = "Count",
             showCategory = 20, 
             color = "p.adjust", 
             title="Enrichment GO Top20 of Turquoise Module Genes")
print(p1)



p1 = dotplot(KEGG_result[order(KEGG_result@result$Count, decreasing = T), asis = T], 
             x = "Count", 
             showCategory = nrow(KEGG_result), 
             color = "p.adjust", 
             title ="Enrichment KEGG All of Turquoise Module Genes")

print(p1)



####### barplot

p1 = barplot(GO_result[order(GO_result@result$Count, decreasing = T), asis = T], 
             showCategory = 20, 
             color = "p.adjust", 
             title="Enrichment GO Top20 of Turquoise Module Genes") 

p1 + scale_x_discrete(labels = function(x) str_wrap(x))



p1 = barplot(KEGG_result[order(KEGG_result@result$Count, decreasing = T), asis = T], 
             showCategory = nrow(KEGG_result), 
             color = "p.adjust", 
             title="Enrichment KEGG All of Turquoise Module Genes")

p1 + scale_x_discrete(labels = function(x) str_wrap(x))




#------------------ genes annotation of TF lncRNA miRNA ------------------------

#--------------------------- data preparing ------------------------------------

animalTFDB_humanTF = read.table( "animalTFDB_Homo_sapiens_TF.txt", sep = "\t", header = T)
dim(animalTFDB_humanTF) 
## 1665 TFs

head(animalTFDB_humanTF)

table(animalTFDB_humanTF$Species) 
## 1665 human TFs

TFcheckpoint_humanTF = read.csv( "TFcheckpoint_humanTF.csv", header = T, row.names = 1)
dim(TFcheckpoint_humanTF) 
## 570  36

table(TFcheckpoint_humanTF$taxon)
TFcheckpoint_humanTF[1:3, 1:3]

intersect_humanTF = intersect(TFcheckpoint_humanTF$gene_symbol, animalTFDB_humanTF$Symbol)
length(intersect_humanTF) 
## 514 same TFs


humanTF = animalTFDB_humanTF


Gencode_lncRNA = read.csv("/Gencode_v33/lncRNA_v33.csv", header = T)
Gencode_miRNA = read.csv("/Gencode_v33/miRNA_v33.csv", header = T)

dim(Gencode_lncRNA) 
## 17952 3

dim(Gencode_miRNA) 
## 1881 3

head(Gencode_lncRNA)
head(Gencode_miRNA)


turquoise_geneinfo = readRDS("turquoise_geneinfo_ID_transform.rds") 
dim(turquoise_geneinfo) 
## 3103 genes

head(turquoise_geneinfo)
fivenum(abs(turquoise_geneinfo$kME))

turquoise_hubgeneinfo = readRDS("turquoise_hubgene_info.rds")
dim(turquoise_hubgeneinfo) 
## 265 genes

head(turquoise_hubgeneinfo)
turquoise_hubgeneinfo = turquoise_geneinfo[match(turquoise_hubgeneinfo$gene_id, turquoise_geneinfo$gene_id),]
dim(turquoise_hubgeneinfo) 
# 265 genes

head(turquoise_hubgeneinfo)
min(abs(turquoise_hubgeneinfo$kME))




CytoscapeInput_edges_turquoise = read.table("/cyt/CytoscapeInput-edges-turquoise.txt", header = T)
CytoscapeInput_nodes_turquoise = read.table("/cyt/CytoscapeInput-nodes-turquoise.txt", sep = "\t", header = T)


dim(CytoscapeInput_edges_turquoise) 
## 3772280  6

head(CytoscapeInput_edges_turquoise)
dim(CytoscapeInput_nodes_turquoise) 
## 3103 5

head(CytoscapeInput_nodes_turquoise)
fivenum(CytoscapeInput_edges_turquoise$weight) 
## 0.02-0.18690548


CYP3A4_info = CytoscapeInput_nodes_turquoise[match("CYP3A4", CytoscapeInput_nodes_turquoise$altName), ] 
CYP3A4_info

CYP3A5_info = CytoscapeInput_nodes_turquoise[match("CYP3A5", CytoscapeInput_nodes_turquoise$altName), ] 
CYP3A5_info

CYP3A7_info = CytoscapeInput_nodes_turquoise[match("CYP3A7", CytoscapeInput_nodes_turquoise$altName), ] 
CYP3A7_info



#----------------- get genes directly co-expressed with CYP3A ------------------

turquoise_hubgene_edegs = CytoscapeInput_edges_turquoise[(CytoscapeInput_edges_turquoise$fromNode %in% turquoise_hubgeneinfo$gene_id) | (CytoscapeInput_edges_turquoise$toNode %in% turquoise_hubgeneinfo$gene_id), ]
dim(turquoise_hubgene_edegs)  
## 787050 edges

head(turquoise_hubgene_edegs)

#saveRDS(turquoise_hubgene_edegs, file = "/turquoise_hubgene_edegs.rds")
#write.csv(turquoise_hubgene_edegs, file ="/turquoise_hubgene_edegs.csv", row.names = F)

CYP3A4_hubgene_edegs = turquoise_hubgene_edegs[(turquoise_hubgene_edegs$fromNode %in% CYP3A4_info$nodeName) | (turquoise_hubgene_edegs$toNode %in% CYP3A4_info$nodeName), ]
dim(CYP3A4_hubgene_edegs) 
# 265 genes

CYP3A4_hubgene_edegs = CYP3A4_hubgene_edegs[order(CYP3A4_hubgene_edegs$weight, decreasing = T), ]
head(CYP3A4_hubgene_edegs)

CYP3A5_hubgene_edegs = turquoise_hubgene_edegs[(turquoise_hubgene_edegs$fromNode %in% CYP3A5_info$nodeName) | (turquoise_hubgene_edegs$toNode %in% CYP3A5_info$nodeName), ]
dim(CYP3A5_hubgene_edegs) 
# 265 genes

CYP3A5_hubgene_edegs = CYP3A5_hubgene_edegs[order(CYP3A5_hubgene_edegs$weight, decreasing = T), ]
head(CYP3A5_hubgene_edegs)

CYP3A7_hubgene_edegs = turquoise_hubgene_edegs[(turquoise_hubgene_edegs$fromNode %in% CYP3A7_info$nodeName) | (turquoise_hubgene_edegs$toNode %in% CYP3A7_info$nodeName), ]
dim(CYP3A7_hubgene_edegs) 
# 265 genes

CYP3A7_hubgene_edegs = CYP3A7_hubgene_edegs[order(CYP3A7_hubgene_edegs$weight, decreasing = T), ]
head(CYP3A7_hubgene_edegs)


if(F)
  {
  saveRDS(CYP3A4_hubgene_edegs, file = "/Allhubgene_CYP3A4.5.7/CYP3A4_hubgene_edegs-turquoise.rds")
  write.csv(CYP3A4_hubgene_edegs, file =  "/Allhubgene_CYP3A4.5.7/CYP3A4_hubgene_edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A5_hubgene_edegs, file =  "/Allhubgene_CYP3A4.5.7/CYP3A5_hubgene_edegs-turquoise.rds")
  write.csv(CYP3A5_hubgene_edegs, file =  "/Allhubgene_CYP3A4.5.7/CYP3A5_hubgene_edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A7_hubgene_edegs, file =  "/Allhubgene_CYP3A4.5.7/CYP3A7_hubgene_edegs-turquoise.rds")
  write.csv(CYP3A7_hubgene_edegs, file = "/Allhubgene_CYP3A4.5.7/CYP3A7_hubgene_edegs-turquoise.csv", row.names = F)
  }






#------------------------------ get all TFs ------------------------------------

table(turquoise_hubgeneinfo$org.Hs.eg.db_Symbol %in% humanTF$Symbol)  
# 15 TFs

TF = intersect(turquoise_hubgeneinfo$org.Hs.eg.db_Symbol, humanTF$Symbol)
TF 
##  ZNF385B    ZFP1     NFIA  ARID3C   ZGPAT 
## SLC2A4RG   KLF12       AR   NR3C2   PATZ1 
##    CREB3  ARNTL2  GATAD2A   TEAD4  ZNF33B

turquoisehubgene_TF = turquoise_hubgeneinfo[match(TF, turquoise_hubgeneinfo$org.Hs.eg.db_Symbol), ]
turquoisehubgene_TF$type = rep("humanTF", nrow(turquoisehubgene_TF))
dim(turquoisehubgene_TF) 
# 15 7

head(turquoisehubgene_TF)

#saveRDS(turquoisehubgene_TF, file = "/TF/turquoisehubgene_TF.rds")
#write.csv(turquoisehubgene_TF, file = "/TF/turquoisehubgene_TF.csv", row.names = F)



TF_edegs = CytoscapeInput_edges_turquoise[(CytoscapeInput_edges_turquoise$fromNode %in% turquoisehubgene_TF$gene_id) | (CytoscapeInput_edges_turquoise$toNode %in% turquoisehubgene_TF$gene_id), ]
dim(TF_edegs)  
## 46425 edges with above TFs

head(TF_edegs)

TF_edegs_gene = c(TF_edegs$fromNode, TF_edegs$toNode)
length(TF_edegs_gene) 
## 92850

head(TF_edegs_gene)
TF_edegs_gene_uniq = unique(TF_edegs_gene)
length(TF_edegs_gene_uniq) 
## 3103 edges with above TFs

CYP3A4_TF_edegs = TF_edegs[(TF_edegs$fromNode %in% CYP3A4_info$nodeName) | (TF_edegs$toNode %in% CYP3A4_info$nodeName), ]
dim(CYP3A4_TF_edegs) 
## 5 TFs

CYP3A4_TF_edegs = CYP3A4_TF_edegs[order(CYP3A4_TF_edegs$weight, decreasing = T), ]
head(CYP3A4_TF_edegs)


CYP3A5_TF_edegs = TF_edegs[(TF_edegs$fromNode %in% CYP3A5_info$nodeName) | (TF_edegs$toNode %in% CYP3A5_info$nodeName), ]
dim(CYP3A5_TF_edegs) 
## 5 TFs 

CYP3A5_TF_edegs = CYP3A5_TF_edegs[order(CYP3A5_TF_edegs$weight, decreasing = T), ]
head(CYP3A5_TF_edegs)

CYP3A7_TF_edegs = TF_edegs[(TF_edegs$fromNode %in% CYP3A7_info$nodeName) | (TF_edegs$toNode %in% CYP3A7_info$nodeName), ]
dim(CYP3A7_TF_edegs) 
## 5 TFs

CYP3A7_TF_edegs = CYP3A7_TF_edegs[order(CYP3A7_TF_edegs$weight, decreasing = T), ]
head(CYP3A7_TF_edegs)

## integrate three CYP3A co-expressed TFs
CYP3A4.5.7_TF_edegs = rbind(CYP3A4_TF_edegs, CYP3A5_TF_edegs, CYP3A7_TF_edegs)
dim(CYP3A4.5.7_TF_edegs) 
## 15 edges

head(CYP3A4.5.7_TF_edegs)


if(F)
  {
  saveRDS(CYP3A4_TF_edegs, file = "/TF/CYP3A4_TF-edegs-turquoise.rds")
  write.csv(CYP3A4_TF_edegs, file = "/TF/CYP3A4_TF-edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A5_TF_edegs, file ="/TF/CYP3A5_TF-edegs-turquoise.rds")
  write.csv(CYP3A5_TF_edegs, file = "/TF/CYP3A5_TF-edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A7_TF_edegs, file =  "/TF/CYP3A7_TF-edegs-turquoise.rds")
  write.csv(CYP3A7_TF_edegs, file = "/TF/CYP3A7_TF-edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A4.5.7_TF_edegs, file = "/TF/CYP3A4.5.7_TF_edegs.rds")
  write.csv(CYP3A4.5.7_TF_edegs, file = "/TF/CYP3A4.5.7_TF_edegs.csv", row.names = F)
  }



#------------------------------ get all lncRNAs --------------------------------

head(Gencode_lncRNA)
table(turquoise_hubgeneinfo$gene_id %in% Gencode_lncRNA$gene_id)  
## 15 lncRNAs

lncRNA = intersect(turquoise_hubgeneinfo$gene_id, Gencode_lncRNA$gene_id)
lncRNA

turquoisehubgene_lncRNA = turquoise_hubgeneinfo[match(lncRNA, turquoise_hubgeneinfo$gene_id), ]
turquoisehubgene_lncRNA$type = rep("lncRNA", nrow(turquoisehubgene_lncRNA))
head(turquoisehubgene_lncRNA)

turquoisehubgene_lncRNA$gene_type = Gencode_lncRNA[match(turquoisehubgene_lncRNA$gene_id, Gencode_lncRNA$gene_id), "gene_type"]
head(turquoisehubgene_lncRNA)

#saveRDS(turquoisehubgene_lncRNA, file = "/lncRNA/turquoisehubgene_lncRNA.rds")
#write.csv(turquoisehubgene_lncRNA, file = "/lncRNA/turquoisehubgene_lncRNA.csv", row.names = F)


lncRNA_edegs = CytoscapeInput_edges_turquoise[(CytoscapeInput_edges_turquoise$fromNode %in% turquoisehubgene_lncRNA$gene_id) | (CytoscapeInput_edges_turquoise$toNode %in% turquoisehubgene_lncRNA$gene_id), ]
dim(lncRNA_edegs)  
## 46425 edges

head(lncRNA_edegs)


lncRNA_edegs_gene = c(lncRNA_edegs$fromNode, lncRNA_edegs$toNode)
length(lncRNA_edegs_gene) 
## 92850

head(lncRNA_edegs_gene)
lncRNA_edegs_gene_uniq = unique(lncRNA_edegs_gene)
length(lncRNA_edegs_gene_uniq) 
## 3103 genes


CYP3A4_lncRNA_edegs = lncRNA_edegs[(lncRNA_edegs$fromNode %in% CYP3A4_info$nodeName) | (lncRNA_edegs$toNode %in% CYP3A4_info$nodeName), ]
dim(CYP3A4_lncRNA_edegs) 
## 15 lncRNAs

CYP3A4_lncRNA_edegs = CYP3A4_lncRNA_edegs[order(CYP3A4_lncRNA_edegs$weight, decreasing = T), ]
head(CYP3A4_lncRNA_edegs)

CYP3A5_lncRNA_edegs = lncRNA_edegs[(lncRNA_edegs$fromNode %in% CYP3A5_info$nodeName) | (lncRNA_edegs$toNode %in% CYP3A5_info$nodeName), ]
dim(CYP3A5_lncRNA_edegs) 
## 15 lncRNAs

CYP3A5_lncRNA_edegs = CYP3A5_lncRNA_edegs[order(CYP3A5_lncRNA_edegs$weight, decreasing = T), ]
head(CYP3A5_lncRNA_edegs)

CYP3A7_lncRNA_edegs = lncRNA_edegs[(lncRNA_edegs$fromNode %in% CYP3A7_info$nodeName) | (lncRNA_edegs$toNode %in% CYP3A7_info$nodeName), ]
dim(CYP3A7_lncRNA_edegs) 
## 15 lncRNAs

CYP3A7_lncRNA_edegs = CYP3A7_lncRNA_edegs[order(CYP3A7_lncRNA_edegs$weight, decreasing = T), ]
head(CYP3A7_lncRNA_edegs)


## integrate
CYP3A4.5.7_lncRNA_edegs = rbind(CYP3A4_lncRNA_edegs, CYP3A5_lncRNA_edegs, CYP3A7_lncRNA_edegs)
dim(CYP3A4.5.7_lncRNA_edegs) 
## 45 edges

head(CYP3A4.5.7_lncRNA_edegs)



if(F)
  {
  saveRDS(CYP3A4_lncRNA_edegs, file = "/lncRNA/CYP3A4_lncRNA-edegs-turquoise.rds")
  write.csv(CYP3A4_lncRNA_edegs, file = "/lncRNA/CYP3A4_lncRNA-edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A5_lncRNA_edegs, file = "/lncRNA/CYP3A5_lncRNA-edegs-turquoise.rds")
  write.csv(CYP3A5_lncRNA_edegs, file =  "/lncRNA/CYP3A5_lncRNA-edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A7_lncRNA_edegs, file = "/lncRNA/CYP3A7_lncRNA-edegs-turquoise.rds")
  write.csv(CYP3A7_lncRNA_edegs, file =  "/lncRNA/CYP3A7_lncRNA-edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A4.5.7_lncRNA_edegs, file = "/lncRNA/CYP3A4.5.7_lncRNA-edegs-turquoise.rds")
  write.csv(CYP3A4.5.7_lncRNA_edegs, file =  "/lncRNA/CYP3A4.5.7_lncRNA-edegs-turquoise.csv", row.names = F)
  }


#------------------------------ get all miRNAs ---------------------------------

table(turquoise_hubgeneinfo$gene_id %in% Gencode_miRNA$gene_id)  
## 1 miRNA

miRNA = intersect(turquoise_hubgeneinfo$gene_id, Gencode_miRNA$gene_id)
miRNA 
## MIR135A1


turquoisehubgene_miRNA = turquoise_hubgeneinfo[match(miRNA, turquoise_hubgeneinfo$gene_id), ]
turquoisehubgene_miRNA$type = rep("miRNA", nrow(turquoisehubgene_miRNA))
head(turquoisehubgene_miRNA)

#saveRDS(turquoisehubgene_miRNA, file = "miRNA/turquoisehubgene_miRNA.rds")
#write.csv(turquoisehubgene_miRNA, file =  "/miRNA/turquoisehubgene_miRNA.csv", row.names = F)


miRNA_edegs = CytoscapeInput_edges_turquoise[(CytoscapeInput_edges_turquoise$fromNode %in% turquoisehubgene_miRNA$gene_id) | (CytoscapeInput_edges_turquoise$toNode %in% turquoisehubgene_miRNA$gene_id), ]
dim(miRNA_edegs)  
## 3102 genes

head(miRNA_edegs)

miRNA_edegs_gene = c(miRNA_edegs$fromNode, miRNA_edegs$toNode)
length(miRNA_edegs_gene) 
## 6204

head(miRNA_edegs_gene)
miRNA_edegs_gene_uniq = unique(miRNA_edegs_gene)
length(miRNA_edegs_gene_uniq) 
##3103 


CYP3A4_miRNA_edegs = miRNA_edegs[(miRNA_edegs$fromNode %in% CYP3A4_info$nodeName) | (miRNA_edegs$toNode %in% CYP3A4_info$nodeName), ]
dim(CYP3A4_miRNA_edegs) 
## 1 miRNA
## MIR135A1

CYP3A4_miRNA_edegs = CYP3A4_miRNA_edegs[order(CYP3A4_miRNA_edegs$weight, decreasing = T), ]
head(CYP3A4_miRNA_edegs)


CYP3A5_miRNA_edegs = miRNA_edegs[(miRNA_edegs$fromNode %in% CYP3A5_info$nodeName) | (miRNA_edegs$toNode %in% CYP3A5_info$nodeName), ]
dim(CYP3A5_miRNA_edegs) 
## 1 miRNA
## MIR135A1

CYP3A5_miRNA_edegs = CYP3A5_miRNA_edegs[order(CYP3A5_miRNA_edegs$weight, decreasing = T), ]
head(CYP3A5_miRNA_edegs)


CYP3A7_miRNA_edegs = miRNA_edegs[(miRNA_edegs$fromNode %in% CYP3A7_info$nodeName) | (miRNA_edegs$toNode %in% CYP3A7_info$nodeName), ]
dim(CYP3A7_miRNA_edegs) 
## 1 miRNA
## MIR135A1

CYP3A7_miRNA_edegs = CYP3A7_miRNA_edegs[order(CYP3A7_miRNA_edegs$weight, decreasing = T), ]
head(CYP3A7_miRNA_edegs)

## integrate
CYP3A4.5.7_miRNA_edegs = rbind(CYP3A4_miRNA_edegs, CYP3A5_miRNA_edegs, CYP3A7_miRNA_edegs)
dim(CYP3A4.5.7_miRNA_edegs) 
## 3 edges

head(CYP3A4.5.7_miRNA_edegs)


if(F)
  {
  saveRDS(CYP3A4_miRNA_edegs, file = "/miRNA/CYP3A4_miRNA-edegs-turquoise.rds")
  write.csv(CYP3A4_miRNA_edegs, file =  "/miRNA/CYP3A4_miRNA-edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A5_miRNA_edegs, file = "/miRNA/CYP3A5_miRNA-edegs-turquoise.rds")
  write.csv(CYP3A5_miRNA_edegs, file =  "/miRNA/CYP3A5_miRNA-edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A7_miRNA_edegs, file = "/miRNA/CYP3A7_miRNA-edegs-turquoise.rds")
  write.csv(CYP3A7_miRNA_edegs, file = "/miRNA/CYP3A7_miRNA-edegs-turquoise.csv", row.names = F)
  saveRDS(CYP3A4.5.7_miRNA_edegs, file =  "/miRNA/CYP3A4.5.7_miRNA-edegs-turquoise.rds")
  write.csv(CYP3A4.5.7_miRNA_edegs, file =  "/miRNA/CYP3A4.5.7_miRNA-edegs-turquoise.csv", row.names = F)
  }




#------------------------- annotation of all genes type ------------------------

turquoisehubgene_lncRNA = readRDS("/lncRNA/turquoisehubgene_lncRNA.rds")
turquoisehubgene_miRNA = readRDS("/miRNA/turquoisehubgene_miRNA.rds")
turquoisehubgene_TF = readRDS("/TF/turquoisehubgene_TF.rds")

dim(turquoisehubgene_lncRNA) 
# 15 lncRNAs

dim(turquoisehubgene_miRNA) 
# 1 miRNAs

dim(turquoisehubgene_TF) 
# 15 TFs

head(turquoisehubgene_lncRNA)
head(turquoisehubgene_miRNA)
head(turquoisehubgene_TF)

turquoisehubgene_lncRNA$gene_type = NULL
head(turquoisehubgene_lncRNA)
CYP3A4.5.7 = turquoise_geneinfo[match(c("CYP3A4", "CYP3A5", "CYP3A7"), turquoise_geneinfo$Description), ]
CYP3A4.5.7$type = "CYP3A"

turquoisehubgene_TF_lncRNA_miRNA = rbind(turquoisehubgene_TF, turquoisehubgene_lncRNA, turquoisehubgene_miRNA)
dim(turquoisehubgene_TF_lncRNA_miRNA) 
## 31 genes
## 15 TFs
## 15 lncRNAs
## 1 miRNA

turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7 = rbind(CYP3A4.5.7, turquoisehubgene_TF_lncRNA_miRNA)
head(turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7)
dim(turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7) 
## 34 genes


#saveRDS(turquoisehubgene_TF_lncRNA_miRNA, file = "turquoisehubgene_TF_lncRNA_miRNA.rds")
#write.csv(turquoisehubgene_TF_lncRNA_miRNA, file = "turquoisehubgene_TF_lncRNA_miRNA.csv", row.names = F)
#saveRDS(turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7, file = "turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7.rds")
#write.csv(turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7, file =  "turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7.csv", row.names = F)



library(clusterProfiler)
library(org.Hs.eg.db)

keytypes(org.Hs.eg.db)
table(!is.na(turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$org.Hs.eg.db_Entrezid)) 
## 23 genes have entrezid


GO = enrichGO(gene = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$org.Hs.eg.db_Entrezid,
              OrgDb = 'org.Hs.eg.db',
              ont = "ALL",   
              pAdjustMethod = "BH",
              pvalueCutoff = 1,  
              qvalueCutoff = 1,
              keyType = 'ENTREZID',
              readable = T)

dim(GO) 
# 289


library(stringr)

GO_geneID = unlist(str_split(GO$geneID, "/"))  
GO_geneID_uniq = GO_geneID[!duplicated(GO_geneID)] 
length(GO_geneID_uniq) 
GO_geneID_uniq


KEGG = enrichKEGG(gene = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$org.Hs.eg.db_Entrezid,
                  organism = "hsa",
                  pvalueCutoff = 1, 
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 1)
dim(KEGG) 
## 47

kk = setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(kk)

KEGG_geneID = unlist(str_split(kk$geneID, "/"))  
KEGG_geneID_uniq = KEGG_geneID[!duplicated(KEGG_geneID)] 
length(KEGG_geneID_uniq) 
## 10 genes 

KEGG_geneID_uniq 
# ZFP1  ZNF33B  TEAD4  AR  CREB3  NR3C2  MIR135A1


#---------------------------------- liver --------------------------------------

hub_TF_lncRNA_miRNA_CYP3A4.5.7 = read.csv(paste0(step5_input_dir, "CYP3A_TF_lncRNA_miRNA.csv"))
turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7 = hub_TF_lncRNA_miRNA_CYP3A4.5.7[hub_TF_lncRNA_miRNA_CYP3A4.5.7$tissue == "liver", ]
dim(turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7) 
## 34  genes 
## 15    TFs
## 15 lncRNs
##  1  miRNA

head(turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7)

Liver_exprs = readRDS(paste0(step1_output_dir, "Liver_adjust.rds"))
dim(Liver_exprs) 
## 22304 193

Liver_exprs[1:3, 1:3]
Liver_exprs = t(Liver_exprs) 
Liver_exprs[1:3, 1:3]


#------------- get all genes exp matrix co-expressed with CYP3As ---------------

TF_gene = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7[turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$type == "humanTF",]
TF_gene 
## 5 TFs

turquoisehubgene_TF_exprs = Liver_exprs[, match(TF_gene$gene_id, colnames(Liver_exprs))]
dim(turquoisehubgene_TF_exprs) 
## 193 15

turquoisehubgene_TF_exprs[1:3, ]
colnames(turquoisehubgene_TF_exprs)
colnames(turquoisehubgene_TF_exprs) = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$Description[match(colnames(turquoisehubgene_TF_exprs),   
                                                                                                    turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$gene_id)] 
colnames(turquoisehubgene_TF_exprs)





lncRNA_gene = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7[turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$type == "lncRNA",]
lncRNA_gene 
## 15 lncRNAs

turquoisehubgene_lncRNA_exprs = Liver_exprs[, match(lncRNA_gene$gene_id, colnames(Liver_exprs))]
dim(turquoisehubgene_lncRNA_exprs) 
## 193 15

turquoisehubgene_lncRNA_exprs[1:3, 1:3]
colnames(turquoisehubgene_lncRNA_exprs)
colnames(turquoisehubgene_lncRNA_exprs) = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$Description[match(colnames(turquoisehubgene_lncRNA_exprs),   
                                                                                                        turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$gene_id)] 
colnames(turquoisehubgene_lncRNA_exprs)




miRNA_gene = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7[turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$type == "miRNA",]
miRNA_gene 
## 1 miRNA: MIR135A1

turquoisehubgene_miRNA_exprs = Liver_exprs[, match(miRNA_gene$gene_id, colnames(Liver_exprs))]
turquoisehubgene_miRNA_exprs = data.frame(turquoisehubgene_miRNA_exprs) 
head(turquoisehubgene_miRNA_exprs)

colnames(turquoisehubgene_miRNA_exprs) = miRNA_gene$Description
head(turquoisehubgene_miRNA_exprs)


CYP3A4_gene = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7[turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$Description == "CYP3A4",]
CYP3A4_gene 

turquoise_CYP3A4_exprs = Liver_exprs[, match(CYP3A4_gene$gene_id, colnames(Liver_exprs))]
length(turquoise_CYP3A4_exprs) 
## 193
head(turquoise_CYP3A4_exprs)

turquoise_CYP3A4_exprs = data.frame(turquoise_CYP3A4_exprs) 
head(turquoise_CYP3A4_exprs)

colnames(turquoise_CYP3A4_exprs) = "CYP3A4"
head(turquoise_CYP3A4_exprs)


CYP3A5_gene = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7[turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$Description == "CYP3A5",]
CYP3A5_gene 

turquoise_CYP3A5_exprs = Liver_exprs[, match(CYP3A5_gene$gene_id, colnames(Liver_exprs))]
length(turquoise_CYP3A5_exprs) 
## 193
turquoise_CYP3A5_exprs = data.frame(turquoise_CYP3A5_exprs) 
head(turquoise_CYP3A5_exprs)

colnames(turquoise_CYP3A5_exprs) = "CYP3A5"
head(turquoise_CYP3A5_exprs)


CYP3A7_gene = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7[turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7$Description == "CYP3A7",]
CYP3A7_gene 

turquoise_CYP3A7_exprs = Liver_exprs[, match(CYP3A7_gene$gene_id, colnames(Liver_exprs))]
length(turquoise_CYP3A7_exprs) 
## 193  

turquoise_CYP3A7_exprs = data.frame(turquoise_CYP3A7_exprs) 
head(turquoise_CYP3A7_exprs)

colnames(turquoise_CYP3A7_exprs) = "CYP3A7"
head(turquoise_CYP3A7_exprs)


### integration of exp matrix
turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7_exprs = cbind(CYP3A4 = turquoise_CYP3A4_exprs,
                                                          CYP3A5 = turquoise_CYP3A5_exprs,
                                                          CYP3A7 = turquoise_CYP3A7_exprs,
                                                          turquoisehubgene_TF_exprs,
                                                          turquoisehubgene_lncRNA_exprs, 
                                                          MIR135A1 = turquoisehubgene_miRNA_exprs)





dim(turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7_exprs)  
## 193 32


colnames(turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7_exprs)



#-------------------------- analysis of relationship ---------------------------

dat = turquoisehubgene_TF_lncRNA_miRNA_CYP3A4.5.7_exprs
colnames(dat)
ncol(dat)

library(WGCNA)

dat_bicor = bicorAndPvalue(dat) 

correlations = dat_bicor$bicor 

p_values = dat_bicor$p 

apply(correlations, 2, function(x) range(abs(x))) 


## correlations and p values
textMatrix = paste(signif(correlations, 2))
dim(textMatrix) = dim(correlations)


## heatmap

par(mar= c(6, 7.5, 3, 3))
labeledHeatmap(Matrix= correlations, 
               xLabels= colnames(dat), 
               yLabels= colnames(dat), 
               ySymbols= colnames(dat), 
               colorLabels= FALSE, 
               colors= blueWhiteRed(50), 
               textMatrix= textMatrix, 
               setStdMargins= FALSE, 
               cex.text= 0.5, 
               zlim= c(-1,1), 
               main= paste("TFs-lncRNAs-miRNA-CYP3A4-CYP3A5-CYP3A7 relationship"))


#-------------------------------- small intestine ------------------------------

hub_TF_lncRNA_miRNA_CYP3A4.5.7 = read.csv("CYP3A_TF_lncRNA_miRNA.csv")
smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7 = hub_TF_lncRNA_miRNA_CYP3A4.5.7[hub_TF_lncRNA_miRNA_CYP3A4.5.7$tissue == "small intestine", ]
dim(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7) 
## 184   genes
## 126     TFs
##  57 lncRNAs
##   1  miRNAs

table(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$type)

SmallIntestine_exprs = readRDS(paste0(step1_output_dir, "SmallIntestine_adjust.rds"))
dim(SmallIntestine_exprs) 
## 25464 175

SmallIntestine_exprs[1:3, 1:3]
SmallIntestine_exprs = t(SmallIntestine_exprs) Òò
SmallIntestine_exprs[1:3, 1:3]



#------------------------------ get exp matrix ---------------------------------

TF_gene = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7[smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$type == "humanTF",]
TF_gene 
## 126 TFs

smallIntestinehub_TF_exprs = SmallIntestine_exprs[, match(TF_gene$gene_id, colnames(SmallIntestine_exprs))]
dim(smallIntestinehub_TF_exprs) 
## 175 126

smallIntestinehub_TF_exprs[1:3, ]
colnames(smallIntestinehub_TF_exprs)
colnames(smallIntestinehub_TF_exprs) = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$Description[match(colnames(smallIntestinehub_TF_exprs),   
                                                                                                       smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$gene_id)] 
colnames(smallIntestinehub_TF_exprs)



lncRNA_gene = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7[smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$type == "lncRNA",]
lncRNA_gene 
## 57 lncRNAs

smallIntestinehub_lncRNA_exprs = SmallIntestine_exprs[, match(lncRNA_gene$gene_id, colnames(SmallIntestine_exprs))]
dim(smallIntestinehub_lncRNA_exprs) 
## 175 57

smallIntestinehub_lncRNA_exprs[1:3, 1:3]
colnames(smallIntestinehub_lncRNA_exprs)
colnames(smallIntestinehub_lncRNA_exprs) = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$Description[match(colnames(smallIntestinehub_lncRNA_exprs),   
                                                                                                           smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$gene_id)] 
colnames(smallIntestinehub_lncRNA_exprs)




miRNA_gene = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7[smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$type == "miRNA",]
miRNA_gene 
## 1¸ömiRNA  MIR621

smallIntestinehub_miRNA_exprs = SmallIntestine_exprs[, match(miRNA_gene$gene_id, colnames(SmallIntestine_exprs))]
smallIntestinehub_miRNA_exprs = data.frame(smallIntestinehub_miRNA_exprs) 
head(smallIntestinehub_miRNA_exprs) 

colnames(smallIntestinehub_miRNA_exprs) = miRNA_gene$Description
head(smallIntestinehub_miRNA_exprs)


CYP3A4_gene = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7[smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$Description == "CYP3A4",]
CYP3A4_gene

smallIntestine_CYP3A4_exprs = SmallIntestine_exprs[, match(CYP3A4_gene$gene_id, colnames(SmallIntestine_exprs))]
length(smallIntestine_CYP3A4_exprs) 
## 175
head(smallIntestine_CYP3A4_exprs)

smallIntestine_CYP3A4_exprs = data.frame(smallIntestine_CYP3A4_exprs) 
head(smallIntestine_CYP3A4_exprs)
colnames(smallIntestine_CYP3A4_exprs) = "CYP3A4"
head(smallIntestine_CYP3A4_exprs)


CYP3A5_gene = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7[smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$Description == "CYP3A5",]
CYP3A5_gene 

smallIntestine_CYP3A5_exprs = SmallIntestine_exprs[, match(CYP3A5_gene$gene_id, colnames(SmallIntestine_exprs))]
length(smallIntestine_CYP3A5_exprs) 
## 193

smallIntestine_CYP3A5_exprs = data.frame(smallIntestine_CYP3A5_exprs) 
head(smallIntestine_CYP3A5_exprs)
colnames(smallIntestine_CYP3A5_exprs) = "CYP3A5"
head(smallIntestine_CYP3A5_exprs)


CYP3A7_gene = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7[smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$Description == "CYP3A7",]
CYP3A7_gene 
## CYP3A7

smallIntestine_CYP3A7_exprs = SmallIntestine_exprs[, match(CYP3A7_gene$gene_id, colnames(SmallIntestine_exprs))]
length(smallIntestine_CYP3A7_exprs) 
## 193 

smallIntestine_CYP3A7_exprs = data.frame(smallIntestine_CYP3A7_exprs) 

head(smallIntestine_CYP3A7_exprs)
colnames(smallIntestine_CYP3A7_exprs) = "CYP3A7"
head(smallIntestine_CYP3A7_exprs)


## integration
smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4_exprs = cbind(smallIntestine_CYP3A4_exprs, 
                                                        smallIntestinehub_TF_exprs,
                                                        smallIntestinehub_lncRNA_exprs, 
                                                        smallIntestinehub_miRNA_exprs,
                                                        MIR621 = smallIntestinehub_miRNA_exprs)

smallIntestine_hub_TF_lncRNA_miRNA_CYP3A5_exprs = cbind(CYP3A5 = smallIntestine_CYP3A5_exprs, 
                                                        smallIntestinehub_TF_exprs,
                                                        smallIntestinehub_lncRNA_exprs, 
                                                        MIR621 = smallIntestinehub_miRNA_exprs)

smallIntestine_hub_TF_lncRNA_miRNA_CYP3A7_exprs = cbind(CYP3A7 = smallIntestine_CYP3A7_exprs, 
                                                        smallIntestinehub_TF_exprs,
                                                        smallIntestinehub_lncRNA_exprs)


smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs = cbind(CYP3A4 = smallIntestine_CYP3A4_exprs,
                                                            CYP3A5 = smallIntestine_CYP3A5_exprs,
                                                            CYP3A7 = smallIntestine_CYP3A7_exprs,
                                                            smallIntestinehub_TF_exprs,
                                                            smallIntestinehub_lncRNA_exprs, 
                                                            MIR621 = smallIntestinehub_miRNA_exprs)




dim(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4_exprs) 
## 175 185

dim(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A5_exprs) 
## 175 185

dim(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A7_exprs) 
## 175 185

dim(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs)  
## 175 187

colnames(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4_exprs)
colnames(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A5_exprs)
colnames(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A7_exprs)
colnames(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs)




dat = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs
colnames(dat)
ncol(dat)

dat_bicor = bicorAndPvalue(dat) 
correlations = dat_bicor$bicor 
p_values = dat_bicor$p 
apply(correlations, 2, function(x) range(abs(x))) 


## correlations and p values

textMatrix = paste(signif(correlations, 2), "\n(", 
                   signif(p_values, 1), ")", sep= "")
dim(textMatrix) = dim(correlations)


## heatmap

par(mar= c(6, 7.5, 3, 3))
labeledHeatmap(Matrix= correlations, 
               xLabels= colnames(dat), 
               yLabels= colnames(dat), 
               ySymbols= colnames(dat), 
               colorLabels= FALSE, 
               colors= blueWhiteRed(50), 
               textMatrix= textMatrix, 
               setStdMargins= FALSE, 
               cex.text= 0.5, 
               zlim= c(-1,1), 
               main= paste("TFs-lncRNAs-miRNA-CYP3A4-CYP3A5-CYP3A7 relationship"))




## module turquoise
turquoise_hubTF_lncRNA_miRNA_CYP3A4.5 = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7[smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$module == "turquoise", ]
dim(turquoise_hubTF_lncRNA_miRNA_CYP3A4.5) 
## 154 genes

head(turquoise_hubTF_lncRNA_miRNA_CYP3A4.5)
smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs[1:3, 1:3]
colnames(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs)[duplicated(colnames(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs))] 
turquiosedat = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs[, 
                                                                   match(turquoise_hubTF_lncRNA_miRNA_CYP3A4.5$Description, colnames(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs))]
colnames(turquiosedat)
ncol(turquiosedat)



turquiosedat_bicor = bicorAndPvalue(turquiosedat) 
correlations = turquiosedat_bicor$bicor 
p_values = turquiosedat_bicor$p 
apply(correlations, 2, function(x) range(abs(x))) 


## correlations and p values

textMatrix = paste(signif(correlations, 2), "\n(", 
                   signif(p_values, 1), ")", sep= "")
dim(textMatrix) = dim(correlations)


## heatmap

par(mar= c(6, 7.5, 3, 3))
labeledHeatmap(Matrix= correlations, 
               xLabels= colnames(turquiosedat), 
               yLabels= colnames(turquiosedat), 
               ySymbols= colnames(turquiosedat), 
               colorLabels= FALSE, 
               colors= blueWhiteRed(50), 
               textMatrix= textMatrix, 
               setStdMargins= FALSE, 
               cex.text= 0.5, 
               zlim= c(-1,1), 
               main= paste("TFs-lncRNAs-miRNA-CYP3A4-CYP3A5 relationship"))


## module red

red_hubTF_lncRNA_miRNA_CYP3A4.5 = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7[smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7$module == "red", ]
dim(red_hubTF_lncRNA_miRNA_CYP3A4.5) 
## 154 genes

head(red_hubTF_lncRNA_miRNA_CYP3A4.5)
red_hubTF_lncRNA_miRNA_CYP3A4.5[1:3, 1:3]
reddat = smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs[, 
                                                             match(red_hubTF_lncRNA_miRNA_CYP3A4.5$Description, colnames(smallIntestine_hub_TF_lncRNA_miRNA_CYP3A4.5.7_exprs))]
colnames(reddat)
ncol(reddat)



reddat_bicor = bicorAndPvalue(reddat) 
correlations = reddat_bicor$bicor 
p_values = reddat_bicor$p 
apply(correlations, 2, function(x) range(abs(x))) 


## correlations and p values

textMatrix = paste(signif(correlations, 2), "\n(", 
                   signif(p_values, 1), ")", sep= "")
dim(textMatrix) = dim(correlations)


## heatmap

par(mar= c(6, 7.5, 3, 3))
labeledHeatmap(Matrix= correlations, 
               xLabels= colnames(reddat), 
               yLabels= colnames(reddat), 
               ySymbols= colnames(reddat), 
               colorLabels= FALSE, 
               colors= blueWhiteRed(50), 
               textMatrix= textMatrix, 
               setStdMargins= FALSE, 
               cex.text= 0.5, 
               zlim= c(-1,1), 
               main= paste("TFs-lncRNAs-CYP3A7 relationship"))


#--------------------------- comparing figures ---------------------------------

#------------------------ with or without ncRNA --------------------------------

tog <- read.table("all_7_si.txt", header = T, sep = "\t")
tf <- read.table("TF_7_si.txt", header = T, sep = "\t")

#which(tf[, 5] == "CYP3A4")

tf_sort <- tf[which(tf[, 5] == "CYP3A7"), ]
tmp <- tf[-which(tf[, 5] == "CYP3A7"), ]
tmp[, 5] -> tmp[, 7]
colnames(tmp)[6] -> colnames(tmp)[7]
colnames(tmp)[5] -> colnames(tmp)[6]
tmp <- tmp[, -5]
tf_sort <- rbind(tf_sort, tmp)
tf <- tf_sort[order(tf_sort[, 3], decreasing = T), ]

tf_sort <- tf[which(tf[, 1] == "ENSG00000160870.12"), ]
tmp <- tf[-which(tf[, 1] == "ENSG00000160870.12"), ]
tmp[, 1] -> tmp[, 7]
tmp[, 2] -> tmp[, 1]
tmp[, 7] -> tmp[, 2]
tmp <- tmp[, -7]
tf_sort <- rbind(tf_sort, tmp)
tf_sort <- tf_sort[order(tf_sort[, 3], decreasing = T), ]

tog <- tog[which(tog[, 7] == "humanTF"), ]     
## " humanTF "

same <- merge(tog, tf_sort, by.x = "toNode", by.y = "toNode", all = F)
data <- rbind(tog[, c(3,6)], tf_sort[, c(3,6)])
data$group <- c(rep("TF and ncRNA", dim(tog)[1]), rep("TF", dim(tf_sort)[1]))
data$group <- as.factor(data$group)
same_sort <- c(same[,3], same[,9])
same_sort <- as.data.frame(same_sort)
same_sort$group <- c(rep("TF and ncRNA", dim(same)[1]), rep("TF", dim(same)[1]))
same_sort$group <- as.factor(same_sort$group)

data_4_si <- data
same_4_si <- same_sort
data_4_si$group_name <- rep("small intestine_CYP3A4", dim(data_4_si)[1])
data_4_si$group_name <- as.factor(data_4_si$group_name)
data <- rbind(data_4_li, data_4_si, data_5_li, data_5_si, data_7_li, data_7_si)
same_4_si$group_name <- rep("small intestine_CYP3A4", dim(same_4_si)[1])
same_4_si$group_name <- as.factor(same_4_si$group_name)
same_sort <- rbind(same_4_li, same_4_si, same_5_li, same_5_si, same_7_li, same_7_si)


ggplot(data = same_sort,
       mapping = aes(x = group_name, y = same_sort) ) +
  geom_violin(data = data,
              mapping = aes(x = group_name, y = weight, fill = group),
              trim = F,
              draw_quantiles = c(0.25, 0.5, 0.75),
              scale = "count") +
  
  geom_quasirandom(mapping = aes(col = group), dodge.width=.8 ) +
  
  scale_color_manual(values = c("red", "blue" )) +
  
  labs(x = "", y= "Weight") +
  
  guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
  
  theme_classic(base_size = 16) 
#theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data = same_sort,
       mapping = aes(x = group_name, y = same_sort) ) +
  
  geom_violin(data = data,
              mapping = aes(x = group_name, y = weight, fill = group),
              trim = F,
              draw_quantiles = c(0.25, 0.5, 0.75),
              scale = "count") +
  
  geom_quasirandom(mapping = aes(col = group), dodge.width=.8 ) +
  scale_color_manual(values = c("#FF3300", "#FF3300")) +
  labs(x = "", y= "Weight") +
  
  theme_classic(base_size = 15) +
  
  guides(fill = guide_legend(title = ""), col = FALSE)


#--------------------- liver vs small intestine --------------------------------

liver <- read.table("liver_3A7.txt", header = T, sep = "\t")
si <- read.table("si_3A7.txt", header = T, sep = "\t")

same <- merge(liver, si, by.x = "toNode", by.y = "toNode", all = F)
data <- rbind(liver[, c(3,6,7)], si[, c(3,6,7)])
data$group <- c(rep("liver_CYP3A7", dim(liver)[1]), 
                rep("small intestine_CYP3A7", dim(si)[1]))
data$group <- as.factor(data$group)

same_sort <- c(same[,3], same[,9] )
same_sort <- as.data.frame(same_sort)
same_sort$group <- c(rep("liver_CYP3A7", dim(same)[1]), 
                     rep("small intestine_CYP3A7", dim(same)[1]))
same_sort$group <- as.factor(same_sort$group)
same_sort$type <- c(same[, 7], same[,13] )
which(data[, 3] == "miRNA")
data[which(data[, 3] == "miRNA"),3] <- "lncRNA"
unique(data[, 3])
data[which(data[, 3] == " humanTF "), 3] <- "humanTF"
data[which(data[, 3] == " lncRNA "), 3] <- "lncRNA"
data[which(data[, 3] == " miRNA "), 3] <- "lncRNA"
which(same_sort[, 3] == "miRNA")
same_sort[which(same_sort[, 3] == "miRNA"),3] <- "lncRNA"
unique(same_sort[, 3])
same_sort[which(same_sort[, 3] == " humanTF "), 3] <- "humanTF"
same_sort[which(same_sort[, 3] == " lncRNA "), 3] <- "lncRNA"
same_sort[which(same_sort[, 3] == " miRNA "), 3] <- "lncRNA"


data_7 <- data
same_7 <- same_sort
#data_5$group_name <- rep("CYP3A5", dim(data_5)[1])
#data_5$group_name <- as.factor(data_5$group_name)
data <- rbind(data_4, data_5, data_7 )
#same_5$group_name <- rep("CYP3A5", dim(same_5)[1])
#same_5$group_name <- as.factor(same_5$group_name)
same_sort <- rbind(same_4, same_5, same_7)

same_sort$type <- as.factor(same_sort$type)
data$Type <- as.factor(data$Type)


ggplot(data = same_sort,
       mapping = aes(x = group_name, y = same_sort) ) +
  geom_violin(data = data,
              mapping = aes(x = group_name, y = weight, fill = group),
              trim = F,
              draw_quantiles = c(0.25, 0.5, 0.75),
              scale = "count") +
  scale_fill_manual(values = c(rep(c("#FF6666", "#00CC99"), 3)) ) +
  
  geom_quasirandom(mapping = aes(col = type), dodge.width=.8 ) +
  scale_color_manual(values = c("red", "blue" )) +
  labs(x = "", y= "Weight") +
  guides(#fill = guide_legend(title = ""), 
    fill = FALSE, 
    color = guide_legend(title = "")
  ) +
  
  theme_classic(base_size = 16)


ggplot() +
  geom_violin(data = data,
              mapping = aes(x = group_name, y = weight, fill = group),
              trim = F,
              draw_quantiles = c(0.25, 0.5, 0.75),
              scale = "count") +
  geom_jitter(data = same_sort,
              mapping = aes(x = group_name, y = same_sort, col = group),
              #col = "#FF3300",
              width = 0.3) +
  labs(x = "", y= "Weight") +
  
  theme_classic(base_size = 15) +
  guides(fill = FALSE) 


#---------------------------- between CYP3A4, 5, 7 -----------------------------

#------------------------------------- liver -----------------------------------

data_4 <- read.table("all_4_liver.txt", header = T, sep = "\t")
data_5 <- read.table("all_5_liver.txt", header = T, sep = "\t")
data_7 <- read.table("all_7_liver.txt", header = T, sep = "\t")

same_45 <- merge(data_4, data_5, by.x = "toNode", by.y = "toNode", all = F)
same_sort45 <- c(same_45[, 3], same_45[, 9])
same_sort45 <- as.data.frame(same_sort45)
same_sort45$type <- c(same_45[, 7], same_45[, 7])
colnames(same_sort45)[1] <- "weight"

same_47 <- merge(data_4, data_7, by.x = "toNode", by.y = "toNode", all = F)
same_sort47 <- c(same_47[, 3], same_47[, 9])
same_sort47 <- as.data.frame(same_sort47)
same_sort47$type <- c(same_47[, 7], same_47[, 7])
colnames(same_sort47)[1] <- "weight"

same_57 <- merge(data_5, data_7, by.x = "toNode", by.y = "toNode", all = F)
same_sort57 <- c(same_57[, 3], same_57[, 9])
same_sort57 <- as.data.frame(same_sort57)
same_sort57$type <- c(same_57[, 7], same_57[, 7])
colnames(same_sort57)[1] <- "weight"

data_liver <- rbind(data_4[, c(3,6,7)], data_5[, c(3,6,7)], data_7[, c(3,6,7)])
data_liver$group_name <- c(rep("liver_CYP3A4", dim(data_4)[1] ), 
                           rep("liver_CYP3A5", dim(data_5)[1]), 
                           rep("liver_CYP3A7", dim(data_7)[1]) )
which(data_liver$Type == "miRNA")
data_liver$Type[which(data_liver$Type == "miRNA")] <- "lncRNA"
unique(data_liver$Type)

same_sort_liver <- rbind(same_sort45, same_sort47, same_sort57)
same_sort_liver$group_name <- c(rep("liver_CYP3A4", dim(same_45)[1] ),
                                rep("liver_CYP3A5-1", dim(same_45)[1] ),
                                rep("liver_CYP3A4-1", dim(same_47)[1] ),
                                rep("liver_CYP3A7", dim(same_47)[1] ), 
                                rep("liver_CYP3A5", dim(same_57)[1] ),
                                rep("liver_CYP3A7-1", dim(same_57)[1] ) )
which(same_sort_liver$type == "miRNA")
same_sort_liver$type[which(same_sort_liver$type == "miRNA")] <- "lncRNA"
unique(same_sort_liver$type)

#--------------------------------- small intestine -----------------------------

data_4 <- read.table("all_4_si.txt", header = T, sep = "\t")
data_5 <- read.table("all_5_si.txt", header = T, sep = "\t")
data_7 <- read.table("all_7_si.txt", header = T, sep = "\t")


same_45 <- merge(data_4, data_5, by.x = "toNode", by.y = "toNode", all = F)
same_sort45 <- c(same_45[, 3], same_45[, 9])
same_sort45 <- as.data.frame(same_sort45)
same_sort45$type <- c(same_45[, 7], same_45[, 7])
colnames(same_sort45)[1] <- "weight"

same_47 <- merge(data_4, data_7, by.x = "toAltName", by.y = "toAltName", all = F)
same_sort47 <- c(same_47[, 3], same_47[, 9])
same_sort47 <- as.data.frame(same_sort47)
same_sort47$type <- c(same_47[, 7], same_47[, 7])
colnames(same_sort47)[1] <- "weight"

same_57 <- merge(data_5, data_7, by.x = "toNode", by.y = "toNode", all = F)
same_sort57 <- c(same_57[, 3], same_57[, 9])
same_sort57 <- as.data.frame(same_sort57)
same_sort57$type <- c(same_57[, 7], same_57[, 7])
colnames(same_sort57)[1] <- "weight"

data_si <- rbind(data_4[, c(3,6,7)], data_5[, c(3,6,7)], data_7[, c(3,6,7)])
data_si$group_name <- c(   rep("small intestine_CYP3A4", dim(data_4)[1] ), 
                           rep("small intestine_CYP3A5", dim(data_5)[1]), 
                           rep("small intestine_CYP3A7", dim(data_7)[1])    )
which(data_si$Type == "miRNA")
data_si$Type[which(data_si$Type == "miRNA")] <- "lncRNA"
unique(data_si$Type)
data_si[which(data_si$Type == " humanTF "), 3] <- "humanTF"
data_si[which(data_si$Type == " lncRNA "), 3] <- "lncRNA"
data_si[which(data_si$Type == " miRNA "), 3] <- "lncRNA"

same_sort_si <- rbind(same_sort45, same_sort47, same_sort57)
same_sort_si$group_name <- c(   rep("small intestine_CYP3A4", dim(same_45)[1] ),
                                rep("small intestine_CYP3A5-1", dim(same_45)[1] ),
                                rep("small intestine_CYP3A4-1", dim(same_47)[1] ),
                                rep("small intestine_CYP3A7", dim(same_47)[1] ), 
                                rep("small intestine_CYP3A5", dim(same_57)[1] ),
                                rep("small intestine_CYP3A7-1", dim(same_57)[1] ) )
which(same_sort_si$type == "miRNA")
same_sort_si$type[which(same_sort_si$type == "miRNA")] <- "lncRNA"
unique(same_sort_si$type)

data <- rbind(data_liver, data_si)
data2 <- rep(data, 2)
data2 <- as.data.frame(data2)
data2$group_name <- paste( data$group_name, "-1", sep = "" )
data2$group_name <- as.factor(data2$group_name)
data_all <- rbind(data, data2[, 1:4])

same_sort <- rbind(same_sort_liver, same_sort_si)

data_all$Type <- as.factor(data$Type)
data_all$group_name <- as.factor(data$group_name)
same_sort$type <- as.factor(same_sort$type)
same_sort$group_name <- as.factor(same_sort$group_name)


ggplot(data = same_sort,
       mapping = aes(x = group_name, y = weight) ) +
  geom_violin(data = data_all,
              mapping = aes(x = group_name, y = weight, fill = group_name),
              trim = F,
              draw_quantiles = c(0.25, 0.5, 0.75),
              scale = "count") +
  
  scale_fill_manual(values = c("#66FFFF", "#FF33FF", "#66CC99", 
                               "#996600", "#336666", "#FFFF00",
                               "#66FFFF", "#FF33FF", "#66CC99", 
                               "#996600", "#336666", "#FFFF00"
  ) ) +
  
  geom_quasirandom(mapping = aes(col = type), dodge.width=1 ) +
  scale_color_manual(values = c("red", "blue" )) +
  labs(x = "", y= "Weight") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = guide_legend(title = ""), 
         #fill = FALSE, 
         color = guide_legend(title = "")  ) +
  
  theme_classic(base_size = 8)



dif_li <- tf_sort[-which(tf_sort$toNode %in% same$toNode), ]
dif_si <- tog[-which(tog$toNode %in% same$toNode), ]

write.csv(same[, c(1:7,9)], file = "CYP3A7_same_TF.csv")
#write.csv(same_7_si[, c(1:7,9)], file = "CYP3A7_same_TF.csv")
#write.csv(same_5_si[, c(1:7,9)], file = "CYP3A5_same_TF.csv")

write.csv(dif_li, file = "CYP3A7_si_ncRNA_TF.csv")
write.csv(dif_si, file = "CYP3A7_si_ncRNA_TF.csv")


