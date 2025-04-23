library(MASS) 
library(glmnet) 
library(multcomp) ## multiT, multinormal
library(poibin)
library(sgof)	## Beta-binomial SGoF procedure, 
              ## for correcting positive correlation

# files
[1] "CYPs.sorted"           ## response data frame                  
[2] "CYPs.sorted.m"         ## response matrix                  
[3] "input.genes"           ## inputs data frame                  
[4] "input.genes.m"         ## inputs matrix                  
[5] "insig.est.gs.with.noise"               ## above average picks from the pool of irrelevant genes and artificial noise 
[6] "scan.all.CYPs"                         ## model fitting without any artificial noise, gene-gene comparison   
[7] "scan.all.CYPs.against.noise"           ## signal genes tested against artificial noise 
[8] "scan.all.CYPs.against.noise.insig.genes"     ## irrelevant genes tested against artificial noise 
[9] "sig.est.gs"                             ## separation of signal genes and irrelavent genes 
[10] "sig.est.gs.with.noise"                  ## above average picks from the pool of signal genes and artificial noise 
 

names(scan.all.CYPs[[1]]) 
[1] "glm.coef.with.intercept" "glm.p.with.intercept"    
[3] "input.se.est"            "input.rho.est"           
[5] "est.gs"                  "est.gs.order"            
[7] "glm.p.order"       


## caculation of the numerator of the sobol main index
t.gs <- function(i, coef, mse, rhom) 
  {
     rho <- rhom[i,-i]
     return ((coef[i] + matrix(coef[-i], nr = 1) %*% (matrix(mse[-i], nc = 1) * matrix(rho, nc = 1)) / mse[i]) ^ 2 * mse[i] ^ 2)
  }

scan.one.CYP <- function(CYP, input.genes) 
  {
     ## fit the sample with GLM via Iteratively Reweighted Least Square (IRLS)
     
     fit <- glm(CYP ~ ., data = data.frame(input.genes))
     glm.coef <- fit$coef
     glm.p <- summary(fit)$coef[, 4]
  
     glm.p.order <- order(glm.p[-1], na.last = TRUE)
  
     ## calculate emprical estimate of sobol ranking using glm.coef
     
     coef.est <- glm.coef[-1]
     
     mse.est <-apply(input.genes, 2, 
                     function(x) {var(x, na.rm = T) ^ 0.5})
  
     rhom.est <-cor(input.genes, use = "pairwise.complete.obs")
  
     est.gs <-apply(matrix(1:length(coef.est), nc = 1), 1, 
                 function(i) {t.gs(i, coef.est, mse.est, rhom.est)})
  
     est.gs.order <-order(est.gs, decreasing = TRUE, na.last = TRUE)
  
     my.list <-list(
       glm.coef.with.intercept = glm.coef,
       glm.p.with.intercept = glm.p,
       input.se.est = mse.est,
       input.rho.est = rhom.est,
       est.gs = est.gs,
       est.gs.order = est.gs.order,
       glm.p.order = glm.p.order )
     
     return(my.list)
}

SI <- function(y, input) 
  {
     fit.1 <-glm(y ~ ., data = data.frame(input))
     
     glm.coef <- fit.1$coef
  
     SI <-var(apply(input, 1, 
                    function(x) {matrix(x, nr = 1) %*% matrix(glm.coef[-1], nc = 1)}), na.rm = TRUE)
     ##return(SI)
     
     return(c( SI,
               fit.1$deviance / fit.1$null.deviance,
               fit.1$null.deviance ))
  
}

## k >=3
comb.m <- function(input, k = 3) 
  {
     d <- dim(input)[2]
     if (d >= 3) 
       {
       combn.index <-apply(matrix(2:(min(d, k)), nr = 1), 2, 
                           function(x) {combn(c(1:d), x)})
    
       combn.list <-lapply(combn.index, 
                           function(x) {cc <-input[, x[1,]]
                                        for (j in 2:dim(x)[1]) 
                                          {
                                          cc <- cc * input[, x[j,]]
                                          }
                           return(list(cc))
                                       })
       combn.input <-do.call(cbind, lapply(combn.list, function(x) {x[[1]]}))
    
       }
     
     if (d == 2) 
       {
       combn.input <- input[, 1] * input[, 2]
       }
     
     return(combn.input)
  
}

## k >=3
comb.m <- function(input, k = 3)
  {
  d <- dim(input)[2]
  if (d >= 3) 
    {
    combn.index <-apply(matrix(2:(min(d, k)), nr = 1), 2, 
                        function(x) { combn(c(1:d), x)  })
    
    combn.list <-lapply(combn.index, 
                        function(x) 
                          { 
                          cc <-input[, x[1,]]
                          for (j in 2:dim(x)[1]) 
                            {
                            cc <- cc * input[, x[j,]]
                            }
                          return(list(cc))
                          })
     combn.input <-do.call(cbind, lapply(combn.list, 
                                         function(x){ x[[1]] }))
    
  }
  if (d == 2) 
    {
    combn.input <- input[, 1] * input[, 2]
    }
   
  return(combn.input)
  
}

SI.comb <- function(y, input, k=3){ 
     main.new.m <- apply(input, 2, 
                         function(x){ apply(matrix(1:k),1, function(i){x^i}) }) 
     
     main.SI <- apply(main.new.m, 2, 
                      function(x){ SI(y,matrix(x, nr=dim(input)[1]))}) 
     
     comb.new.m <- cbind(matrix(main.new.m, nr=dim(input)[1]), comb.m(input, k)) 
     comb.SI <- SI(y, comb.new.m) 
      
     ##return(c(comb.SI, main.SI, comb.SI-sum(main.SI))) 
     return(cbind(comb.SI, main.SI)) 
} 


SI.main.total <- function(i, y, input, k=3){ 
     main.input <- apply(matrix(c(1:k), nc=1), 1, function(j){input[,i]^j}) 
     main.SI <- SI(y, main.input) 
      
     total.input.1 <- apply(input[,-i], 2, function(x){ apply(matrix(1:k),1, function(i){x^i}) }) 
     total.input <- cbind(matrix(total.input.1, nr=dim(input)[1]), comb.m(input[,-i], k)) 
     total.SI <- SI(y, input[,-i]) 
     return(cbind(total.SI, main.SI)) 
} 


interaction3gene <- function(y,gene)
  { 
     xm <- data.frame(gene) 
     m1 <- glm(y ~ .,data=xm); ##print(summary(m1)); 
     m2 <- glm(y ~ .^4, data=xm);##print(summary(m2)); 
     
     return(c(m1$deviance, m2$deviance, (m1$deviance-m2$deviance)/m1$deviance)) 
      
} 

test.sig.against.noise <- function(CYP, input.genes, sig.gene.index)
  { 
  sig.genes <- input.genes[,sig.gene.index]
  sig.genes <- as.matrix(sig.genes)
  noise.mu.est <- quantile(apply(sig.genes,2,function(x){mean(x,na.rm=T)}),c(1:10)/10)[c(1,9)] 
  noise.se.est <- quantile(apply(sig.genes,2,function(x){var(x,na.rm=T)^0.5}),c(1:10)/10)[c(1,9)] 
  noise.rho.est <- quantile( cor(sig.genes, use="pairwise.complete.obs"),c(1:10)/10)[c(1,9)] 
  
  d <- min(length(sig.gene.index),30) 
  noise.rho <- matrix(runif(1, min=0, max=noise.rho.est[2]),nr=d,nc=d) 
  noise.rho[lower.tri(noise.rho)] <- t(noise.rho)[lower.tri(noise.rho)] 
  noise.rho[row(noise.rho)==col(noise.rho)] <- 1 
  noise.se <- matrix(runif(d, min=noise.se.est[1], max=noise.se.est[2]),nc=1) 
  noise.cov <- noise.se %*% t(noise.se) * noise.rho 
  noise.mu <- matrix(runif(d, min=noise.mu.est[1], max=noise.se.est[2]),nc=1) 
   
  noise <- mvrnorm(length(CYP), noise.mu, noise.cov) 
  new.input <- cbind(sig.genes, noise) 
  
  return( list( scan.one.CYP(CYP, new.input), 
                length(sig.gene.index) , 
                names(input.genes)[sig.gene.index] )) 
} 


##----------- start analysis ---------------------------------------

CYPs.sorted <- hub_data[,81]
input.genes <- hub_data[,1:80]
CYPs.sorted.m <- as.matrix(CYPs.sorted) 
input.genes.m <- as.matrix(input.genes) 

##---------------------------------do glm---------------------------------------

scan.all.CYPs <- apply(CYPs.sorted.m, 2, function(y){scan.one.CYP(y,input.genes.m)})
names(scan.all.CYPs[[1]])
[1] "glm.coef.with.intercept" "glm.p.with.intercept"    "input.se.est"           
[4] "input.rho.est"           "est.gs"                  "est.gs.order"           
[7] "glm.p.order" 

glm.p.c <- do.call(rbind, lapply(scan.all.CYPs, 
                                 function(x){x[["glm.p.with.intercept"]]})) 

glm.p.c.fdr <- t(apply(glm.p.c, 1, function(x){p.adjust(x,"fdr")})) 

##-----------number of variables selected by controlling fdr at 0.05----------- 

table(apply(glm.p.c.fdr, 1, function(x){length(which(x<=0.05))})) 

sig.glm.p.fdr <- apply(glm.p.c.fdr, 1, function(x){which(x<=0.05)}) 

## sig.glm.p.fdr.v <- unlist(sig.glm.p.fdr) 
## length(sig.glm.p.fdr.v) 
colnames(glm.p.c.fdr)[sig.glm.p.fdr] 


glm.p.c.fdr.order <- t(apply(glm.p.c.fdr,1,
                             function(x){order(x,decreasing=FALSE,na.last=TRUE)}
                             )) 
dim(glm.p.c.fdr.order) 

input.gene.names <- names(input.genes) 
names <- c("intercept",input.gene.names) 

## apply(glm.p.c.fdr.order[,1:10], 1 , function(x){names[x]})

names[glm.p.c.fdr.order[,1:10]]        
#[1] "AL161668.5"    "LBX2.AS1"      "CTC.490E21.11" "AP000355.2"    "HNF4A.AS1"    
# "RP11.344P13.6" "ARNTL2"        "HLF"           "GATAD2A"       "ARID3C"  

##----------------variable selection via Sobol indices-------------------------- 

est.gs.c <- do.call(rbind, lapply(scan.all.CYPs, function(x){x[["est.gs"]]}))

##----------------------mean cut off of sobol indices--------------------------- 
##-----------------this cuf-off criteria works as long as----------------------- 
##--------the number of important genes is less than unimportant genes---------- 
##------------------------------------------------------------------------------ 

sig.est.gs <- apply(est.gs.c, 1, function(x){which( x > mean(x) )}) 
length(sig.est.gs) 
# [1] 36
sig.est.gs[1:5] 
# [1]  1  6  7  8 10
#unlist(lapply(sig.est.gs, function(x){length(x)})) 
dim(input.genes) 
# [1] 193  80
scan.all.CYPs.against.noise <- apply(matrix(1,nr=1), 2, 
                                    function(i){
                                       test.sig.against.noise(
                                        CYPs.sorted.m[,i],
                                         input.genes,
                                         sig.est.gs)}) 

##---------- 
##---------------------- noise cut off of sobol indices  ----------------------- 

sig.est.gs.with.noise <- lapply(scan.all.CYPs.against.noise, 
                                function(x)
                                  {gs <- x[[1]][["est.gs"]];
                                  return( which(gs[c(1:x[[2]])] > max(gs[-c(1:x[[2]])] ) )); }) 

sig.est.gs.with.noise[[1]] 
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
# [29] 29 30 31 32 33 34 35 36
sig.est.gs 

# sig.est.gs.with.noise[[2]] 
# sig.est.gs[[2]] 


plot(scan.all.CYPs.against.noise[[1]][[1]][["est.gs"]], pch = 19, col = "lightblue", ylab = "est.gs") 
scan.all.CYPs.against.noise[[1]][[3]] 
#[1] "GCFC2"          "CREB3L3"        "TEAD2"          "ESR1"           "HIF1A"         
#[6] "HLF"            "KLF12"          "SLC2A4RG"       "NR1I3"          "ZNF385B"       
#[11] "NR1I2"          "NR3C2"          "NFIA"           "AR"             "AFF1"          
#[16] "ZFP1"           "ZGPAT"          "TEAD4"          "ARID3C"         "RP1.40E16.11"  
#[21] "HNF4A.AS1"      "AC004538.3"     "AC004160.4"     "AC004862.6"     "RP5.881L22.6"  
#[26] "NME2"           "RP11.622A1.2"   "RP11.119D9.1"   "RP11.115C10.1"  "RP11.252E2.2"  
#[31] "RP11.669E14.4"  "RP11.250B2.6"   "KB.68A7.1"      "CTC.490E21.11"  "RP11.1182P23.5"
#[36] "RP11.258F1.2"  

##------------------------ confirm no signal genes ----------------------------- 

index <- c(1:80)    # number of input genes  

scan.all.CYPs.against.noise.insig.genes <- apply(matrix(1,nr=1), 2, 
                                                 function(i)
                                                   {test.sig.against.noise(
                                                     CYPs.sorted.m[,i],
                                                     input.genes,
                                                     index[-sig.est.gs ] )
                                                    }) 


insig.est.gs.with.noise <- lapply(scan.all.CYPs.against.noise.insig.genes, function(x){gs <- x[[1]][["est.gs"]];return( which(gs[c(1:x[[2]])] > max(gs[-c(1:x[[2]])] ) )); }) 


plot(scan.all.CYPs.against.noise.insig.genes[[1]][[1]][["est.gs"]], pch = 19, ylab = "est.gs", col = "lightblue") 
insig.est.gs.with.noise[[1]] 
#[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
#[31] 31 32 33 34 35 36 37 38 39 40 41 42 43 44
index[-sig.est.gs ] 


# insig.est.gs.with.noise[[2]] 
# index[-sig.est.gs[[2]] ] 

# plot(scan.all.CYPs.against.noise[[2]][[1]][["est.gs"]]) 
# dev.new() 
# plot(scan.all.CYPs.against.noise.insig.genes[[2]][[1]][["est.gs"]]) 

##---------------------input normality test-------------------------------------

input.normality.test <- unlist( apply(input.genes.m, 2, 
                                      function(x)
                                        {t <- shapiro.test(x); 
                                        return(t$p.value)}) ) 
length(which(input.normality.test<=0.05)) 
# [1] 60


scan.all.CYPs <- apply(CYPs.sorted.m, 2, 
                       function(y){scan.one.CYP(y,input.genes.m)}) 
glm.p.c <- do.call(rbind, lapply(scan.all.CYPs, 
                                 function(x){x[["glm.p.with.intercept"]]})) 

glm.p.c.fdr <- t(apply(glm.p.c, 1, function(x){p.adjust(x,"fdr")})) 


##-----------------------pairs-----------Example of analyzing one CYP at a time:  
## y1 is the response variable=CYP expression of interest 
## y1 <- CYPs.sorted.m[,30] ## CYP3A4 

y1 <- CYPs.sorted.m[,1]

index.2m.all <- combn(c(1:80), 2)   # number of input gene  
dim(index.2m.all) 
# [1]    2 3160
input.gene.names <- names(input.genes)
gene.names.2m.all <- apply(index.2m.all, 2, function(i){input.gene.names[i]}) 
all.gene.pairs.3 <- apply(index.2m.all, 2, 
                          function(i){SI.comb(y1, input.genes.m[,i], k=3)}) 
order.2m.3 <- order(-all.gene.pairs.3[2,], all.gene.pairs.3[1,],decreasing=TRUE) 
 

gene.names.2m.all[,order.2m.3][,1:10] 
#        [,1]    [,2]           [,3]         [,4]            [,5]    [,6]      [,7]    [,8]       
#[1,] "TEAD2" "TEAD2"        "TEAD2"      "ESR1"          "TEAD2" "TEAD2"   "TEAD2" "ESR1"     
#[2,] "ESR1"  "RP11.250B2.6" "AC004160.4" "CTC.490E21.11" "AR"    "ZNF385B" "NME2"  "HNF4A.AS1"
#        [,9]    [,10]  
#[1,] "TEAD2" "TEAD2"
#[2,] "TEAD4" "HLF"  
all.gene.pairs.3[,order.2m.3][,1:10] 


##-------------------------------------------------------  

index.3m.all <- combn(c(1:80), 3) 
dim(index.3m.all) 


gene.names.3m.all <- apply(index.3m.all, 2, function(i){input.gene.names[i]}) 
all.gene.triplets.3 <- apply(index.3m.all, 2, function(i){SI.comb(y1, input.genes.m[,i], k=3)}) 
order.3m.3 <- order(-all.gene.triplets.3[2,], all.gene.triplets.3[1,],decreasing=TRUE) 
 

gene.names.3m.all[,order.3m.3][,1:30] 
#     [,1]            [,2]            [,3]            [,4]         [,5]           
#[1,] "AC104809.2"    "NR1I2"         "ESR1"          "AR"         "HNF4A.AS1"    
#[2,] "HNF4A.AS1"     "HNF4A.AS1"     "HNF4A.AS1"     "AC104809.2" "RP11.115C10.1"
#[3,] "CTC.490E21.11" "CTC.490E21.11" "CTC.490E21.11" "HNF4A.AS1"  "CTC.490E21.11"
#     [,6]            [,7]            [,8]            [,9]            [,10]          
#[1,] "ESR1"          "AR"            "HNF4A.AS1"     "AC104809.2"    "ESR1"         
#[2,] "RP5.881L22.6"  "AC104809.2"    "RP11.706C16.7" "RP11.70D24.2"  "RP11.115C10.1"
#[3,] "CTC.490E21.11" "CTC.490E21.11" "RP11.115C10.1" "CTC.490E21.11" "CTC.490E21.11"
#     [,11]           [,12]           [,13]           [,14]           [,15]       
#[1,] "HNF4A.AS1"     "RP11.122K13.7" "RP11.115C10.1" "AC104809.2"    "HMGB1"     
#[2,] "AC004538.3"    "HNF4A.AS1"     "LBX2.AS1"      "RP11.252E2.2"  "AC104809.2"
#[3,] "CTC.490E21.11" "CTC.490E21.11" "CTC.490E21.11" "CTC.490E21.11" "HNF4A.AS1" 
#     [,16]           [,17]        [,18]          [,19]           [,20]          
#[1,] "ESR1"          "ESR1"       "AC104809.2"   "ESR1"          "CREB3L3"      
#[2,] "ZGPAT"         "AC104809.2" "HNF4A.AS1"    "AC104809.2"    "HNF4A.AS1"    
#[3,] "CTC.490E21.11" "HNF4A.AS1"  "RP5.881L22.6" "CTC.490E21.11" "CTC.490E21.11"
#     [,21]        [,22]           [,23]           [,24]            [,25]          
#[1,] "CREB3L3"    "ESR1"          "ESR1"          "AC104809.2"     "YBX3"         
#[2,] "AC104809.2" "NR1I3"         "RP11.706C16.7" "HNF4A.AS1"      "AC104809.2"   
#[3,] "HNF4A.AS1"  "CTC.490E21.11" "CTC.490E21.11" "RP11.1182P23.5" "CTC.490E21.11"
#     [,26]           [,27]           [,28]        [,29]           [,30]       
#[1,] "NR1I2"         "ARNTL2"        "NR3C2"      "ESR1"          "AC104809.2"
#[2,] "RP11.115C10.1" "ESR1"          "AC104809.2" "AR"            "HNF4A.AS1" 
#[3,] "CTC.490E21.11" "CTC.490E21.11" "HNF4A.AS1"  "CTC.490E21.11" "AL161668.5"

all.gene.triplets.3[,order.3m.3][,1:30]
 

table(gene.names.2m.all[,order.2m.3][,1:10]) 
# AC004538.3    AC104809.2    AL161668.5            AR         CREB3 CTC.490E21.11 
#          1             3             1             1             1             3 
# ESR1         HMGB1     HNF4A.AS1         PATZ1 
#    1             2             6             1 
table(gene.names.3m.all[,order.3m.3][,1:100])

table(gene.names.2m.all[,order.2m.3][,1:200])

table(gene.names.3m.all[,order.3m.3][,1:200]) 


dev.new();hist(all.gene.triplets.3[2,], breaks=1200) 
dev.new();hist(all.gene.triplets.3[1,order.3m.3][1:200], breaks=50) 

##------------------------------------------------------- 
##----------code for scanning for subset effect using polynomial with degree k=3
##----------and step wise model selection by AIC: 
##----------need to impute all missing values before using the following code 
 
##-------------------------------single-----------------------------------------

all.gene.main.total.SI.3.stepwise <- apply(matrix(1:80, nr=1), 2, 
                                           function(i)
                                             {SI.main.total(i, y1, 
                                                           input.genes.m, 
                                                           k=3)}) 
summary(all.gene.main.total.SI.3.stepwise[1,])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.779   4.891   4.910   4.898   4.915   4.918  

which(all.gene.main.total.SI.3.stepwise[1,]>4.917)
# [1]  35

M.T.SI <- data.frame(names[2:81], t(all.gene.main.total.SI.3.stepwise))
dim(M.T.SI)
# [1] 80  7

names(M.T.SI) <- c("Gene", "T.SI", "T.ResDev.Perc", "T.NullDev", "M.SI", "M.ResDev.Perc", "M.NullDev") 

M.T.SI[1:3,] 
#    Gene     T.SI T.ResDev.Perc T.NullDev     M.SI M.ResDev.Perc M.NullDev
#1  GCFC2 4.758742     0.2071717  1152.429 1.521635     0.7464887  1152.429
#2 ARNTL2 4.751903     0.2083111  1152.429 1.653318     0.7245497  1152.429
#3  HMGB3 4.736358     0.2109010  1152.429 1.318803     0.7802814  1152.429

order.1m.3.stepwise <- order(-M.T.SI$M.ResDev.Perc, M.T.SI$M.SI,decreasing=TRUE) 
M.T.SI.report <- M.T.SI[,-c(2:4)] 
M.T.SI.report <- M.T.SI.report[order.1m.3.stepwise[1:10],] 
#            Gene     M.SI M.ResDev.Perc M.NullDev
#7          TEAD2 2.533869     0.5778457  1152.429
#8           ESR1 2.362026     0.6064755  1152.429
#36         TEAD4 2.127974     0.6454698  1152.429
#21       ZNF385B 2.074340     0.6544055  1152.429
#52  RP11.119D9.1 2.021925     0.6631379  1152.429
#44     HNF4A.AS1 2.019709     0.6635072  1152.429
#25          NFIA 2.005801     0.6658243  1152.429
#75 CTC.490E21.11 1.983347     0.6695652  1152.429
#35         ZGPAT 1.886991     0.6856186  1152.429
#31          ZFP1 1.878830     0.6869782  1152.429
#M.T.SI.report$M.SI.re <- M.T.SI.report$M.SI-1

ggplot(M.T.SI.report, aes(x = reorder(Gene, -M.SI), y = M.SI)) +
  geom_bar(stat = "identity", aes(fill = "blue")) + 
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Genes", y = "Total Effect indices", title = "") +
  guides(fill=FALSE) +
  #coord_cartesian(ylim = c(0.4,0.55)) +
  theme_gray(base_size = 15)


M.SI.single <- M.T.SI.report[order.1m.3.stepwise,] 
M.T.SI.report[order.1m.3.stepwise[1:10],] 
#            Gene     M.SI M.ResDev.Perc M.NullDev
#7          TEAD2 2.533869     0.5778457  1152.429
#8           ESR1 2.362026     0.6064755  1152.429
#36         TEAD4 2.127974     0.6454698  1152.429
#21       ZNF385B 2.074340     0.6544055  1152.429
#52  RP11.119D9.1 2.021925     0.6631379  1152.429
#44     HNF4A.AS1 2.019709     0.6635072  1152.429
#25          NFIA 2.005801     0.6658243  1152.429
#75 CTC.490E21.11 1.983347     0.6695652  1152.429
#35         ZGPAT 1.886991     0.6856186  1152.429
#31          ZFP1 1.878830     0.6869782  1152.429

order.1m.3.stepwise <- order(-M.T.SI$T.ResDev.Perc, M.T.SI$T.SI,decreasing=TRUE) 
M.T.SI.report <- M.T.SI[,-c(5:7)] 
M.T.SI.report[order.1m.3.stepwise[1:30],] 


T.SI.single <- M.T.SI.report[order.1m.3.stepwise,] 
M.T.SI.report[order.1m.3.stepwise[1:10],] 
#           Gene     T.SI T.ResDev.Perc T.NullDev
#35        ZGPAT 4.793267     0.2014196  1152.429
#11          DBP 4.764598     0.2061960  1152.429
#63 RP11.250B2.6 4.764598     0.2061960  1152.429
#77   AC006994.2 4.764598     0.2061960  1152.429
#14       ZBTB47 4.764598     0.2061960  1152.429
#23        NR3C2 4.764598     0.2061960  1152.429
#58 RP11.252E2.2 4.764598     0.2061960  1152.429
#8          ESR1 4.762352     0.2065703  1152.429
#5          YBX3 4.758742     0.2071717  1152.429
#7         TEAD2 4.758742     0.2071717  1152.429

##------------pairs-----Generate csv files for Gephi:---------------------------

# dim(all.gene.pairs.3.stepwise)  # [1]    9 3003 
  dim(all.gene.pairs.3)
# [1]    9 3160

# all.gene.pairs.3.stepwise[,1:5] 
  all.gene.pairs.3[,1:5]
# [,1]         [,2]         [,3]         [,4]         [,5]
# [1,]    1.9950437    1.8796584    1.9106413    1.8795657    1.9770963
# [2,]    0.6676165    0.6868403    0.6816784    0.6868557    0.6706067
# [3,] 1152.4291799 1152.4291799 1152.4291799 1152.4291799 1152.4291799
# [4,]    1.5467488    1.5467488    1.5467488    1.5467488    1.5467488
# [5,]    0.7423045    0.7423045    0.7423045    0.7423045    0.7423045
# [6,] 1152.4291799 1152.4291799 1152.4291799 1152.4291799 1152.4291799
# [7,]    1.6690549    1.3190563    1.3806764    1.1628318    1.6202014
# [8,]    0.7219278    0.7802392    0.7699730    0.8062669    0.7300670
# [9,] 1152.4291799 1152.4291799 1152.4291799 1152.4291799 1152.4291799

dim(gene.names.2m.all)  
# [1]    2 3160
choose(80,2)  
# [1] 3160 
gene.names.2m.all[,1:5] 


Gephi.gene.pairs.edge <- data.frame(t(gene.names.2m.all), 
                                    t(all.gene.pairs.3[1:3,])) 
dim(Gephi.gene.pairs.edge)  
# [1] 3160    5
Gephi.gene.pairs.edge[1:5,] 

Gephi.gene.pairs.edge <- Gephi.gene.pairs.edge[order(Gephi.gene.pairs.edge$X1.1, decreasing = TRUE), ]
write.csv(Gephi.gene.pairs.edge, file="gene.pairs.edge.csv")

Gephi.gene.pairs.edge$Label <- apply(gene.names.2m.all, 2, 
                                     function(x){paste(x[1],"-",x[2])}) 
Gephi.gene.pairs.edge[1:5,]  # file of label

Gephi.gene.pairs.node <- data.frame(c(gene.names.2m.all[1,],
                                      gene.names.2m.all[2,]), 
                                    rbind(t(all.gene.pairs.3[4:6,]), 
                                          t(all.gene.pairs.3[7:9,]))) 
dim(Gephi.gene.pairs.node)  
# [1] 6320    4

Gephi.gene.pairs.node.unique <- Gephi.gene.pairs.node[!duplicated(Gephi.gene.pairs.node[,1]),] 
dim(Gephi.gene.pairs.node.unique)
# [1] 80  4 
write.csv(Gephi.gene.pairs.node.unique, file="gene.pairs.node.csv") 


##--------------------------------triplets-------------------
dim(all.gene.triplets.3.stepwise) 
all.gene.triplets.3.stepwise[,1:5] 
dim(gene.names.3m.all) 
choose(80,3) 

gene.names.3m.all[,1:5] 
Gephi.gene.triplets.edge <- data.frame(t(gene.names.3m.all), t(all.gene.triplets.3.stepwise[1:3,])) 
dim(Gephi.gene.triplets.edge) 


Gephi.gene.triplets.edge[1:5,] 
 

Gephi.gene.triplets.edge$Label <- apply(gene.names.3m.all, 2, function(x){paste(x[1],"-",x[2], "-", x[3])}) 
Gephi.gene.triplets.edge[1:5,] 
#write.csv(Gephi.gene.triplets.edge, file="U:/Data/Gephi_gene_triplets_edge.csv") 
 

Gephi.gene.triplets.edge.1 <- Gephi.gene.triplets.edge 
Gephi.gene.triplets.edge.2 <- Gephi.gene.triplets.edge 
Gephi.gene.triplets.edge.2[,1] <- Gephi.gene.triplets.edge.1[,2] 
Gephi.gene.triplets.edge.2[,2] <- Gephi.gene.triplets.edge.1[,3] 
#write.csv(rbind(Gephi.gene.triplets.edge.1, Gephi.gene.triplets.edge.2), 
#          file="U:/Data/Gephi_gene_triplets_edge_in_pairs.csv") 


##--------------------start analysis of GTEx RNAseq-----------------------------

##----------------------------functions-----------------------------------------
SI <- function(y, input)
  {
  colnames(input) <- NULL 
  m.full <- glm(y~., data=data.frame(input)); 
  fit.step <- stepAIC(m.full, direction="both", trace=F) 
  m.best <- formula(fit.step) 
  fit.1 <- glm(m.best, data=data.frame(input));glm.coef <- fit.1$coef; 
  v.select <- input[, as.integer(sub("X", "", attr(terms(fit.step),"term.labels") ))] 
  SI <- var( matrix(v.select, nr=dim(input)[1]) %*% matrix(glm.coef[-1], nc=1) , na.rm=TRUE ) 
  return(c(SI, fit.1$deviance / fit.1$null.deviance, fit.1$null.deviance)) 
  } 

## code tested up to slection of 4comb with k=3: 
comb.m <- function(input, k=3)
  { 
     d <- dim(input)[2] 
     combn.index <- apply(matrix(2:(min(d,k)), nr=1),2, 
                          function(x){list(combn(c(1:d), x))} )  
      
     combn.list <- lapply(combn.index, 
                          function(x){cc<-input[,x[[1]][1,]];
                                      for(j in 2:dim(x[[1]])[1])
                                        {cc <- cc*input[,x[[1]][j,]]}; 
                                      return(list(cc)); 
                                     }) 
     combn.input <- do.call(cbind, lapply(combn.list, function(x){x[[1]]})) 
      
     return(combn.input) 
  } 

#--------------------------------data----------------------------------------
dim(rpkm)
# [1] 193  81
rpkm[1:3,1:3]
#                               GCFC2     ARNTL2      HMGB3
#GTEX-11DXY-0526-SM-5EGGQ  0.20455075 -0.9764185 -1.1764500
#GTEX-11DXZ-0126-SM-5EGGY  0.38479019 -0.5432659 -0.9081299
#GTEX-11EQ9-0526-SM-5A5JZ -0.04864744  0.1436718 -0.3364388

y.m <- data[,81] 
y.m <- as.matrix(y.m)
gene.tb <- data[,1:80] 

which(names(gene.tb)=="CYP3A4")     # integer(0) 
which(names(gene.tb)=="y")          # integer(0) 
which(names(rpkm)=="y")             # [1] 81

apply(y.m, 2, summary)
#                  [,1]
# Min.    -7.235413e+00
# 1st Qu. -9.788475e-01
# Median   5.020985e-01
# Mean    -8.549223e-09
# 3rd Qu.  1.626908e+00
# Max.     7.721721e+00

#----------------logged CYP3A4 ~ unlogged lnc and TF: triplets------------------

log.CYP3A4 <- y.m[,1]                               # log done before WGCNA
unlog.gene.tb <- 10^(gene.tb)

hist(log.CYP3A4, breaks=20, main = "Frequency of CYP3A4")

index.m <- combn(c(1:80),3) 
dim(index.m) 
# [1]     3 82160
gene.names <- names(unlog.gene.tb) 
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])}) 
dim(unlog.gene.tb)
# [1] 193  80
results <- apply(index.m, 2, 
                 function(i){gene <- unlog.gene.tb[,i];   # CYP3A4 was logged, gene.tb were not
                             return(interaction3gene(log.CYP3A4, gene))}) 
Residual.Deviance.Deduction <- results[3,] 
hist(Residual.Deviance.Deduction,breaks=100, main="log.CYP3A4 Residual Deviance.Deduction")

Residual.Deviance <- results[2,] 
hist(Residual.Deviance,breaks=100,main="log.CYP3A4 Residual Deviance")

summary(Residual.Deviance)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 57.38   82.10   86.28   85.58   89.78   99.65 
summary(Residual.Deviance.Deduction) 
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001495 0.0199949 0.0352311 0.0415804 0.0564617 0.2448966 

#--------------------------threshold changed according to summary---------------

large.interaction.index <- 
  which((Residual.Deviance<=65) & (Residual.Deviance.Deduction>=0.10) )
length(large.interaction.index)
# [1] 9

goodfit.index <- which((Residual.Deviance<=60) )
length(goodfit.index)
# [1] 3

goodfitr <- rbind(index.m[,goodfit.index], 
                  gene.names.m[,goodfit.index], 
                  results[,goodfit.index]) 

goodfitr[,order(goodfitr[8,])] 
#     [,1]                [,2]                [,3]                
#[1,] "27"                "27"                "27"                
#[2,] "36"                "36"                "36"                
#[3,] "45"                "43"                "65"                
#[4,] "AR"                "AR"                "AR"                
#[5,] "TEAD4"             "TEAD4"             "TEAD4"             
#[6,] "AC004538.3"        "AP000355.2"        "RP11.344P13.6"     
#[7,] "66.9962495172548"  "66.8945664030034"  "63.2868754541754"  
#[8,] "57.3790418509942"  "59.2733995349819"  "59.3011045086252"  
#[9,] "0.143548448391633" "0.113928040464574" "0.0629794237264285"

#---------------------logged lncRNA and TFs: triplets---------------------------

log.CYP3A4 <- data[, 81]
log.gene.tb <- data[, 1:80]                              #log has done before WGCNA
log.gene.tb[1:3, 1:5]
#                                GCFC2     ARNTL2      HMGB3      ARID4A       YBX3
# GTEX-11DXY-0526-SM-5EGGQ  0.20455075 -0.9764185 -1.1764500  0.06721183 -2.2118944
# GTEX-11DXZ-0126-SM-5EGGY  0.38479019 -0.5432659 -0.9081299  0.42960775 -0.2467226
# GTEX-11EQ9-0526-SM-5A5JZ -0.04864744  0.1436718 -0.3364388 -0.01874430 -0.6303828

# hist(gene.tb[,21], breaks=20)
# hist(gene.tb[,40], breaks=20) 
# gene.tb <- gene.tb[,-c(21,40)] 
# dim(gene.tb)

index.m <- combn(c(1:80),3)
dim(index.m) 
# [1]     3 82160
gene.names <- names(log.gene.tb) 
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
results <- apply(index.m, 2, 
                 function(i){gene <- log.gene.tb[,i]; 
                             return(interaction3gene(log.CYP3A4, gene))})
Residual.Deviance.Deduction <- results[3,] 
hist(Residual.Deviance.Deduction,breaks=100, col="lightblue", border="white")

Residual.Deviance <- results[2,] 
hist(Residual.Deviance,breaks=100,main="log.CYP3A4 Residual Deviance", 
     col="lightblue", border="white")

summary(Residual.Deviance)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 550.3   724.5   761.7   762.2   802.2   939.1  
summary(Residual.Deviance.Deduction)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003323 0.0260248 0.0408316 0.0451600 0.0596371 0.2043722 

#-------------------threshold changed according to summary----------------------

large.interaction.index <- 
           which((Residual.Deviance<=60) & (Residual.Deviance.Deduction>=0.20) ) 
length(large.interaction.index)
# [1] 12
goodfit.index <- which((Residual.Deviance<=580) ) 
length(goodfit.index)
# [1] 19
goodfitr <- rbind(index.m[,goodfit.index], 
                  gene.names.m[,goodfit.index], 
                  results[,goodfit.index]) 

goodfitr[,order(goodfitr[8,])]
#     [,1]                 [,2]                [,3]               [,4]               
#[1,] "7"                  "7"                 "7"                "7"                
#[2,] "8"                  "8"                 "23"               "63"               
#[3,] "46"                 "69"                "46"               "69"               
#[4,] "TEAD2"              "TEAD2"             "TEAD2"            "TEAD2"            
#[5,] "ESR1"               "ESR1"              "NR3C2"            "RP11.250B2.6"     
#[6,] "AC004160.4"         "RP5.875H18.10"     "AC004160.4"       "RP5.875H18.10"    
#[7,] "599.955872565858"   "619.424675644047"  "652.134630311867" "651.331004208455" 
#[8,] "550.343234935723"   "554.837466675961"  "555.585338992905" "555.80901414404"  
#[9,] "0.0826938111597337" "0.104269673953386" "0.14805116433211" "0.146656599251713"
#     [,5]                [,6]                [,7]                [,8]               
#[1,] "7"                 "7"                 "7"                 "8"                
#[2,] "46"                "24"                "28"                "27"               
#[3,] "69"                "46"                "46"                "44"               
#[4,] "TEAD2"             "TEAD2"             "TEAD2"             "ESR1"             
#[5,] "AC004160.4"        "BATF"              "AFF1"              "AR"               
#[6,] "RP5.875H18.10"     "AC004160.4"        "AC004160.4"        "HNF4A.AS1"        
#[7,] "675.205739839465"  "670.655047635603"  "658.108180279658"  "667.043364114554" 
#[8,] "559.534615750196"  "562.64403348692"   "566.276092490628"  "566.893539800869" 
#[9,] "0.171312412297873" "0.161053010082421" "0.139539502077009" "0.150139900494514"
#     [,9]                [,10]               [,11]                [,12]              
#[1,] "7"                 "7"                 "7"                  "7"                
#[2,] "16"                "27"                "63"                 "46"               
#[3,] "46"                "46"                "75"                 "53"               
#[4,] "TEAD2"             "TEAD2"             "TEAD2"              "TEAD2"            
#[5,] "KLF12"             "AR"                "RP11.250B2.6"       "AC004160.4"       
#[6,] "AC004160.4"        "AC004160.4"        "CTC.490E21.11"      "RP11.706C16.7"    
#[7,] "662.719341541486"  "656.233541741905"  "632.129946300219"   "670.006260426038" 
#[8,] "567.145271378297"  "568.884349099891"  "569.593817892467"   "571.274726006861" 
#[9,] "0.144215000487059" "0.133106869865496" "0.0989292293044627" "0.147359122221331"
#     [,13]               [,14]               [,15]               [,16]              
#[1,] "7"                 "7"                 "7"                 "7"                
#[2,] "46"                "9"                 "46"                "9"                
#[3,] "58"                "46"                "75"                "36"               
#[4,] "TEAD2"             "TEAD2"             "TEAD2"             "TEAD2"            
#[5,] "AC004160.4"        "PATZ1"             "AC004160.4"        "PATZ1"            
#[6,] "RP11.252E2.2"      "AC004160.4"        "CTC.490E21.11"     "TEAD4"            
#[7,] "656.848723193384"  "673.920526408348"  "653.628124074882"  "703.484460252124" 
#[8,] "573.511766910257"  "574.171447667006"  "575.247431700711"  "576.901176659545" 
#[9,] "0.126873895526461" "0.148013118509616" "0.119916340021488" "0.179937568979438"
#     [,17]              [,18]                [,19]              
#[1,] "8"                "7"                  "7"                
#[2,] "44"               "8"                  "25"               
#[3,] "57"               "63"                 "46"               
#[4,] "ESR1"             "TEAD2"              "TEAD2"            
#[5,] "HNF4A.AS1"        "ESR1"               "NFIA"             
#[6,] "RP1.239B22.5"     "RP11.250B2.6"       "AC004160.4"       
#[7,] "656.362952353902" "600.883211107158"   "651.236522916266" 
#[8,] "577.889482924483" "578.191756853339"   "578.677259555163" 
#[9,] "0.11955804200708" "0.0377635018492357" "0.111417681299843"

#------------------------------split peaks--------------------------------------

index.low <- which(Residual.Deviance<=650) 
length(index.low)/(length(Residual.Deviance[-index.low])) # [1] 0.01943073 

length(index.low)
# [1] 1566


par(mfrow=c(1,2)) 
hist(Residual.Deviance[index.low], breaks=30, ylim = c(0,3500));  
hist(Residual.Deviance[-index.low], breaks=70, ylim = c(0,3500)); 
 
 
top <-  table(gene.names.m[,index.low])
top[order(top, decreasing = T)]
#          ESR1          TEAD2  CTC.490E21.11   RP11.250B2.6          TEAD4      HNF4A.AS1 
#          3024           2870           1630           1416           1382           1339 
#          NFIA     AC004160.4          NR1I2             AR            HLF        ZNF385B 
#          1310           1288           1180            958            881            872 
#  RP11.119D9.1          NR3C2   RP11.252E2.2        CREB3L3           ZFP1 RP11.1182P23.5 
#           866            770            714            672            667            666 
# RP11.669E14.4          ZGPAT           NME2     AC004538.3          PATZ1       LBX2.AS1 
#           648            595            582            548            528            510 
# RP11.115C10.1           AFF1           YBX3         ARNTL2          GCFC2   CTD.2194A8.2 
#           500            485            461            448            426            409 
# RP11.122K13.7     AL161668.5        GATAD2A           BATF         PLSCR1         ARID4A 
#           395            394            384            372            370            369 
# RP11.104J23.1          NR1I3   RP5.881L22.6     AC004862.6          MLXIP          KLF12 
#           364            363            362            359            358            357 
#         HIF1A  RP11.94C24.13  CTD.2012K14.8         ZNF33B       SLC2A4RG          HMGB3 
#           344            342            331            331            327            320 
#  RP11.258F1.2         ZNF680     AP000355.2            DBP   RP11.70D24.2  RP5.875H18.10 
#           316            315            313            312            310            310 
#  RP1.40E16.11       MIR135A1  RP13.650J16.1  RP11.263G22.1  RP11.706C16.7           CARF 
#           305            303            298            286            281            278 
#  RP1.239B22.5  CTD.2152M20.2  CTD.2527I21.5     AC006994.2      KB.68A7.1        DBH.AS1 
#           276            274            273            270            268            259 
# RP11.344P13.6           RORC    RP11.96D1.6         ARID3C          HMGB1     AC104809.2 
#           250            248            243            239            239            235 
#          MYCL   RP5.1159O4.2   RP11.622A1.2    RP11.39H3.2  RP11.449P15.2  CTD.2325A15.5 
#           229            219            214            208            208            193 
#         CREB3         ZBTB47 
#           190            161 

top <- as.list(top[order(top, decreasing = T)])
top <- as.data.frame(top)
top <- t(top)
rownames(top) -> names
top <- as.data.frame(top)
top[,2] <- names
colnames(top) <- c("Count", "names")
top[, 3] <- top[,1 ]/length(index.low)
colnames(top)[3] <- c("Freq")

ggplot(top[1:20,], aes(x = reorder(names, -Freq), y = Freq)) +
  geom_bar(stat = "identity", aes(fill = "blue")) + 
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Gene names", y = "Frequency", title = "Interaction Triplets") +
  guides(fill=FALSE) +
  theme_gray(base_size = 15) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))



#-----------------------------logged lnc and TF: pairs--------------------------
interaction2gene <- function(y,gene)
  { 
  xm <- data.frame(gene) 
  m1 <- glm(y ~ .,data=xm); ##print(summary(m1)); 
  m2 <- glm(y ~ .^4, data=xm);##print(summary(m2)); 
  
  return(c(m1$deviance, m2$deviance, (m1$deviance-m2$deviance)/m1$deviance))
  }

log.gene.tb <- data[,1:80]
log.CYP3A4 <- data[,81]

index.m <- combn(c(1:80),2)
dim(index.m) 
# [1]     2 3160
gene.names <- names(log.gene.tb) 
gene.names.m<- apply(index.m, 2, function(x){return(gene.names[x])})
results <- apply(index.m, 2, 
                 function(i){gene <- log.gene.tb[,i]; 
                 return(interaction3gene(log.CYP3A4, gene))})
Residual.Deviance.Deduction <- results[3,] 
hist(Residual.Deviance.Deduction,breaks=100, col="lightblue", border="white")

Residual.Deviance <- results[2,] 
hist(Residual.Deviance,breaks=100,main="log.CYP3A4 Residual Deviance", 
           col="lightblue", border="white")


summary(Residual.Deviance)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 622.0   780.6   825.1   823.9   868.7   962.7 
summary(Residual.Deviance.Deduction)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.001374 0.005862 0.011918 0.016553 0.154229 

#-------------------threshold changed according to summary----------------------

large.interaction.index <- 
  which((Residual.Deviance<=700) & (Residual.Deviance.Deduction>=0.03) ) 
length(large.interaction.index)
# [1] 6
goodfit.index <- which((Residual.Deviance<=680) ) 
length(goodfit.index)
# [1] 11
goodfitr <- rbind(index.m[,goodfit.index], 
                  gene.names.m[,goodfit.index], 
                  results[,goodfit.index]) 

goodfitr[,order(goodfitr[7,])]
#     [,1]                  [,2]                  [,3]                 
#[1,] "8"                   "7"                   "8"                  
#[2,] "75"                  "49"                  "22"                 
#[3,] "ESR1"                "TEAD2"               "ESR1"               
#[4,] "CTC.490E21.11"       "NME2"                "NR1I2"              
#[5,] "642.404544909063"    "666.335623106696"    "658.261618767024"   
#[6,] "641.642607992382"    "665.089190430605"    "656.562692175358"   
#[7,] "0.00118607024610692" "0.00187057787827738" "0.00258092913703201"
#     [,4]                  [,5]                  [,6]                
#[1,] "7"                   "8"                   "8"                 
#[2,] "8"                   "46"                  "44"                
#[3,] "TEAD2"               "ESR1"                "ESR1"              
#[4,] "ESR1"                "AC004160.4"          "HNF4A.AS1"         
#[5,] "626.43134514859"     "669.295424299785"    "672.971941440216"  
#[6,] "622.012061019047"    "664.086661574063"    "665.28783392333"   
#[7,] "0.00705469827423042" "0.00778245679950858" "0.0114181692336859"
#     [,7]                 [,8]                 [,9]                 [,10]               
#[1,] "6"                  "8"                  "7"                  "7"                 
#[2,] "8"                  "63"                 "63"                 "13"                
#[3,] "CREB3L3"            "ESR1"               "TEAD2"              "TEAD2"             
#[4,] "ESR1"               "RP11.250B2.6"       "RP11.250B2.6"       "HLF"               
#[5,] "688.272379491918"   "686.184745556496"   "654.762881215318"   "696.492763066359"  
#[6,] "676.884136652489"   "673.250049063237"   "641.551261433873"   "676.930914368115"  
#[7,] "0.0165461279266149" "0.0188501661936087" "0.0201777164840532" "0.0280862196071081"
#     [,11]               
#[1,] "7"                 
#[2,] "46"                
#[3,] "TEAD2"             
#[4,] "AC004160.4"        
#[5,] "676.237184871062"  
#[6,] "638.444281100829"  
#[7,] "0.0558870535601193"

#-----------------------------------split peak----------------------------------

index.low <- which(Residual.Deviance<=750) 
length(index.low)
# [1] 342
length(index.low)/(length(Residual.Deviance[-index.low])) # [1] 0.1213627

par(mfrow=c(1,2)) 
hist(Residual.Deviance[index.low], breaks=50, ylim = c(0,120));  
hist(Residual.Deviance[-index.low], breaks=50, ylim = c(0,120)); 

top <-  table(gene.names.m[,index.low])
top[order(top, decreasing = T)]
#           ESR1          TEAD2      HNF4A.AS1  CTC.490E21.11           NFIA          TEAD4 
#             79             78             39             37             28             27 
#          NR1I2   RP11.250B2.6     AC004160.4   RP11.119D9.1        CREB3L3           NME2 
#             23             21             19             19             13             13 
# RP11.1182P23.5        ZNF385B             AR            HLF  RP11.669E14.4          ZGPAT 
#             12             12             11             11             11             11 
#   RP11.252E2.2          NR3C2     AC004538.3  RP11.115C10.1           ZFP1           AFF1 
#              9              8              7              7              7              6 
#          KLF12          NR1I3     AL161668.5          GCFC2          HMGB3       LBX2.AS1 
#              6              6              5              5              5              5 
#          PATZ1  RP11.94C24.13   RP5.881L22.6           YBX3         ARID4A         ARNTL2 
#              5              5              5              5              4              4 
#           BATF   CTD.2194A8.2            DBP         PLSCR1  RP11.122K13.7     AC004862.6 
#              4              4              4              4              4              3 
#     AC006994.2     AP000355.2           CARF  CTD.2012K14.8  CTD.2152M20.2  CTD.2527I21.5 
#              3              3              3              3              3              3 
#        DBH.AS1        GATAD2A          HIF1A           RORC   RP1.239B22.5  RP11.104J23.1 
#              3              3              3              3              3              3 
#   RP11.258F1.2  RP11.706C16.7    RP11.96D1.6  RP13.650J16.1       SLC2A4RG         ZNF33B 
#              3              3              3              3              3              3 
#     AC104809.2         ARID3C          CREB3          HMGB1      KB.68A7.1       MIR135A1 
#              2              2              2              2              2              2 
#          MLXIP           MYCL   RP1.40E16.11  RP11.263G22.1  RP11.344P13.6    RP11.39H3.2 
#              2              2              2              2              2              2 
#  RP11.449P15.2   RP11.622A1.2   RP11.70D24.2   RP5.1159O4.2  RP5.875H18.10         ZBTB47 
#              2              2              2              2              2              2 
#         ZNF680  CTD.2325A15.5 
#              2              1 

top <- as.list(top[order(top, decreasing = T)])
top <- as.data.frame(top)
top <- t(top)
rownames(top) -> names
top <- as.data.frame(top)
top[,2] <- names
colnames(top) <- c("Count", "names")
top[, 3] <- top[, 1]/length(index.low)
colnames(top)[3] <- "Freq"


ggplot(top[1:20,], aes(x = reorder(names, -Freq), y = Freq)) +
  geom_bar(stat = "identity", aes(fill = "blue")) + 
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Gene names", y = "Frequency", title = "Interaction Pairs") +
  guides(fill=FALSE) +
  theme_gray(base_size = 15)+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

rownames(results) <- c("RD1", "RD2", "DEDUCTION")
edge <- cbind(t(gene.names.m), t(results))
colnames(edge)[1:2] <- c("from", "to")
edge <- as.data.frame(edge)
edge <- edge[order(edge$RD2), ]
write.csv(edge, "glm_edge_pairs.csv")


#---------------csv of edge----------------------

glm.m <- matrix(nrow = 3160, ncol = 2)
for(i in 1:3160)
{
  paste(glm[i,2],"-", glm[i,3],sep = "") -> glm.m[i,1]
}
glm.m <- as.matrix(glm.m)
glm[,4] -> glm.m[,2]
colnames(glm.m) <- c("gene","RD2")


si.m <- matrix(nrow = 3160, ncol = 2)
for(i in 1:3160)
{
  paste(SI[i,2],"-", SI[i,3], sep = "") -> si.m[i,1]
}
si.m <- as.matrix(si.m)
si.m[,2] <- SI[, 4]
colnames(si.m) <- c("gene", "SI")

edge <- merge(x=si.m, y=glm.m, by.x = "gene", by.y = "gene", all = FALSE)

gene <- strsplit(edge[,1], split = "-")

for (i in 1:3160) 
{
  gene[[i]][1] -> edge[i, 4]
  gene[[i]][2] -> edge[i, 5]
}
write.csv(edge,"merge_edge.csv")


