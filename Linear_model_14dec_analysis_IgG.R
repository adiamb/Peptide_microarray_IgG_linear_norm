require(data.table)
require(readr)
require(plyr)
require(dplyr)
require(samr)
require(glmnet)
require(reshape2)
##function to check unique items in each column
apply(Ignorm6, MARGIN = 2, FUN = function(x) (unique(x) %>% length()))
################################################
Igraw = fread('~/Documents/IgG_Analysis/Raw_aggregate_data.tab', showProgress = T)
IGRAWprobes = select(Igraw, PROBE_ID, PROBE_SEQUENCE, CONTAINER)
Ignorm = fread('~/Documents/IgG_Analysis/linear_normalized/Linear_model_dec1/normalized_dec13_IgG.csv')
sampleIDS = read_csv('~/Documents/PA-00075 EXPERIMENTAL LAYOUT-LING_Nov29.csv')
sampleIDS$sampleID = as.character(interaction(sampleIDS$Array_ID, sampleIDS$Subarray)) ##get unique IDS
Ignorm$V1 = NULL
Ignorm$value = NULL
###read in the peptide sequences
flu_peps = fread('~/Documents/IgG_Analysis/FLU_ARRAYS_PEPTIDE_SEQ.csv', showProgress = T)
flu_peps$V1 = NULL
Ignorm = as.data.frame(Ignorm)
Ignorm = transform(Ignorm,  ID= colsplit(uniqueProbeID, "\\.", names = c("CONTAINER", "PROBE_ID")))
head(Ignorm)

dim(flu_peps)
unique(flu_peps$PROBE_SEQUENCE) %>% length()
uniqueIGrawprobes=unique(IGRAWprobes) ##select only uniqe probes and peptides

###merge the peptide sequences into the main DF
Ignorm4 = merge.data.frame(Ignorm, uniqueIGrawprobes, by.x = "ID.PROBE_ID", by.y = "PROBE_ID") ##merge peptide seq and probe IDS
Ignorm4$sampleID = as.character(interaction(Ignorm4$Array, Ignorm4$Subarray))
Ignorm5 = merge.data.frame(Ignorm4, sampleIDS, by.x = "sampleID", by.y = "sampleID") ##merge the sample IDS into main DF
Ignorm5$Subarray.x = NULL
Ignorm5$Subarray.y = NULL
Ignorm5$Array_ID = NULL
Ignorm5$Array = NULL
setDT(Ignorm5)
Ignorm6=Ignorm5[, list(norm_signal_mean = mean(residuals)), by=list(PROBE_SEQUENCE, Customer_Sample_ID, DBID, Dx, Gender, Age, Interval, sample_group)]
####move the DF into wide format
Ignorm7=dcast(Ignorm6, Customer_Sample_ID+DBID+Dx+Gender +Age +Interval +sample_group~PROBE_SEQUENCE , value.var = "norm_signal_mean")

########remove false positives ###############################################
sd_mean=function(x) {
  y=mean(x)+2*sd(x)
  ifelse(x >= y, return(x), NA)
}

temp2 = as.data.frame(apply(Ignorm7[8:ncol(Ignorm7)], MARGIN = 2, FUN = sd_mean))
temp4=temp2[, colSums(is.na(temp2)) != nrow(temp2)] ###real responses
temp5=temp2[, colSums(is.na(temp2)) == nrow(temp2)] ##all discarded false positives 
###join the IDS into temp4
Ignorm8 = cbind.data.frame(Ignorm7[1:7], temp4)
Ignorm8$Diagnosis = ifelse(Ignorm8$Dx == "K", 0, 1)
require(samr)
require(glmnet)
############################################all samples| all proteins #####################################################
set.seed(2086)
y = ifelse(Ignorm8$Dx == "K", 1, 2)
x = t(Ignorm8[,8:138841])
d <- list(x=x, y=y, geneid=row.names(x), logged2 = TRUE)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=100, assay.type = "array", center.arrays = T)

delta=0.1
samr.plot(samr.obj,delta)
delta.table <- samr.compute.delta.table(samr.obj)
delta.table

### create significant genes table
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
up=as.data.frame(siggenes.table$genes.up)
colnames(up)[3]=c("PROBE_SEQUENCE")
grep_pep2(up$PROBE_SEQUENCE)

### plot FDR as a funcion of score d and add dashed line cutoff at 10% FDR
plot(siggenes.table$genes.lo[,8]~siggenes.table$genes.lo[,4], xlab = "Score(d)", ylab = "FDR")
abline(h=10,col=4,lty=2)
require(ggplot2)
ggplot(recent, aes(factor(y), recent$NKVQLQFSEKQMLTTS, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
ggplot(PX, aes(factor(Dx), PX$ENVETMNSNTLELRSR, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
ggplot(Ignorm8, aes(factor(Dx), KQEELLKNIQNAHQDF, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
######################################LIMMA all samples | all proteins ####################################################
require(Biobase)
require(limma)
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimesnions are right !
design = model.matrix(~Ignorm8$Diagnosis) 
fit = eBayes(lmFit(object, design))
we=topTable(fit, coef="Ignorm8$Diagnosis", adjust="BH", number = 15)
we$PROBE_SEQUENCE = rownames(we)
grep_pep(rownames(we))
################################divide into two groups######################################################################
##recent onset SAMR####################################
xyz=split.data.frame(Ignorm8, Ignorm8$sample_group)
PX = xyz$pandemrix ##recent onset 
recent = xyz$`early onset` ## pandemrix associated 
xyz = NULL
y = ifelse(recent$Dx == "K", 1, 2)
x = t(recent[,8:138841])
d <- list(x=x, y=y, geneid=row.names(x), logged2 = TRUE)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=100, assay.type = "array", center.arrays = T)

delta=0.2
samr.plot(samr.obj,delta)
delta.table <- samr.compute.delta.table(samr.obj)
delta.table

### create significant genes table
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
up=as.data.frame(siggenes.table$genes.up)
colnames(up)[3]=c("PROBE_SEQUENCE")
grep_pep2(up$PROBE_SEQUENCE)
######################################LIMMA recent onset samples | all proteins ####################################################
require(Biobase)
require(limma)
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimesnions are right !
y = as.factor(ifelse(recent$Dx == "K", 0, 1))
design = model.matrix(~y) 
fit = eBayes(lmFit(object, design))
we=topTable(fit, coef="y1", adjust="BH", number = 15)
we$PROBE_SEQUENCE = rownames(we)
grep_pep(rownames(we))

##Pandemrix SAMR #######################################

y = ifelse(PX$Dx == "K", 1, 2)
x = t(PX[,8:138841])
d <- list(x=x, y=y, geneid=row.names(x), logged2 = TRUE)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=100, assay.type = "array", center.arrays = T)
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.2
samr.plot(samr.obj,delta)


### create significant genes table
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
up=as.data.frame(siggenes.table$genes.up)
low = as.data.frame(siggenes.table$genes.lo)
colnames(low)[3]=c("PROBE_SEQUENCE")
colnames(up)[3]=c("PROBE_SEQUENCE")
grep_pep_up(up$PROBE_SEQUENCE)
grep_pep_low(low$PROBE_SEQUENCE)

######################################LIMMA recent onset samples | all proteins ####################################################
require(Biobase)
require(limma)
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimesnions are right !
y = as.factor(ifelse(PX$Dx == "K", 0, 1))
design = model.matrix(~y) 
fit = eBayes(lmFit(object, design))
we=topTable(fit, coef="y1", adjust="BH", number = 30)
we$PROBE_SEQUENCE = rownames(we)
unique(grep_pep(rownames(we)))

###################################### SAMR on differences between recent onset and pandemrix #######################################
abc = split.data.frame(Ignorm8, Ignorm8$Diagnosis)
cases_only = abc$`1`
abc = NULL

y = ifelse(cases_only$sample_group == "early onset", 1, 2)
x = t(cases_only[,8:138841])
d <- list(x=x, y=y, geneid=row.names(x), logged2 = TRUE)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=100, assay.type = "array", center.arrays = T)
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=0.18
samr.plot(samr.obj,delta)


### create significant genes table
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
up=as.data.frame(siggenes.table$genes.up)
low = as.data.frame(siggenes.table$genes.lo)
colnames(low)[3]=c("PROBE_SEQUENCE")
colnames(up)[3]=c("PROBE_SEQUENCE")
grep_pep_up(up$PROBE_SEQUENCE)
grep_pep_low(low$PROBE_SEQUENCE)


######################################LIMMA recent onset samples | all proteins ####################################################
require(Biobase)
require(limma)
object<-new("ExpressionSet", exprs=as.matrix(x)) #this is the x from above in the SAMR
object #see if the dimesnions are right !
y = as.factor(ifelse(cases_only$sample_group == "early onset", 1, 0))
design = model.matrix(~cases_only$sample_group+cases_only$Gender+cases_only$Age) 
fit = eBayes(lmFit(object, design))
we=topTable(fit, coef="cases_only$sample_grouppandemrix", adjust="BH", number = 30)
we$PROBE_SEQUENCE = rownames(we)
unique(grep_pep(rownames(we)))
##################################
ggplot(cases_only, aes(factor(cases_only$sample_group)))
###condense IGrawprobes #####
PROBES1 = as.data.frame(PROBES1)
IGRAWprobes = as.data.frame(IGRAWprobes)
IGRAWprobes1 = transform(IGRAWprobes, CONTAINER = colsplit(CONTAINER, "_", names = c("CONT1", "CONT2")))
IGRAWprobes1$CONT1 = IGRAWprobes1$CONTAINER$CONT1
IGRAWprobes1$CONTAINER=NULL
setDT(IGRAWprobes1)
PROBES1=IGRAWprobes1[, list(PROBE_ID = paste(unique(PROBE_ID), collapse = "/"), CONT1 = paste(unique(CONT1), collapse = "/")), by=list(PROBE_SEQUENCE)]
pH1_seq=PROBES1$PROBE_SEQUENCE[grep("H1N1", PROBES1$PROBE_ID)] #pH1N1
gH1_seq=PROBES1$PROBE_SEQUENCE[grep("gH1N1", PROBES1$PROBE_ID)] #gH1N1
H3_seq=PROBES1$PROBE_SEQUENCE[grep("H3N2", PROBES1$PROBE_ID)]  #H3N2
non_flu = PROBES2$PROBE_SEQUENCE[grepl("gH1N1|pH1N1|H3N2|Vacc", PROBES2$PROBEID2)]
flu_peps1 = PROBES1$PROBE_SEQUENCE[grep("FLU", PROBES1$CONT1)]
gallus=PROBES2$PROBE_SEQUENCE[grep("Gallus*", PROBES2$Protein_fullname)]
##read in text file containng protein IDS from python"
proteinIDS_seqs=read_delim('~/Desktop/ProteinIDS_seq.txt', delim = ";", col_names = F)
colnames(proteinIDS_seqs) = c("Protein_ID", "Protein_fullname", "PROBE_SEQUENCE")

PROBES2 = merge.data.frame(PROBES1, proteinIDS_seqs, by.x = "PROBE_SEQUENCE", by.y = "PROBE_SEQUENCE", all.x = T)
flu1 = PROBES2$PROBE_SEQUENCE[grep("Influenza*", PROBES2$Protein_fullname)]

PROBES2[grep("NKVQLQFSEKQMLTTS", PROBES2$PROBE_SEQUENCE),]
####only for limma ooutput only### query for peptide names#######
grep_pep = function(x) {
  qw= PROBES2[PROBES2$PROBE_SEQUENCE %in% x,]
  we1 = merge.data.frame(we, qw, by.x = "PROBE_SEQUENCE", by.y = "PROBE_SEQUENCE", sort = F)
  return(we1)
}
############only for SAMR oytput query###############  
grep_pep_up = function(x) {
  qw = PROBES2[PROBES2$PROBE_SEQUENCE %in% x,]
  we2 = merge.data.frame(up, qw, by.x = "PROBE_SEQUENCE", by.y = "PROBE_SEQUENCE", sort = F, all.x = T)
  return(we2)
}

grep_pep_low = function(x) {
  qw = PROBES2[PROBES2$PROBE_SEQUENCE %in% x,]
  we2 = merge.data.frame(low, qw, by.x = "PROBE_SEQUENCE", by.y = "PROBE_SEQUENCE", sort = F, all.x = T)
  return(we2)
}
#############read HI titers############
HI_titers = read_csv('~/Dropbox/HI_titers_aus.csv')
sa_HI=HI_titers[HI_titers$DbID %in% Ignorm8$DBID,]
sa=Ignorm8[Ignorm8$DBID %in% HI_titers$DbID,]
sa_HI = select(sa_HI, DbID, Titer)
sa_HI$Titer[sa_HI$Titer == "<10"] = 0
sa_HI$DBID = as.character(sa_HI$DbID)
sa_HI$Titer = as.numeric(sa_HI$Titer)
sa_HI$DbID = NULL

##merge HI titers ###
sa = as.data.frame(sa)
sa = merge.data.frame(sa, sa_HI, by.x = "DBID", by.y = "DBID")
y = ifelse(sa$Titer >= 40, 1, 2)
x = t(Ignorm8[colnames(Ignorm8) %in% gallus])
x = t(cases_only[-c(1:7, 138842)])
d <- list(x=x, y=y, geneid=row.names(x), logged2 = TRUE)
samr.obj<-samr(d,  resp.type="Two class unpaired", nperms=100, assay.type = "array")

delta=0.2
samr.plot(samr.obj,delta)
delta.table <- samr.compute.delta.table(samr.obj)
delta.table

### create significant genes table
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
#################################limma analysis #####################################################
require(Biobase)
require(limma)
require(magrittr)
object<-new("ExpressionSet", exprs=as.matrix(x))
object
design = model.matrix(~Ignorm8$Gender+Ignorm8$Dx)
fit = eBayes(lmFit(object, design))
we=topTable(fit, coef="Ignorm8$DxK", adjust="BH", number = 10)
we$PROBE_SEQUENCE = rownames(we)
grep_pep(rownames(we))
volcanoplot(fit, coef = "y1")



















ggplot(Ignorm8, aes(factor(Diagnosis), IDGVKLESTGIYQILA, color = factor(sample_group)))+geom_point(size = 6)+stat_summary(fun.y = "median", fun.ymax = "median", fun.ymin = "median", geom = "crossbar", color = "blue")
q = ifelse(sa$Titer >= 40, "High", "Low")
