### load packages
#library(DirichletMultinomial)
library(lattice)
library(xtable)
library(parallel)
library(ggplot2)
library(scales)
library(plyr)
library(data.table)
library(dplyr)
library(tidyr)

#cst.path
parameters=paste0('k', k_i)
cst.path=file.path(output.out, 'CST', parameters)
dir.create(cst.path)

#pull out k of model fit from fit object 
(best <- fit[[k_i]])

#save important information
write.csv(mixturewt(best), file.path(cst.path, paste0('k', k_i,'pitheta.csv')))
write.csv(mixture(best), file.path(cst.path, paste0('k', k_i,"bestcomponent.csv")))
### plot contribution of each taxonomic group to community types
pdf(file.path(cst.path, paste0('k', k_i, "fitted.pdf")))
splom(log(fitted(best)))
dev.off()
### poster mean difference - how different are your community types from the population average (i.e. having just one community type)
# prep
p0 <- fitted(fit[[1]], scale=TRUE) # scale by theta
p1 <- fitted(best, scale=TRUE)
colnames(p1) <- paste("m", 1:k_i, sep="") 
# difference
(meandiff <- colSums(abs(p1 - as.vector(p0))))
sum(meandiff)
# Range: [0-200%] where 0% is completely identical and 200% is completely dissimilar to the population average.
write.csv(meandiff, file=file.path(cst.path, paste0('k_', k_i,"meandiff.csv")))

### taxonomic contributions to community types
diff <- rowSums(abs(p1 - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df <- head(cbind(Mean=p0[o], p1[o,], diff=diff[o], cdiff), dim(p1)[1])   # dim(p1)[1] is the number of taxa you have
write.csv(df, file=file.path(cst.path, paste0('k_', k_i, "taxcontribution.csv")))


comp<-read.csv(file.path(cst.path, paste0('k', k_i,'bestcomponent.csv')), row.names=1)
#assign community type
comp$comp<-apply(comp, 1, function(x) which(x==max(x)))
### save the probability of the newly assigned community type - pulling out the max probability among all community types per sample
comp$maxcomp<-NA
comp$maxcomp<-apply(comp[,1:k_i], 1, max, na.rm=TRUE)   
### compare the best community type with the next best community
# transform to matrix
probmatrix<-as.matrix(comp[,1:k_i])
# sort the matrix so the last column has the highest probability
probmatrix_sorted<-as.data.frame(t(apply(probmatrix, 1, sort)))
# subtract probability of best community type by probability of next best community type
probmatrix_sorted$probdiff<-probmatrix_sorted[,k_i]-probmatrix_sorted[,k_i-1]  
#plot probability of DMM component
pdf(file.path(cst.path, paste0('k', k_i, "DMMprobabilitybyCST.pdf")))
plot(maxcomp~as.factor(comp),data=comp, main="Probability of Dirichlet Component by Component",
     xlab="Component", ylab="Probability")
dev.off()

#cutoffs for misclassification
exclude.8<-row.names(subset(probmatrix_sorted, probmatrix_sorted[, k_i]<thresh.8 | probmatrix_sorted$probdiff<thresh2))    # change "2" to the last column
exclude.9<-row.names(subset(probmatrix_sorted, probmatrix_sorted[, k_i]<thresh.9 | probmatrix_sorted$probdiff<thresh2))    # change "2" to the last column
comp$CST_nocut<-comp$comp; comp$CST.8<-comp$comp; comp$CST.9<-comp$comp
comp$maxcomp.9<-comp$maxcomp; comp$maxcomp.8<-comp$maxcomp
comp$CST.8[row.names(comp) %in% exclude.8]<-NA; comp$CST.9[row.names(comp) %in% exclude.9]<-NA
comp$maxcomp.8[row.names(comp) %in% exclude.8]<-NA; comp$maxcomp.9[row.names(comp) %in% exclude.9]<-NA

# histogram
pdf(file.path(cst.path, paste0('k', k_i,'HistogramDMM.pdf')))
par(mfrow=c(3, 1))
hist(comp$maxcomp, breaks=100, main="Density Plot Probability of Dirichlet Component",
     xlab="Probability", ylab="Number of Samples")
hist(comp$maxcomp.8, breaks=100, main="Density Plot Probability of Dirichlet Component",
     xlab="Probability", ylab="Number of Samples")
hist(comp$maxcomp.9, breaks=100, main="Density Plot Probability of Dirichlet Component",
     xlab="Probability", ylab="Number of Samples")
dev.off()

# boxplot by community type
pdf(file.path(cst.path, paste0('k', k_i,'BoxplotDMM.pdf')))
par(mfrow=c(3, 1))
plot(maxcomp~as.factor(comp),data=comp, main="Probability of Dirichlet Component by Component",
     xlab="Component", ylab="Probability")
plot(maxcomp.8~as.factor(comp),data=comp, main="Probability of Dirichlet Component by Component",
     xlab="Component", ylab="Probability")
plot(maxcomp.9~as.factor(comp),data=comp, main="Probability of Dirichlet Component by Component",
     xlab="Component", ylab="Probability")
dev.off()

###### Make a dataset for CST and save #######
### create ID for comp dataset
if(my_data=='cohra'){
        comp$FoxCavID<-row.names(comp)
        CST_df<-comp %>% select(FoxCavID, CST_nocut, maxcomp, CST.8, maxcomp.8, CST.9, maxcomp.9)
}
if(my_data=='PRJ'){
        comp$SampleID<-row.names(comp)
        CST_df<-comp %>% select(SampleID, CST_nocut, maxcomp, CST.8, maxcomp.8, CST.9, maxcomp.9)
}
CST_df<-CST_df %>% mutate_at(vars(contains('CST')), funs(as.factor))
assign(paste0('CST_k', k_i) , CST_df)
save(list=c(paste0('CST_k', k_i)), file=paste(data.out, '/CST/k_', k_i, 'CST.Rdata', sep = ''))
