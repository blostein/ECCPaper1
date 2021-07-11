####### Set Up #######
setwd("/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/DMM")
library(DirichletMultinomial)
library(phyloseq)
library(lattice)
library(xtable)
library(parallel)
options(width=70, digits=2)
.qualitative <- DirichletMultinomial:::.qualitative
dev.off <- function(...) invisible(grDevices::dev.off(...))

### set display parameters
# set width of R output to 70 characters and number of floating point digits displayed to two
options(width=70, digits=2)

# use .qualitative color set
.qualitative <- DirichletMultinomial:::.qualitative

# dev.off redefined to return without displaying values
dev.off <- function(...) invisible(grDevices::dev.off(...))

# WARNING: restart R to set parameters back to default before working on a different script
remove(list= ls())

load("/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/DMM/phyasv_postTax.Rdata")
ps.aggregate = ps.filter
tax = data.frame(tax_table(ps.aggregate)[,1:7])
rowname = rownames(tax)

sciname = tax[,6:7]
rownames(sciname) = rowname
scinames <-paste(sciname$Genus, ifelse(is.na(sciname$Species), rowname, sciname$Species), sep="..")
asvnames<-row.names(ps.aggregate@tax_table)
###### Data #######
### read in your count dataset - option 2: from csv (e.g. a collaborator sends you a count table)
#data<-t(as.matrix(read.csv("DMM_ready.csv", row.names=1)))
otu= ps.aggregate@otu_table
otu= as.data.frame(otu)
data= otu
colnames(data)<-asvnames
data<-as.matrix(data)
#For DMM in flux, we want the rows to be the samples, columns to be the taxa
#cnts <- log10(colSums(data))
#densityplot(cnts, xlim=range(cnts), xlab="Taxon representation (log 10 count)")
#data=as.matrix(data)  ###remember, it should be a matrix before runing mclapply
#data=t(as.matrix(data))

###### Identifying number of community types in study population ######
### fit count data to a series of models with 1 to 5 community types

fit <- mclapply(1:10, dmn, count=data, verbose=TRUE)
# important: the higher the maximum number of community types, the higher the computing time (exponentially increases)
# I've set the maximum number to 5 due to computational time and because the template dataset is tiny.
# You may want to consider going up to at least 10 - WARNING: will take you all day unless you have a supercomputer
# or submit to flux.
fit_asvtrimmed<-fit
### save the models so you don't have to rerun clustering later
save(fit_asvtrimmed, file="DMM_fit_asvtrim.rda")
