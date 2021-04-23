suppressMessages(library(dplyr, quietly=T))
suppressMessages(library(qvalue, quietly=T))

#usage: Rscript positive_selection_LRT_FDR.R <species>
#eg: Rscript positive_selection_LRT_FDR.R Nnot

#get arguments from the command line
args<-commandArgs(TRUE)
stopifnot(length(args)==1)
if(length(args)==0){
  print("No sample specified; cannot make plot")
} else if(length(args)==1){
  sp<-as.character(args[1])
} else {
  print("Too many arguments given from the command line")
}

print(paste0("Working on branch ",sp))
# getting input files
infile_alt<-paste0("alt_",sp,"_all.lnl")
infile_null<-paste0("null_",sp,"_all.lnl")

read.table(infile_alt, header = T, sep = "\t") -> alt
read.table(infile_null, header = T, sep = "\t") -> null

# selecting best lnL from replicates
apply(alt[,2:4], 1, FUN=max) -> best_alt
apply(null[,2:4], 1, FUN=max) -> best_null
lnl <- tibble(alt$GENE, best_alt, best_null)

# calculating LRT statistics and p-value 
lnl$LRT <- 2*(lnl$best_alt - lnl$best_null)
lnl$p <- pchisq(lnl$LRT, df=1, lower.tail = FALSE)
# applying FDR correction
qobj <- qvalue(lnl$p, pi0.meth="bootstrap", fdr.level=0.05)
lnl <- cbind(lnl, qobj$qvalues, qobj$significant)
filter(lnl, qobj$significant==TRUE) -> sig

print(paste0("Number of genes significantly under PS in branch ",sp,":"))
nrow(sig)

#outputting files
outfile_name<-paste0("sig_PS_genes_",sp,".txt")
write.table(sig[,1:6], outfile_name, sep = "\t", row.names = F, 
            col.names = c("GENE", "alt_lnL", "null_lnL", "LRT", "p-value", "q-value"), quote = F)
