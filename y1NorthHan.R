## conditioal GWAS NorthHan 
## sample 2909 from GAB 2509(exHWP) merge 1000GHJD400  -- exclude(geno 0.03 HWP 0.001)
## P < 1e-6
genofile <- " --bfile ConditionalGWAS.GAB2509.1000G400"
phenofile <- " --pheno language8.sample2909"
phename <- " --pheno-name y1"
outpath <- "y1NorthHan/"
outname <- "ancestry.GAB2509.conditional.gwas1e-6.y1"
resultname <- "ancestry.GAB2509.conditional.gwas1e-6.y1.result"
setwd(outpath)
selectSNP <- matrix(NA, 20,4)
for(i in 1: 20){
    if(i ==1){
	   tmp <- paste0("plink",genofile,phenofile,phename," --logistic hide-covar --ci 0.95"," --out ",outpath,outname,"_",i)
	   system(tmp)
	   #filter P < 1e-6 save 
	   tmp <- paste0("awk '$12<1e-6 {print $2,$1,$3,$12}' ",outpath,outname,"_",i,".assoc.logistic"," > ",outpath,resultname,"_",i)
	   system(tmp)
	   # plotManhattan resultfile= resultname_i
	   resultfile <- paste0(resultname,"_",i)
	   dataManhattan <- read.table(resultfile, header = FALSE)
	   names(dataManhattan) <- c("SNP", "CHR", "BP", "P")
	 
	   # select top snpfile = resultname_i.selectSNP_i(snplist).
	   selectSNP[i,] <- as.matrix(dataManhattan[which(dataManhattan$P == min(dataManhattan$P)),])
       rm(dataManhattan)
       gc()
	   covarname <- as.character(selectSNP[i,1])
	   snpfile <- paste0(resultfile,".selectSNP_",i)
	   write.table(covarname, file = snpfile, quote = FALSE, row.names = FALSE)
	   # recode top SNP to raw 012file
	   tmp <- paste0("plink",genofile," --extract ",outpath,snpfile," --recodeA"," --out ",outpath,outname,"_",i)
	   system(tmp)
	   # covar file 
	   covarfile <- read.table(paste0(outname,"_",i,".raw"),header = TRUE)
	   covarname <- names(covarfile)[-c(1:6)]
	}else{
	   # plink covar(top snp)
	   tmp <- paste0("plink",genofile,phenofile,phename," --covar ",outpath,outname,"_",i-1,".raw"," --covar-name ",covarname," --logistic hide-covar --ci 0.95"," --out ",outpath,outname,"_",i)
       system(tmp)
	   #filter P < 1e-6 save 
	   tmp <- paste0("awk '$12<1e-6 {print $2,$1,$3,$12}' ",outpath,outname,"_",i,".assoc.logistic"," > ",outpath,resultname,"_",i)
	   system(tmp)
	   # plotManhattan resultfile= resultname_i
	   resultfile <- paste0(resultname,"_",i)
	   dataManhattan <- read.table(resultfile, header = FALSE)
	   # if no snp P <1e-6 break
	   if(nrow(dataManhattan)==0){
	               break
		               }
	   names(dataManhattan) <- c("SNP", "CHR", "BP", "P")
	   
	   # select top snpfile = resultname_i.selectSNP_i(snplist).
	   selectSNP[i,] <- as.matrix(dataManhattan[which(dataManhattan$P == min(dataManhattan$P)),][1,])
       rm(dataManhattan)
       gc()
	   covarname <- as.character(selectSNP[1:i,1])
	   snpfile <- paste0(resultfile,".selectSNP_",i)
	   write.table(covarname, file = snpfile, quote = FALSE, row.names = FALSE)
	   # recode top SNP to raw 012file
	   tmp <- paste0("plink",genofile," --extract ",outpath,snpfile," --recodeA"," --out ",outpath,outname,"_",i)
	   system(tmp)
	   # covar file 
	   covarfile <- read.table(paste0(outname,"_",i,".raw"),header = TRUE)
	   covarname <- names(covarfile)[-c(1:6)]
	    # replace X to " ", . to ":"(problem from read.table 8: -> X8.)
	   covarname <- gsub("X", " ", covarname)
	   covarname <- gsub("\\.", ":", covarname)
	   covarname <- paste0(covarname, collapse = ",") 
	}
write.table(selectSNP, file = paste0(outpath,"selectSNP.",outname), quote = FALSE, row.names = FALSE)
}