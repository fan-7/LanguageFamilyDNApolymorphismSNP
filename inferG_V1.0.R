# Source code for China province prediction tool 
# version: 1.0

geto <- function(df,newdata,tr=tr){
	# df: internal data consisting of all samples
	# newdata: if missing, infile is read for prediction
	# tr: threshold for AUC loss, predicted probabilities will not be reported if AUC loss > tr
	require(nnet)
	require(pROC)
	require(caret)
	cpginfo <- names(df)[-c(1)]
	dftrain <- df[!is.na(df$y),names(df) %in% cpginfo | names(df) %in% "y"]
	fo <- as.formula(paste("y~",paste(cpginfo,collapse="+")))

	getauc <- function(fo,data,ycol = 1){
		fit <- multinom(fo=fo,data=data)
		yp <- predict(fit,type="probs")
		if(is.null(dim(yp))){yp<-as.matrix(yp)}
		ny <- ncol(yp)
		auc <- rep(NA,ny)
		for (j in 1:ny){
		    yt <- ifelse(data[,ycol]==j,1,0) 
			roc_result <- roc(yt, yp[,j],auc=T, quiet = TRUE)
            auc[j] <- roc_result$auc  # AUC
		}
		return (auc)
	}
	auc <- getauc(fo,dftrain) 
	
	getsenspe <- function(fo, data, ycol = 1){
	    fit <- multinom(fo=fo,data=data)
		yp <- predict(fit)
		if(is.null(length(yp))){yp<-as.matrix(yp)}
		per_fit <-  confusionMatrix(yp,dftrain$y)
		senspe <- matrix(NA, 2, nrow(per_fit$byClass))
		#Sensitivity
		senspe[1,] <- per_fit$byClass[,1] 
		#Specificity
        senspe[2,] <- per_fit$byClass[,2] 
		return (senspe)
	}
    senspe <- getsenspe(fo,dftrain)	 		

	# prediction result
	o <- list() 
	o[[1]] <- matrix(NA,dim(newdata)[1],26)  
	colnames(o[[1]])<-c(
		"XinJiang",
		"NeiMeng",
		"ShanXi",
		"HeNan",
		"ShanDong",
		"JiangXi",
		"SiChuan",
		"GuangXi",
		"Full_AUC_XJ",
		"Full_AUC_NM",
		"Full_AUC_SX",
		"Full_AUC_HN",
		"Full_AUC_SD",
		"Full_AUC_JX",
		"Full_AUC_SC",
		"Full_AUC_GX",
		"Numb_missingCpGs",
		"Name_missingCpGs",
		"AUC_Loss_XJ",
		"AUC_Loss_NM",
		"AUC_Loss_SX",
		"AUC_Loss_HN",
		"AUC_Loss_SD",
		"AUC_Loss_JX",
		"AUC_Loss_SC",
		"AUC_Loss_GX"
	)
	o[[1]] <- as.data.frame(o[[1]])

	for (i in 1:nrow(newdata)){
		ix <- cpginfo %in% names(newdata)[!is.na(newdata[i,])]
		if(sum(ix)>0){
			fo <- as.formula(paste("y~",paste(cpginfo[ix],collapse="+")))
			fit <- multinom(fo,data=df) 
			o[[1]][i,1:8]<-predict(fit,newdata=newdata[i,],type="probs")
			o[[1]][i,9:16]<-auc
			o[[1]][i,17]<-sum(!ix)
			o[[1]][i,18]<-paste(cpginfo[!ix],collapse="/")
			aucnow <- getauc(fo,df)
			o[[1]][i,19:26]<-auc-aucnow
		}	
	}
	if(i==1){
		    if(any(o[[1]][1,19:26] > tr)){o[[1]][1,1:8]<-NA} 
	}else{
		ix <- apply(o[[1]][,19:26]>tr,1,function(x)any(x,na.rm=T))
		o[[1]][ix,1:8] <- NA
	}

	# output
	out <- newdata[,!names(newdata) %in% cpginfo,drop=F]
	out <- cbind(out,o[[1]])
	return(out)
}
load("inferG.RData")
input <- read.csv("input.csv", header = TRUE)
result <- geto(df=inferG,newdata = input, tr=0.5)
write.csv(result,"Result.csv",row.names=F,quote=F)

