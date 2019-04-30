#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-r", "--normalized_read_counts_csv_path"), type = "character", default = NULL,
              help = "Path to normalized read counts table", metavar = "character"),
  make_option(c("-c", "--contrast"), type = "character", default = NULL,
              help = "Contrast, pair of conditions ", metavar = "character"),
  make_option(c("-t", "--sample_file_path"), type = "character", default = NULL,
              help = "Path to sample.tsv", metavar = "character"),
  make_option(c("-x", "--xtail_result_path"), type = "character", default = "NULL",
              help = "Path for writing xtail result file", metavar = "character"),
  make_option(c("-f", "--xtail_fcplot_path"), type = "character", default = "NULL",
              help = "Path for writing xtail fc plot file", metavar = "character"),
  make_option(c("-p", "--xtail_rplot_path"), type = "character", default = "NULL",
              help = "Path for writing xtail rplot file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$normalized_read_counts_csv_path)){
  print_help(option_parser)
  stop("Please supply arguments (-r, -t, -x), see --help \n", call.=FALSE)
}

library(xtail)
library(DESeq2)

#dispersionMatrix.DESeqDataSet <- function(object){
#  if (!"dispersionMatrix" %in% names(assays(object))) return (NULL)
#  disp <- assays(object)[["dispersionMatrix"]]
#  colnames(disp) <- colnames(object)
#  disp
#}

estimateMLEForBeta <- function(object,modelMatrix=NULL,modelMatrixType="standard",
                           maxit=100, useOptim=TRUE, quiet=FALSE,useQR=TRUE) {
  if (is.null(dispersionMatrix(object))) {
    stop("testing requires dispersionMatrix estimates, first call estimateDispersions()")
  }
  stopifnot(length(maxit)==1)
  
  # in case the class of the mcols(mcols(object)) are not character
  object <- sanitizeRowData(object)
  
  if (is.null(mcols(object)$allZero)) {
    object <- getBaseMeansAndVariances(object)
  }

  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]

  if (is.null(modelMatrix)) {
    modelAsFormula <- TRUE
    termsOrder <- attr(terms.formula(design(object)),"order")
    interactionPresent <- any(termsOrder > 1)
    blindDesign <- design(object) == formula(~ 1)

    # store modelMatrixType so it can be accessed by estimateBetaPriorVar
    attr(object, "modelMatrixType") <- modelMatrixType
    hasIntercept <- attr(terms(design(object)),"intercept") == 1
    renameCols <- hasIntercept
  } else {
    message("using supplied model matrix")
    modelAsFormula <- FALSE
    attr(object, "modelMatrixType") <- "user-supplied"
    renameCols <- FALSE
  }


  # fit the negative binomial GLM without a prior
  # (in actuality a very wide prior with standard deviation 1e3 on log2 fold changes)
  fit <- fitNbinomGLMs(objectNZ, maxit=maxit, useOptim=useOptim, useQR=useQR,
                       renameCols=renameCols, modelMatrix=modelMatrix)
  H <- fit$hat_diagonals
  modelMatrix <- fit$modelMatrix
  modelMatrixNames <- fit$modelMatrixNames
  # record the wide prior variance which was used in fitting
  betaPriorVar <- rep(1e6, ncol(fit$modelMatrix))


  # store mu in case the user did not call estimateDispersionsGeneEst
  dimnames(fit$mu) <- NULL
  assays(objectNZ)[["mu"]] <- fit$mu
  assays(object)[["mu"]] <- buildMatrixWithNARows(fit$mu, mcols(object)$allZero)

  # store the prior variance directly as an attribute
  # of the DESeqDataSet object, so it can be pulled later by
  # the results function (necessary for setting max Cook's distance)
  attr(object,"betaPrior") <- FALSE
  attr(object,"betaPriorVar") <- betaPriorVar
  attr(object,"modelMatrix") <- modelMatrix

  # add betas to the object
  modelMatrixNames <- colnames(modelMatrix)
  betaMatrix <- fit$betaMatrix
  colnames(betaMatrix) <- modelMatrixNames
  betaSE <- fit$betaSE
  colnames(betaSE) <- paste0("SE_",modelMatrixNames)
  betaConv <- fit$betaConv

  if (any(!betaConv)) {
    if (!quiet) message(paste(sum(!betaConv),"rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument"))
  }

  
  resultsList <- c(matrixToList(betaMatrix),
                   matrixToList(betaSE),
                   list(betaConv = betaConv,
                        betaIter = fit$betaIter,
                        deviance = -2 * fit$logLike))
  
  Results <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
  
  modelMatrixNamesSpaces <- gsub("_"," ",modelMatrixNames)

  lfcType <- "MLE"
  coefInfo <- paste(paste0("log2 fold change (",lfcType,"):"),modelMatrixNamesSpaces)
  seInfo <- paste("standard error:",modelMatrixNamesSpaces)

  mcols(Results) <- DataFrame(type = rep("results",ncol(Results)),
                                  description = c(coefInfo, seInfo,
                                    "convergence of betas",
                                    "iterations for betas",
                                    "deviance for the fitted model"))
 
  mcols(object) <- cbind(mcols(object),Results)
  return(object)
}

dispersionMatrix <- function(object){
  if (!"dispersionMatrix" %in% names(assays(object))) return (NULL)
  disp <- assays(object)[["dispersionMatrix"]]
  colnames(disp) <- colnames(object)
  disp
}


xtailDispersionFit <- function(rawmeans, rawdisps, quantilePercent=0.35, binnum=50, minDisp=1e-8){
  useForFit <- rawdisps > 100 * minDisp
  if(sum(useForFit) == 0){
    return (rawdisps)
    stop("all gene-wise dispersion estimates are within 2 orders of magnitude
    from the minimum value, and so the gene-wise estimates as final estimates.")
  }
  usedMeans <- rawmeans[useForFit]
  usedDisps <- rawdisps[useForFit]
  sortMeans <- usedMeans[order(usedMeans)]
  sortDisps <- usedDisps[order(usedMeans)]
  genenum <- length(sortMeans)
  quantile_means <- c()
  quantile_disps <- c()
  for(i in seq(1,genenum,binnum)){
    num = min(binnum,genenum-i+1)
    if(num<(binnum)){
      next
    }
    curbin_means <- sortMeans[i:(i+num-1)]
    curbin_disps <- sortDisps[i:(i+num-1)]

    m <- curbin_means[order(curbin_disps)][1:floor(num*quantilePercent)]
    d <- curbin_disps[order(curbin_disps)][1:floor(num*quantilePercent)]
    quantile_means <- c(quantile_means, m)
    quantile_disps <- c(quantile_disps, d)
  }
  dat <- data.frame(logDisps = log(quantile_disps), logMeans=log(quantile_means))
  fit <- locfit(logDisps~logMeans, data=dat[quantile_disps>=minDisp*10,,drop=FALSE])
  dispFit <- function(rawmeans)exp(predict(fit,data.frame(logMeans=log(rawmeans))))
  dispFit
}

xTest <- function(object1, object2,threads,bins,baseLevel, ci){
	intersect.genes <- intersect(rownames(object1), rownames(object2))
	object1 <- object1[intersect.genes,,drop=FALSE]
	object2 <- object2[intersect.genes,,drop=FALSE]
	counts1 <- counts(object1)
	counts1[which(counts1==0)] = 1
	counts2 <- counts(object2)
	counts2[which(counts2==0)] = 1
	intercept1 <- mcols(object1)$Intercept
	intercept2 <- mcols(object2)$Intercept
	dispersion1 <- dispersionMatrix(object1)
	dispersion2 <- dispersionMatrix(object2)
	resName1 <- resultsNames(object1)[2]
	resName2 <- resultsNames(object2)[2]
	log2Ratio1 <- mcols(object1)[[resName1]]
	log2Ratio2 <- mcols(object2)[[resName2]]
	sizefactor1 <- sizeFactors(object1)
	sizefactor2 <- sizeFactors(object2)
	cond1 <- as.numeric(colData(object1)$condition) - 1
	cond2 <- as.numeric(colData(object2)$condition) - 1
	## Estimate the probabilities of log2R and log2FC
	## parallel used for speeding up
	rowNo <- nrow(counts1)
	if (is.na(threads)){
		## automaticly detect the number of CPU cores.
		cluster <- makeCluster(detectCores())
	}else{
		cluster <- makeCluster(threads)
	}
	clusterExport(cluster, c('counts1','counts2','intercept1','intercept2','log2Ratio1','log2Ratio2','dispersion1','dispersion2','sizefactor1','sizefactor2','cond1','cond2','bins','ci'),envir=environment())
	res <- clusterMap(cluster, xTestWrapper, i = c(1:rowNo));
	stopCluster(cluster)
	res <- matrix(unlist(res),ncol=4,byrow=T,dimnames=list(intersect.genes, c("deltaTE","Pval","lowci","highci")))
	res <- data.frame(res)
	res$log2Ratio1 <- log2Ratio1
	res$log2Ratio2 <- log2Ratio2
	res
}


xTestWrapper <- function(i){
	res <- xtail_test(counts1[i,],counts2[i,],intercept1[i],intercept2[i],log2Ratio1[i],log2Ratio2[i],dispersion1[i,],dispersion2[i,],sizefactor1,sizefactor2,cond1,cond2,bins,ci)
}

estimateFun <- function(countData, condition, baseLevel, libsize, dispers){
	## using the DESeqDataSet to store data
	colData <- data.frame(row.names = colnames(countData), condition=condition)
	dataSet <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
	colData(dataSet)$condition <- relevel(colData(dataSet)$condition, baseLevel)

	# normalization
	if (missing(libsize)){
		dataSet <- suppressMessages(estimateSizeFactors(dataSet))
	}else{
		sizeFactors(dataSet) <- libsize
	}

	# estimateDispersion
	if (missing(dispers))
	{
		if (ncol(countData) == 2){
			object = dataSet #creat a copy of a dataSet for estimating dispersion
			design(object) = formula(~1) # take two conditon samples as replicates for estimating dispersion
			object = estimateDispersionsGeneEst(object)
			#fit
			#dispFitFun <- xtailDispersionFit(mcols(object)$baseMean, mcols(object)$dispGeneEst)
                        fittedDisps <- xtailDispersionFit(mcols(object)$baseMean, mcols(object)$dispGeneEst)
			#attr(dispFitFun, "fitType") = "bin quantile"
			#fittedDisps <- dispFitFun(mcols(object)$baseMean)
                        #fittedDisps <- xtailDispersionFit(mcols(object)$baseMean)
			#dispersionMatrix(dataSet) = matrix(rep(fittedDisps, ncol(dataSet)), ncol=ncol(dataSet), byrow=FALSE)
                        dispersionMatrix = matrix(rep(fittedDisps, ncol(dataSet)), ncol=ncol(dataSet), byrow=FALSE)
		}else{
			dataSet <- suppressMessages(estimateDispersions(dataSet))
			#dispersionMatrix(dataSet) <- matrix(rep(dispersions(dataSet), ncol(dataSet)),ncol=ncol(dataSet), byrow=FALSE)
                        dispersionMatrix <- matrix(rep(dispersions(dataSet), ncol(dataSet)),ncol=ncol(dataSet), byrow=FALSE)
		}
	}else{
		#dispersionMatrix(dataSet) <- dispers
                dispersionMatrix <- dispers
	}

	# estimate log2FC or log2R
	dataSet <- suppressMessages(estimateMLEForBeta(dataSet, modelMatrixType="standard"))
	dataSet
}

xRIBOtail <- function(mrna, rpf, condition, baseLevel = NA, minMeanCount = 1, normalize = TRUE, p.adjust.method ="BH", threads=NA,bins=10000,ci = 0)
{
	## default baseLevel is the first coniditon
	if (is.na(baseLevel)) {baseLevel <- condition[1]}
	## if the baseLevel is in conditon
	if (!is.element(baseLevel, condition)) stop("baseLevel is not in condition")
	## if the colnames of mrna and rpf are matched
	if (ncol(mrna) != ncol(rpf)) stop("the mrna and rpf must have same number of columns")
	## if the colnames and condition are match
	if(length(condition)!=ncol(mrna)) stop("condition must have same length as the number of columns of rna")
	## ther must be exactly two different conditon
	if (length(unique(condition))!=2) stop("There must be exactly two different conditions")
	## check confidence interval
	if (ci<0 ) stop("ci must be non-negative")
	## check pvalue.adjust method
	supported.padj <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr")
	if (!is.element(p.adjust.method, supported.padj)) stop(paste0("The adjustment methods must be one of ",paste(supported.padj,collapse=", ")))

	## filter genes with low count number
	if (minMeanCount<1) stop("minMeanCount needs to be at least 1")
	keep.genes <- rowMeans(mrna) >= minMeanCount
	mrna.keep <- rownames(mrna)[keep.genes]
	keep.genes <- rowMeans(rpf) >= minMeanCount
	rpf.keep <- rownames(rpf)[keep.genes]
	## merge genes in mrna and rpf
	keep.genes <- intersect(mrna.keep,rpf.keep)
	mrna <- mrna[keep.genes,,drop=FALSE]
	rpf <- rpf[keep.genes,,drop=FALSE]

	## if the colnames of mrna and rpf are same, add characters to distinguish them.
	if (sum(colnames(mrna) == colnames(rpf)) >= 1){
		colnames(mrna) = paste0("mrna_",colnames(mrna))
		colnames(rpf) = paste0("rpf_",colnames(rpf))
	}
	#normalize together
	if (normalize){
		message("Calculating the library size factors")
		pool_sizeFactor <- estimateSizeFactorsForMatrix(cbind(mrna,rpf))
		mrna_sizeFactor <- pool_sizeFactor[1:ncol(mrna)]
		rpf_sizeFactor <- pool_sizeFactor[(ncol(mrna)+1):(ncol(mrna)+ncol(rpf))]
	}else{
		mrna <- round(mrna)
		rpf <- round(rpf)
		mrna_sizeFactor <- rep(1, ncol(mrna))
		rpf_sizeFactor <- rep(1,ncol(rpf))
	}

	## 1. Estimate the difference of log2FC between mRNA and RPF
	#If no replicate, we assume no more than one-third of genes' expression changed.
	message ("1. Estimate the log2 fold change in mrna")
	mrna_object = estimateFun(mrna,condition,baseLevel,mrna_sizeFactor)
	message ("2. Estimate the log2 fold change in rpf")
	rpf_object = estimateFun(rpf,condition,baseLevel,rpf_sizeFactor)
	message ("3. Estimate the difference between two log2 fold changes")
	result_log2FC = xTest(mrna_object,rpf_object,threads,bins,baseLevel,ci)

	### 2. Estimate the difference of log2R between control and treatment
	condition1_mrna <- mrna[,condition==baseLevel,drop=FALSE]
	condition1_rpf <- rpf[,condition==baseLevel,drop=FALSE]
	condition2_mrna <- mrna[,condition!=baseLevel,drop=FALSE]
	condition2_rpf <- rpf[,condition!=baseLevel,drop=FALSE]
	condition1_counts <- cbind(condition1_mrna,condition1_rpf)
	condition2_counts <- cbind(condition2_mrna,condition2_rpf)

	condition1_sizeFactor <- c(mrna_sizeFactor[condition==baseLevel],rpf_sizeFactor[condition==baseLevel])
	condition2_sizeFactor <- c(mrna_sizeFactor[condition!=baseLevel],rpf_sizeFactor[condition!=baseLevel])
        print("sizefkt1\n")
        print(condition1_sizeFactor)
        print("sizefkt2\n")
        print(condition2_sizeFactor)
	condition1_disper <- cbind(dispersionMatrix(rpf_object)[,condition==baseLevel,drop=FALSE], dispersionMatrix(rpf_object)[,condition==baseLevel,drop=FALSE])
	condition2_disper <- cbind(dispersionMatrix(rpf_object)[,condition!=baseLevel,drop=FALSE], dispersionMatrix(rpf_object)[,condition!=baseLevel,drop=FALSE])
	condition1 <- c(paste0(condition[condition == baseLevel],"_mRNA"), paste0(condition[condition == baseLevel],"_rpf"))
	baseLevel1 <- paste0(baseLevel,"_mRNA")
	condition2 <- c(paste0(condition[condition != baseLevel],"_mRNA"), paste0(condition[condition != baseLevel],"_rpf"))
	baseLevel2 <- paste0(unique(condition)[2],"_mRNA")
	#
	message ("4. Estimate the log2 ratio in first condition")
	condition1_object <- estimateFun(condition1_counts,condition1,baseLevel1,condition1_sizeFactor,condition1_disper)
	message ("5. Estimate the log2 ratio in second condition")
	condition2_object <- estimateFun(condition2_counts,condition2,baseLevel2,condition2_sizeFactor,condition2_disper)
	message ("6. Estimate the difference between two log2 ratios")
	result_log2R = xTest(condition1_object,condition2_object,threads,bins,baseLevel,ci)

	## 3. combine the log2FC and log2R results and report

	intersect.genes <- intersect(rownames(result_log2FC), rownames(result_log2R))
	result_log2R <- result_log2R[intersect.genes,]
	result_log2FC <- result_log2FC[intersect.genes,]

	#result data frame
	condition1_TE <- paste0(baseLevel,"_log2TE")
	condition2_TE <- paste0(unique(condition)[2], "_log2TE")
	final_result <- cbind(result_log2FC[,c("log2Ratio1","log2Ratio2","deltaTE","Pval")],result_log2R[,c("log2Ratio1","log2Ratio2","deltaTE","Pval")])
	colnames(final_result) <- c("mRNA_log2FC","RPF_log2FC","log2FC_TE_v1","pvalue_v1",condition1_TE,condition2_TE,"log2FC_TE_v2","pvalue_v2")
	final_result <- as.data.frame(final_result)
	final_result$log2FC_TE_final <- 0
	final_result$pvalue_final <- 0
	final_result$pvalue.adjust <- 0
	if (ci>0){
		resultCI <- NA
	}
	log2FC_determine_num <- 0
	log2R_determine_num <- 0
	for (i in 1:nrow(final_result)){
		if (is.na(final_result[i,"pvalue_v1"]) || is.na(final_result[i,"pvalue_v2"])){
		  final_result$log2FC_TE_final[i] <- NA
		  final_result$pvalue_final[i] <- NA
			if (ci>0) resultCI[i] <- NA
		}else if(final_result[i,"pvalue_v1"] > final_result[i,"pvalue_v2"]){
		  final_result$log2FC_TE_final[i] <- final_result[i,"log2FC_TE_v1"]
			final_result$pvalue_final[i] <- final_result[i,"pvalue_v1"]
			if (ci>0) resultCI[i] <- paste0("[",round(result_log2FC[i,"lowci"],2),",",round(result_log2FC[i,"highci"],2),"]")
			log2FC_determine_num <- log2FC_determine_num + 1
		}else{
		  final_result$log2FC_TE_final[i] <- final_result[i,"log2FC_TE_v2"]
		  final_result$pvalue_final[i] <- final_result[i,"pvalue_v2"]
			if (ci>0) resultCI[i] <- paste0("[",round(result_log2R[i,"lowci"],2),",",round(result_log2R[i,"highci"],2),"]")
			log2R_determine_num <- log2R_determine_num + 1
		}
	}

	final_result$pvalue.adjust = p.adjust(final_result$pvalue_final,method=p.adjust.method)
	if (ci>0){
		CI_string = paste0("CI(",100*ci,"%)")
		final_result[[CI_string]] = resultCI
	}
	message("Number of the log2FC and log2R used in determining the final p-value")
	message(paste0(" log2FC: ", log2FC_determine_num))
	message(paste0(" log2R: ", log2R_determine_num))
	xtail_results <- list(resultsTable = final_result, log2FC_determine_num = log2FC_determine_num,
						log2R_determine_num=log2R_determine_num,condition1=baseLevel,
						condition2=unique(condition)[2] )
	class(xtail_results) <- c("xtailResults","list")
	xtail_results
}

# read table with mormalized read counts
counts <- read.csv(options$normalized_read_counts_csv_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# get sample sheet
sampleSheet <- read.csv(file=options$sample_file_path ,header=TRUE, sep="\t", stringsAsFactors=FALSE)

#create condition vector
constraststring <- gsub("contrasts/", "", options$contrast)
#print(constraststring)
contrastconditions <- unlist(strsplit(constraststring,"-"))
cond1 <- contrastconditions[1]
cond2 <- contrastconditions[2]

# split data frame into RIBO and RNA
RIBO <- counts[, (sampleSheet$method == "RIBO") & ( sampleSheet$condition == cond1 | sampleSheet$condition == cond2)]
print("RIBO---------------\n")
print(head(RIBO))
RNA <- counts[, (sampleSheet$method == "RNA")  & ( sampleSheet$condition == cond1 | sampleSheet$condition == cond2)]
print("RNA----------------\n")
print(head(RNA))


#contrastconditions <- c( "TruDWT", "TruDWT", "TruDDel", "TruDDel")
numberofreplicates <- max(sampleSheet$replicate)
print("numberofreplicates")
print(numberofreplicates)
contrastconditionsvector <- rep(contrastconditions,each=numberofreplicates)
# run xtail analysis
test_results <- xRIBOtail(RNA, RIBO, contrastconditionsvector, normalize = FALSE)
#test_results <- xRIBOtail(RNA, RIBO, contrastconditionsvector)

# turn results into table
test_tab <- resultsTable(test_results, log2FCs = TRUE)

# write results into file
write.csv(test_tab, options$xtail_result_path, quote = F)

#plot results
pdf(file=options$xtail_fcplot_path, paper = "a4r", height = 10, width = 13)
plotFCs(test_results)
dev.off()
pdf(file=options$xtail_rplot_path, paper = "a4r", height = 10, width = 13)
plotRs(test_results)
dev.off()

