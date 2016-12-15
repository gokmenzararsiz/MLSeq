test_makeDESeqData <- function(){
	library(MLSeq)
	data(cervical)

	# a subset of cervical data with first 150 features.
	data <- cervical[c(1:150), ]  

	# defining sample classes.
	class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))   

	n <- ncol(data)  # number of samples
	p <- nrow(data)  # number of features

	# train set
	data.train <- data
	data.train <- as.matrix(data.train + 1)
	classtr <- data.frame(condition = class)

	# train set in S4 class
	data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
	                                       colData = classtr, formula(~ condition))
	
	dataClass <- class(data.trainS4)                                
	checkTrue(dataClass == "DESeqDataSet")
}