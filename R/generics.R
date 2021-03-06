# setGeneric("MLSeq", function(object) standardGeneric("MLSeq"))

#' Accessors for the 'method' slot of an \code{MLSeq} object
#'
#' This slot stores the name of selected model which is used in \code{classify} function.
#' The trained model is stored in slot \code{trainedModel}.
#' See \code{\link{trained}} for details.
#'
#' \code{method} slot stores the name of the classification method as "svm", support vector machines using
#' radial-based kernel function; "bagsvm", support vector machines with bagging ensemble; "randomForest",
#' random forest algorithm and "cart", classification and regression trees algorithm.
#'
#' @docType methods
#' @name method-methods
#' @rdname method
#' @aliases method
#' @include class.R
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @seealso \code{\link{trained}}
#'
#' @examples
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150),]
#'
#' # defining sample classes.
#' class <- data.frame(condition=factor(rep(c("N","T"), c(29,29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (20% test, 80% train).
#' nTest <- ceiling(n*0.2)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                   colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart",
#'           transformation = "vst", ref = "T", normalize = "deseq",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' method(cart)
#'
#' @export
setGeneric("method", function(object) standardGeneric("method"))


#' Accessors for the 'transformation' slot of an \code{MLSeq} object
#'
#' This slot stores the name of transformation method which is used while transforming the
#' count data (i.e. either "vst" or "voomCPM")
#'
#' @docType methods
#' @name transformation-methods
#' @rdname transformation
#' @aliases transformation
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @examples
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150),]
#'
#' # defining sample classes.
#' class <- data.frame(condition=factor(rep(c("N","T"), c(29,29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (20% test, 80% train).
#' nTest <- ceiling(n*0.2)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                   colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart",
#'           transformation = "vst", ref = "T", normalize = "deseq",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' transformation(cart)
#'
#' @export
setGeneric("transformation", function(object) standardGeneric("transformation"))


#' Accessors for the 'normalization' slot of an \code{MLSeq} object
#'
#' This slot stores the name of normalization method which is used while normalizing the count data such
#' as "deseq", "none" or "tmm"
#'
#' @docType methods
#' @name normalization-methods
#' @rdname normalization
#' @aliases normalization
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @examples
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150),]
#'
#' # defining sample classes.
#' class <- data.frame(condition=factor(rep(c("N", "T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (20% test, 80% train).
#' nTest <- ceiling(n*0.2)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                   colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart",
#'           transformation = "vst", ref = "T", normalize = "deseq",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' normalization(cart)
#'
#' @export
setGeneric("normalization", function(object) standardGeneric("normalization"))



#' Accessors for the 'confusionMat' slot of an \code{MLSeq} object
#'
#' This slot stores the confusion matrix for the trained model using \code{classify} function.
#'
#' \code{confusionMat} slot includes information about cross-tabulation of observed and predicted classes
#' and corresponding statistics such as accuracy rate, sensitivity, specifity, etc. The returned object
#' is in \code{confusionMatrix} class of caret package. See \code{\link[caret]{confusionMatrix}} for details.
#'
#' @docType methods
#' @name confusionMat-methods
#' @rdname confusionMat
#' @aliases confusionMat
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @seealso \code{\link[caret]{confusionMatrix}}
#'
#' @examples
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150),]
#'
#' # defining sample classes.
#' class <- data.frame(condition=factor(rep(c("N", "T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (20% test, 80% train).
#' nTest <- ceiling(n*0.2)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                   colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart",
#'           transformation = "vst", ref = "T", normalize = "deseq",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' confusionMat(cart)
#'
#' @export
setGeneric("confusionMat", function(object) standardGeneric("confusionMat"))



#' Accessors for the 'trainedModel' slot of an \code{MLSeq} object
#'
#' This slot stores the trained model. This object is returned from \code{train.default} function in caret package.
#' Any further request using caret functions is available for \code{trainedModel} since this object is in the
#' same class as the returned object from \code{train.default}. See \code{\link[caret]{train.default}} for details.
#'
#' @docType methods
#' @name trained-methods
#' @rdname trained
#' @aliases trained
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @seealso \code{\link[caret]{train.default}}
#'
#' @examples
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150),]
#'
#' # defining sample classes.
#' class <- data.frame(condition=factor(rep(c("N", "T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (20% test, 80% train).
#' nTest <- ceiling(n*0.2)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                   colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart",
#'           transformation = "vst", ref = "T", normalize = "deseq",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' trained(cart)
#'
#' @export
setGeneric("trained", function(object) standardGeneric("trained"))



#' Accessors for the 'ref' slot of an \code{MLSeq} object
#'
#' This slot stores the information about reference category. Confusion matrix and related statistics are calculated using
#' the user-defined reference category.
#'
#' @docType methods
#' @name ref-methods
#' @rdname ref
#' @aliases ref
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @examples
#' library(DESeq2)
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150),]
#'
#' # defining sample classes.
#' class <- data.frame(condition=factor(rep(c("N", "T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (20% test, 80% train).
#' nTest <- ceiling(n*0.2)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                   colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart",
#'           transformation = "vst", ref = "T", normalize = "deseq",
#'           control = trainControl(method = "repeatedcv", number = 5,
#'                                  repeats = 3, classProbs = TRUE))
#'
#' ref(cart)
#'
#' @export
setGeneric("ref", function(object) standardGeneric("ref"))


#' Show method for MLSeq objects
#'
#' Prints out the information from the trained model using \code{classify} function.
#'
#' @docType methods
#' @name show-methods
#' @rdname show
#' @aliases show
#'
#' @author Gokmen Zararsiz
#'
#' @param object an \code{MLSeq} object returned from \code{classify} function.
#'
#' @seealso \code{\link{classify}}
#'
NULL
# setGeneric("show", function(object) standardGeneric("show"))








