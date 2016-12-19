#' Accessors for the 'method' slot of an \code{MLSeq} object
#'
#' This slot stores the name of selected model which is used in \code{classify} function. The trained model is stored in slot \code{trainedModel}.
#' See \code{\link{trained}} for details.
#'
#' @docType methods
#' @name method-methods
#' @rdname method-methods
#' @aliases method-methods
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @seealso \code{\link{trained}}
#'
#' @examples
#'
#' data(cervical)
#'
#' data <- cervical[c(1:150),]  # a subset of cervical data with first 150 features.
#'
#' class <- data.frame(condition=factor(rep(c("N","T"), c(29,29))))# defining sample classes.
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' nTest <- ceiling(n*0.2)  # number of samples for test set (20% test, 80% train).
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                                        colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq",
#'                  transformation = "vst", ref = "T",
#'                  control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE))
#'
#' method(cart)
#'
#' @export
setMethod("method", "MLSeq", function(object) object@method)



#' Accessors for the 'transformation' slot of an \code{MLSeq} object
#'
#' This slot stores the name of transformation method which is used while transforming the count data (i.e. either "vst" or "voomCPM")
#'
#' @docType methods
#' @name transformation-methods
#' @rdname transformation-methods
#' @aliases transformation-methods
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @examples
#'
#' data(cervical)
#'
#' data <- cervical[c(1:150),]  # a subset of cervical data with first 150 features.
#'
#' class <- data.frame(condition=factor(rep(c("N","T"), c(29,29))))# defining sample classes.
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' nTest <- ceiling(n*0.2)  # number of samples for test set (20% test, 80% train).
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                                        colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq",
#'                  transformation = "vst", ref = "T",
#'                  control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE))
#'
#' transformation(cart)
#'
#' @export
setMethod("transformation", "MLSeq", function(object) object@transformation)


#' Accessors for the 'normalization' slot of an \code{MLSeq} object
#'
#' This slot stores the name of normalization method which is used while normalizing the count data such as "deseq", "none" or "tmm"
#'
#' @docType methods
#' @name normalization-methods
#' @rdname normalization-methods
#' @aliases normalization-methods
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @examples
#'
#' data(cervical)
#'
#' data <- cervical[c(1:150),]  # a subset of cervical data with first 150 features.
#'
#' # defining sample classes.
#' class <- data.frame(condition=factor(rep(c("N", "T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' nTest <- ceiling(n*0.2)  # number of samples for test set (20% test, 80% train).
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                                        colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq",
#'                  transformation = "vst", ref = "T",
#'                  control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE))
#'
#' normalization(cart)
#'
#' @export
setMethod("normalization", "MLSeq", function(object) object@normalization)



#' Accessors for the 'confusionMat' slot of an \code{MLSeq} object
#'
#' This slot stores the confusion matrix for the trained model using \code{classify} function.
#'
#' \code{confusionMat} slot includes information about cross-tabulation of observed and predicted classes
#' and corresponding statistics such as accuracy rate, sensitivity, specifity, etc. The returned object is in \code{confusionMatrix} class of
#' caret package. See \code{\link[caret]{confusionMatrix}} for details.
#'
#' @docType methods
#' @name confusionMat-methods
#' @rdname confusionMat-methods
#' @aliases confusionMat-methods
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @seealso \code{\link[caret]{confusionMatrix}}
#'
#' @examples
#'
#' data(cervical)
#'
#' data <- cervical[c(1:150),]  # a subset of cervical data with first 150 features.
#'
#' # defining sample classes.
#' class <- data.frame(condition=factor(rep(c("N", "T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' nTest <- ceiling(n*0.2)  # number of samples for test set (20% test, 80% train).
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                                        colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq",
#'                  transformation = "vst", ref = "T",
#'                  control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE))
#'
#' confusionMat(cart)
#'
#' @export
setMethod("confusionMat", signature = "MLSeq", function(object) object@confusionMat)


#' Accessors for the 'trainedModel' slot of an \code{MLSeq} object
#'
#' This slot stores the trained model. This object is returned from \code{train.default} function in caret package. Any further request
#' using caret functions is available for \code{trainedModel} since this object is in the same class as the returned object from
#' \code{train.default}. See \code{\link[caret]{train.default}} for details.
#'
#' @docType methods
#' @name trained-methods
#' @rdname trained-methods
#' @aliases trained-methods
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @seealso \code{\link[caret]{train.default}}
#'
#' @examples
#'
#' data(cervical)
#'
#' data <- cervical[c(1:150),]  # a subset of cervical data with first 150 features.
#'
#' # defining sample classes.
#' class <- data.frame(condition=factor(rep(c("N", "T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' nTest <- ceiling(n*0.2)  # number of samples for test set (20% test, 80% train).
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                                        colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq",
#'                  transformation = "vst", ref = "T",
#'                  control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE))
#'
#' trained(cart)
#'
#' @export
setMethod("trained", signature = "MLSeq", function(object) object@trainedModel)



#' Accessors for the 'ref' slot of an \code{MLSeq} object
#'
#' This slot stores the information about reference category. Confusion matrix and related statistics are calculated using
#' the user-defined reference category.
#'
#' @docType methods
#' @name ref-methods
#' @rdname ref-methods
#' @aliases ref-methods
#'
#' @param object an \code{MLSeq} object.
#'
#' @author Gokmen Zararsiz
#'
#' @examples
#'
#' data(cervical)
#'
#' data <- cervical[c(1:150),]  # a subset of cervical data with first 150 features.
#'
#' # defining sample classes.
#' class <- data.frame(condition=factor(rep(c("N", "T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' nTest <- ceiling(n*0.2)  # number of samples for test set (20% test, 80% train).
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                                        colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' # Classification and Regression Trees (CART)
#' cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq",
#'                  transformation = "vst", ref = "T",
#'                  control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE))
#'
#' ref(cart)
#'
#' @export
setMethod("ref", "MLSeq", function(object) object@ref)

