#' Machine learning interface for RNA-Seq data
#'
#' This package applies several machine learning methods, including SVM, bagSVM, Random Forest and CART, to RNA-Seq data.
#'
#' @references
#'
#' @author Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk
#'
#' -----------------
#'
#' Maintainers:
#'
#' Gokmen Zararsiz, \email{gokmenzararsiz@erciyes.edu.tr}
#'
#' Dincer Goksuluk \email{dincer.goksuluk@hacettepe.edu.tr}
#'
#' Selcuk Korkmaz \email{selcukorkmaz@hotmail.com}
#'
#' @docType package
#' @name MLSeq-package
#' @aliases MLSeq-package
#' @keywords package
NULL


#' Cervical cancer data
#'
#' Cervical cancer data measures the expressions of 714 miRNAs of human samples. There are 29 tumor and 29 non-tumor cervical
#' samples and these two groups are treated as two separete classes.
#'
#' @format A data frame with 58 observations on the following 715 variables.
#'
#' @source \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2880020/#supplementary-material-sec}
#'
#' @references Witten, D., et al. (2010) Ultra-high throughput sequencing-based small RNA discovery and discrete statistical biomarker
#' analysis in a collection of cervical tumours and matched controls. BMC Biology, 8:58
#'
#' @author Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk
#'
#' @docType data
#' @name cervical
#' @aliases cervical
#' @keywords cervical data
#'
#' @examples
#' data(cervical)
#'
NULL

#' Fitting classification models to sequencing data
#'
#' This function fits classification algorithms to sequencing data and measures model performances using various statistics
#'
#' In RNA-Seq studies, normalization is used to adjust between-sample differences for further analysis. In this package,
#' "deseq" and "tmm" normalization methods are available. "deseq" estimates the size factors by dividing each sample by
#' the geometric means of the transcript counts. "tmm" trims the lower and upper side of the data by log fold changes to
#' minimize the log-fold changes between the samples and by absolute intensity. After normalization, it is useful to
#' transform the data for classification. MLSeq package has "voomCPM" and "vst" transformation methods. "voomCPM"
#' transformation applies a logarithmic transformation (log-cpm) to normalized count data. Second transformation method is
#' the "vst" transformation and this approach uses an error modeling and the concept of variance stabilizing transformations
#' to estimate the mean-dispersion relationship of data.
#'
#' For model validation, k-fold cross-validation ("cv" option in MLSeq package) is a widely used technique. Using this technique,
#' training data is randomly splitted into k non-overlapping and equally sized subsets. A classification model is trained on (k-1)
#' subsets and tested in the remaining subsets. MLSeq package also has the repeat option as "rpt" to obtain more generalizable models.
#' Giving a number of m repeats, cross validation concept is applied m times.
#'
#' @param data a \code{DESeqDataSet} object, see the constructor functions
#' \code{\link{DESeqDataSet}},
#' \code{\link{DESeqDataSetFromMatrix}},
#' \code{\link{DESeqDataSetFromHTSeqCount}} in DESeq2 package.
#' @param method a character string indicating the name of classification method. There are four methods available to perform classification:
#' \itemize{
#' \code{svm}: support vector machines using radial-based kernel function
#' \code{bagsvm}: support vector machines with bagging ensemble
#' \code{randomForest}: random forest algorithm
#' \code{cart}: classification and regression trees algorithm
#' }
#' @param normalize a character string indicating the name of normalization method for count data.
#' Available options are:
#' \itemize{
#' \item \code{none}: Normalization is not applied. Count data is used for classification.
#' \item \code{deseq}: deseq normalization.
#' \item \code{tmm}: Trimmed mean of \code{M} values.
#' }
#' @param transformation a character string indicating the normalization method. Note that transformation method is applied after normalization.
#' Available options are \code{vst}: variance stabilizing transformation and \code{voomCPM}: voom transformation (log of counts-per-million).
#' @param control a list including all the control parameters passed to model training process. This arguement is a wrapper for the
#' arguement \code{trControl} from caret package. See \bold{?trainControl} for details.
#' @param B an integer. It is the number of bootstrap samples for method \code{bagsvm}. Default is 100.
#' @param ref a character string indicating the user defined reference class. Default is \code{NULL}. If NULL is selected,
#' first category of class labels is used as reference.
#' @param \dots optional arguments for \code{train(...)} function from \code{caret} package.
#'
#' @return an \code{\link{MLSeq}} object for trained model.
#'
#' @author Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk
#'
#' @references
#'
#' Kuhn M. (2008). Building predictive models in R using the caret package. Journal of Statistical Software, (\url{http://www.jstatsoft.org/v28/i05/})
#'
#' Anders S. Huber W. (2010). Differential expression analysis for sequence count data. Genome Biology, 11:R106
#'
#' Witten DM. (2011). Classification and clustering of sequencing data using a poisson model. The Annals of Applied Statistics, 5(4), 2493:2518
#'
#' Charity WL. et al. (2014) Voom: precision weights unlock linear model analysis tools for RNA-Seq read counts, Genome Biology, 15:R29, doi:10.1186/gb-2014-15-2-r29
#'
#' Witten D. et al. (2010) Ultra-high throughput sequencing-based small RNA discovery and discrete statistical biomarker analysis in a collection of cervical tumours and matched controls. BMC Biology, 8:58
#'
#' Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-Seq data. Genome Biology, 11:R25, doi:10.1186/gb-2010-11-3-r25
#'
#' @keywords RNA-seq classification
#'
#' @import BiocGenerics BiocParallel S4Vectors IRanges GenomicRanges SummarizedExperiment Biobase Rcpp methods
#' @importFrom caret train confusionMatrix bagControl predict.train trainControl
#' @importFrom stats xtabs model.matrix predict relevel
#' @exportClass MLSeq
#'
#' @useDynLib DESeq2
#'
#' @seealso \code{\link{predictClassify}}, \code{\link[caret]{train}}, \code{\link[caret]{trainControl}}
#'
#' @examples
#' data(cervical)
#'
#' # a subset of cervical data with first 150 features.
#' data <- cervical[c(1:150), ]
#'
#' # defining sample classes.
#' class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))
#'
#' n <- ncol(data)  # number of samples
#' p <- nrow(data)  # number of features
#'
#' # number of samples for test set (20% test, 80% train).
#' nTest <- ceiling(n*0.2)
#' ind <- sample(n, nTest, FALSE)
#'
#' # train set
#' data.train <- data[ ,-ind]
#' data.train <- as.matrix(data.train + 1)
#' classtr <- data.frame(condition=class[-ind, ])
#'
#' # train set in S4 class
#' data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
#'                                        colData = classtr, formula(~ condition))
#' data.trainS4 <- DESeq(data.trainS4, fitType = "local")
#'
#' ## Number of repeats (repeats) might change model accuracies ##
#' # Classification and Regression Tree (CART) Classification
#' cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq",
#'                  transformation = "vst", ref = "T",
#'                  control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE))
#' cart
#'
#' # Random Forest (RF) Classification
#' # rf <- classify(data = data.trainS4, method = "randomforest", normalize = "deseq",
#' #  transformation = "vst", ref = "T",
#' #  control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE))
#' # rf
#'
#' @export
classify <- function (data, method = c("svm", "bagsvm", "randomforest", "cart"),
                      normalize = c("deseq", "none", "tmm"), transformation = c("vst", "voomCPM"),
                      control = trainControl(method = "repeatedcv", number = 5, repeats = 10),
                      B = 100, ref = NULL, ...){
  if (!is.null(ref)) {
    if (!is.character(ref))
      stop("Reference class should be \"character\"")
  }
  if (is.null(ref)) {
    ref = levels(data$condition)[1]
  }
  if (class(data)[1] != "DESeqDataSet") {
    stop("Data should be a \"DESeqDataSet Object\" of S4 class.")
  }
  if (is.null(method)) {
    stop("Classification method is not specified.")
  }
  method = match.arg(method)
  normalize = match.arg(normalize)
  transformation = match.arg(transformation)
  if ((normalize == "tmm" & transformation == "vst")) {
    transformation = "voomCPM"
    warning("\"vst\" transformation can be applied with \"deseq\" normalization. \"voom-CPM\" transformation is used.")
  }
  if (normalize == "none") {
    transformation = "NULL"
    warning("\"transformation\" method is ignored since normalization is not applied.")
  }
  conditions = as.factor(data$condition)
  conditions = relevel(conditions, which(levels(conditions) == ref))
  counts = counts(data)
  org.classes = conditions
  org.class.levels = levels(org.classes)
  ctrl <- control

  if (normalize == "none") {
    counts = t(counts)
    conditions = conditions
    dataexp = data.frame(counts, conditions)
  }
  if (normalize == "tmm") {
    counts = counts(data)
    y <- DGEList(counts = counts, genes = rownames(counts))
    y <- calcNormFactors(y)
    design <- model.matrix(~conditions)
    v <- voom(y, design, plot = FALSE)
    dataexp = data.frame(t(v$E), conditions)
    counts = dataexp[, -length(dataexp)]
    conditions = as.factor(dataexp[, length(dataexp)])
  }
  if (normalize == "deseq") {
    data = estimateSizeFactors(data)
    if (transformation == "vst") {
      data = estimateDispersions(data, fitType = "local")
      datavst = getVarianceStabilizedData(data)
      datavst <- as.matrix(datavst)
      cond = data.frame(conditions, row.names = colnames(datavst))
      datavst = ExpressionSet(datavst, AnnotatedDataFrame(cond))
      dataexp = data.frame(t(exprs(datavst)), conditions)
      counts = dataexp[, -length(dataexp)]
      conditions = as.factor(dataexp[, length(dataexp)])
    }
    if (transformation == "voomCPM") {
      counts = counts(data, normalized = TRUE)
      y <- DGEList(counts = counts, genes = rownames(counts))
      design <- model.matrix(~conditions)
      v <- voom(y, design, plot = FALSE)
      dataexp = data.frame(t(v$E), conditions)
      counts = dataexp[, -length(dataexp)]
      conditions = as.factor(dataexp[, length(dataexp)])
    }
  }
  if (method == "svm") {
    train <- train(counts, conditions, method = "svmRadial",
                   trControl = ctrl, ...)
  }
  if (method == "bagsvm") {
    train <- train(counts, conditions, method = "bag", B = B,
                   bagControl = bagControl(fit = svmBag$fit, predict = svmBag$pred,
                                           aggregate = svmBag$aggregate), trControl = ctrl,
                   ...)
  }
  if (method == "randomforest") {
    train <- train(counts, conditions, method = "rf", trControl = ctrl,
                   ...)
  }
  if (method == "cart") {
    train <- train(counts, conditions, method = "rpart",
                   trControl = ctrl, ...)
  }
  train.pred = predict(train)
  tbl.trn = table(Predicted = train.pred, Actual = org.classes)
  confM = confusionMatrix(tbl.trn, reference = ref)
  result = new("MLSeq", confusionMat = confM, trainedModel = train,
               method = method, normalization = normalize, transformation = transformation,
               ref = ref)
  result
}
