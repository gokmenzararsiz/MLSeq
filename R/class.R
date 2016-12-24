setOldClass(c("confusionMatrix","train"))

#' \code{MLSeq} object
#'
#' For classification, this is the main class for the \code{MLSeq} package.
#'
#' Objects can be created by calls of the form \code{new("MLSeq", ...)}. This type
#' of objects is created as a result of \code{classify} function of \code{MLSeq} package.
#' It is then used in \code{predictClassify} function for predicting the class labels of new samples.
#'
#' @section Slots:
#'
#' \describe{
#'   \item{\code{method}:}{stores the name of used classification method in the classification model}
#'   \item{\code{transformation}:}{stores the name of used transformation method in the classification model}
#'   \item{\code{normalization}:}{stores the name of used normalization method in the classification model}
#'   \item{\code{confusionMat}:}{stores the information of classification performance results}
#'   \item{\code{trainedModel}:}{stores the information about training process and model parameters that used in the corresponding model}
#'   \item{\code{ref}:}{stores user defined reference class}
#' }
#'
#' @note An \code{MLSeq} class stores the results of \code{classify} function and offers further slots that are populated
#' during the analysis. The slot \code{confusionMat} stores the information of classification performance results. These
#' results contain the classification table and several statistical measures including accuracy rate, sensitivity, specifity,
#' positive and negative predictive rates, etc. \code{method}, \code{normalization} and \code{deseqTransform} slots store
#' the name of used classification method, normalization method and transformation method in the classification model respectively.
#' Lastly, the slot \code{trained} stores the information about training process and model parameters that used in the corresponding model.
#'
#' @author Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk
#'
#' @docType class
#' @name MLSeq-class
#' @rdname MLSeq-class
#' @aliases MLSeq-class
#' @exportClass MLSeq
setClass("MLSeq",
			slots = c(method = "character",
					      transformation = "character",
					      normalization = "character",
					      confusionMat = "confusionMatrix",
					      trainedModel = "train",
					      ref = "character"),
			prototype = prototype(confusionMat = structure(list(), class = "confusionMatrix"),
			                      trainedModel = structure(list(), class = "train")))


setValidity("MLSeq", function( object ) {
    # if (!(method(object)  %in% c("svm", "bagsvm", "randomforest", "cart")))
    # return("Error: 'method' slot must be in one of the following methods: \"svm\", \"bagsvm\", \"randomforest\", \"cart\" ")

    if (!(normalization(object)  %in% c("deseq", "none", "tmm")))
    return("Error: 'normalization' slot must be in one of the following: \"deseq\", \"none\", \"tmm\" ")

    if (!(transformation(object)  %in% c("vst", "voomCPM", "NULL")))
    return("Error: 'transformation' slot must be in one of the following: \"vst\", \"voomCPM\" ")

    if (!is.character(ref(object)))
    return("Error: 'ref' slot must be a character ")

    if ((normalization(object) == "tmm" & transformation(object) == "vst"))
    return("Warning: \"vst\" transformation can be applied only with \"deseq\" normalization. \"voom-CPM\" transformation is used. ")

    TRUE
} )





