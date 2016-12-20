#' @rdname method
#' @export
setMethod("method", signature(object = "MLSeq"), function(object) object@method)


#' @rdname transformation
#' @export
setMethod("transformation", signature(object = "MLSeq"), function(object) object@transformation)


#' @rdname normalization
#' @export
setMethod("normalization", signature(object = "MLSeq"), function(object) object@normalization)


#' @rdname confusionMat
#' @export
setMethod("confusionMat", signature(object = "MLSeq"), function(object) object@confusionMat)


#' @rdname trained
#' @export
setMethod("trained", signature(object = "MLSeq"), function(object) object@trainedModel)


#' @rdname ref
#' @export
setMethod("ref", signature(object = "MLSeq"), function(object) object@ref)


#' @rdname show
#' @export
setMethod("show",
  signature(object = "MLSeq"),
  definition = function(object) {
    if (dim(confusionMat(object)$table)[1] == 2){
      cat("\n", sep = " ")
      cat("  An object of class ", class(object), "\n\n", sep = " ")
      cat("            Method  : ", method(object), "\n\n")
      cat("       Accuracy(%)  : ", round(confusionMat(object)$overall[1],4)*100, "\n")
      cat("    Sensitivity(%)  : ", round(confusionMat(object)$byClass[1],4)*100, "\n")
      cat("    Specificity(%)  : ", round(confusionMat(object)$byClass[2],4)*100, "\n\n")
      cat("  Reference Class   : ", ref(object), "\n\n")
    } else {
      cat("\n", sep = " ")
      cat("  An object of class ", class(object), "\n\n", sep = " ")
      cat("            Method  : ", method(object), "\n\n")
      cat("       Accuracy(%)  : ", round(confusionMat(object)$overall[1],4)*100, "\n")
      cat("    Sensitivity(%)  : ", round(confusionMat(object)$byClass[1,1],4)*100, "\n")
      cat("    Specificity(%)  : ", round(confusionMat(object)$byClass[1,2],4)*100, "\n\n")
      cat("  Reference Class   : ", ref(object), "\n\n")
      invisible(NULL)
    }
  })
