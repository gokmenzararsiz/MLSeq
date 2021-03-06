% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictClassify.R
\name{predictClassify}
\alias{predictClassify}
\title{Extract predictions from \code{classify()} object}
\usage{
predictClassify(model, test.data, ...)
}
\arguments{
\item{model}{a model of \code{MLSeq} class}

\item{test.data}{a \code{DESeqDataSet} instance of new observations.}

\item{\dots}{further arguments to be passed to or from methods. These arguements are used in
\code{\link[caret]{predict.train}} from caret package.}
}
\value{
\code{MLSeqObject} an MLSeq object returned from \code{classify}. See details.

\code{Predictions} a data frame or vector including either the predicted class
probabilities or class labels of given test data.
}
\description{
This function predicts the class labels of test data for a given model.
}
\details{
\code{predictClassify} function returns the predicted class information along with trained model.
Predicted values are given either as class labels or estimated probabilities of each class for
each sample. If \code{type = "raw"}, as can be seen in the example below, the predictions are
extracted as raw class labels.In order to extract estimated class probabilities, one should follow the steps below:
\itemize{
\item set \code{classProbs = TRUE} within \code{control} arguement in \code{\link{classify}}
\item set \code{type = "prob"} within \code{predictClassify}
}
}
\examples{
library(DESeq2)
data(cervical)

# a subset of cervical data with first 150 features.
data <- cervical[c(1:150), ]

# defining sample classes.
class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))

n <- ncol(data)  # number of samples
p <- nrow(data)  # number of features

# number of samples for test set (20\% test, 80\% train).
nTest <- ceiling(n*0.2)
ind <- sample(n, nTest, FALSE)

# train set
data.train <- data[ ,-ind]
data.train <- as.matrix(data.train + 1)
classtr <- data.frame(condition = class[-ind, ])

# train set in S4 class
data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
                  colData = classtr, formula(~ condition))
data.trainS4 <- DESeq(data.trainS4, fitType = "local")

# test set
data.test <- data[ ,ind]
data.test <- as.matrix(data.test + 1)
classts <- data.frame(condition=class[ind, ])

# test set in S4
data.testS4 <- DESeqDataSetFromMatrix(countData = data.test,
                 colData = classts, formula(~ condition))
data.testS4 <- DESeq(data.testS4, fitType = "local")

## Number of repeats (repeats) might change model accuracies ##
# Classification and Regression Tree (CART) Classification
cart <- classify(data = data.trainS4, method = "cart",
          transformation = "vst", ref = "T", normalize = "deseq",
          control = trainControl(method = "repeatedcv", number = 5,
                                 repeats = 3, classProbs = TRUE))
cart

# Random Forest (RF) Classification
rf <- classify(data = data.trainS4, method = "randomforest",
        transformation = "vst", ref = "T", normalize = "deseq",
        control = trainControl(method = "repeatedcv", number = 5,
                               repeats = 3, classProbs = TRUE))
rf

# predicted classes of test samples for CART method (class probabilities)
pred.cart = predictClassify(cart, data.testS4, type = "prob")
pred.cart

# predicted classes of test samples for RF method (class labels)
pred.rf = predictClassify(rf, data.testS4, type = "raw")
pred.rf

}
\author{
Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk
}
\references{
Kuhn M. (2008). Building predictive models in R using the caret package. Journal of Statistical Software,
(\url{http://www.jstatsoft.org/v28/i05/})

Anders S. Huber W. (2010). Differential expression analysis for sequence count data. Genome Biology, 11:R106

Witten DM. (2011). Classification and clustering of sequencing data using a poisson model. The Annals of Applied Statistics, 5(4), 2493:2518

Charity WL. et al. (2014) Voom: precision weights unlock linear model analysis tools for RNA-Seq read counts,
Genome Biology, 15:R29, doi:10.1186/gb-2014-15-2-r29

Witten D. et al. (2010) Ultra-high throughput sequencing-based small RNA discovery and discrete statistical
biomarker analysis in a collection of cervical tumours and matched controls. BMC Biology, 8:58

Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-Seq data.
Genome Biology, 11:R25, doi:10.1186/gb-2010-11-3-r25
}
\seealso{
\code{\link{classify}}, \code{\link[caret]{train}}, \code{\link[caret]{trainControl}}
}
\keyword{RNA-seq}
\keyword{classification}

