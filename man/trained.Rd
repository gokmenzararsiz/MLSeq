\name{trained-methods}
\docType{methods}
\alias{trained}
\alias{trained,MLSeq-method}
\title{Accessors for the 'trained' slot of an MLSeq object}
\usage{
\S4method{trained}{MLSeq}(object)
}
\arguments{
  \item{object}{an \code{MLSeq} object}
}

\description{
Details about the trained model which is returned by \code{classify} function. This object is a list returned of class \code{train} from CRAN package "caret".}
\details{
\code{trained} slot stores information about the training process such as optimum model parameters and resampling properties on the fitted classification model. For detailed information on returned values, see the help manual of \code{train(...)} function in caret.}

\author{
Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk
}

\seealso{
\code{\link[caret]{train}}, \code{\link[MLSeq]{classify}}
}

\examples{
data(cervical)

# a subset of cervical data with first 150 features.
data <- cervical[c(1:150), ]

# defining sample classes.
class <- data.frame(condition = factor(rep(c("N","T"), c(29, 29))))   

n <- ncol(data)  # number of samples
p <- nrow(data)  # number of features

# number of samples for test set (20% test, 80% train).
nTest <- ceiling(n*0.2)  
ind <- sample(n, nTest, FALSE)

# train set
data.train <- data[ ,-ind]
data.train <- as.matrix(data.train + 1)
classtr <- data.frame(condition=class[-ind, ])

# train set in S4 class
data.trainS4 <- DESeqDataSetFromMatrix(countData = data.train,
                                       colData = classtr, formula(~ condition))
data.trainS4 <- DESeq(data.trainS4, fitType = "local")

## Number of repeats (repeats) might change model accuracies ##
# Classification and Regression Tree (CART) Classification
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", 
	deseqTransform = "vst", ref = "T",
	control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE))

trained(cart)
}
