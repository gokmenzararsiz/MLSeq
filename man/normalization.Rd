\name{normalization-methods}
\docType{methods}
\alias{normalization}
\alias{normalization,MLSeq-method}
\title{Accessors for the 'normalization' slot of an MLSeq object}

\usage{
\S4method{normalization}{MLSeq}(object)
}

\arguments{
  \item{object}{an \code{MLSeq} object}
}

\description{
Used normalization method for the trained model using \code{classify} function.  
}

\details{
\code{normalization} slot stores the name of the normalization method "deseq", "none" or "tmm"		
}

\author{
Gokmen Zararsiz, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Izzet Parug Duru, Turgay Unver, Ahmet Ozturk
}


\examples{
data(cervical)

# a subset of cervical data with first 150 features.
data <- cervical[c(1:150),]  

# defining sample classes.
class <- data.frame(condition = factor(rep(c("N", "T"),c(29, 29))))

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

# Classification and Regression Trees (CART)
cart <- classify(data = data.trainS4, method = "cart", normalize = "deseq", 
	deseqTransform = "vst", ref = "T",
	control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE))

normalization(cart)	
}