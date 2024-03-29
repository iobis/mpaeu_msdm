% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdm_addons.R
\name{eval_metrics}
\alias{eval_metrics}
\title{Evaluate models with common metrics}
\usage{
eval_metrics(
  original,
  predicted,
  metrics = c("auc", "cbi", "pr", "prg", "tss", "spec", "sens", "kap", "fmeas", "opr",
    "upr"),
  thresholds = c("maxsss", "mtp", "p10")
)
}
\arguments{
\item{original}{a vector with the original values}

\item{predicted}{a vector with predicted values}

\item{metrics}{a character vector with one or more metrics (see details)}

\item{thresholds}{a character vector with one or more thresholds to be used
for those metrics that need thresholding (TSS, Specificity, Sensitivity and
Kappa; see details). Accepted values are "p10", "maxsss" and "mtp".}
}
\value{
a data.frame with the chosen metrics
}
\description{
Evaluate models with common metrics
}
\details{
Eight metrics are available:
\itemize{
\item \code{auc}: Area Under the Curve
\item \code{cbi}: Continuous Boyce Index
\item \code{pr}: Area under the Precision-Recall curve
\item \code{prg}: Precision-Recall-Gain curve
\item \code{tss}: True Skill Statistics
\item \code{spec}: Specificity (True Negative Rate)
\item \code{sens}: Sensitivity (True Positive Rate)
\item \code{kap}: Kappa statistics
\item \code{fmeas}: F-measure (Sorensen's similarity index)
\item \code{opr}: Overprediction rate
\item \code{upr}: (underprediction rate)
}

For the thresholding methods, three options are available:
\itemize{
\item \code{p10}: The value at which 90\% of presence points are included
\item \code{mtp}: The minimum training presence, i.e., the value at which all
training presences are included.
\item \code{maxsss}: The value that maximizes the sum of specificity plus sensitivity.
}
}
\examples{
\dontrun{

}
}
