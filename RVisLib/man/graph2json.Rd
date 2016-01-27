\name{graph2json}
\alias{graph2json}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
graph2json(g, zeroBasedIndex = TRUE, attributes = c("name"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{zeroBasedIndex}{
%%     ~~Describe \code{zeroBasedIndex} here~~
}
  \item{attributes}{
%%     ~~Describe \code{attributes} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (g, zeroBasedIndex = TRUE, attributes = c("name")) 
{
    gjso = list(nodes = c(), links = c())
    for (i in 1:length(V(g))) {
        gjso$nodes[[i]] = c(list(id = if (zeroBasedIndex) (i - 
            1) else i), sapply(attributes, FUN = function(an) {
            return(vertex_attr(g, an, index = i))
        }))
    }
    for (i in 1:length(E(g))) {
        if (!is.na(E(g)[[i]]$weight)) {
            gjso$links[[i]] = list(source = as.integer(head_of(g, 
                i)) - (if (zeroBasedIndex) 1 else 0), target = as.integer(tail_of(g, 
                i)) - (if (zeroBasedIndex) 1 else 0), value = E(g)[[i]]$weight)
        }
    }
    return(gjso)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
