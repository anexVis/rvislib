\name{matrix2json}
\alias{matrix2json}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
matrix2json(m, mode = "undirected", weighted = TRUE, zeroBasedIndex = TRUE, node.attributes = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{mode}{
%%     ~~Describe \code{mode} here~~
}
  \item{weighted}{
%%     ~~Describe \code{weighted} here~~
}
  \item{zeroBasedIndex}{
%%     ~~Describe \code{zeroBasedIndex} here~~
}
  \item{node.attributes}{
%%     ~~Describe \code{node.attributes} here~~
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
function (m, mode = "undirected", weighted = TRUE, zeroBasedIndex = TRUE, 
    node.attributes = NA) 
{
    g = graph.adjacency(m, mode = mode, weighted = weighted)
    for (attr in attributes) {
        vertex.attributes(g, attr.name) = attr.value
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
