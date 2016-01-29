
#' gettree
#'
#' @param merge_data 'merge' object resulted from hclust function
#' @param root_index the largest index of the node in hclust tree
gettree <- function(merge_data, root_index ) {
    children = merge_data[root_index,]
    subtree = list("name"= root_index)
#     names(subtree$children) = 1:length(children)
    for (i in 1:length(children)) {
        child = children[i]
        if (child < 0) {
            subtree$children[[i]] = list(name=as.character(child))
        } else {
            subtree$children[[i]] = gettree(merge_data, child)
        }
    }
    return(subtree)
}

#' gettree.hclust
#' 
#' This function takes a hclust result and returns a tree structure
#' @param hc result of hclust function
#' @param root_index the index of the node to start as root, generally not necessary to specify if one wants to get the full tree from hclust
gettree.hclust <- function(hc, root_index=NULL) {
    if (is.null(root_index)) { root_index = dim(hc$merge)[1] }
    children = hc$merge[root_index,]
    subtree = list("name"= root_index, height=hc$height[root_index])
#     names(subtree$children) = 1:length(children)
    for (i in 1:length(children)) {
        child = children[i]
        if (child < 0) {
            subtree$children[[i]] = list(name=hc$labels[-child], height=0)
        } else {
            subtree$children[[i]] = gettree.hclust(hc, child)
        }
    }
    return(subtree)
}

library(igraph)

#' graph2json
#'
#' Convert an 'igraph' object to json-compatible object
graph2json <- function(g, zeroBasedIndex=TRUE, attributes= c('name') ) {
    gjso = list(nodes=c(), links=c())
    
    for (i in 1:length(V(g)) ) {
#         z <- function(an) {return(V(g)[[i]]$an)}
        gjso$nodes[[i]] = c( list(id = if (zeroBasedIndex) (i-1) else i), sapply(attributes, FUN= function(an) {return( vertex_attr(g,an, index=i) )}) )
                               
    }

    for (i in 1:length(E(g)) ) {
        if ( ! is.na( E(g)[[i]]$weight ) ) {
            gjso$links[[i]] = list(source=as.integer(head_of(g, i)) - (if (zeroBasedIndex) 1 else 0), 
                               target=as.integer(tail_of(g,i))  - (if (zeroBasedIndex) 1 else 0), 
                               value=E(g)[[i]]$weight )    
            
        }
        
    }
    return(gjso)
}

#' matrix2json
#' 
#' Convert an adjacency matrix to json-compatible object
matrix2json <- function(m, mode='undirected', weighted=TRUE, zeroBasedIndex=TRUE, node.attributes= NA ) {
    g = graph.adjacency(m, mode =mode, weighted=weighted)
    for (attr in attributes) {
        vertex.attributes(g, attr.name) = attr.value    
    }

    for (i in 1:length(E(g)) ) {
        if (! is.na( E(g)[[i]]$weight) ) {
            gjso$links[[i]] = list(source=as.integer(head_of(g, i)) - (if (zeroBasedIndex) 1 else 0), 
                               target=as.integer(tail_of(g,i))  - (if (zeroBasedIndex) 1 else 0), 
                               value=E(g)[[i]]$weight )    
            
        }
        
    }
    return(gjso)
}


library(jsonlite) 

#' heatmap
#'
#' Return 3 json files for visualizing a heatmap: matrix, row dendrogram and column dendrogram
heatmap.js <- function(m) {
    mat = matrix2json(m, mode='undirected', weighted=TRUE)
    rdist = dist(m)
    cdist = dist(t(m))
    rowdend = gettree.hclust(hclust(dist(m)))
    coldend = gettree.hclust(hclust(dist(t(m)))) 
    return(toJSON(mat,auto_unbox=T),
            toJSON(rowdend, auto_unbox=T),
            toJSON(coldend, auto_unbox=T))
}

