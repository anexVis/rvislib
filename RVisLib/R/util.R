#' gettree
#'
#' @param merge_data 'merge' object resulted from hclust function
#' @param root_index the largest index of the node in hclust tree
#' @export
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
#' @export
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


#' graph2json
#'
#' Convert an 'igraph' object to json-compatible object
#' @export
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
#' @export
#' @import igraph
matrix2json <- function(m, mode='undirected', weighted=TRUE, zeroBasedIndex=TRUE, node.attributes= NA ) {
    library(igraph)
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


#' Generate heatmap from an adjacency matrix (symmetric)
#'
#' Return 3 json files for visualizing a heatmap: matrix, row dendrogram and column dendrogram
#' @importFrom jsonlite toJSON
#' @export
heatmap.adjacency <- function(m) {
    mat = matrix2json(m, mode='undirected', weighted=TRUE,zeroBasedIndex=zeroBasedIndex)
    rdist = dist(t(m))
    cdist = dist(m)
    rowdend = gettree.hclust(hclust(dist(m)))
    coldend = gettree.hclust(hclust(dist(t(m)))) 
    return(c(toJSON(mat,auto_unbox=T),
             toJSON(rowdend, auto_unbox=T),
             toJSON(coldend, auto_unbox=T)))
}
#' Generate heatmap from a generic mxn matrix
#'
#' @importFrom jsonlite toJSON
#' @import igraph
#' @export
heatmap.generic <- function(m, rowNodeType='row', colNodeType='col') {
    cdend = gettree.hclust(hclust(dist(m)))
    rdend = gettree.hclust(hclust(dist(t(m))))
    # Extend the matrix to make it adjacency-like
    esize = sum(dim(m))
    extended_mat = matrix(rep(NA, esize^2), nrow=esize, ncol=esize)
    extended_mat[1:dim(m)[1], 1:dim(m)[2]] = m[,]
    rownames(extended_mat) = c(rownames(m), colnames(m))
    colnames(extended_mat) = c(colnames(m), rownames(m))    
    gr = graph.adjacency(extended_mat, mode='directed', weighted=T)
    # vertex_attr(gr,'nodeType') = c(rep("tissue", times=dim(m)[2]), rep("gene", times=dim(m)[1]))
    vertex_attr(gr,'nodeType') = c(rep(colNodeType, times=dim(m)[2]), rep(rowNodeType, times=dim(m)[1]))
    gr = delete_edges(gr, E(gr)[is.na(weight)])
    return(c(toJSON(graph2json(gr,attributes=c('name','nodeType'))),
            toJSON(rdend,auto_unbox=TRUE),
            toJSON(cdend,auto_unbox=TRUE)))
}

#' Generate heatmap from uploaded matrix
#'
#' @export
heatmap.files <- function(matrixFile) {
    m = read.csv(matrixFile,sep=' ', comment.char="#")
    return(heatmap.generic(m) )
}
