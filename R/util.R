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

# #' matrix2json
# #'
# #' Convert an adjacency matrix to json-compatible object
# #' @export
# #' @import igraph
# matrix2json <- function(m, mode='undirected', weighted=TRUE, zeroBasedIndex=TRUE, node.attributes="name") {
#     library(igraph)
#     g = graph.adjacency(m, mode =mode, weighted=weighted)
#     for (attr in node.attributes) {
#         vertex.attributes(g, attr.name) = attr.value
#     }
#     gjso = list(nodes=c(), links=c())
#     for (i in 1:length(E(g)) ) {
#         if (! is.na( E(g)[[i]]$weight) ) {
#             gjso$links[[i]] = list(source=as.integer(head_of(g, i)) - (if (zeroBasedIndex) 1 else 0),
#                                target=as.integer(tail_of(g,i))  - (if (zeroBasedIndex) 1 else 0),
#                                value=E(g)[[i]]$weight )
#
#         }
#
#     }
#     return(gjso)
# }

#' tsv2json
#'
#' Convert a tsv-formatted directed network to json-compatible object
#' @export
#' @import igraph
tsv2json <- function(tsv_file, mode='signed',zeroBasedIndex=TRUE) {
    library(igraph)
    tsv = read.table(tsv_file,header=F,sep='\t',stringsAsFactors = F)
    g = graph_from_edgelist(as.matrix(tsv[,1:2]),directed = TRUE)
    if (mode == 'signed') {
        tsv[tsv[,3] == '+',3] = 1
        tsv[tsv[,3] == '-',3] = -1
        E(g)$weight = as.numeric(tsv[,3])
    } else {
        E(g)$weight = tsv[,3]   # might not work if the tsv file also includes 0-weighted edges
    }
    gjso = list(nodes=c(), links=c())
    for (i in 1:length(V(g))) {
        gjso$nodes[[i]] = list(id = i - (if (zeroBasedIndex) 1 else 0),
                               name=V(g)[[i]]$name)
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

#' tsv2jsoncy
#'
#' Convert a tsv-formatted directed network to json-compatible object that could be imported directly by cytoscape
#' @export
#' @import igraph
tsv2jsoncy <- function(tsv_file, mode='signed',zeroBasedIndex=TRUE) {
    library(igraph)
    tsv = read.table(tsv_file,header=F,sep='\t',stringsAsFactors = F)
    g = graph_from_edgelist(as.matrix(tsv[,1:2]),directed = TRUE)
    if (mode == 'signed') {
        tsv[tsv[,3] == '+',3] = 1
        tsv[tsv[,3] == '-',3] = -1
        E(g)$weight = as.numeric(tsv[,3])
    } else {
        E(g)$weight = tsv[,3]   # might not work if the tsv file also includes 0-weighted edges
    }
    gjso = list();
    nNodes = length(V(g))
    for (i in 1:nNodes) {
        gjso[[i]] = list(group= 'nodes',
                         data = list(
                            id = as.character(i - (if (zeroBasedIndex) 1 else 0)),
                            name=V(g)[[i]]$name))
    }
    for (i in 1:length(E(g)) ) {
        if (! is.na( E(g)[[i]]$weight) ) {
            gjso[[i + nNodes]] = list(group = 'edges',
                                      data = list(
                                        id = paste('e',i,sep=''),
                                        source=as.character(as.integer(head_of(g, i)) - (if (zeroBasedIndex) 1 else 0)),
                                        target=as.character(as.integer(tail_of(g,i))  - (if (zeroBasedIndex) 1 else 0)),
                                        value=E(g)[[i]]$weight))
        }
    }
    return(gjso)
}


#' Generate heatmap from an adjacency matrix (symmetric)
#'
#' Return an array of 2 json objects for visualizing a heatmap: matrix and row dendrogram.
#' Since we return the row dendrogram, rownames in m would be used as leave names in the dendrogram.
#' To keep the attribute consistent, be sure to specify add.rownames='name' in calling this function
#' @importFrom jsonlite toJSON
#' @export
heatmap.adjacency <- function(m, zeroBasedIndex=TRUE,add.colnames='label', add.rownames='name') {
    gr = graph.adjacency(m, mode='directed', weighted=T, add.colnames = add.colnames, add.rownames = add.rownames)
    gr = delete_edges(gr, E(gr)[is.na(weight)])

    rowdend = gettree.hclust(hclust(dist(m)))
    return(c(toJSON(graph2json(gr, attributes= c(add.colnames,add.rownames)),auto_unbox=T),
             toJSON(rowdend, auto_unbox=T)))
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
