library(futile.logger)
flog.appender(appender.file("test_util.log"), name="log")

context("Utility functions")
data(iris)
somedata = iris[1:50,1:4]
corrMatrix = cor(somedata)

test_that("Generate json data for heatmap plotting. Please verify output in log file.", {
    # confirm that corrMatrix has both rownames and colnames
    rownames(corrMatrix) = paste("V", c(1:4), sep="")
    hm_data = heatmap.adjacency(corrMatrix)
    flog.info("Input adjacency matrix", corrMatrix,name="log", capture=TRUE)
    expect_equal(length(hm_data),2)
    
    flog.info("Matrix data",
              jsonlite::prettify(hm_data[1]),
              name="log", 
              capture=TRUE)
})

# This function is from igraph, and should work. Test is just to demonstrate the proper use
test_that("graph.adjacency incorporates row and column names in vertex attribute", {
    rownames(corrMatrix) = paste("V", c(1:4), sep="")
    gr = graph.adjacency(corrMatrix, mode="directed",add.colnames = 'name', add.rownames = 'label')
    expect_equal(vertex_attr(gr,'name'), colnames(corrMatrix))
    expect_equal(vertex_attr(gr,'label'), rownames(corrMatrix))
})