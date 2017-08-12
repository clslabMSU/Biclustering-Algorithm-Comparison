# PERFORM ISA BICLUSTERING ALGORITHM
ISA.cluster = function(input.matrix, ...) {
    isa.result = isa(data.matrix(input.matrix[-1, -1]))
    bicluster.rows = apply(isa.result$rows, MARGIN=2, function(b) which(b != 0, arr.ind=TRUE))
    bicluster.cols = apply(isa.result$columns, MARGIN=2, function(b) which(b != 0, arr.ind=TRUE))
    biclusters = list()
    for (i in 1:length(bicluster.rows)) {
        biclusters = c(biclusters, list(list(genes=bicluster.rows[[i]], samples=bicluster.cols[[i]])))
    }
    return(biclusters)
}


ISA.plot = function(input.data, isa.result, ...) {
    if (interactive()) {
        image(data.matrix(input.data[-1, -1]), main="Input Data")
        for (bicluster.index in 1:dim(isa.result$rows)[[2]]) {
            png(filename=paste("Module ", bicluster.index, ".png"))
            image(
                outer(
                    isa.result$rows[,bicluster.index], isa.result$columns[,bicluster.index]
                ), main=paste("Module ", bicluster.index)
            )
            dev.off()
        }
    }
}