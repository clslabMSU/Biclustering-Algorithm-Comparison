Misc.ConvertBiclust = function(biclust.object) {
    biclusters = list()
    for (i in 1:dim(biclust.object@RowxNumber)[[2]]) {
        biclusters = c(biclusters, list(list(
            genes=which(biclust.object@RowxNumber[,i] == TRUE),
            samples=which(biclust.object@NumberxCol[i,] == TRUE)
        )))
    }
    return(biclusters)
}


Misc.GatherResults = function(result, outfile, genuine_biclusters) {
    recovery = Metrics.MatchScore(genuine_biclusters, result)
    relevance = Metrics.MatchScore(result, genuine_biclusters)
    write(paste("\t\tRecovery:  ", recovery, "\n\t\tRelevance: ", relevance), file=outfile, append=TRUE, sep="\n")
    return(list(recovery=recovery, relevance=relevance))
}