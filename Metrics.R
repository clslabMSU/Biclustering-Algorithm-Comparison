# compute the match score of two sets of biclusters
Metrics.MatchScore = function(bicluster.set1, bicluster.set2) {
    score = 0L
    jaccards = list()
    for (bicluster in 1:length(bicluster.set1)) {
        score = score + max(sapply(bicluster.set2, FUN=Metrics.Jaccard, bicluster2=bicluster.set1[[bicluster]]))
    }
    return(score / length(bicluster.set1))
}


Metrics.Jaccard = function(bicluster1, bicluster2) {
    n = length(intersect(bicluster1$genes, bicluster2$genes))
    u = length(union(bicluster1$genes, bicluster2$genes))
    return (n/u)
}