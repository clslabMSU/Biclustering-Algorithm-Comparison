FABIA.parseBiclusters = function(FABIA.biclusters) {
    biclusters = list()
    for (i in 1:dim(FABIA.result$bic)[[1]]) {
        genes = as.numeric(unlist(lapply(FABIA.biclusters$bic[i,]$bixn, FUN=function(x) return(substring(x, 2)))))
        samples = as.numeric(unlist(lapply(FABIA.biclusters$bic[i,]$biypn, FUN=function(x) return(substring(x, 2)))))
        biclusters = c(biclusters, list(list(genes=genes, samples=samples)))
    }
    return(biclusters)
}