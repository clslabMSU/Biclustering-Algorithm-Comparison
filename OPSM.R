OPSM.cluster = function(input.matrix, opsm.path, ...) {
    A = as.matrix(input.matrix[-1, -1])
    colnames(A) = NULL
    write.matrix(A, "../bin/opsm_tmp.txt", sep="\t")
    command = paste("java -jar", opsm.path, "../bin/opsm_tmp.txt", dim(A)[[1]], dim(A)[[2]], "../bin/opsm_tmp.out", 10)
    system(command)
    result = OPSM.ParseOutput("../bin/opsm_tmp.out")
    file.remove("../bin/opsm_tmp.txt")
    file.remove("../bin/opsm_tmp.out")
    return(result)
}

OPSM.ParseOutput = function(filename) {
    opsm.output = readLines(filename)
    opsm.output = opsm.output[which(opsm.output != "")]
    biclusters = list()
    for (i in seq(from=1, to=length(opsm.output), by=2)) {
        genes = as.numeric(unlist(strsplit(opsm.output[[i]], " ")))
        samples = as.numeric(unlist(strsplit(opsm.output[[i+1]], " ")))
        biclusters = c(biclusters, list(list(genes=genes, samples=samples)))
    }
    return(biclusters)
}