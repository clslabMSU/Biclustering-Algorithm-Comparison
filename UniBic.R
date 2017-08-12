UniBic.cluster = function(input.filename) {
    system(paste("./", UniBic_Binary, " -i ", input.filename, sep=""))
    biclusters = UniBic.parseBiclusters()
    return(biclusters)
}


UniBic.parseBiclusters = function() {
    
    lines = readLines("../data/out.blocks")
    biclusters = list()
    for (line in 1:length(lines)) {
        if (!is.na(str_match(lines[[line]], "BC[0-9]+")[1,1])) {
            genes_line = lines[[line+1]]
            samples_line = lines[[line+2]]
            genes = as.numeric(sapply(strsplit(substring(genes_line, str_locate(genes_line, "Genes \\[[0-9]+\\]: ")[1,2]+1, nchar(genes_line)-1), " "), FUN=function(x) substring(x, 2)))
            samples = as.numeric(sapply(strsplit(substring(genes_line, str_locate(samples_line, "Conds \\[[0-9]+\\]: ")[1,2]+1, nchar(samples_line)-1), " "), FUN=function(x) substring(x, 2)))
            if (is.na(genes[[1]])) {
                genes = genes[-1]
            }
            if (is.na(samples[[1]])) {
                samples = samples[-1]
            }
            biclusters = c(biclusters, list(list(genes=genes, samples=samples)))
        }
    }
    return(biclusters)
}