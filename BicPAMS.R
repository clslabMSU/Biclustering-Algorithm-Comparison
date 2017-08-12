BicPAMS.cluster = function(input.filename, BicPAMS.path, is.windows=FALSE, ...) {
    if (is.windows) {
        command = paste("java -jar ", BicPAMS.path, " --input=\\", input.filename, " --output=\\..\\bin\\bicpams_tmp.out", sep="")
    } else {
        command = paste("java -jar ", BicPAMS.path, " --input=/", input.filename, " --output=/../bin/bicpams_tmp.out", sep="")
    }
    print(command)
    system(command)
    print(input.filename)
    return(BicPAMS.ParseOutput("../bin/bicpams_tmp.out"))
}

BicPAMS.ParseOutput = function(filename) {
    html = suppressWarnings(readLines(filename))
    bracketed_expressions = str_locate_all(html, "\\[(.*?)\\]")[[1]]
    biclusters = list()
    for (i in seq(from=1, to=dim(bracketed_expressions)[[1]]/2, by=2)) {
        genes = as.numeric(unlist(strsplit(substr(html, bracketed_expressions[i, 1]+1, bracketed_expressions[i, 2]-1), ", ")))
        samples = as.numeric(unlist(strsplit(substr(html, bracketed_expressions[i+1, 1]+1, bracketed_expressions[i+1, 2]-1), ", ")))
        biclusters = c(biclusters, list(list(genes=genes, samples=samples)))
    }
    return(biclusters)
}