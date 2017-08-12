# must have working directory set to data files location
IO.GatherDatasets = function() {
    datafiles = list.files(pattern="data.\\.txt$")
    resultfiles = list.files(pattern="_hiddenBics\\.txt")
    datasets = list()
    for (i in seq(from=1, to=length(datafiles), by=1)) {
        datasets = c(datasets, list(list(data=datafiles[[i]], ground_truth=resultfiles[[i]])))
    }
    return(datasets)
}


IO.Read.UniBic.InputData = function(filename) {
    input_data = read.table(current_dataset$data, header=FALSE, sep="\t")
    colnames(input_data) = as.character(unlist(input_data[1, ]))
    rownames(input_data) = input_data[, 1]
    return(input_data)
}


IO.Read.UniBic.Biclusters = function(filename) {
    
    # read file
    data_file = readLines(filename)
    biclusters = list()
    
    # parse each line in bicluster file, ignoring first line since it only specifies the number of biclusters
    for (line in 2:length(data_file)) {
        genes_stop_index = str_locate(data_file[[line]], "\\]")
        genes = as.numeric(unlist(as.list(strsplit(substr(data_file[[line]], 12, genes_stop_index - 1), ", ")[[1]])))
        samples = as.numeric(unlist(as.list(strsplit(substr(data_file[[line]], genes_stop_index + 4, nchar(data_file[[line]]) - 2), ", ")[[1]])))
        biclusters = c(biclusters, list(list(genes=genes, samples=samples)))
    }
    
    return(biclusters)
}





