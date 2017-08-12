#######################################################
#              Install required packages:             #
#######################################################
#                                                     #
# ISA:                                                #
#     install.packages("isa2")                        #
#                                                     #
# QUBIC:                                              #
#     source("https://bioconductor.org/biocLite.R")   #
#     biocLite("QUBIC")                               #
#                                                     #
# FABIA:                                              #
#     source("https://bioconductor.org/biocLite.R")   #
#     biocLite("fabia")                               #
#                                                     #
#######################################################
#                                                     #
# BICLUSTER FORMAT:                                   #
#     list(                                           #
#         list(genes=[], samples=[]),                 #
#         list(genes=[], samples=[]),                 #
#         ...                                         #
#     )                                               #
#                                                     #
#######################################################


# PARAMETERS
DATA_PATH = "/path/to/code/Reporting/data"          # must be set to data folder where synthetic data is stored.
OUTFILE = "../out.txt"                              # path to output file, it's handy to tail this file as the program runs (tail -f out.txt)
RESULTS_FILE = "../results.txt"                     # results file, formatted to import easily into Excel or other tool.
OPSM_JAR = "../bin/OPSM.jar"                        # acquire the JAR for OPSM (e.g. from BicAT) and put path here.
BicPAMS_JAR = "../bin/bicpams_4.0.3_console.jar"    # acquire the JAR file for BicPAMS (console version) and put the path here. Website: bicpams.com
UniBic_Binary = "../bin/unibic"                     # use our slightly modified binary file for UniBic, the only change is the output filename.


# prepare environment
library(stringr)
library(isa2)
library(QUBIC)
library(fabia)
setwd(DATA_PATH)
rm(DATA_PATH)


# load biclustering algorithms helpers
source("../ISA.R")
source("../OPSM.R")
source("../BicPAMS.R")
source("../FABIA.R")
source("../UniBic.R")
source("../Metrics.R")
source("../IO.R")
source("../Misc.R")


# LOAD DATA FILES
datasets = IO.GatherDatasets()


# RUN EACH ALGORITHM ON EACH DATASET
iteration = 1

ISA.performance     = list(dataset=list(), relevance=list(), recovery=list())
QUBIC.performance   = list(dataset=list(), relevance=list(), recovery=list())
OPSM.performance    = list(dataset=list(), relevance=list(), recovery=list())
BicPAMS.performance = list(dataset=list(), relevance=list(), recovery=list())
FABIA.performance   = list(dataset=list(), relevance=list(), recovery=list())
UniBic.performance  = list(dataset=list(), relevance=list(), recovery=list())


cat(paste("Executed ", Sys.time()), file=OUTFILE, sep="\n")
cat(paste("Algorithm", "Dataset", "Relevance", "Recovery", sep="\t"), file=RESULTS_FILE, sep="\n")
for (current_dataset in datasets) {
    
    # log progress
    cat(paste("( ", round(100*iteration/length(datasets), digits=2), "% ) Running dataset: ", current_dataset$data, sep=""), file=OUTFILE, append=TRUE, sep="\n")
    
    # load dataset
    input_data = IO.Read.UniBic.InputData(current_dataset$data)
    genuine_biclusters = IO.Read.UniBic.Biclusters(current_dataset$ground_truth)
    
    # UniBic
    write("\tUniBic (UNIversal BIClustering Algorithm)", file=OUTFILE, append=TRUE, sep="\n")
    UniBic.result = UniBic.cluster(current_dataset$data)
    UniBic.metrics = Misc.GatherResults(UniBic.result, OUTFILE, genuine_biclusters)
    UniBic.performance$relevance = c(UniBic.performance$relevance, UniBic.metrics$relevance)
    UniBic.performance$recovery = c(UniBic.performance$recovery, UniBic.metrics$recovery)
    UniBic.performance$dataset = c(UniBic.performance$dataset, current_dataset$data)
    write(paste("UniBic", current_dataset$data, UniBic.metrics$relevance, UniBic.metrics$recovery), file=RESULTS_FILE, append=TRUE, sep="\t")
    
    # ISA
    write("\tISA (Iterative Signature Algorithm):", file=OUTFILE, append=TRUE, sep="\n")
    ISA.result = ISA.cluster(input_data)
    ISA.metrics = Misc.GatherResults(ISA.result, OUTFILE, genuine_biclusters)
    ISA.performance$relevance = c(ISA.performance$relevance, ISA.metrics$relevance)
    ISA.performance$recovery = c(ISA.performance$recovery, ISA.metrics$recovery)
    ISA.performance$dataset = c(ISA.performance$dataset, current_dataset$data)
    write(paste("ISA", current_dataset$data, ISA.metrics$relevance, ISA.metrics$recovery), file=RESULTS_FILE, append=TRUE, sep="\t")
    
    # QUBIC
    write("\tQUBIC (QUalitative BIClustering Algorithm):", file=OUTFILE, append=TRUE, sep="\n")
    QUBIC.result = biclust(data.matrix(input_data[-1, -1]), method=BCQU())
    QUBIC.result = Misc.ConvertBiclust(QUBIC.result)
    QUBIC.metrics = Misc.GatherResults(QUBIC.result, OUTFILE, genuine_biclusters)
    QUBIC.performance$relevance = c(QUBIC.performance$relevance, QUBIC.metrics$relevance)
    QUBIC.performance$recovery = c(QUBIC.performance$recovery, QUBIC.metrics$recovery)
    QUBIC.performance$dataset = c(QUBIC.performance$dataset, current_dataset$data)
    write(paste("QUBIC", current_dataset$data, QUBIC.metrics$relevance, QUBIC.metrics$recovery), file=RESULTS_FILE, append=TRUE, sep="\t")
    
    # OPSM
    write("\tOPSM (Order-Preserving SubMatrix):", file=OUTFILE, append=TRUE, sep="\n")
    OPSM.result = OPSM.cluster(input_data, OPSM_JAR)
    OPSM.metrics = Misc.GatherResults(OPSM.result, OUTFILE, genuine_biclusters)
    OPSM.performance$relevance = c(OPSM.performance$relevance, OPSM.metrics$relevance)
    OPSM.performance$recovery = c(OPSM.performance$recovery, OPSM.metrics$recovery)
    OPSM.performance$dataset = c(OPSM.performance$dataset, current_dataset$data)
    write(paste("OPSM", current_dataset$data, OPSM.metrics$relevance, OPSM.metrics$recovery), file=RESULTS_FILE, append=TRUE, sep="\t")
    
    # BicPAMS
    write("\tBicPAMS (Biclustering PAttern Mining Software):", file=OUTFILE, append=TRUE, sep="\n")
    BicPAMS.result = BicPAMS.cluster(current_dataset$data, BicPAMS_JAR)
    BicPAMS.metrics = Misc.GatherResults(BicPAMS.result, OUTFILE, genuine_biclusters)
    BicPAMS.performance$relevance = c(BicPAMS.performance$relevance, BicPAMS.metrics$relevance)
    BicPAMS.performance$recovery = c(BicPAMS.performance$recovery, BicPAMS.metrics$recovery)
    BicPAMS.performance$dataset = c(BicPAMS.performance$dataset, current_dataset$data)
    write(paste("BicPAMS", current_dataset$data, BicPAMS.metrics$relevance, BicPAMS.metrics$recovery), file=RESULTS_FILE, append=TRUE, sep="\t")
    
    # FABIA
    write("\tFABIA (Factor Analysis for Bicluster Acquisition)", file=OUTFILE, append=TRUE, sep="\n")
    FABIA.result = extractBic(fabia(data.matrix(input_data[-1, -1]), 5, 0.01, 500))
    FABIA.metrics = Misc.GatherResults(FABIA.parseBiclusters(FABIA.result), OUTFILE, genuine_biclusters)
    FABIA.performance$relevance = c(FABIA.performance$relevance, FABIA.metrics$relevance)
    FABIA.performance$recovery = c(FABIA.performance$recovery, FABIA.metrics$recovery)
    FABIA.performance$dataset = c(FABIA.performance$dataset, current_dataset$data)
    write(paste("FABIA", current_dataset$data, FABIA.metrics$relevance, FABIA.metrics$recovery), file=RESULTS_FILE, append=TRUE, sep="\t")

    iteration = iteration + 1
}


# clean up workspace
rm(iteration, current_dataset, input_data, genuine_biclusters)
rm(ISA.result, ISA.metrics)
rm(QUBIC.result, QUBIC.metrics)
rm(BicPAMS.result, BicPAMS.metrics)
rm(FABIA.result, FABIA.metrics)
rm(OPSM.result, OPSM.metrics)
rm(UniBic.result, UniBic.metrics)
