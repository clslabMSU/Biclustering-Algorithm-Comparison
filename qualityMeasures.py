import numpy as np
import os; os.environ["R_HOME"] = "/path/to/R/R-3.4.0"; os.environ["R_USER"] = "current_use_account"
from pickle import dump, load
import rpy2.robjects as robjects


# RData file saved from R containing biclusters of each algorithm
SYNTHETIC_RESULTS = r"D:\Work\Biclustering\Reporting\SyntheticDataWithPValue.RData"
robjects.r["load"](SYNTHETIC_RESULTS)

# pickle file to use for storing results dictionary
PICKLE_FILE = "Synthetic_Results.pickle"

# Some algorithms call the first bicluster by index 1...
ONE_INDEXED_ALGORITHMS = ["QUBIC", "ISA"]


def parseSyntheticData(data_file: str, pickle_file: str = None, calculateValidation: bool = False, algorithms: list = None, useCache: bool = False) -> dict:
    """
    Read RData file and parse result into dictionary. Optionally calculate validation measures on biclusters.

    :param data_file: RData file to load.
    :param pickle_file: filename to use when working with pickled data.
    :param calculateValidation: whether or not to calculate bicluster quality measures.
    :param algorithms: list of functions to use for validation. Functions must take a bicluster (np.ndarray) as a parameter and return a float.
    :param useCache: whether or not to use the already present pickle file to save time.
    :return: - dictionary of biclusters, with keys as biclusters[algorithm][dataset][bicluster_number][genes | samples | ASR | MSR | ...]
    """

    if useCache and pickle_file is not None:
        print("Reading synthetic data results from cache... ", end="")
        with open(PICKLE_FILE, "rb") as f:
            results_ = load(f)
        print("done.")

    else:

        print("Parsing synthetic data from %s... " % data_file, end="")
        results_ = {
            "ISA": {},
            "OPSM": {},
            "QUBIC": {},
            "UniBic": {},
            "FABIA": {},
            "BicPAMS": {},
        }

        for algorithm in results_:
            algorithm_results = robjects.r[algorithm + ".results"]
            for dataset_number in range(len(algorithm_results)):
                dataset = str(list(algorithm_results[dataset_number])[0]).split("\"")[1]
                biclusters = algorithm_results[dataset_number][1]
                results_[algorithm][dataset] = []
                for bic_number in range(len(biclusters)):
                    genes = list(map(int, biclusters[bic_number][0]))
                    samples = list(map(int, biclusters[bic_number][1]))
                    if algorithm in ONE_INDEXED_ALGORITHMS:
                        genes = list(map(lambda x:x-1, genes))
                        samples = list(map(lambda x:x-1, samples))
                    if algorithm == "UniBic":
                        results_[algorithm][dataset].append({"genes": genes, "samples": samples, "pvalue": biclusters[bic_number][2]})
                    else:
                        results_[algorithm][dataset].append({"genes": genes, "samples": samples})

        print("done.")


    if calculateValidation and algorithms is not None:
        results_ = externalValidation(results_, algorithms)

    if pickle_file is not None:
        with open(pickle_file, "wb") as f:
            dump(results_, f)

    return results_


def getBicluster(genes, samples, expression_data):
    return expression_data[genes, samples]


def externalValidation(results: dict, validation_functions: list, datasets_directory="synthetic_datasets/") -> dict:
    print("Calculating validation measures: " + ",".join(list(map(lambda x:x.__name__, validation_functions))))
    # TODO - only load each dataset file once
    for algorithm in results:
        print("\t" + algorithm)
        for dataset in results[algorithm]:
            expression_data = np.loadtxt(datasets_directory + dataset, skiprows=1, dtype=str)[:, 1:].astype(np.float)
            for bicluster_number in range(len(results[algorithm][dataset])):
                genes = results[algorithm][dataset][bicluster_number]["genes"]
                samples = results[algorithm][dataset][bicluster_number]["samples"]
                bicluster = expression_data[genes][:, samples]
                for f in validation_functions:
                    results[algorithm][dataset][bicluster_number][f.__name__] = f(bicluster) if bicluster.size != 0 else None
    print("Done.")
    return results


# ====================================== BICLUSTER QUALITY MEASURES ====================================== #

def MSR(bicluster: np.ndarray) -> float:
    """
    Mean Squared Residue Score
    Cheng, Y., & Church, G. M. (2000, August). Biclustering of expression data. In Ismb (Vol. 8, No. 2000, pp. 93-103).

    :param bicluster: - np.ndarray of expression levels of bicluster
    :return: - mean squared residue score, lower is better
    """
    column_means = np.mean(bicluster, axis=0)
    row_means = np.mean(bicluster, axis=1)
    bicluster_mean = np.mean(bicluster.flatten())
    msr = 0
    for i in range(bicluster.shape[0]):
        for j in range(bicluster.shape[1]):
            msr += (bicluster[i, j] - row_means[i] - column_means[j] + bicluster_mean)**2
    return msr / (bicluster.shape[0] * bicluster.shape[1])


def SMSR(bicluster: np.ndarray) -> float:
    """
    Scaled Mean Squared Residue Score
    Mukhopadhyay, A., Maulik, U., & Bandyopadhyay, S. (2009). A novel coherence measure for discovering scaling biclusters from gene expression data. Journal of Bioinformatics and Computational Biology, 7(05), 853-868.

    :param bicluster: - np.ndarray of expression levels of bicluster
    :return: - scaled mean squared residue score, lower is better
    """
    column_means = np.mean(bicluster, axis=0)
    row_means = np.mean(bicluster, axis=1)
    bicluster_mean = np.mean(bicluster.flatten())
    smsr = 0
    for i in range(bicluster.shape[0]):
        for j in range(bicluster.shape[1]):
            smsr += (row_means[i] * column_means[j] - bicluster[i, j] * bicluster_mean)**2 / (row_means[i]**2 * column_means[j]**2)
    return smsr / (bicluster.shape[0] * bicluster.shape[1])


def VE(bicluster: np.ndarray) -> float:
    """
    Virtual Error of a bicluster
    Divina, F., Pontes, B., Giráldez, R., & Aguilar-Ruiz, J. S. (2012). An effective measure for assessing the quality of biclusters. Computers in biology and medicine, 42(2), 245-256.

    :param bicluster: - np.ndarray of expression levels of bicluster
    :return: virtual error score, lower is better
    """
    rho = np.mean(bicluster, axis=0)
    rho_std = np.std(rho)
    if rho_std != 0:
        rho_hat = (rho - np.mean(rho)) / np.std(rho)
    else:
        rho_hat = (rho - np.mean(rho))
    bic_hat = _standardize_bicluster_(bicluster)
    ve = 0
    for i in range(bicluster.shape[0]):
        for j in range(bicluster.shape[1]):
            ve += abs(bic_hat[i, j] - rho_hat[j])
    ve /= (bicluster.shape[0] * bicluster.shape[1])
    return ve


def VEt(bicluster: np.ndarray) -> float:
    """
    Transposed virtual error of a bicluster
    Pontes, B., Giráldez, R., & Aguilar-Ruiz, J. S. (2010, September). Measuring the Quality of Shifting and Scaling Patterns in Biclusters. In PRIB (pp. 242-252). Chicago


    :param bicluster: - np.ndarray of expression levels of bicluster
    :return: transposed virtual error, lower is better
    """
    return VE(np.transpose(bicluster))


def ASR(bicluster: np.ndarray) -> float or None:
    """
    Average Spearman's Rho for a bicluster.
    Ayadi, W., Elloumi, M., & Hao, J. K. (2009). A biclustering algorithm based on a Bicluster Enumeration Tree: application to DNA microarray data. BioData mining, 2(1), 9.

    :param bicluster: - np.ndarray of expression levels of bicluster
    :return: average spearman's rho, close to 1 or -1 is better
    """
    if bicluster.shape[0] <= 1 or bicluster.shape[1] <= 1:
        return None
    spearman_genes = 0
    spearman_samples = 0
    for i in range(bicluster.shape[0] - 1):
        for j in range(i+1, bicluster.shape[0]):
            spearman_genes += spearman(bicluster[i, :], bicluster[j, :])
    for k in range(bicluster.shape[1] - 1):
        for l in range(k+1, bicluster.shape[1]):
            spearman_samples += spearman(bicluster[:, k], bicluster[:, l])
    spearman_genes /= bicluster.shape[0]*(bicluster.shape[0]-1)
    spearman_samples /= bicluster.shape[1]*(bicluster.shape[1]-1)
    asr = 2*max(spearman_genes, spearman_samples)
    return asr


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    """
    Spearman's Rho for two vectors.
    Ayadi, W., Elloumi, M., & Hao, J. K. (2009). A biclustering algorithm based on a Bicluster Enumeration Tree: application to DNA microarray data. BioData mining, 2(1), 9.

    :param x: first vector
    :param y: second vector
    :return: spearman's rho between vectors, close to 1 or -1 is better
    """
    assert(x.shape == y.shape)
    rx = rankdata(x)
    ry = rankdata(y)
    m = len(x)
    coef = 6.0/(m**3-m)
    ans = 0
    for k in range(m):
        ans += (rx[k] - ry[k])**2
    return 1 - coef*ans


def _standardize_bicluster_(bicluster: np.ndarray) -> np.ndarray:
    """
    Standardize a bicluster by subtracting the mean and dividing by standard deviation.
    Pontes, B., Girldez, R., & Aguilar-Ruiz, J. S. (2015). Quality measures for gene expression biclusters. PloS one, 10(3), e0115497.

    Note that UniBic synthetic data was generated with mean 0 and standard deviation 1, so it is already standardized.

    :param bicluster: np.ndarray of expression levels of bicluster
    :return: standardized bicluster
    """
    bic = np.copy(bicluster)
    for i in range(bic.shape[0]):
        gene = bic[i, :]
        std = np.std(gene)
        if std != 0:
            bic[i, :] = (gene - np.mean(gene)) / std
        else:
            bic[i, :] = (gene - np.mean(gene))
    return bic

# ======================================================================================================== #


results = parseSyntheticData(SYNTHETIC_RESULTS, pickle_file="SyntheticDataWithPValue.pickle", calculateValidation=True, algorithms=[MSR, SMSR, VEt, ASR], useCache=False)


# include pvalue for UniBic
for algorithm in results:
    for dataset in results[algorithm]:
        for i in range(len(results[algorithm][dataset])):
            if "pvalue" in results[algorithm][dataset][i]:
                results[algorithm][dataset][i]["pvalue"] = str(results[algorithm][dataset][i]["pvalue"][0])

# write results to tab separated files
print("Writing results to file...")
for algorithm in results:
    print("\t" + algorithm)
    lines = [["Dataset", "Bicluster Number", "Number of Genes", "Number of Samples", "ASR", "MSR", "SMSR", "VEt", "PValue"]]
    for dataset in results[algorithm]:
        for bic in range(len(results[algorithm][dataset])):
            bicluster = results[algorithm][dataset][bic]
            lines.append(list(map(str, [
                dataset, bic, len(bicluster["genes"]), len(bicluster["samples"]), bicluster["ASR"], bicluster["MSR"], bicluster["SMSR"], bicluster["VEt"], str(bicluster["pvalue"][0]) if "pvalue" in bicluster else ""
            ])))
    with open("SyntheticDataSummary\\" + algorithm + ".tsv", "w") as f:
        f.write("\n".join(list(map(lambda x:"\t".join(x), lines))))
pass