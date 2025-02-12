# **Read Me:**

Repository for implementing kmeans to classify NSCLC types.

# **Project Overview:**

The goal of this repo is to evaluate the kmeans mathematical model from scikit-learn to discriminate between two cancer subtypes - adenocarcinoma and squamous cell carcinoma - using gene expression data from GEO. 

In conjunction with evaluating the accuracy, precision, recall, and f1 score of correctly identifying lung cancer subtypes, further biological relevance is examined by calculating the statistical significance of gene expression value differences and cross-validating potential biomarkers with published literature.

Additional information on the biological relevance of this analysis and in-depth step-by-step details of the workflow can be found in the *analysis_report.ipynb* notebook.

# **Data:**

The data used in this project - a SOFT formatted family file *GSE10245_RAW.tar* - can be found from the Gene Expression Omnibus (GEO) link [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10245). 

The file contains 58 samples total (40 AC and 18 SCC) and 54,675 genes with corresponding expression values for each respective sample.

The zipped dataset mentioned in the __Data__ section above should be downloaded and moved to this directory. The dataset should be located at the same level as this _README.md_ file and the *run_model.py* script. See below for an example:
```
--NSCLCEval
    |--analysis
        |--...
    |--notebooks
        |--...
    |--parameters
        |--params.yaml
    |--dataset (GSE10245_family.soft.gz)
    |--README.md
    |--requirements.txt
    |--run_model.py
```
# **Running the Model:**

To run the mode, navigate to this directory and enter the following:
`python run_model.py`

Parameters for adjusting the dataset file path, number of folds for cross-validation, number of significant genes to be displayed, and more can be adjusted in the _params.yaml_ file.

# **Output:**

Running the script *run_model.py* using the directions listed in the __Running the Model__ section above will display the following outputs:

1) Number and IDs of samples contained in dataset
2) Metric performance averages across all folds of kmeans classification
3) Barplot of metric performance for individual folds of kmeans classification
4) Mean gene expression differences between NSCLC subtypes
5) Statistically significant genes
6) Gene search engine prompt

# **Dependencies:**

Prerequisites for running this model are included in the _requirements.txt_ file.

# **License:**

This project is licensed under the MIT license. See the LICENSE file/tab for more info.
