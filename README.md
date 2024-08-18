# **Read Me:**

Repository for BIOINF 575 final project.

# **Project Overview:**

The goal of this project (described and executed in the file `BIOINF575_FinalProject.ipynb`) is to evaluate a mathematical model from scikit-learn to discriminate between two cancer subtypes: adenocarcinoma and squamous cell carcinoma. 

In conjunction with evaluating the accuracy score of correctly identifying lung cancer types, further biological relevance is examined by calculating the statistical significance of gene expression value differences between the two subtypes and cross-validating potential biomarkers with published literature.

# **Data:**

The data used in this project - a SOFT formatted family file `GSE10245_RAW.tar` - can be found from the Gene Expression Omnibus (GEO) link [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10245). 

The file contains 58 samples total (40 AD and 18 SCC) and 54,675 genes with corresponding expression values for each respective sample.
