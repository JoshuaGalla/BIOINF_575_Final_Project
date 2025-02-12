import pandas as pd
from sklearn.cluster import KMeans
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

def cross_validate(gene_expr_df, cancer_type_df, k, n_folds, seed):
    """
    Assigns kmeans cluster value to each cross-validated GSM training sample subset and initializes performance metrics

    Args:
        gene_expr_df (df): pandas df containing sample ID and gene expression values
        cancer_type_df (df): pandas df containing sample ID and cancer type label
        n_folds (int): number of folds/data splits to use for cross-validation
        k (int): number of clusters to be assigned during kmeans
        seed (int): set value for kmeans cluster initialization

    Returns:
        acc_list (list): list of model predicition accuracies for each dataset fold after kmeans
        prec_list (list): list of model predicition precision values for each dataset fold after kmeans
        recall_list (list): list of model predicition recall values for each dataset fold after kmeans
        f1_list (list): list of model predicition f1 scores for each dataset fold after kmeans
    """

    #initialize metrics
    acc_list, prec_list, recall_list, f1_list = [], [], [], []

    #define cross-validation split
    cross_val = StratifiedKFold(n_splits = n_folds, shuffle=True, random_state=seed)

    #edit data to only contain expression values and labels
    gene_expr_df = gene_expr_df.values
    cancer_type_df = cancer_type_df['Type']

    #calculate performance metrics for each cross-validation fold (80/20 split)
    for train_idx, val_idx in cross_val.split(gene_expr_df, cancer_type_df):
        gene_expr_train, gene_expr_val = gene_expr_df[train_idx], gene_expr_df[val_idx]
        cancer_type_train, cancer_type_val = cancer_type_df.iloc[train_idx], cancer_type_df.iloc[val_idx]

        #train kmeans on expression data from training data subset
        kmeans = KMeans(n_clusters = k, random_state=seed)
        kmeans.fit(gene_expr_train)

        #predict on validation subset
        cancer_type_pred = kmeans.predict(gene_expr_val)
        cancer_type_val = cancer_type_val.astype(int)

        #calculate metrics
        acc, prec, recall, f1 = eval_metrics(cancer_type_val, cancer_type_pred)

        #add metric from fold to list of values
        acc_list.append(acc)
        prec_list.append(prec)
        recall_list.append(recall)
        f1_list.append(f1)

    return acc_list, prec_list, recall_list, f1_list

def eval_metrics(label_true, label_pred):
    """
    Calculates performance metrics for each data subset label predicted by kmeans

    Args:
        label_true (list): list of ints representing true cancer type labels for given sample
        label_pred (list): list of ints representing kmeans predicted cancer type labels for given sample

    Returns:
        acc (float): model predicition accuracy for given dataset fold after kmeans
        prec (float): model predicition precision value for given dataset fold after kmeans
        recall (float): model predicition recall values for given dataset fold after kmeans
        f1 (float): model predicition f1 scores for given dataset fold after kmeans
    """

    #compute individual metrics for given fold
    acc = accuracy_score(label_true, label_pred)

    #map correctly predicted clusters to cancer type label (AC=0, SCC=1) using majority voting to avoid arbitrary cluster label flipping (<0.50)
    if acc <= 0.50: 
        acc = 1-acc #overall correctness of label
        prec = 1 - precision_score(label_true, label_pred) #correct labeling of specific class
        recall = 1 - recall_score(label_true, label_pred) #correct labeling of true positives
        f1 = 1 - f1_score(label_true, label_pred) #balance between prec and recall - measures false positives and false negatives
    else:
        prec = precision_score(label_true, label_pred)
        recall = recall_score(label_true, label_pred)
        f1 = f1_score(label_true, label_pred)

    return acc, prec, recall, f1