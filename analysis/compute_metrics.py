import pandas as pd
from sklearn.cluster import KMeans
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

def cross_validate(gene_expr_df, cancer_type_df, k, n_folds, seed):
    """
    Function for assigning cluster value to each GSM sample

    Args:
        k (int): # of clusters to be assigned (two by default)

    Returns:
        cluster_sample_df (dataframe): dataframe containing GSM sample, corresponding cancer subtype, and custer value assigned by scikit kmeans clustering
        cluster_val_df (dataframe): dataframe containing GSM sample and the assigned clustering value from scikit kmeans clustering
    """

    #initialize metrics
    acc_list, prec_list, recall_list, f1_list, auc_list = [], [], [], [], []

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
        acc, prec, recall, f1, auc = eval_metrics(cancer_type_val, cancer_type_pred)

        #add metric from fold to list of values
        acc_list.append(acc)
        prec_list.append(prec)
        recall_list.append(recall)
        f1_list.append(f1)
        auc_list.append(auc)

    return acc_list, prec_list, recall_list, f1_list, auc_list

def eval_metrics(label_true, label_pred):
    """
    """

    #compute individual metrics for given fold
    acc = accuracy_score(label_true, label_pred)

    if acc <= 0.50: #map correctly predicted clusters to cancer type label (AC=0, SCC=1) using majority voting to avoid arbitrary cluster label flipping (<0.50)
        acc = 1-acc
        prec = 1 - precision_score(label_true, label_pred)
        recall = 1 - recall_score(label_true, label_pred)
        f1 = 1 - f1_score(label_true, label_pred)
        auc = 1 - roc_auc_score(label_true, label_pred)
    else:
        prec = precision_score(label_true, label_pred)
        recall = recall_score(label_true, label_pred)
        f1 = f1_score(label_true, label_pred)
        auc = roc_auc_score(label_true, label_pred)

    return acc, prec, recall, f1, auc