import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

def eval_metrics(cluster_sample_df, cluster_val_df):
    """
    Function for manually computing accuracy score of kmeans clustering for cancer subtypes of each GSM sample

    Args:
         method (string): "manual" for manual calculation or "scikit" for scikit's accuracy_score calculation of clustering accuracy

    Returns:
        cluster_acc (float): value as percentage representing clustering accuracy using respectve calculation method
    """

    #assigning cluster val to respective cancer subtype
    cluster_sample_df['cluster'] = cluster_sample_df['Type'].apply(lambda x: 1 if "Adenocarcinoma" in x else 0 if "Squamous Cell Carcinoma" in x else -1)
    #print(cluster_sample_df)

    #merging dataframes to compare actual cluster of cancer type (correct cluster) to predicted cluster group (predicted cluster)
    cluster_comparison = cluster_val_df.merge(cluster_sample_df, left_on = "sample", right_on = "sample")
    cluster_comparison = cluster_comparison.rename(columns = {"cluster_x":"predicted cluster", "cluster_y":"correct cluster"})

    #compute eval metrics using scikit-learn
    acc = float(accuracy_score(cluster_comparison["predicted cluster"], cluster_comparison["correct cluster"]))*100
    prec = float(precision_score(cluster_comparison["predicted cluster"], cluster_comparison["correct cluster"]))*100
    recall = float(recall_score(cluster_comparison["predicted cluster"], cluster_comparison["correct cluster"]))*100
    f1 = float(f1_score(cluster_comparison["predicted cluster"], cluster_comparison["correct cluster"]))*100
    auc = float(roc_auc_score(cluster_comparison["predicted cluster"], cluster_comparison["correct cluster"]))*100
        
    return acc, prec, recall, f1, auc, cluster_sample_df