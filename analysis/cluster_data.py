import pandas as pd
from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split

def assign_cluster(gse_file, gene_expression_df, cancer_type_df, k, seed):
    """
    Function for assigning cluster value to each GSM sample

    Args:
        k (int): # of clusters to be assigned (two by default)

    Returns:
        cluster_sample_df (dataframe): dataframe containing GSM sample, corresponding cancer subtype, and custer value assigned by scikit kmeans clustering
        cluster_val_df (dataframe): dataframe containing GSM sample and the assigned clustering value from scikit kmeans clustering
    """

    #assign random cluster (k = 2) to each GSM sample
    kmeans = KMeans(n_clusters = k, random_state = seed)
    gene_expression_df = gene_expression_df.values
    cancer_type_df = cancer_type_df['Type']
    labels_prediction = kmeans.fit_predict(gene_expression_df)

    #split data into training (80%) and testing (20%) sets
    x_train, x_test, y_train, y_test = train_test_split(gene_expression_df, cancer_type_df, test_size=0.2, random_state=seed)



    #create dataframe that contains GSM sample ID and cluster value (0 or 1)
    #cluster_val_df = pd.DataFrame(gse_file.gsms.keys(), columns = ["sample"])
    #cluster_val_df['cluster'] = labels_prediction.tolist()
    
    #merge dfs to produce new df with GMS sample ID, cancer subtype, and cluster value
    #cluster_sample_df = cancer_type_df.merge(cluster_val_df, left_on = "sample", right_on = "sample")
    #cluster_sample_df = cluster_sample_df.loc[:,["sample", "Type", "cluster"]]
            
    return #cluster_sample_df, cluster_val_df