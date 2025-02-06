import yaml
from analysis.process_data import load_data
from analysis.process_data import load_GSE
from analysis.cluster_data import assign_cluster
from analysis.compute_metrics import eval_metrics

def main():
    #load data path
    with open('parameters/params.yaml', 'r') as file:
        params = yaml.safe_load(file)

    raw_data_path = params['data']['raw_data_path']
    seed = params['data']['seed']
    
    #organize sample and metadata data
    gse_file, metadata = load_data(raw_data_path)

    #assign gene expression data to samples
    gene_expression_df, cancer_type_df = load_GSE(gse_file, metadata)
    print("Total samples loaded:", len(gene_expression_df))
    #print(cancer_type_df)
    #print(gene_expression_df)

    #load k (number of clusters) and perform clustering on samples based on subtype label
    k = params['clusters']['k']
    cluster_sample_df, cluster_val_df = assign_cluster(gse_file, gene_expression_df, cancer_type_df, k, seed)
    #print(cluster_val_df)

    #computing clustering metrics
    #acc, prec, recall, f1, auc = eval_metrics(cluster_sample_df, cluster_val_df)
    #print("Clustering metrics:", acc, prec, recall, f1, auc)

if __name__ == '__main__':
    main()