import yaml
from analysis.process_data import load_data
from analysis.process_data import load_GSE
from analysis.compute_metrics import cross_validate
from analysis.plot_metrics import metrics_plot

def main():
    #load data path
    with open('parameters/params.yaml', 'r') as file:
        params = yaml.safe_load(file)

    raw_data_path = params['data']['raw_data_path']
    seed = params['data']['seed']
    
    #organize sample and metadata data
    gse_file, metadata = load_data(raw_data_path)

    #assign gene expression data to samples
    gene_expr_df, cancer_type_df = load_GSE(gse_file, metadata)
    print("Total samples loaded:", len(gene_expr_df))

    #load k (number of clusters) and perform clustering on samples based on subtype label
    k = params['clusters']['k']
    n_folds = params['clusters']['n_folds']
    acc_list, prec_list, recall_list, f1_list, auc_list = cross_validate(gene_expr_df, cancer_type_df, k, n_folds, seed)

    #plot metrics for each fold and average
    avgs = metrics_plot(acc_list, prec_list, recall_list, f1_list, auc_list, n_folds)
    
    print(f'Metric avgs across {n_folds} folds:')
    print("Accuracy:", round(avgs[0], 2))
    print("Precision:", round(avgs[1], 2))
    print("Recall:", round(avgs[2], 2))
    print("f1 score:", round(avgs[3], 2))
    print("ROC_AUC score:", round(avgs[4], 2))

if __name__ == '__main__':
    main()