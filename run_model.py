import yaml
from analysis.process_data import load_data
from analysis.process_data import load_GSE
from analysis.compute_metrics import cross_validate
from analysis.plot_metrics import metrics_plot
from analysis.expr_analysis import expr_dif, calc_pval

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
    print('')

    #load k (number of clusters) and perform clustering on samples based on subtype label
    k = params['clusters']['k']
    n_folds = params['clusters']['n_folds']
    acc_list, prec_list, recall_list, f1_list, auc_list = cross_validate(gene_expr_df, cancer_type_df, k, n_folds, seed)

    #plot metrics for each fold and average
    #avgs = metrics_plot(acc_list, prec_list, recall_list, f1_list, auc_list, n_folds)
    
    #print(f'Metric avgs across {n_folds} folds:')
    #print("Accuracy:", round(avgs[0], 2)) #overall correctness of predicitions
    #print("Precision:", round(avgs[1], 2)) #
    #print("Recall:", round(avgs[2], 2))
    #print("f1 score:", round(avgs[3], 2))
    #print("ROC_AUC score:", round(avgs[4], 2))

    #calculate mean gene expression differences between cancer types
    sort = params['genes']['sort']
    mean_expr_all_df, AC_expr_df, SCC_expr_df = expr_dif(gene_expr_df, cancer_type_df, sort)
    n_genes = params['genes']['n_genes']

    print('')
    if sort == 'ascending':
        print(f'Mean gene expression differences between AC and SCC cancer types (Bottom {n_genes} genes):')
    elif sort == 'descending':
        print(f'Mean gene expression differences between AC and SCC cancer types (Top {n_genes} genes):')
    print(mean_expr_all_df.head(n_genes))

    alpha = params['genes']['alpha']
    sort_by = params['genes']['sort_by']
    print('')

    #determine significance of gene expression differences between cancer types
    print('Calculating gene expression difference significance...')
    sig_test_df = calc_pval(AC_expr_df, SCC_expr_df, alpha, sort, sort_by)

    if sort == 'ascending':
        print(f'Statistically significant genes (Sorted by {sort_by}; Bottom {n_genes} genes):')
    elif sort == 'descending':
        print(f'Statistically significant genes (Sorted by {sort_by}; Top {n_genes} genes):') 
    print(sig_test_df.head(n_genes))
    print('')

    #search specific gene
    while True:
        gene_name = input('Enter Gene ID to search (or type "exit" to quit):')

        if gene_name.lower() == 'exit':
            break

        #search values for gene
        gene_searched_expr = mean_expr_all_df.loc[mean_expr_all_df['Gene'] == gene_name]
        gene_searched_pval = sig_test_df.loc[sig_test_df['Gene'] == gene_name]

        #check if gene/input exists
        if not gene_searched_expr.empty:
            print(gene_searched_expr)
            print(gene_searched_pval)
        else:
            print(f"Gene {gene_name} not found in dataset")
        print('')

if __name__ == '__main__':
    main()