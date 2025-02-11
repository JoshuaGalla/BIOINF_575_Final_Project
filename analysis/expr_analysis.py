import pandas as pd
import numpy as np
import scipy.stats as stats

def expr_dif(gene_expr_df, cancer_type_df, sort):
    """
    """

    #reformatting gene expression df index column
    gene_expr_df = gene_expr_df.reset_index()

    #splitting data into gene expression dfs for respective cancer type
    expr_type_df = pd.merge(cancer_type_df, gene_expr_df, left_on='sample', right_on='index', how = "inner")
    AC_expr_df = expr_type_df[expr_type_df['Label'] == 'AC']
    SCC_expr_df = expr_type_df[expr_type_df['Label'] == 'SCC']
    
    #calculate mean expression for each gene across both cancer types
    mean_AC_expr_df = AC_expr_df.drop(columns=['sample', 'Label', 'Type', 'index']).mean(axis=0)
    mean_SCC_expr_df = SCC_expr_df.drop(columns=['sample', 'Label', 'Type', 'index']).mean(axis=0)
    
    #save results to new df
    mean_AC_expr_df = mean_AC_expr_df.reset_index()
    mean_AC_expr_df.columns = ['Gene', 'Mean Expression (AC)']
    mean_SCC_expr_df = mean_SCC_expr_df.reset_index()
    mean_SCC_expr_df.columns = ['Gene', 'Mean Expression (SCC)']

    #merge AC and SCC dfs into one
    mean_expr_all_df = pd.merge(mean_AC_expr_df, mean_SCC_expr_df, on='Gene')

    #calculate mean expression difference of genes between cancer types
    mean_expr_all_df["Mean Expression Diff (abs)"] = abs(mean_expr_all_df['Mean Expression (AC)'] - mean_expr_all_df['Mean Expression (SCC)'])
    if sort == 'ascending':
        mean_expr_all_df = mean_expr_all_df.sort_values(by = "Mean Expression Diff (abs)", ascending = True)
    elif sort == 'descending':
        mean_expr_all_df = mean_expr_all_df.sort_values(by = "Mean Expression Diff (abs)", ascending = False)

    return mean_expr_all_df, AC_expr_df, SCC_expr_df

def calc_pval(AC_expr_df, SCC_expr_df, alpha, sort, sort_by):
    """
    """

    #initialize list for gene expression significance values
    results = []

    #iterate through each gene and calculate significant statistics
    for gene in AC_expr_df.drop(columns=['sample', 'Label', 'Type', 'index']).columns:
        AC_gene = AC_expr_df[gene]
        SCC_gene = SCC_expr_df[gene]

        #perform t-test (Welch's)
        t_stat, pval = stats.ttest_ind(AC_gene, SCC_gene, equal_var=False)

        #calculate log pval
        if pval > 0:
            log_pval = -np.log10(pval)
        else:
            log_pval = np.inf

        #save stats to results dict
        results.append({'Gene': gene, 'p-value': pval, 'Log p-value': log_pval})

    #convert results to df
    sig_test_df = pd.DataFrame(results)

    #perform Bonferroni correction
    num_tests = len(sig_test_df) #54675 genes
    bonferroni_val = alpha/num_tests

    #append Bonferroni values
    sig_test_df['Bonferroni p-value'] = sig_test_df['p-value'] * num_tests
    sig_test_df['Significant (Bonferroni)'] = sig_test_df['Bonferroni p-value'] < alpha

    #edit display of significant genes
    if sort == 'ascending':
        sig_test_df = sig_test_df.sort_values(by = sort_by, ascending = True)
    elif sort == 'descending':
        sig_test_df = sig_test_df.sort_values(by = sort_by, ascending = False)

    return sig_test_df