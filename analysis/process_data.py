import pandas as pd
import GEOparse
import warnings

def load_data(raw_data_path):
    """
    Loads in raw GEO data and metadata from provided datapath

    Args:
        raw_data_path (str): path to GEO dataset defined in params.yaml file - zipped dataset should be located in same directory as run_model.py

    Returns:
        gse_file (GEOparse.GSE): series of all 58 GSM samples
        metadata (dataframe): pandas df containing sample ID, GEO accession ID, total gene #, etc.
    """

    #load sample and metadata
    gse_file = GEOparse.get_GEO(filepath=raw_data_path)
    metadata = gse_file.phenotype_data;

    return gse_file, metadata

def load_GSE(gse_file, metadata):
    """
    Function for matching and labeling gene expression values for each respective GSM sample loaded from raw GSE data and metadata

    Args:
        gse_file (GEOparse.GSE): series of all 58 GSM samples
        metadata (dataframe): pandas df containing sample ID, GEO accession ID, total gene #, etc.

    Returns:
        gene_expr_df (dataframe): dataframe of each GSM sample within the GSE file and its respective expression value per gene
        cancer_type_df (dataframe): dataframe of each GSM sample and its corresponding cancer subtype
    """
    
    #ignore output of warnings from reading in GSE data and df manipulation
    with warnings.catch_warnings():
        
        #suppress output of warnings/error pertaining to df manipulation or future package deprecation
        warnings.simplefilter("ignore", FutureWarning)
        warnings.simplefilter("ignore", pd.errors.SettingWithCopyWarning) 
        
        #create empty dataframe with range of total gene IDs and column headers
        col_id = gse_file.gsms['GSM258562'].table["ID_REF"]
        gene_expr_df = pd.DataFrame(index=range(0),columns=range(len(col_id)))
        gene_expr_df.columns = col_id
        
        #combine gene expression values from each GSM sample to a single dataframe
        for sample in gse_file.gsms.keys():
            gse_table = gse_file.gsms[sample].table
            gse_pivot_table = gse_table.pivot_table(columns = "ID_REF", values = "VALUE")
            gene_expr_df = pd.concat([gene_expr_df, gse_pivot_table])
        
        #add GSM sample ID as row header
        gene_expr_df.index = gse_file.gsms.keys()
        
        #create supplementary df containing cancer subtype corresponding to each sample ID
        phenotypes = metadata.iloc[:,[0,1]]
        phenotypes["Label"] = ['AC' if "AC" in x else 'SCC' for x in phenotypes['title']]
        phenotypes["Type"] = ['0' if "AC" in x else '1' for x in phenotypes['title']]
        cancer_type_df = pd.DataFrame(data = phenotypes.iloc[:,:])
        cancer_type_df = cancer_type_df.reset_index().drop(columns = ["geo_accession", "title"]).rename(columns = {"index":"sample"})
                    
    return gene_expr_df, cancer_type_df
