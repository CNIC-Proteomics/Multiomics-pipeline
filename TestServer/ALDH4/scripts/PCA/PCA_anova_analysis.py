import argparse
import os
import sys
import pandas as pd
from sklearn.decomposition import PCA
import json
import logging
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import tukey_hsd

#
# Constants
#
n_components = 10


#
# Local Functions
#

def create_folder_if_not_exists(folder_path):
    if not os.path.exists(folder_path):
        try:
            os.makedirs(folder_path)
            print(f"Folder '{folder_path}' created successfully.")
        except OSError as e:
            print(f"Error creating folder '{folder_path}': {e}")
    else:
        print(f"Folder '{folder_path}' already exists.")


# Read the JSON formatted DataFrame
def read_dataframe(file_name):
    try:
        df = pd.read_json(file_name, orient='records')
        return df
    except FileNotFoundError:
        logging.error(f"The file '{file_name}' does not exist.")
        return None

# Perform PCA and save the results in JSON files
def perform_pca_and_save(x, n_components, outfolder_path):
    # Create a PCA object with the desired number of components
    pca = PCA(n_components=n_components)
    
    # Perform PCA on the DataFrame's data
    principal_components = pca.fit_transform(x)
    
    # Create a DataFrame for the projections
    projections_df = pd.DataFrame(principal_components, columns=[f'PCA_{i+1}' for i in range(n_components)], index=x.index)
    
    # Create a DataFrame for the loadings
    loadings_df = pd.DataFrame(pca.components_.T, columns=[f'PCA_{i+1}' for i in range(n_components)], index=x.columns)
    
    # Calculate the fraction of variability explained by each component
    explained_variance = pca.explained_variance_ratio_
    explained_variance_df = pd.DataFrame({'Explained_Variance': explained_variance}, index=projections_df.columns)
    
    # Save the results in JSON files
    projections_df.to_json(os.path.join(outfolder_path, "projections.json"), orient='records')
    loadings_df.to_json(os.path.join(outfolder_path, "loadings.json"), orient='records')
    explained_variance_df.to_json(os.path.join(outfolder_path, "explained_variance.json"), orient='records')
    logging.info('PCA results written')
    
    return projections_df, loadings_df, explained_variance_df


# Calculate ANOVA and Tukey
def get_anova_tukey(projections_df, mdata, mdataType, outfolder_path):

    anova_res = {}
    
    for ycol in projections_df.columns.tolist():
        anova_res[ycol] = {}
        
        for xcol in mdata.columns.tolist():    
            anova_res[ycol][xcol] = {}
            
            df_lm = projections_df[[ycol]].join(mdata[[xcol]], how='inner').dropna(axis=0).rename(columns={ycol:'Y', xcol:'X'})
            
            
            try: 
                logging.info(f'Calculating ANOVA: {ycol} - {xcol}')
                if mdataType[xcol]['type']=='categorical':                
                    lm = ols('Y ~ C(X)', df_lm).fit()
                    anova = sm.stats.anova_lm(lm, typ=2).rename(columns={'PR(>F)': 'pvalue'})
                    anova_res[ycol][xcol] = dict(anova.loc['C(X)', :])
                    
                elif mdataType[xcol]['type']=='numeric':
                    lm = ols('Y ~ X', df_lm).fit()
                    anova = sm.stats.anova_lm(lm, typ=2).rename(columns={'PR(>F)': 'pvalue'})
                    anova_res[ycol][xcol] = dict(anova.loc['X', :])            
            
            except Exception as e:
                #anova_res[(ycol, xcol)] = None
                logging.error(f'ANOVA error in {(ycol, xcol)}')
                logging.error(e)
                continue
            
            
            if mdataType[xcol]['type']=='categorical' and mdataType[xcol]['nlevels']>2:
                try:
                    logging.info(f'Calculating Tukey-HSD: {ycol} - {xcol}')
                    anova_res[ycol][xcol]['tukey_hsd'] = {}
                    
                    groups, values = list(zip(*list(df_lm.groupby('X'))))    
                    values = [i.Y.values for i in values]
                    tukey = tukey_hsd(*values)    
                
                    anova_res[ycol][xcol]['tukey_hsd']['groups'] = groups
                    anova_res[ycol][xcol]['tukey_hsd']['pvalues'] = tukey.pvalue.tolist()
                
                except Exception as e:
                    logging.error(f'Tukey-HSD error in {(ycol, xcol)}')
                    logging.error(e)
                    continue
            
    with open(os.path.join(outfolder_path, 'anova_res.json'), 'w') as f:
        json.dump(anova_res, f)
                
    return anova_res


#
# Main
#

def main(args):
    
    omic, x_filepath, mdata_filepath, mdataType_filepath, index_filepath, outfolder_path =\
        args.omic, args.x_filepath, args.mdata_filepath, args.mdataType_filepath, args.index_filepath, args.outfolder_path
    
    #
    # Read tables
    #
    logging.info('Reading tables...')
    
    with open(index_filepath, 'r') as f:
        myIndex = json.load(f)
    logging.info(f'index file read: {index_filepath}')
        
    with open(mdataType_filepath, 'r') as f:
        mdataType = json.load(f)
    logging.info(f'mdataType file read: {mdataType_filepath}')
    
    x = read_dataframe(x_filepath)
    logging.info(f'x file read: {x_filepath} read')
    
    mdata = read_dataframe(mdata_filepath)
    logging.info(f'mdataType file read: {mdataType_filepath}')
    

    #
    # Set indexes
    #
    x.index = myIndex[f'x{omic}']
    mdata.index = myIndex['mdata']
    
    
    #
    # Calculate PCA
    #
    logging.info('Calculating PCA...')
    projections_df, loadings_df, explained_variance_df = perform_pca_and_save(
        x, 
        n_components,
        outfolder_path
        )
    logging.info("PCA completed")
    
    logging.info('Calculating ANOVA and Tukey-HSD...')
    anova_res = get_anova_tukey(projections_df, mdata, mdataType, outfolder_path)
    logging.info("ANOVA and Tukey-HSD completed")
    
    return anova_res


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
        Data preprocessing tool. 
        Exec example: 
            python pca_anova_analysis.py omicType path/to/x path/to/mdata path/to/mdataType path/to/index path/to/outfolder
    """)
    
    parser.add_argument("omic", help="Omic type {q, m}")
    parser.add_argument("x_filepath", help="Path to file with quantitative data")
    parser.add_argument("mdata_filepath", help="Path to file with mdata")
    parser.add_argument("mdataType_filepath", help="Path to file with mdataType")
    parser.add_argument("index_filepath", help="Path to file with index")
    parser.add_argument("outfolder_path", help="Path to output folder where files will be written")
    args = parser.parse_args()
    
    create_folder_if_not_exists(args.outfolder_path)
    
    logFile=f'{args.outfolder_path}/preprocessing.log'
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(stream=sys.stdout),
            logging.FileHandler(logFile)
        ]
    )
    
    logging.info('Start main')
    anova_res = main(args)
    logging.info('End main')