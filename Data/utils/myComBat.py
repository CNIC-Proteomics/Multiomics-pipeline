#
# Import modules
#

import os
import pandas as pd
import subprocess

#
# Define Function
#

def myComBat(data, metadata, batch, catVars, conVars, Rpath=None, 
             Rengine=r"C:\Users\rbarreror\AppData\Local\Programs\R\R-4.2.0\bin\Rscript.exe"):

    if not Rpath: Rpath = r'S:\U_Proteomica\UNIDAD\software\MacrosRafa\data\Metabolomics\PESA_Integromics\Data\utils\RDataExchange'
    if not os.path.exists(Rpath):
        os.makedirs(Rpath)

    df = data.T

    metadata = metadata.set_index('Seqn').loc[df.columns, [*catVars, *conVars, batch]]
    covar = pd.get_dummies(metadata[catVars], drop_first=True).join(metadata.loc[:,conVars])
    covar['Intercept'] = 1
    covar = covar.loc[:, [covar.columns[-1], *covar.columns[:-1]]]

    df.to_csv(os.path.join(Rpath, 'dat.tsv'), sep='\t')
    covar.to_csv(os.path.join(Rpath, 'mod.tsv'), sep='\t')
    metadata.loc[:,[batch]].rename(columns={batch:'batch'}).to_csv(os.path.join(Rpath, 'batch.tsv'), sep='\t')

    res = subprocess.check_output([
        Rengine,
        "--vanilla",
        r"S:\U_Proteomica\UNIDAD\software\MacrosRafa\data\Metabolomics\PESA_Integromics\Data\utils\ComBat.R",
        Rpath
    ],stderr=subprocess.STDOUT)

    print(res.decode('utf-8'))

    dfa = pd.read_csv(os.path.join(Rpath, 'dfa.tsv'), sep='\t')
    dfa.columns = df.columns
    dfa.index = df.index

    #_ = [os.remove(os.path.join(Rpath, i)) for i in ['dat.tsv', 'batch.tsv', 'mod.tsv', 'dfa.tsv']]

    return dfa.T