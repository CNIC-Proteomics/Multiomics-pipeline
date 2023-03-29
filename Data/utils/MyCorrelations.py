#
# Import libraries
#

import os
import pandas as pd
import subprocess
import numpy as np
from dotmap import DotMap
from functools import reduce
from sklearn.covariance import GraphicalLassoCV
from sklearn.cross_decomposition import CCA

from time import time

#
# Define Class
#

class Xsplit:
    def __init__(self, xd, mdata, factor, c, d) -> None:
        self.mdata = mdata.loc[xd.index,:]
        self.x = DotMap()
        self.x.a = xd
        self.x.c = xd.loc[self.mdata[factor] == c, :]
        self.x.d = xd.loc[self.mdata[factor] == d, :]


class MyCorrelations:

    def __init__(self, xq, xm, mdata, factor='Group', c='C', d='D') -> None:
        
        self.xqo = Xsplit(xq, mdata, factor, c, d)
        self.xmo = Xsplit(xm, mdata, factor, c, d)
        self.xqmo = Xsplit(xq.join(xm, how='inner'), mdata, factor, c, d)

        self.psk = DotMap()
        self.gl = DotMap()
        self.rpc = DotMap()
        self.pc = DotMap()
        self.rcca = DotMap()
        self.cca = DotMap()
    
    
    def PSK(self, methods=['pearson', 'spearman', 'kendall']):
        print(f'Start PSK'); starttime = time()

        self.psk = DotMap()

        for xo, omic in zip([self.xqo, self.xmo, self.xqmo], ['qq', 'mm', 'qm']):

            for ss in ['a', 'c', 'd']:

                self.psk[omic][ss] = reduce(
                    lambda df1, df2: df1.add(df2, fill_value=0),
                    [xo.x[ss].corr(m) for m in methods]
                )/len(methods)

            self.psk[omic].dc = self.psk[omic].d - self.psk[omic].c

        for ss in ['a', 'c', 'd', 'dc']:
            self.psk.qm[ss] = self.psk.qm[ss].loc[self.xqo.x.a.columns, self.xmo.x.a.columns]

        print(f'End PSK | {round(time()-starttime, 2)}s')

    def gLasso(self, test=False):
        print(f'Start gLasso'); starttime = time()

        self.gl = DotMap()

        for xo, omic in zip([self.xqo, self.xmo], ['qq', 'mm']):
            print(f'omic: {omic}')
            for ss in ['a', 'c', 'd']:
                print(f'ss: {ss}')
                self.gl[omic][ss] = GraphicalLassoCV(
                    #alphas=[0.2,0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
                    n_refinements= 4 if not test else 1,
                    cv=5 if not test else 2,
                    n_jobs=-1, 
                    mode='cd',
                    assume_centered=True,
                    verbose=True
                    ).fit(xo.x[ss])


        self.pc = DotMap()

        for xo, omic in zip([self.xqo, self.xmo], ['qq', 'mm']):
            for ss in ['a', 'c', 'd']:
                print(f"Omic: {omic} | Subset: {ss} | Alpha: {round(self.gl[omic][ss].alpha_, 4)}")
                mp = self.gl[omic][ss].precision_
                self.rpc[omic][ss] = - mp / np.sqrt(
                    np.outer(
                        np.diagonal(mp), 
                        np.diagonal(mp)
                    )
                )
                self.rpc[omic][ss] = pd.DataFrame(
                    self.rpc[omic][ss],
                    columns=xo.x[ss].columns,
                    index=xo.x[ss].columns,
                )
            
            # Add differential partial correlation
            self.rpc[omic].dc = self.rpc[omic].d - self.rpc[omic].c

        print(f'End gLasso | {round(time()-starttime, 2)}s')

    def PartialCorrelation(self):
        #
        # Calculate non-regularised Partial Correlation
        # We are using Moore-Penrose Pseudo Inverse as in 
        # https://search.r-project.org/CRAN/refmans/ppcor/html/pcor.html
        # https://pubmed.ncbi.nlm.nih.gov/26688802/
        #
        print(f'Start Partial Correlation'); starttime = time()
        self.pc = DotMap()

        for ss in ['a', 'c', 'd']:

            # qq
            Mpr = np.linalg.pinv(
                np.cov(self.xqo.x[ss], rowvar=False)
            )
            Mpc = -Mpr/np.sqrt(np.outer(
                np.diag(Mpr),
                np.diag(Mpr)
            ))
            self.pc.qq[ss] = pd.DataFrame(
                Mpc, 
                columns=self.xqo.x[ss].columns, 
                index=self.xqo.x[ss].columns
            )

            # mm
            Mpr = np.linalg.pinv(
                np.cov(self.xmo.x[ss], rowvar=False)
            )
            Mpc = -Mpr/np.sqrt(np.outer(
                np.diag(Mpr),
                np.diag(Mpr)
            ))
            self.pc.mm[ss] = pd.DataFrame(
                Mpc, 
                columns=self.xmo.x[ss].columns, 
                index=self.xmo.x[ss].columns
            )

        self.pc.qq.dc = self.pc.qq.d - self.pc.qq.c
        self.pc.mm.dc = self.pc.mm.d - self.pc.mm.c

        print(f'End Partial Correlation | {round(time()-starttime, 2)}s')

    def rCCA(
        self, 
        basePath, 
        Rengine=r"C:\Users\rbarreror\AppData\Local\Programs\R\R-4.2.0\bin\Rscript.exe"
        ):

        print(f'Start rCCA'); starttime = time()
        
        if not os.path.exists(basePath):
            os.makedirs(basePath)

        self.rcca.qm = DotMap()

        for ss in ['a', 'c', 'd']:
            qqpath = os.path.join(basePath, f'qq_{ss}.tsv')
            mmpath = os.path.join(basePath, f'mm_{ss}.tsv')
            qmpath = os.path.join(basePath, f'qm_{ss}.tsv')
            
            cidx = np.intersect1d(self.xqo.x[ss].index, self.xmo.x[ss].index)
            
            self.xqo.x[ss].loc[cidx].to_csv(qqpath, sep="\t", index=False)
            self.xmo.x[ss].loc[cidx].to_csv(mmpath, sep="\t", index=False)

            res = subprocess.check_output([
                Rengine,
                "--vanilla",
                r"S:\U_Proteomica\UNIDAD\software\MacrosRafa\data\Metabolomics\PESA_Integromics\Data\utils\rCCA.R",
                qqpath, mmpath, qmpath
            ],stderr=subprocess.STDOUT)

            print(res.decode('utf-8'))

            self.rcca.qm[ss] = pd.read_csv(qmpath, sep='\t')
        
        self.rcca.qm.dc = self.rcca.qm.d - self.rcca.qm.c

        print(f'End rCCA | {round(time()-starttime, 2)}s')

    def CCA(self):
        print(f'Start CCA'); starttime = time()
        
        for ss in ['a', 'c', 'd']:

            X, Y = self.xqo.x[ss], self.xmo.x[ss]
            cidx = np.intersect1d(X.index, Y.index)
            X, Y = X.loc[cidx], Y.loc[cidx]

            cca = CCA(n_components=2)
            variates = cca.fit_transform(X, Y)

            # Equiangular vectors
            Z = variates[0] + variates[1]

            # Project feature over canonical variates
            Xproj = np.array([
                (
                    np.corrcoef(Z[:,0], i, rowvar=False)[0,1],
                    np.corrcoef(Z[:,1], i, rowvar=False)[0,1]
                )
                for i in X.to_numpy().T
            ])

            Yproj = np.array([
                (
                    np.corrcoef(Z[:,0], i, rowvar=False)[0,1],
                    np.corrcoef(Z[:,1], i, rowvar=False)[0,1]
                )
                for i in Y.to_numpy().T
            ])

            # Calculate correlation through projected inner product
            ccaCor = Xproj @ Yproj.T

            self.cca.qm[ss] = pd.DataFrame(
                ccaCor,
                columns=self.xmo.x[ss].columns,
                index=self.xqo.x[ss].columns
            )

        self.cca.qm.dc = self.cca.qm.d - self.cca.qm.c

        print(f'End CCA | {round(time()-starttime, 2)}s')