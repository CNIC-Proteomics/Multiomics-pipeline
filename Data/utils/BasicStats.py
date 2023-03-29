#
# Import libraries
#

import numpy as np
import pandas as pd
import itertools

from statsmodels.discrete.discrete_model import Logit
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.tools import add_constant
import scipy.stats
import scikit_posthocs as sp

from statsmodels.formula.api import ols
import statsmodels.api as sm

#
# Define functions & classes
#


class Xobj:
    """
    Class used to generate Y and X numpy arrays without missing values.
    Observations with missing values in Y are deleted (X should not contain nan) 
    """
    def __init__(self, mdata, xd) -> None:
        self.mdata = mdata
        self.xd = xd

    def getYX(self, md_col):
        df = self.mdata.loc[:, [md_col]].join(self.xd).dropna(axis=0, how='any')
        return (
            df[md_col].to_numpy(),
            df.drop(md_col, axis=1).to_numpy()
        )


class Xstats:
    """
    """

    def __init__(self, mdata, xd, qualCols, quanCols) -> None:
        
        # Object to obtain pairs Y, X
        self.xdo = Xobj(mdata, xd)
        
        # Binary columns in metadata (1s and 0s)
        self.qualCols = qualCols

        # Quantitative columns in metadata
        self.quanCols = quanCols

        # Dictionary with results
        self.res = {i:{} for i in qualCols+quanCols}

    
    #
    # Logistic Regression
    #

    def LogR(self):
        '''
        Simple Logistic Regression
        '''

        for i in self.qualCols:
            self.res[i]['LogR'] = {}

            Y, X = self.xdo.getYX(i)

            LR = [Logit(Y, add_constant(j)).fit() for j in X.T]
            self.res[i]['LogR']['params'] = np.array([j.params[1] for j in LR])
            self.res[i]['LogR']['pvalues'] = np.array([j.pvalues[1] for j in LR])

    def MLogR(self, topN=None, pvalueThr=None):
        '''
        Multiple Logistic Regression
        Features are selected based on Simple Logistic Regression:
            - topN: Select best N features in LogR
            - pvalueThr: Select features below
            - If both are None, use all
        '''

        for i in self.qualCols:
            self.res[i]['MLogR'] = {}
            
            # Get X and Y
            Y, X = self.xdo.getYX(i)

            # feature index
            if pvalueThr:
                fidx = self.res[i]['LogR']['pvalues']<pvalueThr
            
            elif topN:
                fidx = self.res[i]['LogR']['pvalues'].argsort().argsort() < topN
            
            else:
                fidx = np.ones(X.shape[1], dtype=bool)


            # Initialize values
            self.res[i]['MLogR']['params'] = np.zeros_like(fidx, dtype=np.float32)
            self.res[i]['MLogR']['params'][:] = np.nan
            self.res[i]['MLogR']['pvalues'] = np.ones_like(fidx, dtype=np.float32)
            self.res[i]['MLogR']['pvalues'][:] = np.nan

            
            try:
                myModel = Logit(Y, add_constant(X[:,fidx])).fit()
                self.res[i]['MLogR']['params'][fidx] = myModel.params[1:]
                self.res[i]['MLogR']['pvalues'][fidx] = myModel.pvalues[1:]
            except Exception as e:
                print(f'MLogR Exception with {i}:')
                print(e)

    def RMLogR(self, topN=None, pvalueThr=None):
        '''
        Regularised Multiple Logistic Regression
        Features are selected based on Simple Logistic Regression:
            - topN: Select best N features in LogR
            - pvalueThr: Select features below
            - If both are None, use all
        '''

        for i in self.qualCols:
            
            self.res[i]['RMLogR'] = {}
            
            # Get X & Y
            Y, X = self.xdo.getYX(i)
            
            # Get feature index
            if pvalueThr:
                fidx = self.res[i]['LogR']['pvalues']<pvalueThr
            
            elif topN:
                fidx = self.res[i]['LogR']['pvalues'].argsort().argsort() < topN
            
            else:
                fidx = np.ones(X.shape[1], dtype=bool)

            # Initialize values
            self.res[i]['RMLogR']['params'] = np.zeros_like(fidx, dtype=np.float32)
            self.res[i]['RMLogR']['params'][:] = np.nan
            self.res[i]['RMLogR']['pvalues'] = np.ones_like(fidx, dtype=np.float32)
            self.res[i]['RMLogR']['pvalues'][:] = np.nan

            try:
                myModel = Logit(Y, add_constant(X[:, fidx])).fit_regularized(method='l1', alpha=1)
                self.res[i]['RMLogR']['params'][fidx] = myModel.params[1:]
                self.res[i]['RMLogR']['pvalues'][fidx] = myModel.pvalues[1:]

            except Exception as e:
                print(f'RMLogR Exception with {i}:')
                print(e)

    #
    # Linear Regression
    #

    def LinR(self):
        '''
        Simple Linear Regression
        Features are selected based on Simple Linear Regression:
            - topN: Select best N features in LogR
            - pvalueThr: Select features below
            - If both are None, use all
        '''

        for i in self.quanCols:

            self.res[i]['LinR'] = {}

            Y, X = self.xdo.getYX(i)

            myModels = [OLS(Y, add_constant(j)).fit() for j in X.T]
            self.res[i]['LinR']['params'] = np.array([j.params[1] for j in myModels])
            self.res[i]['LinR']['pvalues'] = np.array([j.pvalues[1] for j in myModels])
    
    def MLinR(self, topN=None, pvalueThr=None):
        '''
        Multiple Linear Regression
        '''

        for i in self.quanCols:

            self.res[i]['MLinR'] = {}
            Y, X = self.xdo.getYX(i)

            # Get feature index
            if pvalueThr:
                fidx = self.res[i]['LinR']['pvalues']<pvalueThr
            
            elif topN:
                fidx = self.res[i]['LinR']['pvalues'].argsort().argsort() < topN
            
            else:
                fidx = np.ones(X.shape[1], dtype=bool)        


            self.res[i]['MLinR']['params'] = np.zeros_like(fidx, dtype=np.float32)
            self.res[i]['MLinR']['params'][:] = np.nan
            self.res[i]['MLinR']['pvalues'] = np.ones_like(fidx, dtype=np.float32)
            self.res[i]['MLinR']['pvalues'][:] = np.nan


            # Do not regularize as it does not generates pvalue
            #myModel = OLS(Y, add_constant(X[:, fidx])).fit_regularized(method='elastic_net', alpha=1, L1_wt=1)
            myModel = OLS(Y, add_constant(X[:, fidx])).fit()
            self.res[i]['MLinR']['params'][fidx] = myModel.params[1:]
            self.res[i]['MLinR']['pvalues'][fidx] = myModel.pvalues[1:]


    #
    # Correlations (Pearson, Spearman, Kendall)
    #

    def correlations(self):

        corrTypes = ['Pearson', 'Spearman', 'Kendall']

        for i in self.quanCols:

            _ = [self.res[i].update({j:{}}) for j in corrTypes]

            Y, X = self.xdo.getYX(i)

            corrs = list(zip(*[
                (
                    list(scipy.stats.pearsonr(j, Y)),
                    list(scipy.stats.spearmanr(j, Y)),#._asdict().values()),
                    list(scipy.stats.kendalltau(j, Y))#._asdict().values()),
                )
                for j in X.T
            ]))

            for n,j in enumerate(corrTypes):
                self.res[i][j]['correlations'], self.res[i][j]['pvalues'] = [np.array(j) for j in list(zip(*corrs[n]))]


    #
    # T-test and U-test (Mann-Whitney)
    #

    def TUtest(self):
        '''
        Calculate T-test and U-test
        '''
        for i in self.qualCols:
            self.res[i]['ttest'], self.res[i]['utest'] = {}, {}

            Y, X = self.xdo.getYX(i)

            # Get both groups
            values = list(set(Y))
            g1 = X[Y == values[0]]
            g2 = X[Y == values[1]]

            # T-test
            myTests = [
                scipy.stats.ttest_ind(
                    a=g1[:,j], 
                    b=g2[:,j], 
                    equal_var=True, 
                    alternative='two-sided'
                ) 
                for j in range(X.shape[1])
            ]

            self.res[i]['ttest']['t'], self.res[i]['ttest']['pvalues'] = [
                np.array(arr) 
                for arr in zip(*[
                    (j.statistic, j.pvalue) for j in myTests
                ])
            ]

            # U-Test
            myTests = [
                scipy.stats.mannwhitneyu(
                    x=g1[:,j], 
                    y=g2[:,j], 
                    alternative='two-sided'
                )
                for j in range(X.shape[1])
            ]

            self.res[i]['utest']['u'], self.res[i]['utest']['pvalues'] = [
                np.array(arr) 
                for arr in zip(*[
                    (j.statistic, j.pvalue) for j in myTests
                ])
            ]
    
    #
    # FDR
    #

    def FDR(self):
        '''
        Calculate FDR for all pvalues using Benjamini/Hochberg method
        '''

        for i in self.res:
            for j in self.res[i]:
                pva = self.res[i][j]['pvalues']
                fidx = ~np.isnan(pva)
                fdr = np.zeros_like(pva)
                fdr[:] = np.nan
                
                if fidx.sum() > 0:
                    fdr[fidx] = multipletests(
                        pva[fidx], 
                        method='fdr_bh', 
                        is_sorted=False
                    )[1] # extract pvalues
                
                self.res[i][j]['fdr'] = fdr


    #
    # Export results as a pandas DataFrame
    #

    def export(self):
        '''
        Export results as a pandas DataFrame
        '''

        resdf = pd.DataFrame(
            {
                (i,j,k): self.res[i][j][k]
                for i in self.res for j in self.res[i] for k in self.res[i][j]
            }, 
            index=self.xdo.xd.columns
        )

        return resdf


#
# ANOVA & Kruskal-Wallis
#

def anova_calc(xd, mdata, factor):
    res = {}

    for i in xd.columns:

        ii = i.replace('-', '_')

        tmp = mdata.loc[:, [factor]].join(xd.loc[:, [i]], how='inner').rename(columns={i:ii})

        aov = {'pvalue': sm.stats.anova_lm(
            ols(f"{ii} ~ {factor}", data=tmp).fit()
        ).loc[factor, 'PR(>F)']}

        hoc = sp.posthoc_tukey(tmp, group_col=factor, val_col=ii)
        [
            aov.update({f'{j[0]}_vs_{j[1]}': hoc.loc[j[0], j[1]]}) 
            for j in itertools.combinations(set(tmp[factor]), 2)
        ]

        res.update({i: aov})

    xds = pd.DataFrame(res).T
    xds['fdr'] = multipletests(
        xds['pvalue'], 
        method='fdr_bh', 
        is_sorted=False
    )[1]
    xds.columns = pd.MultiIndex.from_tuples([(factor, 'anova', i) for i in xds.columns])
    return xds

def kruskal_calc(xd, mdata, factor):
    res = {}

    for i in xd.columns:
        tmp = mdata[[factor]].join(xd[[i]], how='inner')#.groupby(factor).groups.values()
        aov = {'pvalue': scipy.stats.kruskal(*tmp.groupby(factor).agg(list)[i].to_list()).pvalue}

        hoc = sp.posthoc_conover(tmp, val_col=i, group_col=factor, p_adjust = 'fdr_bh')
        [
            aov.update({f'{j[0]}_vs_{j[1]}': hoc.loc[j[0], j[1]]}) 
            for j in itertools.combinations(set(tmp[factor]), 2)
        ]

        res.update({i: aov})

    xds = pd.DataFrame(res).T
    xds['fdr'] = multipletests(
        xds['pvalue'], 
        method='fdr_bh', 
        is_sorted=False
    )[1]
    xds.columns = pd.MultiIndex.from_tuples([(factor, 'kruskal', i) for i in xds.columns])
    return xds


sign_calc = lambda xd, mdata, factor, g, ref: pd.DataFrame(
    mdata.loc[mdata[factor] == g, []].join(xd, how='inner').mean() -\
        mdata.loc[mdata[factor] == ref, []].join(xd, how='inner').mean(),
    columns=pd.MultiIndex.from_tuples([(factor, 'sign', f'{g}-{ref}')])
)