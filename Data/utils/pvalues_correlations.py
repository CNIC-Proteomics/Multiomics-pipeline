#
# Import modules
#

from dotmap import DotMap
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from statsmodels.stats.multitest import multipletests


#
# Define Functions
#

def getH0(corrSL, corrType, omics, subsets):
    '''
    From a list of corr objects, get null hypothesis
    distribution
    '''

    h0c = DotMap()

    for n, corrS in enumerate(corrSL):
        for omic in omics:
            for ss in subsets:
                if n==0:
                    h0c[omic][ss] = np.array([])
                
                arr = getattr(corrS, corrType)[omic][ss].to_numpy()

                if omic != 'qm':
                    arr = arr[np.triu_indices(
                        n=arr.shape[0],
                        k=1
                    )]
                else:
                    arr = arr.flatten()

                h0c[omic][ss] = np.concatenate((
                    h0c[omic][ss],
                    arr
                ))
    
    return h0c



def saveToFile(file, fig):
        with open(file, 'a') as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', default_height='50%', default_width='80%'))


def plotNullAlt(corr, h0c, ctype_omic, absolute=False, size=5*10**5, prefix=''):
    '''
    '''
    palette = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']
    subsets = ['a', 'c', 'd', 'dc']
    subsetTitle = ['All', 'Control', 'Disease', 'Disease - Control']

    for ctype in ctype_omic:

        for omic in ctype_omic[ctype]:

            fig = make_subplots(rows=1, cols=len(subsetTitle), subplot_titles=subsetTitle)

            for n,ss in enumerate(subsets):

                nullX = getattr(h0c, ctype)[omic][ss]

                if nullX.shape[0] > size: nullX = np.random.choice(nullX, size=size, replace=False)

                if absolute: nullX = np.abs(nullX)

                nullX = np.sort(nullX)
                nullY = np.arange(nullX.shape[0])/nullX.shape[0]

                fig.add_trace(go.Scatter(
                    x= nullX,
                    y= nullY,
                    name='Null',
                    marker_color=palette[0],
                    showlegend=True if n==0 else False
                ), row=1, col=n+1)

                realX = getattr(corr, ctype)[omic][ss].to_numpy()

                if omic != 'qm':
                    realX = realX[np.triu_indices(realX.shape[0], 1)]
                else:
                    realX = realX.flatten()

                if absolute: realX = np.abs(realX)

                realX = np.sort(realX)
                realY = np.arange(realX.shape[0])/realX.shape[0]
                
                fig.add_trace(go.Scatter(
                    x= realX,
                    y= realY,
                    name='Real',
                    marker_color=palette[1],
                    showlegend=True if n==0 else False
                ), row=1, col=n+1)

            fig.update_layout(title=f'{ctype} | {omic}')
            saveToFile(f'plots/{prefix}_{omic}.html', fig)


def get_pvalues(corr, h0c, ctype_omic, size=1000):
    pvals = DotMap()
    for ctype in ctype_omic:

        for omic in ctype_omic[ctype]:

            for ss in ['a', 'c', 'd', 'dc']:

                nullD = np.abs(h0c[ctype][omic][ss])
                if nullD.shape[0] > size: nullD = np.random.choice(nullD, size=size, replace=False)

                pvals[ctype][omic][ss] = getattr(corr, ctype)[omic][ss].applymap(lambda x: ((nullD>=x).sum()+1)/(nullD.shape[0]+1))
    
    return pvals


def adjust_pvalues(pvals, ctype_omic):
    adpvals = DotMap()
    for ctype in ctype_omic:
        for omic in ctype_omic[ctype]:
            for ss in ['a', 'c', 'd', 'dc']:
                arr = pvals[ctype][omic][ss].to_numpy().flatten()
                arr = multipletests(arr, method='fdr_bh')[1]
                adpvals[ctype][omic][ss] = pd.DataFrame(
                    arr.reshape(pvals[ctype][omic][ss].shape),
                    columns=pvals[ctype][omic][ss].columns,
                    index=pvals[ctype][omic][ss].index
                )
    
    return adpvals


def graph_size(ctype_omic, pvals, pvalueThr=0.01):
    edges = {}
    nodes = {}
    for omic in ['qq', 'mm', 'qm']:
        for ctype in ctype_omic:
            if omic not in ctype_omic[ctype]: continue
            edges[(omic, ctype)] = {}
            nodes[(omic, ctype)] = {}

            for ss in ['a', 'c', 'd', 'dc']:

                df = pvals[ctype][omic][ss] < pvalueThr
                
                
                if omic == 'qm': 
                    edges[(omic, ctype)][ss] = df.sum().sum()
                    nodes[(omic, ctype)][ss] = df.any(axis=0).sum()+df.any(axis=1).sum()
                else: 
                    
                    edges[(omic, ctype)][ss] = int((df.sum().sum()-np.diag(df).sum())/2)
                    nodes[(omic, ctype)][ss] = df.any().sum()
    
    return (pd.DataFrame(edges), pd.DataFrame(nodes))

def get_corrThr(h0c, pvalueThr, ctype_omic):
    corrThr = {}

    for ctype in ctype_omic:
        for omic in ctype_omic[ctype]:
            corrThr[(omic, ctype)] = {}
            for ss in ['a', 'c', 'd', 'dc']:
                corrThr[(omic, ctype)][ss] = np.quantile(
                    np.abs(h0c[ctype][omic][ss]), 
                    1-pvalueThr
                )

    return pd.DataFrame(corrThr).loc[:, ['qq', 'mm', 'qm']]