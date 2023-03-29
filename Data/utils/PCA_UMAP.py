#
# Import Libraries
#

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import umap
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import statsmodels.api as sm
from statsmodels.formula.api import ols


#
# Define Class
#

class PCA_UMAP:

    def __init__(self, xfnv, mdata, file=False):

        #self.palette = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52', '#FD3216' '#00FE35', '#6A76FC', '#FED4C4', '#FE00CE', '#0DF9FF']
        self.palette = px.colors.qualitative.Alphabet

        self.xfnv = xfnv
        self.mdata = mdata[np.isin(mdata['Seqn'], xfnv.index)]

        self.pca = PCA()
        self.pca.fit(xfnv)

        self.xpca = pd.DataFrame(
            self.pca.transform(xfnv),
            index=xfnv.index
        )
        self.xumap = pd.DataFrame(
            umap.UMAP(n_neighbors=5, min_dist=0.3, metric='correlation').fit_transform(xfnv),
            index=xfnv.index
        )

        self.file = file

    def saveToFile(self, fig):
        with open(self.file, 'a') as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', default_height='50%', default_width='80%'))


    def plotReduction(self, column, pcacomp=[0,1], titleLabel=''):

        fig = make_subplots(
            rows=1, cols=2, 
            vertical_spacing=0.05, shared_xaxes=True, 
            subplot_titles=[f'PCA', f'UMAP']
        )

        for n,colval in enumerate(sorted(set(self.mdata[column]))):

            sn_colVal = self.mdata['Seqn'][self.mdata[column]==colval]

            fig.add_trace(go.Scatter(
                x=self.xpca.loc[sn_colVal, pcacomp[0]],
                y=self.xpca.loc[sn_colVal, pcacomp[1]],
                name=colval,
                mode='markers',
                marker_color=self.palette[n]
            ), row=1, col=1)

            fig.add_trace(go.Scatter(
                x=self.xumap.loc[sn_colVal, 0],
                y=self.xumap.loc[sn_colVal, 1],
                name=colval,
                mode='markers',
                marker_color=self.palette[n],
                showlegend=False,
            ), row=1, col=2)

        fig.update_xaxes(title='UMAP1', row=1, col=2)
        fig.update_yaxes(title='UMAP2', row=1, col=2)
        fig.update_xaxes(title=f'PCA{pcacomp[0]+1}: {round(self.pca.explained_variance_ratio_[pcacomp[0]],4)}', row=1, col=1)
        fig.update_yaxes(title=f'PCA{pcacomp[1]+1}: {round(self.pca.explained_variance_ratio_[pcacomp[1]],4)}', row=1, col=1)

        fig.update_layout(
            title=f'Dimensionality Reduction coloured by {column} {titleLabel}',
            # height=800
        )
        fig.show() if not self.file else self.saveToFile(fig)


def PCA_Var(df, mdata, conVars, catVars, n_comp=10):
    # PCA
    pca = PCA(n_components=n_comp).fit(df)
    X = pca.transform(df)

    pv = {'%Var PCA': pca.explained_variance_ratio_*100}; rs = {}

    # Cont vs Cont
    for i in conVars:
        pv[i], rs[i] = [], []
        for j in range(n_comp):
            model = sm.OLS(
                X[:, j],
                sm.add_constant(mdata.set_index('Seqn').loc[df.index, i])
            ).fit()
            pv[i].append(round(model.pvalues[-1], 4))
            rs[i].append(model.rsquared)

    for i in catVars:
        pv[i] = []
        for j in range(n_comp):
            pv[i].append(round(sm.stats.anova_lm(ols(
                'y ~ C(x)', 
                pd.DataFrame({
                    'y': X[:, j],
                    'x': mdata.set_index('Seqn').loc[df.index, i]
                })
            ).fit())['PR(>F)']['C(x)'],4))

    return pd.DataFrame(pv, index=range(1, n_comp+1))