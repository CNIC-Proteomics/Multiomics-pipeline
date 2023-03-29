#
# Import libraries
#

import numpy as np
import os
import pandas as pd

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from sklearn.preprocessing import StandardScaler, RobustScaler
from scipy.stats import kruskal
from scikit_posthocs import posthoc_mannwhitney


#
# Define Class and Functions
#

class PlotEDA:

    def __init__(self, xf, mdata, file=False, scaler='s'):
        
        self.xf = xf
        self.xfn = pd.DataFrame(
            StandardScaler().fit_transform(self.xf) if scaler=='s' else RobustScaler().fit_transform(self.xf),
            columns=self.xf.columns,
            index=self.xf.index
        )

        self.mdata = mdata[np.isin(mdata['Seqn'], self.xf.index)]

        self.palette = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']

        self.file = file
        #if self.file and os.path.exists(self.file): os.remove(self.file)

    def saveToFile(self, fig):
        with open(self.file, 'a') as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', default_height='50%', default_width='80%'))

    def plotSummary(
        self, 
        r11=(-4,4), r12=(-4,4), r21=(-4,4), r22=(-4,4), r3=(-4,4), ry3=(0,2), vl3=[],
        plots = [1,2,3],
        binsize=0.02,
        titleLabel = ''
        ):
        
        if 1 in plots: self.plot_mean_std(self.xf, title=f'Pre-standardization distribution {titleLabel}', r1=r11, r2=r12, binsize=binsize)
        if 2 in plots: self.plot_mean_std(self.xfn, title=f'Post-standardization distribution {titleLabel}', r1=r21, r2=r22, binsize=binsize)
        if 3 in plots: self.plot_prepost_histogram(r=r3, ry=ry3, vl=vl3, binsize=binsize)
    
    
    def plot_mean_std(self, df, title='', r1=(-4,4), r2=(-4,4), binsize=0.02):

        fig = make_subplots(
            rows=1, cols=2, 
            subplot_titles=['Mean Feature Distribution', 'Std Feature Distribution']
        )

        fig.add_trace(go.Histogram(
            x=df.mean(),
            name='Mean',
            opacity=0.75,
            hoverinfo='skip',
            histnorm='probability density',
            xbins=dict(size=binsize)
        ), row=1, col=1)

        fig.add_trace(go.Histogram(
            x=df.std(),
            name='Std',
            opacity=0.75,
            hoverinfo='skip',
            histnorm='probability density',
            xbins=dict(size=binsize)
        ), row=1, col=2)

        fig.update_layout(
            title=title
        )

        fig.update_xaxes(range=r1, row=1, col=1)
        fig.update_xaxes(range=r2, row=1, col=2)

        fig.show() if not self.file else self.saveToFile(fig)
    

    def plot_prepost_histogram(self, r=[-4,4], ry=(0,2), vl=[], binsize=0.02):
        fig = go.Figure()
        fig.add_trace(
            go.Histogram(
                x=self.xf.to_numpy().flatten(),
                name='Pre-Standardized',
                hoverinfo='skip',
                opacity=0.75,
                histnorm='probability density',
                xbins=dict(size=binsize)
            )
        )
        fig.add_trace(
            go.Histogram(
                x=self.xfn.to_numpy().flatten(),
                name='Post-Standardized',
                hoverinfo='skip',
                opacity=0.75,
                histnorm='probability density',
                xbins=dict(size=binsize)
            )
        )
        [fig.add_vline(i, line_width=0.5, line_dash='dash') for i in vl]
        fig.update_layout(
            title='Feature values distribution',
            barmode='overlay'
        )
        fig.update_xaxes(range=r)
        fig.update_yaxes(range=ry)
        fig.show() if not self.file else self.saveToFile(fig)
    
    def _kruskal(self, df, column, showTest=False):
        samples = [
            df.loc[i, :].to_numpy().flatten()
            for i in self.mdata.loc[:, ['Seqn', column]].groupby(column).agg(list)['Seqn'].tolist()
        ]

        test = kruskal(
            *samples,
            nan_policy='omit'#, axis=1
        )

        if showTest:
            print('Kruskal-Wallis:')
            print(test)
            print('PostHoc-MannWhitney')
            print(posthoc_mannwhitney(samples,p_adjust='bonferroni'))

        return test

    def plotByGroup(
        self, 
        column, 
        f1=lambda fig: fig, f2=lambda fig:fig, 
        r1=(-4,4), r2=(-4,4), vl1=[], vl2=[], binsize=0.05, 
        plotN=True, titleLabel=''
        ):

        ncols = 2
        dfs = [self.xf, self.xfn]
        tests = [self._kruskal(df, column) for df in dfs]
        subplot_titles= [f'{n} | Kruskal-Wallis - pvalue={round(t.pvalue,4)}; statistic={round(t.statistic,4)}' for t,n in zip(tests, ['Pre-Standardized', 'Post-Standardized'])]
        if not plotN:
            ncols = 1
            dfs = [dfs[0]]
            subplot_titles=[subplot_titles[0]]

        fig=make_subplots(rows=2, cols=ncols, subplot_titles=subplot_titles, vertical_spacing=0.02, horizontal_spacing=0.1,shared_xaxes=True)
        for n,df in enumerate(dfs):
            for color,i in enumerate(sorted(set(self.mdata[column].tolist()))):
                data = df.loc[self.mdata[self.mdata[column]==i]['Seqn'], :].to_numpy().flatten()
                fig.add_trace(
                    go.Histogram(
                        x=data,
                        name=i,
                        opacity=0.7,
                        marker_color=self.palette[color],
                        hoverinfo='skip',
                        showlegend=True if n==0 else False,
                        histnorm='probability density',
                        xbins=dict(size=binsize)
                    ), row=1, col=n+1
                )
                fig.add_trace(
                    go.Box(
                        x=data,
                        name=i,
                        opacity=0.7,
                        marker_color=self.palette[color],
                        hoverinfo='skip',
                        showlegend=False
                    ), row=2, col=n+1
                )

        [fig.update_xaxes(range=r1, row=i, col=1) for i in [1,2]]
        [fig.add_vline(i, line_width=0.5, line_dash='dash', row=1, col=1) for i in vl1]
        
        if plotN:
            [fig.update_xaxes(range=r2, row=i, col=2) for i in [1,2]]
            [fig.add_vline(i, line_width=0.5, line_dash='dash', row=1, col=2) for i in vl2]

        fig.update_layout(
            barmode='overlay',
            title=f"Feature values distribution per {column} {titleLabel}"
        )
        f1(fig)
        f2(fig)
        fig.show() if not self.file else self.saveToFile(fig)