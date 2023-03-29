#
# Import libraries
#

import numpy as np
import os
import pandas as pd

import plotly.graph_objects as go
from plotly.subplots import make_subplots


#
# Define Class and Functions
#

class PlotMV:

    def __init__(self, x, mdata, file=False):
        self.x = x
        self.mdata = mdata[np.isin(mdata['Seqn'], x.index)]

        self.file = file
        #if self.file and os.path.exists(self.file): os.remove(self.file)

    def saveToFile(self, fig):
        with open(self.file, 'a') as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', default_height='50%', default_width='80%'))


    def plotSummary(self, titleLabel='', col='Group', g1='C', g2='D'):

        fig = make_subplots(
            rows=1, cols=2, 
            subplot_titles=['MVR in Control vs Disease', 'Accepted features vs MVR threshold']
        )

        datax = self.x.loc[self.mdata['Seqn'][self.mdata[col]==g1], :].isna().sum()/(self.mdata[col]==g1).sum()
        datay = self.x.loc[self.mdata['Seqn'][self.mdata[col]==g2], :].isna().sum()/(self.mdata[col]==g2).sum()
        fig.add_trace(
            go.Scatter(
                x=datax,
                y=datay,
                mode='markers',
            ), row=1, col=1
        )

        tmp = np.arange(0,1,0.01)
        fig.add_trace(
            go.Scatter(
                x=tmp,
                y=[((self.x.isna().sum()/self.x.shape[0])<=i).sum() for i in tmp]
            ), row=1, col=2
        )
 
        fig.add_trace(go.Scatter(
            x=[0,1], y=[0,1], line=dict(dash='dash', width=0.5, color='black'), opacity=0.7, showlegend=False, mode='lines'
        ), row=1, col=1)
        

        r = (-0.01, 1.05*max(datax.max(), datay.max()))
        fig.update_xaxes(title=f'MVR in {g1}', row=1, col=1, range=r)
        fig.update_yaxes(title=f'MVR in {g2}', row=1, col=1, range=r)
        fig.update_xaxes(title='MVR Threshold', row=1, col=2)
        #fig.update_yaxes(title='#features', row=1, col=2)

        fig.update_layout(
            title=f'Missing Value Ratio (MVR) Analysis {titleLabel}',
            showlegend=False
        )

        fig.show() if not self.file else self.saveToFile(fig)

    def plotSummaryObs(self, vline=0.1):
        data = np.arange(0,1, 0.01)

        fig = make_subplots(rows=1, cols=2, subplot_titles=['Observations accepted vs MV Threshold', 'N. Observations by MV ratio'])

        fig.add_trace(go.Scatter(
            x = data,
            y = [
                (self.x.isna().sum(axis=1)/self.x.shape[1]<i).sum()
                for i in data
            ],
            showlegend=False
        ), row=1, col=1)

        fig.add_trace(go.Histogram(
            x = self.x.isna().sum(axis=1)/self.x.shape[1],
            showlegend=False
        ), row=1, col=2)

        fig.update_xaxes(title='MV Threshold', col=1)
        fig.update_xaxes(title='MV Ratio', col=2)
        fig.update_yaxes(title='Number of observations')
        fig.update_layout()

        fig.add_vline(vline, line_width=0.5, line_dash='dash')

        fig.show() if not self.file else self.saveToFile(fig)