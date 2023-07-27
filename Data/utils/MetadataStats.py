#
# Import Libraries
#

import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from scipy.stats import mannwhitneyu, chi2_contingency

#
# Define Classes
#

class MetadataStats:

    def __init__(self, mdata, groupCol='Group', C='C', D='D', file=False) -> None:
        
        
        self.mdata = mdata
        self.groupCol = groupCol; self.c = C; self.d = D
        
        self.cbool = mdata[self.groupCol] == self.c
        self.dbool = mdata[self.groupCol] == self.d

        self.file = file

        self.palette = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']


    def saveToFile(self, fig):
        with open(self.file, 'a') as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', default_height='50%', default_width='90%'))


    def plotQuanCols(self, quanCols):
        fig = make_subplots(rows=1, cols=len(quanCols), subplot_titles=quanCols, horizontal_spacing=0.4/(len(quanCols)-1))
        for ncol, col in enumerate(quanCols):
            for n,g in enumerate([self.c, self.d]):
                fig.add_trace(go.Box(
                    y=self.mdata.loc[self.mdata[self.groupCol] == g, col],
                    marker_color=self.palette[n],
                    name=g,
                    hoverinfo='skip',
                    showlegend=False,
                    boxpoints='all',
                    pointpos=-1.8,
                    jitter=0.3, width=0.5, marker_size=2
                ), row=1, col=ncol+1)

        fig.update_layout(
            title=f"Quantitative Metadata Distribution between Groups"
        )
        
        fig.show() if not self.file else self.saveToFile(fig)
        self._myUtest(quanCols)


    def _myUtest(self, quanCols):
        mdata = self.mdata
        print('Mann-Whitney U Test\n')
        for i in quanCols:
            print(i)
            tmp = mannwhitneyu(
                mdata.loc[mdata[self.groupCol]==self.d, i],
                mdata.loc[mdata[self.groupCol]==self.c, i],
                use_continuity=False, nan_policy='omit'
            )
            print(f"Statistic = {tmp.statistic} | p-value = {tmp.pvalue}")
            print()


    def plotQualCols(self, qualCols):

        fig = make_subplots(rows=1, cols=len(qualCols), subplot_titles=qualCols)

        #col = qualCols[1]
        for n,col in enumerate(qualCols):
            gc = self.mdata.loc[self.mdata[self.groupCol]==self.c,col].value_counts()
            gd = self.mdata.loc[self.mdata[self.groupCol]==self.d,col].value_counts()


            fig.add_trace(go.Bar(
                name=self.c, x=['No', 'Yes'], y=[gc[0], gc[1]], marker_color=self.palette[0],
                showlegend=True if n==0 else False,
                width=0.6
            ), row=1, col=n+1)

            fig.add_trace(go.Bar(
                name=self.d, x=['No', 'Yes'], y=[gd[0], gd[1]], marker_color=self.palette[1],
                showlegend=True if n==0 else False,
                width=0.6
            ), row=1, col=n+1)

        fig.update_layout(barmode='stack', title='Qualitative Metadata Distribution between Groups')
        fig.show() if not self.file else self.saveToFile(fig)
        self._myChiTest(qualCols)

    def _myChiTest(self, qualCols):
        mdata = self.mdata
        print('Chi-Square Homogeneity Sample Test\n')
        for i in qualCols:
            print(i)
            tmp = chi2_contingency(
                pd.crosstab(mdata[self.groupCol], mdata[i]),
                correction=False
            )

            print(f"Statistic = {tmp.statistic} | p-value = {tmp.pvalue}")
            print()