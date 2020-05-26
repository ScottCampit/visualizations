"""
This is an application that visualizes microRNA data from the Cancer Cell Line Encyclopedia.

@author: Scott Campit
"""

# Load data
from sklearn import datasets
iris = datasets.load_iris()

import pandas as pd
import numpy as np
import graphobjs
data = pd.DataFrame(data=np.c_[iris['data'], iris['target']],
                    columns=iris['feature_names']+['Species'])

# Create pairwise coordinate plot object
pxplt = graphobjs.parallel(data)

# Reformat dataset for other plots
targets = {0.0:'setosa',
           1.0:'versicolor',
           2.0:'virginica'}
data['Species'] = data['Species'].replace(targets)
melted_data = pd.melt(data,
                      id_vars=['Species'],
                      value_vars=iris.feature_names)

# Create other plots
hist = graphobjs.histograms(melted_data)
boxplots = graphobjs.boxplots(melted_data)
splom = graphobjs.splom(data)

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, external_stylesheets=dbc.themes.BOOTSTRAP)
app.layout = [html.Div(
    children=[
        html.H1(
            children="Iris Dataset Dashboard"
        )
    ]),
    html.Div(
        dbc.Row(
            html.Div(
                dcc.Graph(
                    id='pairwisecoordinate',
                    figure=pxplt
                )
            )
        ),
        dbc.Row(
            html.Div(
                dbc.Col(
                    html.Div(
                        dbc.Row(
                            html.Div(
                                dcc.Graph(
                                    id='histogram',
                                    figure=hist
                                )
                            )
                        ),
                    ),
                    html.Div(
                        dbc.Row(
                            html.Div(
                                dcc.Graph(
                                    id='boxplt',
                                    figure=boxplots
                                )
                            )
                        )
                    )
                )
            )
        ),
        html.Div(
            dbc.Col(
                html.Div(
                    dcc.Graph(
                        id='scattermatrix',
                        figure=splom
                    )
                )
            )
        )
    )
]

if __name__ == "__main__":
    app.run_server(debug=True)
