"""

@author: Scott Campit
"""

import plotly
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

def histograms(df):
    """

    :return:
    """
    return go.Figure(
        px.histogram(df,
                     x="value",
                     color="variable",
                     nbins=50,
                     opacity=0.6,
                     labels={'value':'Feature distribution'},
                     title="Distribution of the Iris dataset")
    )


def boxplots(df):
    """

    :param df:
    :return:
    """
    return go.Figure(
        px.box(
            df,
            x='Species',
            y='value',
            color='Species',
            points='all')
    )

def splom(df):
    """

    :param df:
    :return:
    """
    splom = px.scatter_matrix(
        df,
        dimensions=["sepal width (cm)", "sepal length (cm)",
                    "petal width (cm)", "petal length (cm)"],
        color="Species")
    splom.update_traces(showupperhalf=False)
    return go.Figure(splom)

def parallel(df):
    """

    :param df:
    :return:
    """
    return go.Figure(
        px.parallel_coordinates(df,
                                color='Species')
    )
