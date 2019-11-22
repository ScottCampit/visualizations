# -*- coding: utf-8 -*-
"""
<<<<<<< HEAD
Visualization script for exploratory data analysis and final figure preparation
"""
# Custom scripts
import preprocess

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec

import scipy
import numpy as np
import pandas as pd
import sys

import plotly
import plotly.graph_objs as go
import plotly.plotly as py
import plotly.tools as tls
import cufflinks as cf

plotly.tools.set_credentials_file(
    username='ScottCampit', api_key='Ngt9S2Je5CbBjBSkNIei')

matplotlib.use('Agg')


def snapshot(fil='filename', fig_name='figurename', ext='file_ext'):
    """
    snapshot does some basic exploratory data analysis of metabolomics data that has been compiled into a single file.
    """





    #df = df.apply(zscore, axis=1)

    #df = pd.DataFrame(df)

    #print(df)
    #df['Value'] = df['Value'].fillna(value=0)

    # Interactive implementation using plotly and dash

    # Heatmap
    #data = [go.Heatmap(z=df.values.tolist(), colorscale='RdBu')]
    #py.plot(data, filename=fig_name)
    # Static output
    sns.heatmap(data=df['Value'], cmap='RdYlBu')
    plt.savefig(fig_name+ext, dpi=300)


snapshot(sys.argv[1], sys.argv[2], sys.argv[3])


def hist(df, filename):
    """
    hist creates a histogram of a single time point.
    """
    #color_lab = df.columns
    #rgb_values = sns.color_palette("muted", len(df.columns))
    #colormap = dict(zip(color_lab, rgb_values))
    sns.distplot(df)  # , color=df.columns.map(colormap))
    plt.xlabel("Scaled Metabolite Intensity")
    #plt.ylabel("")
    plt.savefig(filename+'.png')


def reg(data, filename):
    """
    reg shows a dot plot with a linear regression trend of a single metabolite.
    """
    sns.set()

    x = np.arange(len(data))
    y = data
    sns.regplot(x, y, truncate=True)
    plt.savefig(filename+'.png')


def pairplot(df, filename):
    """
    pairplot will show a pairplot, with the axes corresponding to the time points in the corresponding dataframe.
    """
    sns.pairplot(df)
    plt.savefig(filename+'.png')


def pca(df, filename):
    """
    pca will show the principle component plots, where each point is a specific metabolite and the Scree plot, which compares the size of each eigenvalue.
    """
    from sklearn.decomposition import PCA

    U, S, V = np.linalg.svd(df)
    eigvals = S**2/np.cumsum(S)[-1]

    fig = plt.figure()

    # Create Scree Plot
    pcs = np.arange(len(eigenvals)) + 1
    plt.plot(pcs, eigvals, 'ro-', linewidth=2)
    plt.title('Scree plot')
    plt.xlabel('PC Scores')
    plt.ylabel('Eigenvalue of PC score')
    #pca = PCA(svd_solver='full')
    #principle_components = pca.fit_transform(X)


def volcano_plot_day(cond1_df, cond2_df, filename):
    """
    volcano_plot_day will take post-processed metabolomics data to identify signficant differentially expressed metabolites that are shared between the two cell lines per day.

    Steps:
        1. Get dataframe
        2. Calculate average fold change for each metabolite. This will be the fill value.
        3. Calculate p-value and T-statistic for each metabolite compared to the distribution given.
        4. Create the volcano plots
    """

    # Get rid of positions in both dataframes
    if 'e' in cond1_df.columns:
        pos = cond1_df.pop('e')
        pos = cond2_df.pop('e')
    else:
        pos = pd.concat([cond1_df.pop(x) for x in ['n', 'c', 'm']], 1)
        pos = pd.concat([cond2_df.pop(x) for x in ['n', 'c', 'm']], 1)
    idx = cond2_df.index.intersection(cond1_df.index)
    cond1_df = cond1_df.loc[idx].drop_duplicates(keep='first')
    cond2_df = cond2_df.loc[idx].drop_duplicates(keep='first')

    # I will create 6 dataframes: one for each day and each condition
    cond1_day0 = pd.concat([cond1_df.pop(x) for x in ["0", "0.1", "0.2"]], 1)
    cond1_day2 = pd.concat([cond1_df.pop(x) for x in ["2", "2.1", "2.2"]], 1)
    cond1_day4 = pd.concat([cond1_df.pop(x) for x in ["4", "4.1", "4.2"]], 1)

    cond2_day0 = pd.concat([cond2_df.pop(x) for x in ["0", "0.1", "0.2"]], 1)
    cond2_day2 = pd.concat([cond2_df.pop(x) for x in ["2", "2.1", "2.2"]], 1)
    cond2_day4 = pd.concat([cond2_df.pop(x) for x in ["4", "4.1", "4.2"]], 1)

    # normalized day 2 and day 4:
    cond1_day2n = cond1_day2.div(cond1_day0.values)
    cond1_day4n = cond1_day4.div(cond1_day0.values)

    cond2_day2n = cond2_day2.div(cond2_day0.values)
    cond2_day4n = cond2_day4.div(cond2_day0.values)

    # Get the p value and T-statistic between the two conditions for each day
    try:
        day0_t, day0_pval = scipy.stats.ttest_ind(
            cond1_day0, cond2_day0, axis=1, equal_var=False)
        day2_t, day2_pval = scipy.stats.ttest_ind(
            cond1_day2, cond2_day2, axis=1, equal_var=False)
        day4_t, day4_pval = scipy.stats.ttest_ind(
            cond1_day4, cond2_day4, axis=1, equal_var=False)

    except ValueError:
        day0_t, day0_pval = "NA"
        day2_t, day2_pval = "NA"
        day4_t, day4_pval = "NA"
        #day2_tn, day2_pvaln = "NA"
        #day4_tn, day4_pvaln = "NA"

    idx2 = cond1_day2n.index.intersection(cond2_day2n.index)
    cond1_day2n = cond1_day2n.loc[idx2].drop_duplicates(keep='first')
    cond2_day2n = cond2_day2n.loc[idx2].drop_duplicates(keep='first')

    idx4 = cond1_day4n.index.intersection(cond2_day4n.index)
    cond1_day4n = cond1_day4n.loc[idx4].drop_duplicates(keep='first')
    cond2_day4n = cond2_day4n.loc[idx4].drop_duplicates(keep='first')

    day2_tn, day2_pvaln = scipy.stats.ttest_ind(
        cond1_day2n, cond2_day2n, axis=1, equal_var=False)
    day4_tn, day4_pvaln = scipy.stats.ttest_ind(
        cond1_day4n, cond2_day4n, axis=1, equal_var=False)
    # Get the p value and T-statistic between the days for condition 1
    try:
        cond1_t40, cond1_pval40 = scipy.stats.ttest_ind(
            cond1_day4, cond1_day0, axis=1, equal_var=False)
        cond1_t20, cond1_pval20 = scipy.stats.ttest_ind(
            cond1_day2, cond1_day0, axis=1, equal_var=False)
        cond1_t42, cond1_pval42 = scipy.stats.ttest_ind(
            cond1_day4, cond1_day2, axis=1, equal_var=False)
    except (ValueError, TypeError):
        cond1_t40, cond1_pval40 = "NA"
        cond1_t20, cond1_pval20 = "NA"
        cond1_t42, cond1_pval42 = "NA"

    # Get the p value and T-statistic between the days for condition 2
    try:
        cond2_t40, cond2_pval40 = scipy.stats.ttest_ind(
            cond2_day4, cond2_day0, axis=1, equal_var=False)
        cond2_t20, cond2_pval20 = scipy.stats.ttest_ind(
            cond2_day2, cond2_day0, axis=1, equal_var=False)
        cond2_t42, cond2_pval42 = scipy.stats.ttest_ind(
            cond2_day4, cond2_day2, axis=1, equal_var=False)
    except ValueError:
        cond2_t40, cond2_pval40 = "NA"
        cond2_t20, cond2_pval20 = "NA"
        cond2_t42, cond2_pval42 = "NA"

    # -log p-values
    #day0_pval = -1*np.log10(day0_pval)
    #day2_pval = -1*np.log10(day2_pval)
    #day4_pval = -1*np.log10(day4_pval)
    day2_pvaln = -1*np.log10(day2_pvaln)
    day4_pvaln = -1*np.log10(day4_pvaln)

    c1day40_pval = -1*np.log10(cond1_pval40)
    c1day20_pval = -1*np.log10(cond1_pval20)
    c1day42_pval = -1*np.log10(cond1_pval42)

    c2day40_pval = -1*np.log10(cond2_pval40)
    c2day20_pval = -1*np.log10(cond2_pval20)
    c2day42_pval = -1*np.log10(cond2_pval42)

    # Get the average fold change for each metabolite, scaled between [-1,1] .
    day0_fc = (np.log2(cond1_day0)-np.log2(cond2_day0))/np.log2(cond2_day0)
    day2_fc = (np.log2(cond1_day2)-np.log2(cond2_day2))/np.log2(cond2_day2)
    day4_fc = (np.log2(cond1_day4)-np.log2(cond2_day4))/np.log2(cond2_day4)

    c1day40_fc = (np.log2(cond1_day4)-np.log2(cond1_day0))/np.log2(cond1_day0)
    c1day20_fc = (np.log2(cond1_day2)-np.log2(cond1_day0))/np.log2(cond1_day0)
    c1day42_fc = (np.log2(cond1_day4)-np.log2(cond1_day0))/np.log2(cond1_day0)

    c2day40_fc = (np.log2(cond2_day4)-np.log2(cond2_day0))/np.log2(cond2_day0)
    c2day20_fc = (np.log2(cond2_day2)-np.log2(cond2_day0))/np.log2(cond2_day0)
    c2day42_fc = (np.log2(cond2_day4)-np.log2(cond2_day2))/np.log2(cond2_day2)

    # Params for significant metabolites
    fdr_lower = -5
    fdr_upper = 5
    origin = 0
    pval_cutoff = -1*(np.log10(0.01))

    # These variables are for condition comparisons on the same day
    x1 = day4_t
    y1 = day4_pval

    x2 = day2_t
    y2 = day2_pval

    x3 = day0_t
    y3 = day0_pval

    # These variables are for time changes in condition 1
    x4 = cond1_t40
    y4 = c1day40_pval

    x5 = cond1_t20
    y5 = c1day20_pval

    x6 = cond1_t42
    y6 = c1day42_pval

    # These variables are for time changes in condition 2
    x7 = cond2_t40
    y7 = c2day40_pval

    x8 = cond2_t20
    y8 = c2day20_pval

    x9 = cond2_t42
    y9 = c2day42_pval

    # Normalized values:
    x10 = day2_tn
    y10 = day2_pvaln
    day2 = pd.DataFrame(
        {"Day 2 T-score": x10, "Day 2 p-value": y10}, index=cond1_day2n.index)

    x11 = day4_tn
    y11 = day4_pvaln
    day4 = pd.DataFrame(
        {"Day 4 T-score": x11, "Day 4 p-value": y11}, index=cond1_day4n.index)

    stats = day2.merge(day4, how='inner', left_index=True, right_index=True)

    def col_fill(x, y, fc, idx):
        """
        Create the annotations for specific metabolites, including the color of the metabolite point, the fill, and the name of the metabolite.
        """
        col = []
        face = []
        name = []

        for i in range(0, len(x)):
            if ((x[i] <= fdr_upper or x[i] >= fdr_lower) and y[i] <= pval_cutoff):
                col.append('gray')
                face.append('None')
                name.append(' ')
            elif (x[i] < fdr_lower and y[i] > pval_cutoff):
                col.append('b')
                face.append('b')
                name.append(str(idx[i]))
            elif (x[i] > fdr_upper and y[i] > pval_cutoff):
                col.append('r')
                face.append('r')
                name.append(str(idx[i]))
            else:
                col.append('gray')
                face.append('None')
                name.append(' ')

        return col, face, name

    # original
    #col1, face1, name1 = col_fill(x1, y1, day0_fc, cond1_df.index)
    #col2, face2, name2 = col_fill(x2, y2, day2_fc, cond1_df.index)
    #col3, face3, name3 = col_fill(x3, y3, day4_fc, cond1_df.index)

    # normalized
    col4, face4, name4 = col_fill(x10, y10, day0_fc, cond1_df.index)
    col5, face5, name5 = col_fill(x11, y11, day2_fc, cond1_df.index)

    # Figure generation -- I'm thinking of stacking the volcano plots on top of each other
    """
    fig, axarr = plt.subplots(
        3, sharex=True, sharey=True, figsize=(6, 6), tight_layout=True)

    axarr[0].scatter(x1, y1, facecolors=face1, color=col1, s=10)
    axarr[1].scatter(x2, y2, facecolors=face2, color=col2, s=10)
    axarr[2].scatter(x3, y3, facecolors=face3, color=col3, s=10)

    axarr[0].set_title("Day 4", loc="center")
    axarr[1].set_title("Day 2")
    axarr[1].set_ylabel("-log10 p-value")
    axarr[2].set_title("Day 0")
    axarr[2].set_xlabel("Z score (with T-statistic)")

    for i, txt1 in enumerate(name1):
        axarr[0].annotate(txt1, (x1[i], y1[i]), fontsize=5,
                          horizontalalignment='left', verticalalignment='top')
    for j, txt2 in enumerate(name2):
        axarr[1].annotate(txt2, (x2[j], y2[j]), fontsize=5,
                          horizontalalignment='left', verticalalignment='top')
    for k, txt3 in enumerate(name3):
        axarr[2].annotate(txt3, (x3[k], y3[k]), fontsize=5,
                          horizontalalignment='left', verticalalignment='top')

    for ax in axarr.flat:
        ax.set_xlim(-10, 10)
        ax.set_ylim(-1.0, 5)
        ax.axvline(x=origin, linestyle='-', color='k', linewidth=1)

    fig.savefig(r"./volcano/"+filename+".svg", format="svg", dpi=600)
    """
    fig, axarr = plt.subplots(
        2, sharex=True, sharey=True, figsize=(6, 6), tight_layout=True)

    axarr[0].scatter(x11, y11, facecolors=face5, color=col5, s=10)
    axarr[0].set_title("Day 4", loc="center")

    axarr[1].scatter(x10, y10, facecolors=face4, color=col4, s=10)
    axarr[1].set_title("Day 2")
    axarr[1].set_xlabel("Z score (with T-statistic)")

    #for i, txt1 in enumerate(name4):
    #    axarr[1].annotate(txt1, (x10[i], y10[i]), fontsize=5, horizontalalignment='left', verticalalignment='top')
    #for j, txt2 in enumerate(name5):
    #    axarr[0].annotate(txt2, (x11[j], y11[j]), fontsize=5, horizontalalignment='left', verticalalignment='top')

    for ax in axarr.flat:
        #ax.set_xlim(-10, 10)
        ax.set_ylim(-1.0, 5)
        ax.axvline(x=origin, linestyle='-', color='k', linewidth=1)
        ax.axhline(2, linestyle='--', color='k', linewidth=1)
        ax.axvline(5, linestyle='--', color='k', linewidth=1)
        ax.axvline(-5, linestyle='--', color='k', linewidth=1)

    fig.savefig(r"./volcano/"+filename+".svg", format="svg", dpi=600)

    return stats


def DiffRxn(cond1, cond2, filename):
    """
    DFE will show a bar plot that shows differentially expressed reactions
    similar to the plot in Chandrasekaran et al., 2017.
    """

    # Read in growth rates from rxn knockout
    df1 = pd.read_excel('./DFA/rxnko.xlsx', sheet_name=cond1,
                        names=["Growth Rate"], index_col=0).fillna(0)
    df2 = pd.read_excel('./DFA/rxnko.xlsx', sheet_name=cond2,
                        names=["Growth Rate"], index_col=0).fillna(0)

    def twoSampZ(X1, X2, mudiff, sd1, sd2, n1, n2):
        """
        T test for two samples
        """
        from numpy import sqrt, abs, round
        from scipy.stats import norm
        pooledSE = sqrt(sd1**2/n1 + sd2**2/n2)
        z = ((X1 - X2) - mudiff)/pooledSE
        pval = 2*(1 - norm.cdf(abs(z)))
        return z, pval

    m1 = df1.mean()
    m2 = df2.mean()
    mu_diff = m1-m2
    sd1 = df1.std()
    sd2 = df2.std()
    n1 = len(df1.index)
    n2 = len(df2.index)

    zscr, _ = twoSampZ(df1, df2, mu_diff, sd1, sd2, n1, n2)

    fig = plt.figure()
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=2)
    ax1.hist(zscr.values, bins=15, edgecolor='k', color="#3498db")
    ax1.set_yscale('log')

    title_hist = str("Total reactions screened: 3744")
    ax1.set_title(title_hist, loc="center", fontsize=12)
    ax1.set_title(title_hist, loc="center", fontsize=12)
    ax1.grid(True, which="both", axis="y")
    ax1.set_axisbelow(True)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.set_ylim([0, 10000])

    zscr = pd.DataFrame(zscr)
    ax2 = plt.subplot2grid((3, 3), (2, 0), colspan=3, rowspan=1, sharex=ax1)
    ax2 = sns.stripplot(data=zscr, color="#3498db", s=5,
                        orient="h", jitter=True, edgecolor='black', linewidth=0.5)
    #ax2.hist(zscr, bins=100, edgecolor='k', color="#3498db")
    ax2.set_xlabel("Z-score")
    ax2.get_yaxis().set_ticks([])
    ax2.get_yaxis().set_ticklabels([])
    plt.setp(ax1.get_xticklabels(), visible=False)

    fig.text(0.05, 0.5, "Number of metabolic reactions",
             ha='center', va='center', rotation='vertical')

    plt.savefig(r'./diffrxn/'+str(filename)+'.png', dpi=600)


def DiffGene(cond1, cond2, filename):
    """
    DFE will show a bar plot that shows differentially expressed reactions
    similar to the plot in Chandrasekaran et al., 2017.
    """

    # Read in growth rates from rxn knockout
    df1 = pd.read_excel('./DFA/geneko.xlsx', sheet_name=cond1,
                        names=["Growth Rate"], index_col=0).fillna(0)
    df2 = pd.read_excel('./DFA/geneko.xlsx', sheet_name=cond2,
                        names=["Growth Rate"], index_col=0).fillna(0)

    def twoSampZ(X1, X2, mudiff, sd1, sd2, n1, n2):
        """
        T test for two samples
        """
        from numpy import sqrt, abs, round
        from scipy.stats import norm
        pooledSE = sqrt(sd1**2/n1 + sd2**2/n2)
        z = ((X1 - X2) - mudiff)/pooledSE
        pval = 2*(1 - norm.cdf(abs(z)))
        return z, pval

    m1 = df1.mean()
    m2 = df2.mean()
    mu_diff = m1-m2
    sd1 = df1.std()
    sd2 = df2.std()
    n1 = len(df1.index)
    n2 = len(df2.index)

    zscr, _ = twoSampZ(df1, df2, mu_diff, sd1, sd2, n1, n2)

    fig = plt.figure()
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=2)
    ax1.hist(zscr.values, bins=15, edgecolor='k', color="#3498db")
    ax1.set_yscale('log')
    title_hist = str("Total genes screened: 1488")
    ax1.set_title(title_hist, loc="center", fontsize=12)
    ax1.grid(True, which="both", axis="y")
    ax1.set_axisbelow(True)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.set_ylim([0, 10000])

    ax2 = plt.subplot2grid((3, 3), (2, 0), colspan=3, rowspan=1, sharex=ax1)
    ax2 = sns.stripplot(data=zscr, color="#3498db", s=5,
                        orient="h", jitter=True, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel("Z-score")
    ax2.get_yaxis().set_ticks([])
    ax2.get_yaxis().set_ticklabels([])

    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.setp(ax1.get_xticklabels(), visible=False)

    fig.text(0.05, 0.5, "Number of metabolic genes",
             ha='center', va='center', rotation='vertical')

    plt.savefig(r'./diffexp/'+str(filename)+'.png', dpi=600)


def sensitivity_plots(fil=sys.argv[1], type=sys.argv[2], filename=sys.argv[3]):
    """
    """

    # Argument that will control several parameters in the figure
    if type == "Genes":
        type = "genes"
        number_of_data_points = '1488'
        ylim_max = 1500
    elif type == "Reactions":
        type = "reactions"
        number_of_data_points = '3744'
        ylim_max = 4000

    # Read in data from MATLAB
    df = pd.read_csv(fil)

    # Create figure
    fig = plt.figure()
    # Make the histogram containing the distribution of gene / reactions
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=2)
    ax1.hist(df['Tstatistic'], bins=30, edgecolor='k', color="#0000EE")

    # Set histogram parameters
    ax1.set_yscale('log')
    title_hist = str("Total "+type+" screened: "+number_of_data_points)
    ax1.set_title(title_hist, loc="center", fontsize=8)
    ax1.grid(True, which="both", axis="y")
    ax1.set_axisbelow(True)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.set_ylim([0, ylim_max])
    plt.setp(ax1.get_xticklabels(), visible=False)

    # Create a strpplot of the data to show individual data points within the distribution from the histogram.
    ax2 = plt.subplot2grid((3, 3), (2, 0), colspan=3, rowspan=1, sharex=ax1)
    ax2 = sns.stripplot(data=df['Tstatistic'], color="#ff0000", s=5,
                        orient="h", jitter=True, edgecolor='black', linewidth=0.5)
    #if sys.argv[1] == './../data/GOT1/data/GOT1-rxn.csv':
    # Reactions I want to color specifically
    #    list_of_reactions = ["hypoxanthine phosphoribosyltransferase (Hypoxanthine)", "purine-nucleoside phosphorylase (Guanosine)",
    #                         "purine-nucleoside phosphorylase (Inosine)", "Guanine exchange",  "Guanine transport"]
    #    df['Special'] = df['ID'].isin(list_of_reactions)

    # Get colors corresponding to reactions of interest
    #    df['Color'] = df['Special'].replace(True, "#ff0000")
    #    df['Color'] = df['Color'].replace(False, #0000EE)

    #    colors = dict(zip(df['ID'], df['Color']))
    #    ax2 = sns.stripplot(x='ID', y='Tstatistic', data=df, hue='Color', color=colors,
    #                        s=5, orient="h", jitter=True, edgecolor='black', linewidth=0.5)
    #else:
    #    ax2 = sns.stripplot(data=df['Tstatistic'], color="#0000EE", s=5,
    #                        orient="h", jitter=True, edgecolor='black', linewidth=0.5)

    # Strip plot parameters
    ax2.set_xlabel("Z-score (T-statistic)")
    ax2.get_yaxis().set_ticks([])
    ax2.get_yaxis().set_ticklabels([])

    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    #plt.setp(ax1.get_xticklabels(), visible=False)

    #fig.text(0.05, 0.5, "Number of metabolic"+type+':', ha='center', va='center', rotation='vertical')

    plt.savefig(r'./../figures/sensitivity_plots/'
                + str(filename)+'.svg', dpi=300)


#sensitivity_plots(fil=sys.argv[1], type=sys.argv[2], filename=sys.argv[3])
