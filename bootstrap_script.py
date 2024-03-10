# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 17:01:12 2023

The script is to run bootstrap analysis on a dataset, the tutuorial follows: 
https://thenode.biologists.com/quantification-of-differences-as-alternative-for-p-values/research/

It will take an excel, the excel should have following format:
    First row of column would be the title of the data
    Columns of data    
The script works by changing the position of the excel sheets to load different data.
The reference column is the first column of your data

@author: Chalmers Chau @ University of Leeds
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tkinter as tk
from tkinter import filedialog
from palettable.cartocolors.qualitative import Vivid_10
import matplotlib.ticker as ticker
import math
import ptitprince as pt
from cmcrameri import cm
import cmasher as cmr
import os

#%%
"""
change the sheet number and axis title here
The first tab on in an excel file is default to 0, second tab is 1, third is 2
change the num_sheet parameter to change the tab of the excel
"""
num_sheet = 0 # the excel sheet position
Axis_title = 'Peak Area (nA·ms)' # change title here to your data title.
bootstrap_rounds = 1000 # how many bootstraps to perform
#%% Figure export function
def file_path(): 
    root = tk.Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    file_full_path = filedialog.askopenfilename()
    root.destroy()
    folder_path = os.path.dirname(file_full_path)
    return file_full_path, folder_path
file_path, folder_path = file_path()
File_name = str(os.path.basename(file_path).split('.')[0])

def figure_export_png(Path = folder_path, Name = File_name + 'figure' + '_png'):
    return plt.savefig(os.path.join(Path, Name)+'.png', format = 'png', dpi = 300)
    
def figure_export_svg(Path = folder_path, Name = File_name + 'figure' + '_svg'):
    return plt.savefig(os.path.join(Path, Name)+'.svg', format = 'svg')

#%% df dataframe construction
df = pd.read_excel(file_path, sheet_name=num_sheet)

def df_melt(dataframe):
    column_names = dataframe.columns.tolist()
    # Melt the dataframe to merge columns 'A' and 'B' into a single column
    melted_df = pd.melt(dataframe,
                        value_vars = column_names,
                        var_name = 'Condition',
                        value_name = 'Value')
    return melted_df

melt_df = df_melt(dataframe = df).dropna()

def dfs_dict(dataframe):
    dfs_dict = {}
    for condition in dataframe['Condition'].unique():
        condition_df = dataframe[dataframe['Condition'] == condition]
        dfs_dict[condition] = condition_df

    medians = {}
    for condition, df in dfs_dict.items():
        medians[condition] = df['Value'].median()
    return dfs_dict, medians

dfs_dict, df_medians = dfs_dict(dataframe = melt_df)
#%% Bootstrap sampling
def bootstrap_sample(dataframe):
    sample_indices = np.random.choice(len(dataframe), size=(len(dataframe)), replace=True)
    bootstrap_data = dataframe.iloc[sample_indices.flatten()]
    return bootstrap_data
def generate_bootstrap_dataframes(df_dictionary, N):
    bootstrap_dict = {}
    for key, df in df_dictionary.items():
        bootstrap_dict[key] = {}
        for i in range(1, N+1):
            bootstrap_data, df_name = bootstrap_sample(dataframe = df), f"df_bs_{i}"
            bootstrap_dict[key][df_name] = bootstrap_data
    return bootstrap_dict

bootstrap_dfs_dict = generate_bootstrap_dataframes(df_dictionary = dfs_dict, N = bootstrap_rounds)
#%% Bootstrap stats
def bootstrap_stat(df_dictionary):
    Stat_bootstrap_dfs_dict = {}
    for key, inner_dict in df_dictionary.items():
        summary_data = []
        for df_name, df in inner_dict.items():
            median_value = df['Value'].median()
            mean_value = df['Value'].mean()
            summary_data.append([key, df_name, median_value, mean_value])
        summary_df = pd.DataFrame(summary_data, columns = ['Condition', 'DataFrame', 'Median', 'Mean'])
        Stat_bootstrap_dfs_dict[key] = summary_df
    return Stat_bootstrap_dfs_dict

Stat_bootstrap_dfs_dict = bootstrap_stat(df_dictionary = bootstrap_dfs_dict)

def bootstrap_stat_effectsize(df_dictionary,
                              reference_df = 0,
                              Data_type = 'Median' # Mean or Median
                              ):
    dfs_list = list(df_dictionary.values())
    reference_df_Median = dfs_list[reference_df][Data_type]
    percentile_data = []
    median_dict = {}
    for i in range(0, len(dfs_list)):
        df = dfs_list[i]
        df_name = list(df_dictionary.keys())[i]
        df['Effect Size'] = df[Data_type] - reference_df_Median

        percentiles = np.percentile(df['Effect Size'], [2.5, 97.5])
        percentile_data.append([df_name, percentiles[0], percentiles[1]])

        median = np.median(df[Data_type])
        median_dict[df_name] = median

    percentile_df = pd.DataFrame(percentile_data, columns=['df_name', '2.5', '97.5'])
    return percentile_df, median_dict

effect_size_percentile, bootstrap_median = bootstrap_stat_effectsize(df_dictionary = Stat_bootstrap_dfs_dict)

def bootstrap_stat_effectsize_median(df_dictionary,
                                     Data_type = 'Effect Size'):
    dfs_list = list(df_dictionary.values())
    effect_size_median_dict = {}
    for i in range(0, len(dfs_list)):
        df = dfs_list[i]
        df_name = list(df_dictionary.keys())[i]
        es_median = np.median(df[Data_type])
        effect_size_median_dict[df_name] = es_median
    return effect_size_median_dict
effect_size_median = bootstrap_stat_effectsize_median(df_dictionary = Stat_bootstrap_dfs_dict)
#%%
Concatenated_Stat_bootstrap_df = pd.concat(Stat_bootstrap_dfs_dict.values(), ignore_index=True)

#%% All plot
colourmap_extract = cmr.take_cmap_colors(cm.batlowS,
                                         100,
                                         return_fmt='hex')

sns.set(style = 'ticks',
        rc = {'axes.facecolor': 'white'},
        palette = colourmap_extract
        )
#%% plot raw strip data
def Plot_data_median(data,
                     median_data):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax = sns.stripplot(data = data,
                       x = 'Condition',
                       y = 'Value',
                       size = 8,
                       linewidth = 0.5,
                       zorder = 1,
                       alpha = 0.2,
                       jitter = True
                       )

    ax = sns.boxplot(medianprops = {'ls': '--',
                                    'lw': 2,
                                    'color': '#D3D3D3'},
                     whiskerprops = {'visible': False},
                     zorder = 10,
                     data = data,
                     x = 'Condition',
                     y = 'Value',
                     showfliers = False,
                     showbox = False,
                     showcaps = False,
                     width = 0.5,
                     )

    x_ticks = ax.get_xticks()
    for i, (condition, median_value) in enumerate(median_data.items()):
        x_coordinate = x_ticks[i]
        ax.scatter(x_coordinate,
                   median_value,
                   s = 200,
                   facecolors = 'none',
                   edgecolors = Vivid_10.hex_colors[0],
                   linewidth = 2,
                   zorder = 100)

    ax.spines['top'].set_visible(False) #despine
    ax.spines['right'].set_visible(False) #despine
    ax.set_ylabel('Integral', fontsize = 18)
    ax.set_xlabel('', fontsize = 18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plt.tight_layout()
    return plt.show()

#%% Plot medians strip
def Plot_Bootstrapped_medians(data,
                              median_data):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax = sns.stripplot(data = data,
                       x = 'Condition',
                       y = 'Median',
                       size = 8,
                       linewidth = 0.5,
                       zorder = 1,
                       alpha = 0.2,
                       jitter=True)

    ax = sns.boxplot(medianprops = {'ls': '--',
                                    'lw': 2,
                                    'color': '#D3D3D3'},
                     whiskerprops = {'visible': False},
                     data = data,
                     x = 'Condition',
                     y = 'Median',
                     showfliers = False,
                     showbox = False,
                     showcaps = False,
                     width = 0.5,
                     zorder = 10
                     )

    x_ticks = ax.get_xticks()
    for i, (condition, median_value) in enumerate(median_data.items()):
        x_coordinate = x_ticks[i]
        ax.scatter(x_coordinate,
                   median_value,
                   s = 200,
                   facecolors = 'none',
                   edgecolors = Vivid_10.hex_colors[0],
                   linewidth = 2,
                   zorder = 100
                   )

    ax.spines['top'].set_visible(False) #despine
    ax.spines['right'].set_visible(False) #despine
    ax.set_ylabel('Bootstrapped Median of Integral', fontsize = 18)
    ax.set_xlabel('', fontsize = 18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plt.tight_layout()
    return plt.show()

#%% effect size difference
def Plot_Effect_size_difference(data, median_data):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax = sns.stripplot(data = data,
                       x = 'Condition',
                       y = 'Effect Size',
                       size = 8,
                       linewidth = 0.5,
                       zorder = 1,
                       alpha = 0.2,
                       jitter=True,
                       )

    ax = sns.boxplot(medianprops = {'ls': '--',
                                    'lw': 2,
                                    'color': '#D3D3D3'},
                     whiskerprops = {'visible': False},
                     zorder = 10,
                     data = data,
                     x = 'Condition',
                     y = 'Effect Size',
                     showfliers = False,
                     showbox = False,
                     showcaps = False,
                     width = 0.5,
                     )

    x_ticks = ax.get_xticks()
    for i, (condition, median_value) in enumerate(median_data.items()):
        x_coordinate = x_ticks[i]
        ax.scatter(x_coordinate,
                   median_value,
                   s = 200,
                   facecolors = 'none',
                   edgecolors = Vivid_10.hex_colors[0],
                   linewidth = 2,
                   zorder = 100)

    ax.spines['top'].set_visible(False) #despine
    ax.spines['right'].set_visible(False) #despine
    ax.set_ylabel('Difference in Effect Size', fontsize = 18)
    ax.set_xlabel('', fontsize = 18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plt.tight_layout()
    return plt.show()
#%% Cloudplot
def Plot_Difference_95CI(data, median_data, percentile):
    fig, ax = plt.subplots(figsize=(7, 7))

    ax = pt.half_violinplot(
                            data = data,
                            y = 'Effect Size', # for horizontal graph, y-axis must be categorical value
                            x = 'Condition', #for vertical graph, x-axis must be categorical value
                            bw = 0.15,
                            width = 0.8, #width of the violin plot (or height),
                            inner = None, #box, quartile, point, stick or None.
                            orient = 'v',
                            linewidth = 0.0, #turn off the edge of the cloud
                            cut = 0.0,
                            scale = 'width', #either area, count, or width
                            offset = 0.05,
                            zorder = 0,
                            )

    ax = sns.boxplot(medianprops = {'ls': '--',
                                    'lw': 0.01,
                                    'color': '#D3D3D3'},
                     whiskerprops = {'visible': False},
                     data = data,
                     x = 'Condition',
                     y = 'Effect Size',
                     showfliers = False,
                     showbox = False,
                     showcaps = False,
                     width = 0.5,
                     zorder = 1,
                     )

    x_ticks = ax.get_xticks()
    for i in range(len(percentile)):
        df_name = percentile['df_name'].iloc[i]
        lower_bound = percentile['2.5'].iloc[i]
        upper_bound = percentile['97.5'].iloc[i]
        x_coordinate = x_ticks[i]

        # Draw the vertical line
        ax.vlines(x_coordinate,
                  lower_bound,
                  upper_bound,
                  color = Vivid_10.hex_colors[1],
                  linestyle = '-',
                  linewidth = 4.0,
                  zorder = 10)

    for i, (condition, median_value) in enumerate(median_data.items()):
        x_coordinate = x_ticks[i]
        ax.scatter(x_coordinate,
                   median_value,
                   s = 200,
                   facecolors = 'none',
                   edgecolors = Vivid_10.hex_colors[0],
                   linewidth = 2.0,
                   zorder = 100)

    ax.axhline(y=0, color='#B8B8B8', linestyle='-')

    ax.spines['top'].set_visible(False) #despine
    ax.spines['right'].set_visible(False) #despine
    ax.set_ylabel('Difference and 95% CI', fontsize = 18)
    ax.set_xlabel('', fontsize = 18)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    plt.tight_layout()
    return plt.show()
#%% grid plot
def Plot_full_process(ax1_raw_data = melt_df,
                      ax1_raw_median = df_medians,
                      ax234_data = Concatenated_Stat_bootstrap_df,
                      ax2_median = bootstrap_median,
                      ax34_median = effect_size_median,
                      percentile = effect_size_percentile,
                      title = 'Current Peak (pA)'
                      ):
    """Initiate grid"""
    fig = plt.figure(figsize=(18, 7), constrained_layout=True)
    gs = fig.add_gridspec(nrows = 1, ncols = 4)

    """ax1 - raw data stripplot"""
    ax1 = fig.add_subplot(gs[0, 0])
    ax1 = sns.stripplot(data = ax1_raw_data,
                       x = 'Condition',
                       y = 'Value',
                       size = 8,
                       linewidth = 0.5,
                       zorder = 1,
                       alpha = 0.2,
                       jitter = True
                       )

    ax1 = sns.boxplot(medianprops = {'ls': '--',
                                    'lw': 2,
                                    'color': '#D3D3D3'},
                     whiskerprops = {'visible': False},
                     zorder = 10,
                     data = ax1_raw_data,
                     x = 'Condition',
                     y = 'Value',
                     showfliers = False,
                     showbox = False,
                     showcaps = False,
                     width = 0.5,
                     )

    x_ticks = ax1.get_xticks()
    for i, (condition, median_value) in enumerate(ax1_raw_median.items()):
        x_coordinate = x_ticks[i]
        ax1.scatter(x_coordinate,
                   median_value,
                   s = 200,
                   facecolors = 'none',
                   edgecolors = Vivid_10.hex_colors[0],
                   linewidth = 2,
                   zorder = 100)

    ax1.spines['top'].set_visible(False) #despine
    ax1.spines['right'].set_visible(False) #despine
    ax1.set_ylabel(title, fontsize = 18)
    ax1.set_xlabel('', fontsize = 18)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax1.set_title('Data and Median', fontsize = 18, loc = 'left')

    """ax2 - Median strip"""
    ax2 = fig.add_subplot(gs[0, 1])
    ax2 = sns.stripplot(data = ax234_data,
                       x = 'Condition',
                       y = 'Median',
                       size = 8,
                       linewidth = 0.5,
                       zorder = 1,
                       alpha = 0.2,
                       jitter=True)

    ax2 = sns.boxplot(medianprops = {'ls': '--',
                                    'lw': 2,
                                    'color': '#D3D3D3'},
                     whiskerprops = {'visible': False},
                     data = ax234_data,
                     x = 'Condition',
                     y = 'Median',
                     showfliers = False,
                     showbox = False,
                     showcaps = False,
                     width = 0.5,
                     zorder = 10
                     )

    x_ticks = ax2.get_xticks()
    for i, (condition, median_value) in enumerate(ax2_median.items()):
        x_coordinate = x_ticks[i]
        ax2.scatter(x_coordinate,
                   median_value,
                   s = 200,
                   facecolors = 'none',
                   edgecolors = Vivid_10.hex_colors[0],
                   linewidth = 2,
                   zorder = 100
                   )

    ax2.spines['top'].set_visible(False) #despine
    ax2.spines['right'].set_visible(False) #despine
    ax2.set_ylabel('Bootstrapped Median of ' + title, fontsize = 18)
    ax2.set_xlabel('', fontsize = 18)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax2.set_title('Bootstrapped Medians', fontsize = 18, loc = 'left')

    """ax3 - Effect size difference"""
    ax3 = fig.add_subplot(gs[0, 2])
    ax3 = sns.stripplot(data = ax234_data,
                       x = 'Condition',
                       y = 'Effect Size',
                       size = 8,
                       linewidth = 0.5,
                       zorder = 1,
                       alpha = 0.2,
                       jitter=True,
                       )

    ax3 = sns.boxplot(medianprops = {'ls': '--',
                                    'lw': 2,
                                    'color': '#D3D3D3'},
                     whiskerprops = {'visible': False},
                     zorder = 10,
                     data = ax234_data,
                     x = 'Condition',
                     y = 'Effect Size',
                     showfliers = False,
                     showbox = False,
                     showcaps = False,
                     width = 0.5,
                     )

    x_ticks = ax3.get_xticks()
    for i, (condition, median_value) in enumerate(ax34_median.items()):
        x_coordinate = x_ticks[i]
        ax3.scatter(x_coordinate,
                   median_value,
                   s = 200,
                   facecolors = 'none',
                   edgecolors = Vivid_10.hex_colors[0],
                   linewidth = 2,
                   zorder = 100)

    ax3.spines['top'].set_visible(False) #despine
    ax3.spines['right'].set_visible(False) #despine
    ax3.set_ylabel('Difference in Effect Size', fontsize = 18)
    ax3.set_xlabel('', fontsize = 18)
    ax3.tick_params(axis='x', labelsize=16)
    ax3.tick_params(axis='y', labelsize=16)
    ax3.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax3.set_title('Differences between medians', fontsize = 18, loc = 'left')

    """ax4 - difference and CI"""
    ax4 = fig.add_subplot(gs[0, 3])
    ax4 = pt.half_violinplot(
                            data = ax234_data,
                            y = 'Effect Size', # for horizontal graph, y-axis must be categorical value
                            x = 'Condition', #for vertical graph, x-axis must be categorical value
                            bw = 0.15,
                            width = 0.8, #width of the violin plot (or height),
                            inner = None, #box, quartile, point, stick or None.
                            orient = 'v',
                            linewidth = 0.0, #turn off the edge of the cloud
                            cut = 0.0,
                            scale = 'width', #either area, count, or width
                            offset = 0.05,
                            zorder = 0,
                            )

    ax4 = sns.boxplot(medianprops = {'ls': '--',
                                    'lw': 0.01,
                                    'color': '#D3D3D3'},
                     whiskerprops = {'visible': False},
                     data = ax234_data,
                     x = 'Condition',
                     y = 'Effect Size',
                     showfliers = False,
                     showbox = False,
                     showcaps = False,
                     width = 0.5,
                     zorder = 1,
                     )

    x_ticks = ax4.get_xticks()
    for i in range(len(percentile)):
        df_name = percentile['df_name'].iloc[i]
        lower_bound = percentile['2.5'].iloc[i]
        upper_bound = percentile['97.5'].iloc[i]
        x_coordinate = x_ticks[i]

        # Draw the vertical line
        ax4.vlines(x_coordinate,
                  lower_bound,
                  upper_bound,
                  color = Vivid_10.hex_colors[1],
                  linestyle = '-',
                  linewidth = 4.0,
                  zorder = 10)

    for i, (condition, median_value) in enumerate(ax34_median.items()):
        x_coordinate = x_ticks[i]
        ax4.scatter(x_coordinate,
                   median_value,
                   s = 200,
                   facecolors = 'none',
                   edgecolors = Vivid_10.hex_colors[0],
                   linewidth = 2.0,
                   zorder = 100)

    ax4.axhline(y=0, color='#B8B8B8', linestyle='-')

    ax4.spines['top'].set_visible(False) #despine
    ax4.spines['right'].set_visible(False) #despine
    ax4.set_ylabel('Difference and 95% CI', fontsize = 18)
    ax4.set_xlabel('', fontsize = 18)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    ax4.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax4.set_title('Difference and 95%CI', fontsize = 18, loc = 'left')

    plt.tight_layout()

    return plt.show()
#%% raw to difference and 95%CI
def Plot_of_differences(ax1_raw_data = melt_df,
                      ax1_raw_median = df_medians,
                      ax234_data = Concatenated_Stat_bootstrap_df,
                      ax34_median = effect_size_median,
                      percentile = effect_size_percentile,
                      title = 'Current Peak (pA)',
                      ):
    """Initiate grid"""
    fig = plt.figure(figsize=(10, 7), constrained_layout=True)
    gs = fig.add_gridspec(nrows = 1, ncols = 2)

    """ax1 - raw data stripplot"""
    ax1 = fig.add_subplot(gs[0, 0])
    ax1 = sns.stripplot(data = ax1_raw_data,
                       x = 'Condition',
                       y = 'Value',
                       size = 8,
                       linewidth = 0.5,
                       zorder = 1,
                       alpha = 0.2,
                       jitter = True
                       )

    ax1 = sns.boxplot(medianprops = {'ls': '--',
                                    'lw': 2,
                                    'color': '#D3D3D3'},
                     whiskerprops = {'visible': False},
                     zorder = 10,
                     data = ax1_raw_data,
                     x = 'Condition',
                     y = 'Value',
                     showfliers = False,
                     showbox = False,
                     showcaps = False,
                     width = 0.5,
                     )

    x_ticks = ax1.get_xticks()
    for i, (condition, median_value) in enumerate(ax1_raw_median.items()):
        x_coordinate = x_ticks[i]
        ax1.scatter(x_coordinate,
                   median_value,
                   s = 200,
                   facecolors = 'none',
                   edgecolors = Vivid_10.hex_colors[0],
                   linewidth = 2,
                   zorder = 100)

    ax1.spines['top'].set_visible(False) #despine
    ax1.spines['right'].set_visible(False) #despine
    ax1.set_ylabel(title, fontsize = 18)
    ax1.set_xlabel('', fontsize = 18)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax1.set_title('Data and Median', fontsize = 18, loc = 'left')

    """ax4 - difference and CI"""
    ax4 = fig.add_subplot(gs[0, 1])
    ax4 = pt.half_violinplot(
                            data = ax234_data,
                            y = 'Effect Size', # for horizontal graph, y-axis must be categorical value
                            x = 'Condition', #for vertical graph, x-axis must be categorical value
                            bw = 0.15,
                            width = 0.8, #width of the violin plot (or height),
                            inner = None, #box, quartile, point, stick or None.
                            orient = 'v',
                            linewidth = 0.0, #turn off the edge of the cloud
                            cut = 0.0,
                            scale = 'width', #either area, count, or width
                            offset = 0.05,
                            zorder = 0,
                            )

    ax4 = sns.boxplot(medianprops = {'ls': '--',
                                    'lw': 0.01,
                                    'color': '#D3D3D3'},
                     whiskerprops = {'visible': False},
                     data = ax234_data,
                     x = 'Condition',
                     y = 'Effect Size',
                     showfliers = False,
                     showbox = False,
                     showcaps = False,
                     width = 0.5,
                     zorder = 1,
                     )

    x_ticks = ax4.get_xticks()
    for i in range(len(percentile)):
        df_name = percentile['df_name'].iloc[i]
        lower_bound = percentile['2.5'].iloc[i]
        upper_bound = percentile['97.5'].iloc[i]
        x_coordinate = x_ticks[i]

        # Draw the vertical line
        ax4.vlines(x_coordinate,
                  lower_bound,
                  upper_bound,
                  color = Vivid_10.hex_colors[1],
                  linestyle = '-',
                  linewidth = 4.0,
                  zorder = 10)

    for i, (condition, median_value) in enumerate(ax34_median.items()):
        x_coordinate = x_ticks[i]
        ax4.scatter(x_coordinate,
                   median_value,
                   s = 200,
                   facecolors = 'none',
                   edgecolors = Vivid_10.hex_colors[0],
                   linewidth = 2.0,
                   zorder = 100)

    ax4.axhline(y=0, color='#B8B8B8', linestyle='-')

    ax4.spines['top'].set_visible(False) #despine
    ax4.spines['right'].set_visible(False) #despine
    ax4.set_ylabel('Difference and 95% CI', fontsize = 18)
    ax4.set_xlabel('', fontsize = 18)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    ax4.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax4.set_title('Difference and 95%CI', fontsize = 18, loc = 'left')

    plt.tight_layout()
    return plt.show()

#%% plot function summary
"""Plot raw data"""
Plot_data_median(data = melt_df, median_data = df_medians)
figure_export_svg(Path = folder_path, Name = File_name + '_figure' + '_01_raw_data')
figure_export_png(Path = folder_path, Name = File_name + '_figure' + '_01_raw_data')
plt.close('all')

"""Plot boostrapped medians"""
Plot_Bootstrapped_medians(data = Concatenated_Stat_bootstrap_df, median_data = bootstrap_median)
figure_export_svg(Path = folder_path, Name = File_name + '_figure' + '_02_bootstrapped_medians')
figure_export_png(Path = folder_path, Name = File_name + '_figure' + '_02_bootstrapped_medians')
plt.close('all')

"""Plot effect_size_differnece"""
Plot_Effect_size_difference(data = Concatenated_Stat_bootstrap_df, median_data = effect_size_median)
figure_export_svg(Path = folder_path, Name = File_name + '_figure' + '_03_effect_size_differences')
figure_export_png(Path = folder_path, Name = File_name + '_figure' + '_03_effect_size_differences')
plt.close('all')

"""Plot difference and 95% CI"""
Plot_Difference_95CI(data = Concatenated_Stat_bootstrap_df, median_data = effect_size_median, percentile = effect_size_percentile)
figure_export_svg(Path = folder_path, Name = File_name + '_figure' + '_04_difference_95pc_CI')
figure_export_png(Path = folder_path, Name = File_name + '_figure' + '_04_difference_95pc_CI')
plt.close('all')

"""Plot full transform process"""
Plot_full_process(ax1_raw_data = melt_df,
                  ax1_raw_median = df_medians,
                  ax234_data = Concatenated_Stat_bootstrap_df,
                  ax2_median = bootstrap_median,
                  ax34_median = effect_size_median,
                  percentile = effect_size_percentile,
                  title = Axis_title
                  )
figure_export_svg(Path = folder_path, Name = File_name + '_figure' + '_05_full_process')
figure_export_png(Path = folder_path, Name = File_name + '_figure' + '_05_full_process')

"""Plot raw data"""
Plot_of_differences(ax1_raw_data = melt_df,
                    ax1_raw_median = df_medians,
                    ax234_data = Concatenated_Stat_bootstrap_df,
                    ax34_median = effect_size_median,
                    percentile = effect_size_percentile,
                    title = Axis_title)
figure_export_svg(Path = folder_path, Name = File_name + '_figure' + '_06_plot_of_differences')
figure_export_png(Path = folder_path, Name = File_name + '_figure' + '_06_plot_of_differences')
plt.close('all')
