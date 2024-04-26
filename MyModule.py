import scipy.io
import os
import subprocess
import sys
from IPython.display import display, FileLink
import warnings
import h5py
import numpy as np
from scipy.signal import find_peaks, savgol_filter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# from scipy.stats import sem
from scipy.stats import wilcoxon
import pandas as pd
import json
import jdata as jd
import re
import pickle
import math
import time

from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert import NotebookExporter
import nbformat
from tqdm import tqdm  # Importer tqdm pour la barre de progression

import shutil

import random

from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.image as mpimg

from matplotlib.gridspec import GridSpec
from PIL import Image

# from aquarel import load_theme

from operator import itemgetter

from prettytable import PrettyTable 

import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

import plotly.graph_objects as go

import seaborn as sns


# # qt for popup window (savable as pdf, svg...), inline for inline plot, notebook for interactive plot, widget for interactive plot
# %matplotlib widget 
# plt.ioff()



























































































































def scatter3D(x,y,z,colors,colorlabel,xlabel,ylabel,zlabel,title,filename,filepath,save=True,show=True,anim=True,s=35,alpha=1):

    # theme = load_theme("arctic_light").set_overrides({
    #     "ytick.minor.visible": False,
    #     "xtick.minor.visible": False,
    # })

    # with theme:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    scatter = ax.scatter(x, y, z, c=colors, cmap='coolwarm', marker='o', s=s, alpha=alpha)

    cbar = plt.colorbar(scatter, pad=0.2)
    cbar.set_label(colorlabel)


    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.invert_zaxis()

    plt.title(title)

    if anim:
        def update(frame):
            ax.view_init(30, frame)
            return scatter,

        ani = FuncAnimation(fig, update, frames=np.arange(0, 360, 2), interval=100)

    if save and anim:
        os.makedirs(filepath, exist_ok=True)
        ani.save(os.path.join(filepath, filename), writer='pillow')
    elif save:
        os.makedirs(filepath, exist_ok=True)
        plt.savefig(os.path.join(filepath, filename), format='pdf')

    if show:
        plt.show()


































def getPSTHparameters(StudiedSpikeTimes, timeObject, binResolution):
    local_trial_number = len(StudiedSpikeTimes)

    spike_number_per_trial = [[] for _ in range(local_trial_number)]
    edges = []

    for trial in range(local_trial_number):
        spike_number_per_trial[trial], edges = np.histogram(StudiedSpikeTimes[trial], bins=np.arange(timeObject[0], round(timeObject[-1])+binResolution, binResolution))


    frequency_per_trial = [[spike_number_per_trial[trial][bin]/binResolution for bin in range(len(edges)-1)] for trial in range(local_trial_number)]
    mean_frequency = [np.mean([frequency_per_trial[trial][bin] for trial in range(local_trial_number)]) for bin in range(len(edges)-1)]
    baseline_mean_frequency = [np.mean([frequency_per_trial[trial][bin] for trial in range(local_trial_number)]) for bin in range(len(edges)-1) if edges[bin] < 0]


    Zscore = (mean_frequency - np.mean(baseline_mean_frequency)) / np.std(baseline_mean_frequency) if np.std(baseline_mean_frequency) != 0 else np.zeros(len(mean_frequency))
    Zscore[-1]=Zscore[-2]
    Zunitary = (frequency_per_trial - np.mean(mean_frequency)) / np.std(mean_frequency) if np.std(mean_frequency) != 0 else np.zeros(len(frequency_per_trial))
    SEM = np.std(Zunitary)/np.sqrt(len(Zunitary)) if np.std(mean_frequency) != 0 else np.zeros(len(mean_frequency))

    return edges, Zscore, SEM















def plotPSTH(AllData, animal, condition, unit, plotvelocity=False, color='k',shadedcolor='c',binResolution = 0.03,xlabel=True,ylabel=True, title='', save=False, filename='PSTH.png', show=True, velocitycolor='red', velocityalpha=0.13, extra=None, xlim=None, ylim=None, smooth=True):

    if type(unit)==list:
        if condition == 'phototagging':
            StudiedSpikeTimes = np.concatenate([AllData[animal]['SpikeTimes'][unit_i] for unit_i in unit])
        else:
            StudiedSpikeTimes = np.concatenate([AllData[animal]['SpikeTimes'][condition][unit_i] for unit_i in unit])
    else:
        if condition == 'phototagging':
            StudiedSpikeTimes = AllData[animal]['SpikeTimes'][unit]
        else:
            StudiedSpikeTimes = AllData[animal]['SpikeTimes'][condition][unit]

    timeBef, timeAft = AllData[animal]['timeBef'], AllData[animal]['timeAft']
    if plotvelocity:
        rotationSpeed = AllData[animal]['rotationSpeed']
        duration = AllData[animal]['duration']
        MeanRotation = AllData[animal]['MeanRotation'][condition]

    edges, Zscore, SEM = getPSTHparameters(StudiedSpikeTimes, AllData[animal]['duration'], binResolution)

    if smooth:
        Zscore = savgol_filter(Zscore, 5, 3)

    if plotvelocity:
        normalization = 1/rotationSpeed*max(abs(Zscore))
        
        if min(MeanRotation) < -rotationSpeed/2:
            MeanRotation = MeanRotation * normalization
            haxis = np.zeros(len(duration))
        else:
            MeanRotation = MeanRotation * normalization
            haxis = np.zeros(len(duration))
        plt.plot(duration, MeanRotation, color=velocitycolor, alpha=velocityalpha)
        plt.fill_between(duration, MeanRotation, haxis, color=velocitycolor, alpha=velocityalpha*0.8)

    # plt.figure(figsize=(15,6))
    plt.plot(edges[:-1], Zscore, color=color)
    plt.fill_between(edges[:-1], Zscore-SEM, Zscore+SEM, alpha=0.1, color=shadedcolor)

    if xlim!=None:
        plt.xlim(xlim)
    else:
        plt.xlim(-timeBef,timeAft)
    if ylim!=None:
        plt.ylim(ylim)
    if ylabel:
        plt.ylabel('Z-Score FR')
    if xlabel:
        plt.xlabel('Time (s)')
    plt.grid(False)

    plt.title(title)
            
    if save:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()



















def plotRaster(AllData, animal, condition, unit, color='black', xlabel='Time (s)', ylabel='# Trial', extra=None, title='', plotvelocity=False, velocitycolor='k', velocityalpha=0.13, save=False, filename='Raster.png', show=True, xlim=None, ylim=None, psth=False, smooth=True, binResolution=0.03, psthcolor='b', shadedcolor='b'):
    if type(unit)==list:
        if condition == 'phototagging':
            spikeTimesObject = np.concatenate([AllData[animal]['SpikeTimes'][unit_i] for unit_i in unit])
        else:
            spikeTimesObject = np.concatenate([AllData[animal]['SpikeTimes'][condition][unit_i] for unit_i in unit])
    else:
        if condition == 'phototagging':
            spikeTimesObject = AllData[animal]['SpikeTimes'][unit]
        else:
            spikeTimesObject = AllData[animal]['SpikeTimes'][condition][unit]

    timeBef, timeAft = AllData[animal]['timeBef'], AllData[animal]['timeAft']

    if plotvelocity:
        rotationSpeed = AllData[animal]['rotationSpeed']
        duration = AllData[animal]['duration']
        MeanRotation = AllData[animal]['MeanRotation'][condition]

    linelengths = 1
    
    if ylabel:
        plt.ylabel(ylabel)
    if xlabel:
        plt.xlabel(xlabel)
    
    if xlim!=None:
        plt.xlim(xlim)
    else:
        plt.xlim(-timeBef,timeAft)
    

    ### annulate the offset due to python indexation
    def custom_formatter(x, pos):
        return f"{int(x) + 1}"
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(custom_formatter))

    if ylim!=None:
        plt.ylim(ylim)
    else:
        plt.ylim(0-linelengths/2,len(spikeTimesObject)-1+linelengths/2)
    # plt.grid(False)

    if plotvelocity:
        normalization = 1/rotationSpeed*len(spikeTimesObject)
        
        if min(MeanRotation) < -rotationSpeed/2:
            MeanRotation = MeanRotation * normalization - min(MeanRotation * normalization) - linelengths/2
            haxis = len(spikeTimesObject)*np.ones(len(duration))
        else:
            MeanRotation = MeanRotation * normalization - linelengths/2
            haxis = np.zeros(len(duration)) - linelengths/2
        plt.plot(duration, MeanRotation, color=velocitycolor, alpha=velocityalpha)
        plt.fill_between(duration, MeanRotation, haxis, color=velocitycolor, alpha=velocityalpha*0.8)
        
    plt.eventplot(spikeTimesObject, linelengths=linelengths, colors=color)
    
    if extra is not None:
        if type(extra)==list:
            for ex in extra:
                ex
        else:
            extra

    plt.title(title)

    if psth:
        duration = AllData[animal]['duration']
        edges, Zscore, SEM = getPSTHparameters(spikeTimesObject, duration, binResolution)
        if smooth:
            Zscore = savgol_filter(Zscore, 9, 3)
        plt.twinx()
        plt.plot(edges[:-1], Zscore, color=psthcolor)
        plt.fill_between(edges[:-1], Zscore-SEM, Zscore+SEM, alpha=0.1, color=shadedcolor)
        if ylabel:
            plt.ylabel('Z-Score FR')
        if xlabel:
            plt.xlabel('Time (s)')

    if save:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()




























def dirMI_function(AllData, filename='', path='', save=False, show=True, s=100,alpha=0.6, scale=0.008, ML=True, AP=True, stats=False, hist=False, interest=None, color='c', binsNumber=30):
    # INITIATION
    stereotaxic_title=[]
    stereotaxic_label=[]
    if ML:
        stereotaxic_title.append('Mediolateral')
        stereotaxic_label.append('ML_pos')
    if AP:
        stereotaxic_title.append('Anteroposterior')
        stereotaxic_label.append('AP_pos')

    preference_quantity = dict() if (stats and hist) else None


    # LOOPS
    for pos_title, posOrientation in zip(stereotaxic_title, stereotaxic_label):

        # IMPORTATION OF DATA
        NtotClust = np.sum([sum(AllData[animal]['good_baseline']) for animal in AllData])
        AllDepth = np.concatenate([AllData[animal]['AllDepth'][AllData[animal]['good_baseline']] for animal in AllData])
        dirMI = np.concatenate([AllData[animal]['dirMI'][AllData[animal]['good_baseline']] for animal in AllData])
        pos = [AllData[animal][posOrientation] for animal in AllData]
        position = np.concatenate([np.random.normal(loc=pos[i], scale=scale, size=(sum(AllData[animal]['good_baseline']), 1)) for animal, i in zip(AllData, range(len(pos)))])


        maximum_value = max([max(AllData[animal]['dirMI'][AllData[animal]['good_baseline']]) for animal in AllData])
        minimum_value = min([min(AllData[animal]['dirMI'][AllData[animal]['good_baseline']]) for animal in AllData])
        nul_value = min([min(AllData[animal]['dirMI'][AllData[animal]['good_baseline']], key=lambda x: abs(x)) for animal in AllData])


        

        # CREATING DATA FOR PLOT
        if stats:
            AllPreferences = np.concatenate([[value for i, value in enumerate(AllData[animal]['preference']) if AllData[animal]['good_baseline'][i]] for animal in AllData])
            # AllPreferences = np.concatenate([[AllData[animal]['preference'][unit] for unit in range(AllData[animal]['Nclust'])] for animal in AllData])

            CW_preferrers_bool = AllPreferences == 'CW'
            CCW_preferrers_bool = AllPreferences == 'CCW'
            Non_preferrers_bool = AllPreferences == None

            
            positionCW = [position[i] for i in range(NtotClust) if CW_preferrers_bool[i]]
            depthCW = [AllDepth[i] for i in range(NtotClust) if CW_preferrers_bool[i]]
            dirMICW = [dirMI[i] for i in range(NtotClust) if CW_preferrers_bool[i]]

            positionCCW = [position[i] for i in range(NtotClust) if CCW_preferrers_bool[i]]
            depthCCW = [AllDepth[i] for i in range(NtotClust) if CCW_preferrers_bool[i]]
            dirMICCW = [dirMI[i] for i in range(NtotClust) if CCW_preferrers_bool[i]]

            positionNone = [position[i] for i in range(NtotClust) if Non_preferrers_bool[i]]
            depthNone = [AllDepth[i] for i in range(NtotClust) if Non_preferrers_bool[i]]
            dirMINone = [dirMI[i] for i in range(NtotClust) if Non_preferrers_bool[i]]

            positionData = [positionCW, positionCCW, positionNone]
            depthData = [depthCW, depthCCW, depthNone]
            dirMIData = [dirMICW, dirMICCW, dirMINone]

            if hist:
                histdata2 = [] ; bins2 = [] ; histdataDist2 = []
                histdata3 = [] ; bins3 = [] ; histdataDist3 = []

                for data in [depthCW, depthCCW]:
                    histdatafoo, binsfoo = np.histogram(data, bins=binsNumber)
                    histdata2.append(histdatafoo) ; bins2.append(binsfoo)
                    histdataDist2.append(histdatafoo / NtotClust * 100)

                for data in [positionCW, positionCCW]:
                    histdatafoo, binsfoo = np.histogram(data, bins=binsNumber)
                    histdata3.append(histdatafoo) ; bins3.append(binsfoo)
                    histdataDist3.append(histdatafoo / NtotClust * 100)

                modulation_label = ['CW-preferring', 'CCW-preferring', 'No preference']
                preference_quantity = [len(positionCW), len(positionCCW), len(positionNone)]
        elif not stats:
            if hist:
                histdata3, bins3 = np.histogram(position, bins=binsNumber)
                histdataDist3 = histdata3 / NtotClust * 100

                histdata2, bins2 = np.histogram(AllDepth, bins=binsNumber)
                histdataDist2 = histdata2 / NtotClust * 100

            positionData = position
            depthData = AllDepth
            dirMIData = dirMI
        ############################################################





        # mettre en évidence les neurones d'intérêt
        if interest!=None:
            wanted_array = []
            for value in interest:
                if value == 'max':
                    wanted_array.append(maximum_value)
                elif value == 'min':
                    wanted_array.append(minimum_value)
                elif value == 'nul':
                    wanted_array.append(nul_value)
                else:
                    wanted_array.append(value)
        ############################################################





        # PLOTTING DATA
        fig = plt.figure(figsize=(13,10))

        ## Creating axis
        if hist:
            gs = GridSpec(nrows=4, ncols=4)
            ax1 = fig.add_subplot(gs[1:4,0:3]) ; # scatter plot on the left
            ax2 = fig.add_subplot(gs[1:4,3], sharey=ax1) ;  # histogram on the right
            ax3 = fig.add_subplot(gs[0,0:3], sharex=ax1) ;  # histogram on the top
            # ax4 = fig.add_subplot(gs[0,3]) if stats else None
        ############################################################




        ## Highlighting units of interest
        if interest!=None:
            for value in wanted_array:
                ax = ax1 if hist else plt
                if type(value) == list:
                    values_of_interest = [dirMI_value for dirMI_value in dirMI if value[0] <= dirMI_value <= value[1]]
                    for value in values_of_interest:
                        unit_index = np.where(dirMI == value)[0]
                        ax.scatter(position[unit_index], AllDepth[unit_index], marker='s', s=s*4, edgecolors='black', facecolors='none', linewidths=2) if len(unit_index) > 0 else None
                else:
                    unit_index = np.where(dirMI == value)[0]
                    if len(unit_index) > 0:
                        ax.scatter(position[unit_index], AllDepth[unit_index], marker='s', s=s*4, edgecolors='black', facecolors='none', linewidths=2)
        ############################################################




        ## Plots of data
        if stats:                    
            if hist:
                for position, depth, transparence, color,  label in zip(positionData, depthData, [alpha, alpha, alpha/3], ['red', 'blue', 'grey'], ['CW-preferring', 'CCW-preferring', 'No preference']):
                    ax1.scatter(position, depth, s=s, alpha=transparence, color=color, label=label)
                ax1.legend()
                for histdata, bins, color in zip(histdataDist2, bins2, ['red', 'blue']):
                    ax2.plot(histdata, bins[:-1], color=color)
                for histdata, bins, color in zip(histdataDist3, bins3, ['red', 'blue']):
                    ax3.plot(bins[:-1], histdata, color=color)
                # ax4.bar(modulation_label, preference_quantity, color=['red', 'blue', 'grey'], alpha=0.6)
            else:
                for position, depth, transparence, color,  label in zip(positionData, depthData, [alpha, alpha, alpha/3], ['red', 'blue', 'grey'], ['CW-preferring', 'CCW-preferring', 'No preference']):
                    plt.scatter(position, depth, s=s, alpha=transparence, color=color, label=label)
                plt.legend()
        else:
            if hist:
                ax1.scatter(positionData, depthData, s=s, alpha=alpha, color=color)
                ax2.plot(histdataDist2, bins2[:-1], color=color)
                ax3.plot(bins3[:-1], histdataDist3, color=color)
            else:
                plt.scatter(positionData, depthData, c=dirMI, cmap='coolwarm', s=s, alpha=alpha, clim=(minimum_value, maximum_value))                
                plt.colorbar(label=f"direction modulation index"+r" : $(n_{CW} - n_{CCW})/(n_{CW} + n_{CCW})$") if ((not hist) and (not stats)) else None


        ## MEP
        xlabel = f"{pos_title} position (mm)"
        ylabel = r"Depth ($\mu$m)"
        histlabel = 'Density (% of all neurons)' if hist else None
        # statsylabel = 'Number of units' if stats else None #ax4
        suptitle = f"Direction preference in {pos_title} axis"

        xfontsize = 12
        yfontsize = 12
        suptitlefontsize = 16
        ############################################################

        ## Aesthetics
        if hist:
            ax1.invert_yaxis()
            ax1.invert_xaxis() if posOrientation == 'ML_pos' else None
            ax1.set_xlabel(xlabel, fontsize=xfontsize)
            ax1.set_ylabel(ylabel, fontsize=yfontsize)

            ax2.set_xlabel(histlabel, fontsize=xfontsize)
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False) 
            ax2.spines['bottom'].set_visible(True) 
            ax2.spines['left'].set_visible(False)
            ax2.tick_params(axis='both', which='both', bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)

            ax3.set_ylabel(histlabel, fontsize=yfontsize)
            ax3.spines['top'].set_visible(False)
            ax3.spines['right'].set_visible(False)
            ax3.spines['bottom'].set_visible(False)
            ax3.spines['left'].set_visible(True)
            ax3.tick_params(axis='both', which='both', bottom=False, top=False, left=True, right=False, labelleft=True, labelbottom=False)

            # # if stats:
            # #     ax4.set_ylabel(statsylabel, fontsize=statsyfontsize)
            # #     ax4.set_xticks(np.arange(len(modulation_label)))
            # #     ax4.set_xticklabels(modulation_label, rotation=20)
            # img = mpimg.imread(r'C:\Users\ayazici\Documents\Analyses\Vestibular_experiments\brain_sagittal.jpg')
            # ax4.imshow(img)
            # ax4.axis('off')
        else:
            plt.gca().invert_yaxis()
            plt.gca().invert_xaxis() if posOrientation == 'ML_pos' else None
            plt.xlabel(xlabel, fontsize=xfontsize)
            plt.ylabel(ylabel, fontsize=yfontsize)
                
        plt.suptitle(suptitle, fontsize=suptitlefontsize)
        ############################################################


        # Saving and showing
        if save:
            if path=='':
                KeyError('You must specify a path to save the figure')
            else:
                os.makedirs(path, exist_ok=True)
                if filename:
                    plt.savefig(os.path.join(path , f"{direction}_modulation_{pos_title}_{filename}.png"))
                else:
                    plt.savefig(os.path.join(path , f"{direction}_modulation_{pos_title}.png"))
        
        plt.show() if show else plt.close()
        ############################################################

        foo=dict()
        for animal in AllData:
            foo[animal] = AllData[animal][posOrientation]
        sorted_keys = [key for key, value in sorted(foo.items(), key=itemgetter(1))] 

        myTable = PrettyTable(["Animal", pos_title]) 

        for animal in sorted_keys:
            myTable.add_row([animal, foo[animal]])

        print(myTable)
        
    if stats and hist:
        plt.figure()
        plt.bar(modulation_label, preference_quantity, color=['red', 'blue', 'grey'])
        plt.ylabel('Number of units')
        plt.title('Preference of units')
        plt.show()






















def vMI_function(AllData,
                s=100, alpha=0.6, color='c', scale=0.008, binsNumber=30,
                ML=True, AP=True, CW=True, CCW=True, stats=False, hist=False, interest=None,
                size=(13,10),
                xfontsize=10, yfontsize=10, suptitlefontsize=14,
                save=False, filename='',
                show=True):
    # INITIATION
    stereotaxic_title, stereotaxic_label, direction_selected = [], [], []

    if ML:
        stereotaxic_title.append('Mediolateral')
        stereotaxic_label.append('ML_pos')
    if AP:
        stereotaxic_title.append('Anteroposterior')
        stereotaxic_label.append('AP_pos')
    if CW:
        direction_selected.append('CW')
    if CCW:
        direction_selected.append('CCW')

    modulation_quantity = dict() if (stats and hist) else None


    # LOOPS
    for pos_title, posOrientation in zip(stereotaxic_title, stereotaxic_label):
        for direction in direction_selected:

            # IMPORTATION OF DATA
            NtotClust = np.sum([sum(AllData[animal]['good_baseline']) for animal in AllData])
            AllDepth = np.concatenate([AllData[animal]['AllDepth'][AllData[animal]['good_baseline']] for animal in AllData])
            vMI = np.concatenate([np.array(AllData[animal]['vMI'][direction])[AllData[animal]['good_baseline']] for animal in AllData])
            pos = [AllData[animal][posOrientation] for animal in AllData]
            position = np.concatenate([np.random.normal(loc=pos[i], scale=scale, size=(sum(AllData[animal]['good_baseline']), 1)) for animal, i in zip(AllData, range(len(pos)))])


            maximum_value = max([max(np.array(AllData[animal]['vMI'][direction])[AllData[animal]['good_baseline']]) for animal in AllData])
            minimum_value = min([min(np.array(AllData[animal]['vMI'][direction])[AllData[animal]['good_baseline']]) for animal in AllData])
            nul_value = min([min(np.array(AllData[animal]['vMI'][direction])[AllData[animal]['good_baseline']], key=lambda x: abs(x)) for animal in AllData])


            

            # CREATING DATA FOR PLOT
            if stats:
                AllHow = np.concatenate([AllData[animal]['modulation']['type'][AllData[animal]['good_baseline']] for animal in AllData])
                AllWho = np.concatenate([AllData[animal]['modulation']['selectivity'][AllData[animal]['good_baseline']] for animal in AllData])

                condition_positive, condition_negative, condition_unmodulated = [], [], []

                if direction == 'CW':
                    for i in range(NtotClust):
                        condition_positive.append((AllWho[i] == 'CW' or AllWho[i] == 'both') and (AllHow[i] == '+' or AllHow[i] == '+/+' or AllHow[i] == '+/-'))
                        condition_negative.append((AllWho[i] == 'CW' or AllWho[i] == 'both') and (AllHow[i] == '-' or AllHow[i] == '-/+' or AllHow[i] == '-/-'))
                        condition_unmodulated.append(AllWho[i] == 'unmodulated' or AllWho[i] == 'CCW')
                elif direction == 'CCW':
                    for i in range(NtotClust):
                        condition_positive.append((AllWho[i] == 'CCW' or AllWho[i] == 'both') and (AllHow[i] == '+' or AllHow[i] == '+/+' or AllHow[i] == '-/+'))
                        condition_negative.append((AllWho[i] == 'CCW' or AllWho[i] == 'both') and (AllHow[i] == '-' or AllHow[i] == '+/-' or AllHow[i] == '-/-'))
                        condition_unmodulated.append(AllWho[i] == 'unmodulated' or AllWho[i] == 'CW')


                positionSigni_P = [position[i] for i in range(NtotClust) if condition_positive[i]]
                depthSigni_P = [AllDepth[i] for i in range(NtotClust) if condition_positive[i]]
                vMISigni_P = [vMI[i] for i in range(NtotClust) if condition_positive[i]]

                positionSigni_N = [position[i] for i in range(NtotClust) if condition_negative[i]]
                depthSigni_N = [AllDepth[i] for i in range(NtotClust) if condition_negative[i]]
                vMISigni_N = [vMI[i] for i in range(NtotClust) if condition_negative[i]]

                positionNot = [position[i] for i in range(NtotClust) if condition_unmodulated[i]]
                depthNot = [AllDepth[i] for i in range(NtotClust) if condition_unmodulated[i]]
                vMINot = [vMI[i] for i in range(NtotClust) if condition_unmodulated[i]]

                positionData = [positionSigni_P, positionSigni_N, positionNot]
                depthData = [depthSigni_P, depthSigni_N, depthNot]
                vMIData = [vMISigni_P, vMISigni_N, vMINot]

                if hist:
                    histdata2 = [] ; bins2 = [] ; histdataDist2 = []
                    histdata3 = [] ; bins3 = [] ; histdataDist3 = []

                    for data in [depthSigni_P, depthSigni_N]:
                        histdatafoo, binsfoo = np.histogram(data, bins=binsNumber)
                        histdata2.append(histdatafoo) ; bins2.append(binsfoo)
                        histdataDist2.append(histdatafoo / NtotClust * 100)

                    for data in [positionSigni_P, positionSigni_N]:
                        histdatafoo, binsfoo = np.histogram(data, bins=binsNumber)
                        histdata3.append(histdatafoo) ; bins3.append(binsfoo)
                        histdataDist3.append(histdatafoo / NtotClust * 100)

                    modulation_label = ['Excitation', 'Suppression', 'No modulation']
                    modulation_quantity[direction] = [len(positionSigni_P), len(positionSigni_N), len(positionNot)]
            elif not stats:
                if hist:
                    histdata3, bins3 = np.histogram(position, bins=binsNumber)
                    histdataDist3 = histdata3 / NtotClust * 100

                    histdata2, bins2 = np.histogram(AllDepth, bins=binsNumber)
                    histdataDist2 = histdata2 / NtotClust * 100

                positionData = position
                depthData = AllDepth
                vMIData = vMI
            ############################################################





            # mettre en évidence les neurones d'intérêt
            if interest!=None:
                wanted_array = []
                for value in interest:
                    if value == 'max':
                        wanted_array.append(maximum_value)
                    elif value == 'min':
                        wanted_array.append(minimum_value)
                    elif value == 'nul':
                        wanted_array.append(nul_value)
                    else:
                        wanted_array.append(value)
            ############################################################





            # PLOTTING DATA
            fig = plt.figure(figsize=size)

            ## Creating axis
            if hist:
                gs = GridSpec(nrows=4, ncols=4)
                ax1 = fig.add_subplot(gs[1:4,0:3]) ; # scatter plot on the left
                ax2 = fig.add_subplot(gs[1:4,3], sharey=ax1) ;  # histogram on the right
                ax3 = fig.add_subplot(gs[0,0:3], sharex=ax1) ;  # histogram on the top
                # ax4 = fig.add_subplot(gs[0,3]) if stats else None
            ############################################################




            ## Highlighting units of interest
            if interest!=None:
                for value in wanted_array:
                    ax = ax1 if hist else plt
                    if type(value) == list:
                        values_of_interest = [vMI_value for vMI_value in vMI if value[0] <= vMI_value <= value[1]]
                        for value in values_of_interest:
                            unit_index = np.where(vMI == value)[0]
                            ax.scatter(position[unit_index], AllDepth[unit_index], marker='s', s=s*4, edgecolors='black', facecolors='none', linewidths=2) if len(unit_index) > 0 else None
                    else:
                        unit_index = np.where(vMI == value)[0]
                        if len(unit_index) > 0:
                            ax.scatter(position[unit_index], AllDepth[unit_index], marker='s', s=s*4, edgecolors='black', facecolors='none', linewidths=2)
            ############################################################




            ## Plots of data
            if stats:                    
                if hist:
                    for position, depth, transparence, color,  label in zip(positionData, depthData, [alpha, alpha, alpha/3], ['red', 'blue', 'grey'], ['Excited units', 'Suppressed units', 'Non-significant modulation']):
                        ax1.scatter(position, depth, s=s, alpha=transparence, color=color, label=label)
                    ax1.legend()
                    for histdata, bins, color in zip(histdataDist2, bins2, ['red', 'blue']):
                        ax2.plot(histdata, bins[:-1], color=color)
                    for histdata, bins, color in zip(histdataDist3, bins3, ['red', 'blue']):
                        ax3.plot(bins[:-1], histdata, color=color)
                    # ax4.bar(modulation_label, modulation_quantity[direction], color=['red', 'blue', 'grey'], alpha=1)
                else:
                    for position, depth, transparence, color,  label in zip(positionData, depthData, [alpha, alpha, alpha/3], ['red', 'blue', 'grey'], ['Excited units', 'Suppressed units', 'Non-significant modulation']):
                        plt.scatter(position, depth, s=s, alpha=transparence, color=color, label=label)
                    plt.legend()
            else:
                if hist:
                    ax1.scatter(positionData, depthData, s=s, alpha=alpha, color=color)
                    ax2.plot(histdataDist2, bins2[:-1], color=color)
                    ax3.plot(bins3[:-1], histdataDist3, color=color)
                else:
                    plt.scatter(positionData, depthData, c=vMI, cmap='coolwarm', s=s, alpha=alpha, clim=(minimum_value, maximum_value))                
                    plt.colorbar(label=f"{direction} modulation index"+r" : $(n_{during} - n_{before})/(n_{during} + n_{before})$") if ((not hist) and (not stats)) else None


            ## MEP
            xlabel = f"{pos_title} position (mm)"
            ylabel = r"Depth ($\mu$m)"
            histlabel = 'Density (%)' if hist else None
            # statsylabel = 'Number of units' if stats else None #(ax4)
            suptitle = f"{direction} modulation in {pos_title} axis"

            xfontsize = xfontsize
            yfontsize = yfontsize
            # statsyfontsize = 12 if stats else None #(ax4)
            suptitlefontsize = suptitlefontsize
            ############################################################

            ## Aesthetics
            if hist:
                ax1.invert_yaxis()
                ax1.invert_xaxis() if posOrientation == 'ML_pos' else None
                ax1.set_xlabel(xlabel, fontsize=xfontsize)
                ax1.set_ylabel(ylabel, fontsize=yfontsize)

                ax2.set_xlabel(histlabel, fontsize=xfontsize)
                ax2.spines['top'].set_visible(False)
                ax2.spines['right'].set_visible(False) 
                ax2.spines['bottom'].set_visible(True) 
                ax2.spines['left'].set_visible(False)
                ax2.tick_params(axis='both', which='both', bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)

                ax3.set_ylabel(histlabel, fontsize=yfontsize)
                ax3.spines['top'].set_visible(False)
                ax3.spines['right'].set_visible(False)
                ax3.spines['bottom'].set_visible(False)
                ax3.spines['left'].set_visible(True)
                ax3.tick_params(axis='both', which='both', bottom=False, top=False, left=True, right=False, labelleft=True, labelbottom=False)

                # if stats:
                #     ax4.set_ylabel(statsylabel, fontsize=statsyfontsize)
                #     ax4.set_xticks(np.arange(len(modulation_label)))
                #     ax4.set_xticklabels(modulation_label, rotation=20)
            else:
                plt.gca().invert_yaxis()
                plt.gca().invert_xaxis() if posOrientation == 'ML_pos' else None
                plt.xlabel(xlabel, fontsize=xfontsize)
                plt.ylabel(ylabel, fontsize=yfontsize)
                    
            plt.suptitle(suptitle, fontsize=suptitlefontsize)
            ############################################################


            # Saving and showing
            if save:
                os.makedirs(os.path.dirname(filename), exist_ok=True)
                if filename=='':
                    plt.savefig(os.path.join(saving_path , f"{direction}_modulation_{pos_title}.png"))
                else:
                    plt.savefig(f"{filename}.png")
            
            plt.show() if show else plt.close()
            ############################################################

            foo=dict()
            for animal in AllData:
                foo[animal] = AllData[animal][posOrientation]
            sorted_keys = [key for key, value in sorted(foo.items(), key=itemgetter(1))] 

            myTable = PrettyTable(["Animal", pos_title]) 

            for animal in sorted_keys:
                myTable.add_row([animal, foo[animal]])

            print(myTable)

        if stats and hist:
            plotdata = pd.DataFrame({'CW':modulation_quantity['CW'], 'CCW':modulation_quantity['CCW']}, index=modulation_label)
            plotdata.plot(kind="bar",figsize=(10, 5), color=['red', 'blue'], rot=25)
            plt.ylabel('Number of units')
            plt.title('Rotation modulation')
            if save:
                os.makedirs(os.path.dirname(filename), exist_ok=True)
                if filename=='':
                    plt.savefig(os.path.join(saving_path , f"modulation_{pos_title}_barplot.png"))
                else:
                    plt.savefig(f"{filename}.png")
            if show:
                plt.show()
            else:
                plt.close()
            print('\n')




























#load_data_boolean = False
#SUA_analysis_boolean = False
#analysis_script_path = ''
#selected_path = ''
#saving_path = ''
#
#def toggle_sua_analysis():
#    global SUA_analysis_boolean
#    if sua_checkbutton_var.get():
#        SUA_analysis_boolean = True
#        script_button.pack()
#        files_button.pack()
#    else:
#        SUA_analysis_boolean = False
#        script_button.pack_forget()
#        files_button.pack_forget()
#
#def load_data():
#    global load_data_boolean
#    if load_checkbutton_var.get():
#        load_data_boolean = True
#    else:
#        load_data_boolean = False
#
#def choose_script_path():
#    global analysis_script_path
#    analysis_script_path = filedialog.askopenfilename(initialdir=r'C:\Users\ayazici\BOUVIER')
#    script_path_label.config(text="Script d'analyse : " + analysis_script_path)
#
#def choose_files_path():
#    global selected_path
#    selected_path = filedialog.askdirectory(initialdir=r'P:\SharedFiles\Abdussamed\Pulvinar_rec_dark_80degs')
#    files_path_label.config(text="Files to analyze : " + selected_path)
#
#def choose_data_path():
#    global saving_path
#    saving_path = filedialog.askdirectory(initialdir=r'C:\Users\ayazici\Documents\Analyses')
#    data_path_label.config(text="Path for output : " + saving_path)
#
#root = Tk()
#root.geometry("700x250")
#
#close_button = Button(root, text="Close", command=root.destroy)
#close_button.pack(side=TOP, anchor='ne', padx=5, pady=5)
#
#sua_checkbutton_var = BooleanVar()
#sua_checkbutton = Checkbutton(root, text="Single-to-Multi Animal Analysis", variable=sua_checkbutton_var, command=toggle_sua_analysis)
#sua_checkbutton.pack()
#
#script_button = Button(root, text="Pipeline for Single-Animal", command=choose_script_path)
#files_button = Button(root, text="Files to analyze", command=choose_files_path)
#
#load_checkbutton_var = BooleanVar()
#load_checkbutton = Checkbutton(root, text="Load Data", variable=load_checkbutton_var, command=load_data)
#load_checkbutton.pack()
#
#path_button = Button(root, text="Path for output", command=choose_data_path)
#path_button.pack()
#
## Ajout des étiquettes pour afficher les chemins sélectionnés
#script_path_label = Label(root, text="Script d'analyse : ")
#script_path_label.pack(side=BOTTOM, anchor='w')
#
#files_path_label = Label(root, text="Files to analyze : ")
#files_path_label.pack(side=BOTTOM, anchor='w')
#
#data_path_label = Label(root, text="Path for output : ")
#data_path_label.pack(side=BOTTOM, anchor='w')
#
#root.mainloop()
#
#
#print("Script Path:", analysis_script_path) if analysis_script_path else None
#print("Files Path:", selected_path) if selected_path else None
#print("Output Path:", saving_path) if saving_path else None
#
#selected_folders = select_folders(selected_path) if selected_path else None
#
#MultiAnimal_function(selected_folders, saving_path, analysis_script_path) if SUA_analysis_boolean else None
#
#if load_data_boolean:
#    AllData = load(saving_path, all=True)
#
#    variables = create_variables(AllData)
#
#    for nom_variable, valeur in variables.items():
#        globals()[nom_variable] = valeur