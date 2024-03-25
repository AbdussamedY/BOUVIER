import scipy.io
import os
import h5py
import numpy as np
from scipy.signal import find_peaks, savgol_filter
import matplotlib.pyplot as plt
# from scipy.stats import sem
from scipy.stats import wilcoxon
import pandas as pd
import json
import jdata as jd
import re
import pickle

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

from matplotlib.gridspec import GridSpec
from PIL import Image

from aquarel import load_theme

from operator import itemgetter

from prettytable import PrettyTable 


# qt for popup window (savable as pdf, svg...), inline for inline plot, notebook for interactive plot, widget for interactive plot
#%matplotlib widget 
#plt.ioff()












# SPLIT PATH FUNCTION

def split_path(path):
    drive, path_without_drive = os.path.splitdrive(path)
    path_segments = []
    while True:
        head, tail = os.path.split(path_without_drive)
        if tail:
            path_segments.insert(0, tail)
            path_without_drive = head
        else:
            if head:
                path_segments.insert(0, head)
            break
    path_segments.insert(0, drive)
    return path_segments











# GENERAL BROWSER FUNCTION

import tkinter as tk
from tkinter import filedialog

def browse(type,extra,multiple=False, initialdir=None):
    if type == 'file':
        if multiple:
            browsed = filedialog.askopenfilenames(title='Select your files '+extra, initialdir=initialdir)
        else:
            browsed = filedialog.askopenfilename(title='Select a file '+extra, initialdir=initialdir)
    elif type == 'folder':
        if multiple:
            dirselect = filedialog.Directory(title='Select your folders '+extra, initialdir=initialdir)
            browsed = []
            while True:
                d = dirselect.show()
                if not d: break
                browsed.append(d)
        else:
            browsed = filedialog.askdirectory(title='Select a folder '+extra, initialdir=initialdir)
    
    return browsed














# MULTI ANIMAL FUNCTION

def MultiAnimal_function(selected_folders, analyse_path, analysis_script):
    # Nombre total de dossiers à analyser
    total_folders = len(selected_folders)

    # tqdm pour afficher une barre de progression
    for i, selected_path in enumerate(tqdm(selected_folders, desc="Analyzing folders", unit="folder")):

        exp_id = "_".join([element for element in split_path(selected_path) if 'animal' in element][0].split('_')[0:2])

        executed_notebook = analysis_script
        result_notebook = os.path.join(analyse_path, f"Analysis_{exp_id}.ipynb")  

        jd.save(selected_path, os.path.join(analyse_path, 'path.json'))

        # Charger le notebook à exécuter
        with open(executed_notebook, 'r', encoding='utf-8') as f:
            notebook = nbformat.read(f, as_version=4)

        # Créer un préprocesseur pour exécuter le notebook
        preprocessor = ExecutePreprocessor(timeout=None)

        # Exécuter le notebook
        preprocessor.preprocess(notebook)

        # Exporter le notebook exécuté
        exporter = NotebookExporter()
        body, resources = exporter.from_notebook_node(notebook)

        # Écrire le notebook exécuté dans un fichier
        with open(result_notebook, 'w', encoding='utf-8') as f:
            f.write(body)

        os.makedirs(os.path.join(analyse_path, 'Path'), exist_ok=True)
        shutil.move(os.path.join(analyse_path, 'path.json'), os.path.join(analyse_path, 'Path', f"{exp_id}.json"))











# vMI FUNCTION

def vMI_function(AllData, filename='', path='', save=False, show=True, s=100,alpha=0.6, scale=0.008, ML=True, AP=True, CW=True, CCW=True, stats=False, hist=False, interest=None, color='c', binsNumber=30):
    stereotaxic_title=[]
    stereotaxic_label=[]
    direction_selected=[]

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


    for pos_title, posOrientation in zip(stereotaxic_title, stereotaxic_label):
        for condition in ['second']:
            for direction in direction_selected:

                NtotClust = np.sum([AllData[animal]['SUA_data']['Nclust'] for animal in AllData])
                AllDepth = np.concatenate([AllData[animal]['MUA_data']['AllDepth'] for animal in AllData])
                vMI = np.concatenate([AllData[animal]['Statistics_data']['vMI'][condition][direction] for animal in AllData])
                pos = [AllData[animal]['informative_data'][posOrientation] for animal in AllData]
                position = np.concatenate([np.random.normal(loc=pos[i], scale=scale, size=(AllData[animal]['SUA_data']['Nclust'], 1)) for animal, i in zip(AllData, range(len(pos)))])


                maximum_value = max([max(AllData[animal]['Statistics_data']['vMI'][condition][direction]) for animal in AllData])
                minimum_value = min([min(AllData[animal]['Statistics_data']['vMI'][condition][direction]) for animal in AllData])
                nul_value = min([min(AllData[animal]['Statistics_data']['vMI'][condition][direction], key=lambda x: abs(x)) for animal in AllData])


                fig = plt.figure(figsize=(17, 9))

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

                    


                if hist:
                    gs = GridSpec(nrows=4, ncols=4)
                    ax1 = fig.add_subplot(gs[1:4,0:3]) ; # scatter plot on the left
                    ax2 = fig.add_subplot(gs[1:4,3], sharey=ax1) ; ax2.set_xlabel('Density') # histogram on the right
                    ax3 = fig.add_subplot(gs[0,0:3], sharex=ax1) ; ax3.set_ylabel('Density') # histogram on the top
                    ax4 = fig.add_subplot(gs[0,3]) if stats else None
                    

                    ax2.spines['top'].set_visible(False) ; ax3.spines['top'].set_visible(False)
                    ax2.spines['right'].set_visible(False) ; ax3.spines['right'].set_visible(False)
                    ax2.spines['bottom'].set_visible(True) ; ax3.spines['bottom'].set_visible(False)
                    ax2.spines['left'].set_visible(False) ; ax3.spines['left'].set_visible(True)
                    ax2.tick_params(axis='both', which='both', bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
                    ax3.tick_params(axis='both', which='both', bottom=False, top=False, left=True, right=False, labelleft=True, labelbottom=False)

                    if not stats:
                        histdata, bins = np.histogram(position, bins=binsNumber)
                        histdataDist = histdata / NtotClust
                        ax3.plot(bins[:-1], histdataDist, color=color)

                        histdata, bins = np.histogram(AllDepth, bins=binsNumber)
                        histdataDist = histdata / NtotClust
                        ax2.plot(histdataDist, bins[:-1], color=color)

                        ax1.scatter(position, AllDepth, s=s, alpha=alpha, color=color)

                    # ax=ax1
                # else:
                    # ax=plt


                # mettre en évidence les neurones d'intérêt
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




                if stats:                     
                    AllHow = np.concatenate([[AllData[animal]['Statistics_data']['modulation']['second'][unit]['type'] for unit in range(AllData[animal]['SUA_data']['Nclust'])] for animal in AllData])
                    AllWho = np.concatenate([[AllData[animal]['Statistics_data']['modulation']['second'][unit]['selectivity'] for unit in range(AllData[animal]['SUA_data']['Nclust'])] for animal in AllData])

                    condition_positive = []
                    condition_negative = []
                    condition_unmodulated = []

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
                    
                    if hist:
                        ax2.spines['top'].set_visible(False) ; ax3.spines['top'].set_visible(False)
                        ax2.spines['right'].set_visible(False) ; ax3.spines['right'].set_visible(False)
                        ax2.spines['bottom'].set_visible(True) ; ax3.spines['bottom'].set_visible(False)
                        ax2.spines['left'].set_visible(False) ; ax3.spines['left'].set_visible(True)
                        ax2.tick_params(axis='both', which='both', bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
                        ax3.tick_params(axis='both', which='both', bottom=False, top=False, left=True, right=False, labelleft=True, labelbottom=False)

                        histdata1, bins1 = np.histogram(positionSigni_P, bins=binsNumber)
                        histdataDist1 = histdata1 / NtotClust
                        histdata2, bins2 = np.histogram(positionSigni_N, bins=binsNumber)
                        histdataDist2 = histdata2 / NtotClust
                        ax3.plot(bins1[:-1], histdataDist1, color='red')
                        ax3.plot(bins2[:-1], histdataDist2, color='blue')

                        histdata1, bins1 = np.histogram(depthSigni_P, bins=binsNumber)
                        histdataDist1 = histdata1 / NtotClust
                        histdata2, bins2 = np.histogram(depthSigni_N, bins=binsNumber)
                        histdataDist2 = histdata2 / NtotClust
                        ax2.plot(histdataDist1, bins1[:-1], color='red')
                        ax2.plot(histdataDist2, bins2[:-1], color='blue')

                        modulation_label = ['Excited units', 'Suppressed units', 'Non-significant modulation']
                        modulation_quantity = [len(positionSigni_P), len(positionSigni_N), len(positionNot)]
                        ax4.bar(modulation_label, modulation_quantity, color=['red', 'blue', 'grey'], alpha=0.6)


                        ax=ax1
                    else:
                        ax=plt
                    ax.scatter(positionSigni_P, depthSigni_P, s=s, alpha=alpha, color='red', label='Excited units')
                    ax.scatter(positionSigni_N, depthSigni_N, s=s, alpha=alpha, color='blue', label='Suppressed units')
                    ax.scatter(positionNot, depthNot, s=s, alpha=alpha/3, color='grey', label='Non-significant modulation')
                    ax.legend()
                    
                if not hist and not stats:
                    plt.scatter(position, AllDepth, c=vMI, cmap='coolwarm', s=s, alpha=alpha, clim=(minimum_value, maximum_value))
                
                plt.colorbar(label=f"{direction} modulation index"+r" : $\frac{n_{during} - n_{before}}{n_{during} + n_{before}}$") if ((not hist) and (not stats)) else None

                plt.gca().invert_yaxis() if not hist else ax1.invert_yaxis()
                (plt.gca().invert_xaxis() if posOrientation == 'ML_pos' else None) if not hist else (ax1.invert_xaxis() if posOrientation == 'ML_pos' else None)

                xlabel = f"{pos_title} position (mm)" ; plt.xlabel(xlabel, fontsize=14) if not hist else ax1.set_xlabel(xlabel)
                ylabel = r"Depth ($\mu$m)" ; plt.ylabel(ylabel, fontsize=14) if not hist else ax1.set_ylabel(ylabel)
                plt.suptitle(f"{direction} modulation of units in {pos_title} axis", fontsize=16)
                
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

                foo=dict()
                for animal in AllData:
                    foo[animal] = AllData[animal]['informative_data'][posOrientation]
                sorted_keys = [key for key, value in sorted(foo.items(), key=itemgetter(1))] 

                myTable = PrettyTable(["Animal", pos_title]) 

                for animal in sorted_keys:
                    myTable.add_row([animal, foo[animal]])

                print(myTable)
        
















def dirMI_function(AllData, path='', filename='', save=False, show=True, s=100,alpha=0.6, scale=0.008, ML=True, AP=True, stats=False, hist=False, interest=None, color='c', binsNumber=30, All=False):
    stereotaxic_title=[]
    stereotaxic_label=[]

    if ML:
        stereotaxic_title.append('Mediolateral')
        stereotaxic_label.append('ML_pos')
    if AP:
        stereotaxic_title.append('Anteroposterior')
        stereotaxic_label.append('AP_pos')


    for pos_title, posOrientation in zip(stereotaxic_title, stereotaxic_label):
        for condition in ['second']:

            NtotClust = np.sum([AllData[animal]['SUA_data']['Nclust'] for animal in AllData])
            AllDepth = np.concatenate([AllData[animal]['MUA_data']['AllDepth'] for animal in AllData])
            dirMI = np.concatenate([AllData[animal]['Statistics_data']['dirMI'][condition] for animal in AllData])
            pos = [AllData[animal]['informative_data'][posOrientation] for animal in AllData]
            position = np.concatenate([np.random.normal(loc=pos[i], scale=scale, size=(AllData[animal]['SUA_data']['Nclust'], 1)) for animal, i in zip(AllData, range(len(pos)))])


            maximum_value = max([max(AllData[animal]['Statistics_data']['dirMI'][condition]) for animal in AllData])
            minimum_value = min([min(AllData[animal]['Statistics_data']['dirMI'][condition]) for animal in AllData])
            nul_value = min([min(AllData[animal]['Statistics_data']['dirMI'][condition], key=lambda x: abs(x)) for animal in AllData])


            fig = plt.figure(figsize=(17, 9))

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

                


            if hist:
                gs = GridSpec(nrows=4, ncols=4)
                ax1 = fig.add_subplot(gs[1:4,0:3]) ; # scatter plot on the left
                ax2 = fig.add_subplot(gs[1:4,3], sharey=ax1) ; ax2.set_xlabel('Density') # histogram on the right
                ax3 = fig.add_subplot(gs[0,0:3], sharex=ax1) ; ax3.set_ylabel('Density') # histogram on the top
                # ax4 = fig.add_subplot(gs[0,3]) ; ax4.set_axis_off()# for legend
                

                ax2.spines['top'].set_visible(False) ; ax3.spines['top'].set_visible(False)
                ax2.spines['right'].set_visible(False) ; ax3.spines['right'].set_visible(False)
                ax2.spines['bottom'].set_visible(True) ; ax3.spines['bottom'].set_visible(False)
                ax2.spines['left'].set_visible(False) ; ax3.spines['left'].set_visible(True)
                ax2.tick_params(axis='both', which='both', bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
                ax3.tick_params(axis='both', which='both', bottom=False, top=False, left=True, right=False, labelleft=True, labelbottom=False)

                if not stats:
                    histdata, bins = np.histogram(position, bins=binsNumber)
                    histdataDist = histdata / NtotClust
                    ax3.plot(bins[:-1], histdataDist, color=color)

                    histdata, bins = np.histogram(AllDepth, bins=binsNumber)
                    histdataDist = histdata / NtotClust
                    ax2.plot(histdataDist, bins[:-1], color=color)

                    ax1.scatter(position, AllDepth, s=s, alpha=alpha, color=color)


            # mettre en évidence les neurones d'intérêt
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




            if stats:                     
                AllPreferences = np.concatenate([[AllData[animal]['Statistics_data']['preference']['second'][unit] for unit in range(AllData[animal]['SUA_data']['Nclust'])] for animal in AllData])

                CW_preferrers_bool = AllPreferences == 'CW'
                CCW_preferrers_bool = AllPreferences == 'CCW'
                Non_preferrers_bool = AllPreferences == 'None'

                
                positionCW = [position[i] for i in range(NtotClust) if CW_preferrers_bool[i]]
                depthCW = [AllDepth[i] for i in range(NtotClust) if CW_preferrers_bool[i]]
                dirMICW = [dirMI[i] for i in range(NtotClust) if CW_preferrers_bool[i]]

                positionCCW = [position[i] for i in range(NtotClust) if CCW_preferrers_bool[i]]
                depthCCW = [AllDepth[i] for i in range(NtotClust) if CCW_preferrers_bool[i]]
                dirMICCW = [dirMI[i] for i in range(NtotClust) if CCW_preferrers_bool[i]]

                positionNone = [position[i] for i in range(NtotClust) if Non_preferrers_bool[i]]
                depthNone = [AllDepth[i] for i in range(NtotClust) if Non_preferrers_bool[i]]
                dirMINone = [dirMI[i] for i in range(NtotClust) if Non_preferrers_bool[i]]
                
                if hist:
                    ax2.spines['top'].set_visible(False) ; ax3.spines['top'].set_visible(False)
                    ax2.spines['right'].set_visible(False) ; ax3.spines['right'].set_visible(False)
                    ax2.spines['bottom'].set_visible(True) ; ax3.spines['bottom'].set_visible(False)
                    ax2.spines['left'].set_visible(False) ; ax3.spines['left'].set_visible(True)
                    ax2.tick_params(axis='both', which='both', bottom=True, top=False, left=False, right=False, labelleft=False, labelbottom=True)
                    ax3.tick_params(axis='both', which='both', bottom=False, top=False, left=True, right=False, labelleft=True, labelbottom=False)

                    histdata1, bins1 = np.histogram(positionCW, bins=binsNumber)
                    histdataDist1 = histdata1 / NtotClust
                    histdata2, bins2 = np.histogram(positionCCW, bins=binsNumber)
                    histdataDist2 = histdata2 / NtotClust
                    ax3.plot(bins1[:-1], histdataDist1, color='red')
                    ax3.plot(bins2[:-1], histdataDist2, color='blue')

                    if All:
                        histdata3, bins3 = np.histogram(positionNone, bins=binsNumber)
                        histdataDist3 = histdata3 / NtotClust
                        ax3.plot(bins3[:-1], histdataDist3, color='grey', alpha=0.6)

                    histdata1, bins1 = np.histogram(depthCW, bins=binsNumber)
                    histdataDist1 = histdata1 / NtotClust
                    histdata2, bins2 = np.histogram(depthCCW, bins=binsNumber)
                    histdataDist2 = histdata2 / NtotClust
                    ax2.plot(histdataDist1, bins1[:-1], color='red')
                    ax2.plot(histdataDist2, bins2[:-1], color='blue')

                    if All:
                        histdata3, bins3 = np.histogram(depthNone, bins=binsNumber)
                        histdataDist3 = histdata3 / NtotClust
                        ax2.plot(histdataDist3, bins3[:-1], color='grey', alpha=0.6)

                    ax=ax1
                else:
                    ax=plt
                ax.scatter(positionCW, depthCW, s=s, alpha=alpha, color='red', label='CW-preferring units')
                ax.scatter(positionCCW, depthCCW, s=s, alpha=alpha, color='blue', label='CCW-preferring units')
                ax.scatter(positionNone, depthNone, s=s, alpha=alpha/3, color='grey', label='No directional preference')
                ax.legend()
                
            if not hist and not stats:
                plt.scatter(position, AllDepth, c=dirMI, cmap='coolwarm', s=s, alpha=alpha, clim=(minimum_value, maximum_value))
            
            plt.colorbar(label=r"Direction Modulation Index : $\frac{CW - CCW}{CW + CCW}$") if ((not hist) and (not stats)) else None

            plt.gca().invert_yaxis() if not hist else ax1.invert_yaxis()
            (plt.gca().invert_xaxis() if posOrientation == 'ML_pos' else None) if not hist else (ax1.invert_xaxis() if posOrientation == 'ML_pos' else None)

            xlabel = f"{pos_title} position (mm)" ; plt.xlabel(xlabel, fontsize=14) if not hist else ax1.set_xlabel(xlabel)
            ylabel = r"Depth ($\mu$m)" ; plt.ylabel(ylabel, fontsize=14) if not hist else ax1.set_ylabel(ylabel)
            plt.suptitle(f"CW vs CCW preference in {pos_title} axis", fontsize=16)
            

            if save:
                if path=='':
                    KeyError('You must specify a path to save the figure')
                else:
                    os.makedirs(path, exist_ok=True)
                    if filename:
                        plt.savefig(os.path.join(path , f"Direction_preference_{pos_title}_{filename}.png"))
                    else:
                        plt.savefig(os.path.join(path , f"Direction_preference_{pos_title}.png"))


            plt.show() if show else plt.close()

            foo=dict()
            for animal in AllData:
                foo[animal] = AllData[animal]['informative_data'][posOrientation]
            sorted_keys = [key for key, value in sorted(foo.items(), key=itemgetter(1))] 

            myTable = PrettyTable(["Animal", pos_title]) 

            for animal in sorted_keys:
                myTable.add_row([animal, foo[animal]])

            print(myTable)

















































# LOAD FUNCTION

def load(analyse_path, all):
    AllData = {}

    data_folder = os.path.join(analyse_path, 'Data')
    


    if all:
        files_to_load = sorted(os.listdir(data_folder))
        
        for i, file in enumerate(tqdm(files_to_load, desc="Data loading", unit="folder")):
            # AllData[file.split('_')[1]] = jd.load(os.path.join(data_folder, file))
            with open(os.path.join(data_folder, file), 'rb') as FILE_READER:
                AllData['_'.join(file.split('_')[0:2])] = pickle.load(FILE_READER)
    else:
        files_to_load = browse('file',multiple=True, extra='to load')

        for i, file in enumerate(tqdm(files_to_load, desc="Data loading", unit="folder")):
            # AllData[file.split('_')[1]] = jd.load(file)
            with open(file, 'rb') as FILE_READER:
                AllData['_'.join(file.split('_')[0:2])] = pickle.load(FILE_READER)

    return AllData










def create_variables(AllData):
    variables = {}


    for animal in AllData:
        for keys in AllData[animal]:
            for keys2 in AllData[animal][keys]:
                if keys2 not in variables:
                    variables[keys2] = {}
                variables[keys2][animal] = AllData[animal][keys][keys2]

    return variables









def select_folders(selected_path):    
    selected_folders = []

    for chemin in os.listdir(selected_path):
        foo = []
        for element in os.listdir(os.path.join(selected_path, chemin)):
            # foo = []
            pattern = r'\d{2}[a-zA-Z]\d{2}[a-zA-Z]\d{1}[a-zA-Z]\d{1}_\d{6}_\d{6}'
            condition1 = os.path.isdir(os.path.join(selected_path,chemin, element))
            condition2 = re.match(pattern, element)
            foo.append(element) if (condition1) and (condition2) else None
        
        if len(foo)>0:
            the_good_one = foo[0]
            
            if len(foo)>1:
                for foo_element in foo:
                    the_good_one = foo_element if len(foo_element) < len(the_good_one) else the_good_one

            selected_folders.append(os.path.join(selected_path, chemin, the_good_one))
            selected_folders = sorted(selected_folders)

    print('Selected folders:')
    for folder in selected_folders:
        print(folder)

    return selected_folders










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













def getPSTH(Zscore, SEM, edges, ax='', xlabel='', ylabel='', title='', show=True):
    if ax:
        ax.plot(edges[:-1], Zscore)
        ax.fill_between(edges[:-1], Zscore-SEM, Zscore+SEM, alpha=0.5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    else:
        # with load_theme("arctic_light"):
        plt.plot(edges[:-1], Zscore)
        plt.fill_between(edges[:-1], Zscore-SEM, Zscore+SEM, alpha=0.5)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)

        plt.show() if show else None

def PSTH(StudiedSpikeTimes, timeBef=1, timeAft=5, binResolution=0.03, xlabel='', ylabel='', title='', ax='', show=True):
    local_trial_number = len(StudiedSpikeTimes)

    spike_number_per_trial = [[] for _ in range(local_trial_number)]
    edges = []
    unitary_firing_rate = [[] for _ in range(local_trial_number)]

    for trial in range(local_trial_number):
        spike_number_per_trial[trial], edges = np.histogram(StudiedSpikeTimes[trial], bins=np.arange(-timeBef, timeAft + binResolution, binResolution))

    frequency_per_trial = [[spike_number_per_trial[trial][bin]/binResolution for bin in range(len(edges)-1)] for trial in range(local_trial_number)]
    mean_frequency = [np.mean([frequency_per_trial[trial][bin] for trial in range(local_trial_number)]) for bin in range(len(edges)-1)]

    Zscore = (mean_frequency - np.mean(mean_frequency)) / np.std(mean_frequency) if np.std(mean_frequency) != 0 else np.zeros(len(mean_frequency))
    Zunitary = (frequency_per_trial - np.mean(mean_frequency)) / np.std(mean_frequency) if np.std(mean_frequency) != 0 else np.zeros(len(frequency_per_trial))
    SEM = np.std(Zunitary)/np.sqrt(len(Zunitary)) if np.std(mean_frequency) != 0 else np.zeros(len(mean_frequency))

    getPSTH(Zscore, SEM, edges, xlabel=xlabel, ylabel=ylabel, title=title, ax=ax, show=show)


























































































































# def vMI_function(AllData, analyse_path, save = False,s=100,alpha=0.6, show=True, scale=0.008, ML=True, AP=True, CW=True, CCW=True, seuil_min=None, seuil_max=None, outcolor='green', outshow=True, outscoeff=0.1, incolor=None, colormap=True, interest=None, hist=False):
#     stereotaxic_title=[]
#     stereotaxic_label=[]
#     direction_selected=[]

#     if hist:
#         outscoeff = 1
#         if seuil_min!=None:
#             outcolor='blue'
#             incolor='red'
#         elif seuil_max!=None:
#             outcolor='red'
#             incolor='blue'


#     if ML:
#         stereotaxic_title.append('Mediolateral')
#         stereotaxic_label.append('ML_pos')

#     if AP:
#         stereotaxic_title.append('Anteroposterior')
#         stereotaxic_label.append('AP_pos')

#     if CW:
#         direction_selected.append('CW')

#     if CCW:
#         direction_selected.append('CCW')



#     for pos_title, posOrientation in zip(stereotaxic_title, stereotaxic_label):
#         for condition in ['second']:
#             for direction in direction_selected:
#                 AllPositionHist = []
#                 AllBooleanHist = [] if seuil_min!=None or seuil_max!=None else None
#                 with load_theme("arctic_light"):
#                     fig = plt.figure(figsize=(17, 9))
#                     plt.rcParams['figure.constrained_layout.use'] = False
#                     if hist:
#                         gs = GridSpec(nrows=4, ncols=4)
#                         ax1 = fig.add_subplot(gs[1:4,0:3]) # scatter plot on the left
#                         ax2 = fig.add_subplot(gs[1:4,3]) # histogram on the right
#                         ax3 = fig.add_subplot(gs[0,0:3]) # histogram on the top
#                         ax4 = fig.add_subplot(gs[0,3]) # for legend

#                     for animal in AllData:
#                         Nclust = AllData[animal]['SUA_data']['Nclust']
#                         AllDepth = AllData[animal]['MUA_data']['AllDepth']
#                         vMI = AllData[animal]['Statistics_data']['vMI'][condition][direction]
#                         pos = AllData[animal]['informative_data'][posOrientation]
#                         position = np.random.normal(loc=pos, scale=scale, size=(Nclust, 1))
#                         AllPositionHist.extend(position)
                       
#                         maximum_value = max([max(AllData[animal]['Statistics_data']['vMI'][condition][direction]) for animal in AllData])
#                         minimum_value = min([min(AllData[animal]['Statistics_data']['vMI'][condition][direction]) for animal in AllData])
#                         nul_value = min([min(AllData[animal]['Statistics_data']['vMI'][condition][direction], key=lambda x: abs(x)) for animal in AllData])


#                         # mettre en évidence les neurones d'intérêt
#                         if interest!=None:
#                             wanted_array = []
#                             for value in interest:
#                                 if value == 'max':
#                                     wanted_array.append(maximum_value)
#                                 elif value == 'min':
#                                     wanted_array.append(minimum_value)
#                                 elif value == 'nul':
#                                     wanted_array.append(nul_value)
#                                 else:
#                                     wanted_array.append(value)

#                             for value in wanted_array:
#                                 ax = ax1 if hist else plt
#                                 if type(value) == list:
#                                     values_of_interest = [vMI_value for vMI_value in vMI if value[0] <= vMI_value <= value[1]]
#                                     for value in values_of_interest:
#                                         unit_index = np.where(vMI == value)[0]
#                                         if len(unit_index) > 0:
#                                             ax.scatter(position[unit_index], AllDepth[unit_index], marker='s', s=s*4, edgecolors='black', facecolors='none', linewidths=2)
#                                 else:
#                                     unit_index = np.where(vMI == value)[0]
#                                     if len(unit_index) > 0:
#                                         ax.scatter(position[unit_index], AllDepth[unit_index], marker='s', s=s*4, edgecolors='black', facecolors='none', linewidths=2)


#                         # Scatter des neurones en fonction de leur position
#                         ax = ax1 if hist else plt
#                         if seuil_min!=None or seuil_max!=None:
#                             seuil_pos = np.where(np.array(vMI)>=seuil_min)[0] if seuil_min!=None else np.where(np.array(vMI)<seuil_max)[0]
#                             not_seuil_pos = np.where(np.array(vMI)<seuil_min)[0] if seuil_min!=None else np.where(np.array(vMI)>=seuil_max)[0]

#                             AllBooleanHist.extend([i in seuil_pos for i in range(len(vMI))])

#                             if outshow:
#                                 ax.scatter([position[i] for i in not_seuil_pos], [AllDepth[i] for i in not_seuil_pos], color=outcolor, s=s, alpha=alpha*outscoeff)
                            
#                             if incolor==None:
#                                 ax.scatter([position[i] for i in seuil_pos], [AllDepth[i] for i in seuil_pos], c=[vMI[i] for i in seuil_pos], cmap='coolwarm', s=s, alpha=alpha, clim=(minimum_value, maximum_value))
#                             else:
#                                 ax.scatter([position[i] for i in seuil_pos], [AllDepth[i] for i in seuil_pos], s=s, alpha=alpha, color=incolor)
#                         else:
#                             ax.scatter(position, AllDepth, c=vMI, cmap='coolwarm', s=s, alpha=alpha, clim=(minimum_value, maximum_value))

                    
#                     if hist:
#                         AllDepthHist = np.concatenate([AllData[animal]['MUA_data']['AllDepth'] for animal in AllData])
#                         NtotClust = np.sum([AllData[animal]['SUA_data']['Nclust'] for animal in AllData])

#                         if seuil_min!=None or seuil_max!=None:
#                             histdata1, bins1 = np.histogram([v for v, b in zip(AllDepthHist, AllBooleanHist) if b], bins=50)
#                             histdataDist1 = histdata1 / NtotClust
#                             histdata2, bins2 = np.histogram([v for v, b in zip(AllDepthHist, AllBooleanHist) if not b], bins=50)
#                             histdataDist2 = histdata2 / NtotClust
#                             ax2.plot(histdataDist1, bins1[:-1], color=incolor if incolor!=None else 'green')
#                             ax2.plot(histdataDist2, bins2[:-1], color=outcolor if outcolor!=None else 'green')

#                             histdata1, bins1 = np.histogram([v for v, b in zip(AllPositionHist, AllBooleanHist) if b], bins=50)
#                             histdataDist1 = histdata1 / NtotClust
#                             histdata2, bins2 = np.histogram([v for v, b in zip(AllPositionHist, AllBooleanHist) if not b], bins=50)
#                             histdataDist2 = histdata2 / NtotClust
#                             ax3.plot(bins1[:-1], histdataDist1, color=incolor if incolor!=None else 'green')
#                             ax3.plot(bins2[:-1], histdataDist2, color=outcolor if outcolor!=None else 'green')
#                         else:
#                             histdata, bins = np.histogram(AllDepthHist, bins=50)
#                             histdataDist = histdata / NtotClust
#                             ax2.plot(histdataDist, bins[:-1])
#                             histdata, bins = np.histogram(AllPositionHist, bins=100)
#                             histdataDist = histdata / NtotClust
#                             ax3.plot(bins[:-1], histdataDist)
                        


                    






#                         if seuil_min!=None or seuil_max!=None:
#                             ax4.plot([0, 0.5], [0, 0], color=incolor, linewidth=1)
#                             seuil = seuil_min if seuil_min!=None else seuil_max
#                             ax4.text(0.6, 0, fr'$vMI>${seuil}', ha='left', va='center', fontsize=12, color='k')
#                             ax4.plot([0, 0.5], [-0.5, -0.5], color=outcolor, linewidth=1)
#                             ax4.text(0.6, -0.5, fr'$vMI<${seuil}', ha='left', va='center', fontsize=12, color='k')
#                             ax4.set_xlim(-1,3)
#                             ax4.set_ylim(-2,2)
#                         ax4.set_axis_off()

#                         for ax in [ax2, ax3]:
#                             ax.spines['top'].set_visible(False)
#                             ax.spines['right'].set_visible(False)
#                             ax.spines['bottom'].set_visible(False)
#                             ax.spines['left'].set_visible(False)
#                         ax2.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelleft=False, labelbottom=True)
#                         ax3.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelleft=True, labelbottom=False)
                    
#                         ax1.invert_yaxis()
#                         ax2.invert_yaxis()
#                         ax3.invert_xaxis() if posOrientation == 'ML_pos' else None
#                         ax1.invert_xaxis() if posOrientation == 'ML_pos' else None

#                         ax1.set_xlabel(f"{pos_title} position (mm)", fontsize=14)
#                         ax1.set_ylabel("Depth (µm)", fontsize=14)

#                         ax2.set_xlabel("Density", fontsize=14)
#                         ax3.set_ylabel("Density", fontsize=14)

#                         fig.suptitle(f"{direction} modulation of units in {pos_title} axis", fontsize=16)

#                     else:
#                         plt.gca().invert_yaxis()
#                         plt.gca().invert_xaxis() if posOrientation == 'ML_pos' else None

#                         plt.xlabel(f"{pos_title} position (mm)", fontsize=14)
#                         plt.ylabel("Depth (µm)", fontsize=14)
#                         plt.title(f"{direction} modulation of units in {pos_title} axis", fontsize=16)

#                 # if incolor==None:
#                 #     ax = ax1 if hist else plt
#                 #     try:
#                 #         ax.colorbar(label=f"{direction} modulation index"+r" : $\frac{n_{during} - n_{before}}{n_{during} + n_{before}}$")
#                 #     except:
#                 #         print('Colorbar with histogram ? Pourquoi faire ?')

#                 if not hist:
#                     plt.colorbar(label=f"{direction} modulation index"+r" : $\frac{n_{during} - n_{before}}{n_{during} + n_{before}}$")
                

#                 if save:
#                     direction_modulation_folder = os.path.join(analyse_path, 'Direction_modulation')
#                     os.makedirs(direction_modulation_folder, exist_ok=True)
#                     if hist:
#                         plt.savefig(os.path.join(direction_modulation_folder , f"{direction}_modulation_{pos_title}_hist.png"))
#                     else:
#                         plt.savefig(os.path.join(direction_modulation_folder , f"{direction}_modulation_{pos_title}.png"))
#                 if show:
#                     plt.show()
#                 else:
#                     plt.close()
                
#                 foo=dict()
#                 for animal in AllData:
#                     foo[animal] = AllData[animal]['informative_data'][posOrientation]
#                 sorted_keys = [key for key, value in sorted(foo.items(), key=itemgetter(1))] 
 
#                 myTable = PrettyTable(["Animal", pos_title]) 

#                 for animal in sorted_keys:
#                     myTable.add_row([animal, foo[animal]])

#                 print(myTable)

#     # if interest!=None:
#     #     MyTable = PrettyTable(["Wanted vMI", "Value"])

#     #     for vMI_of_interest_element, wanted_array_element in zip(interest, wanted_array):
#     #         MyTable.add_row([vMI_of_interest_element, wanted_array_element])
#     #     print(MyTable)
    

























# def vMI_function(AllData, analyse_path, save = False,s=35,alpha=0.6, show=True, scale=0.008, ML=True, AP=True, CW=True, CCW=True, seuil=None):
#     foo=[]
#     foo2=[]
#     foo3=[]

#     if ML:
#         foo.append('Mediolateral')
#         foo2.append('ML_pos')
#     if AP:
#         foo.append('Anteroposterior')
#         foo2.append('AP_pos')
#     if CW:
#         foo3.append('CW')
#     if CCW:
#         foo3.append('CCW')
#     for pos_title, posOrientation in zip(foo, foo2):
#         for condition in ['second']:
#             for direction in foo3:

#                 with load_theme("arctic_light"):
#                     plt.figure(figsize=(17, 5))

#                     for animal in AllData:
#                         Nclust = AllData[animal]['SUA_data']['Nclust']
#                         AllDepth = AllData[animal]['MUA_data']['AllDepth']
#                         vMI = AllData[animal]['Statistics_data']['vMI']
#                         pos = AllData[animal]['informative_data'][posOrientation]
#                         position = np.random.normal(loc=pos, scale=scale, size=(Nclust, 1))
                        
#                         if seuil:
#                             seuil_pos = np.where(np.array(vMI)>=seuil)[0]
#                             not_seuil_pos = np.where(np.array(vMI)<seuil)[0]
#                             plt.scatter([position[i] for i in seuil_pos], [AllDepth[i] for i in seuil_pos], c=[vMI[condition][direction][i] for i in seuil_pos], cmap='coolwarm', s=s, alpha=alpha)
#                             plt.scatter([position[i] for i in not_seuil_pos], [AllDepth[i] for i in not_seuil_pos], color='gray', s=s, alpha=alpha)
#                         else:
#                             plt.scatter(position, AllDepth, c=vMI[condition][direction], cmap='coolwarm', s=s, alpha=alpha)
#                         # plt.scatter(np.random.normal(loc=pos, scale=1, size=(Nclust, 1)), AllDepth, c=vMI[condition][direction], cmap='coolwarm', s=s, alpha=alpha)

#                     plt.gca().invert_yaxis()
#                     if posOrientation == 'ML_pos':
#                         plt.gca().invert_xaxis()

#                     plt.xlabel(f"{pos_title} position (mm)", fontsize=14)
#                     plt.ylabel("Depth (µm)", fontsize=14)
#                     plt.title(f"{direction} modulation of units in {pos_title} axis", fontsize=16)

#                 plt.colorbar(label=f"{direction} modulation index")

#                 if save:
#                     direction_modulation_folder = os.path.join(analyse_path, 'Direction_modulation')
#                     os.makedirs(direction_modulation_folder, exist_ok=True)
#                     plt.savefig(os.path.join(direction_modulation_folder , f"{direction}_modulation_{pos_title}.png"))
#                 if show:
#                     plt.show()
#                 else:
#                     plt.close()
#                 print('\n')