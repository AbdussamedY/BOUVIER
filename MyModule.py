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

    Zscore = (mean_frequency - np.mean(mean_frequency)) / np.std(mean_frequency) if np.std(mean_frequency) != 0 else np.zeros(len(mean_frequency))
    Zscore[-1]=Zscore[-2]
    Zunitary = (frequency_per_trial - np.mean(mean_frequency)) / np.std(mean_frequency) if np.std(mean_frequency) != 0 else np.zeros(len(frequency_per_trial))
    SEM = np.std(Zunitary)/np.sqrt(len(Zunitary)) if np.std(mean_frequency) != 0 else np.zeros(len(mean_frequency))

    return edges, Zscore, SEM
















def plotPSTH(AllData, animal, condition, unit, plotvelocity=False, color='k',shadedcolor='c',binResolution = 0.03,xlabel=True,ylabel=True, title='', save=False, filename='PSTH.png', show=True, velocitycolor='red', velocityalpha=0.13, extra=None, xlim=None, ylim=None, smooth=True):

    if type(unit)==list:
        StudiedSpikeTimes = np.concatenate([AllData[animal]['SpikeTimes'][condition][unit_i] for unit_i in unit])
    else:
        StudiedSpikeTimes = AllData[animal]['SpikeTimes'][condition][unit]

    timeBef, timeAft = AllData[animal]['timeBef'], AllData[animal]['timeAft']
    if plotvelocity:
        rotationSpeed = AllData[animal]['rotationSpeed']
        duration = AllData[animal]['duration']
        MeanRotation = AllData[animal]['MeanRotation'][condition]

    edges, Zscore, SEM = getPSTHparameters(StudiedSpikeTimes, duration, binResolution)

    if smooth:
        Zscore = savgol_filter(Zscore, 9, 3)


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

    if extra is not None:
        if type(extra)==list:
            for ex in extra:
                ex
        else:
            extra
            
    if save:
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
    if show:
        plt.show()



















def plotRaster(AllData, animal, condition, unit, color='black', xlabel='Time (s)', ylabel='# Trial', extra=None, title='', plotvelocity=False, velocitycolor='k', velocityalpha=0.13, save=False, filename='Raster.png', show=True, xlim=None, ylim=None, psth=False, smooth=True, binResolution=0.03, psthcolor='b', shadedcolor='b'):
    if type(unit)==list:
        spikeTimesObject = np.concatenate([AllData[animal]['SpikeTimes'][condition][unit_i] for unit_i in unit])
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