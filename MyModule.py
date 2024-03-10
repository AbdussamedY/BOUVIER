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

def vMI_function(AllData, analyse_path, save = False,s=35,alpha=0.6, show=True, scale=0.008):
    for pos_title, posOrientation in zip(['Mediolateral', 'Anteroposterior'], ['ML_pos', 'AP_pos']):
        for condition in ['second']:
            for direction in ['CW','CCW']:

                with load_theme("arctic_light"):
                    plt.figure(figsize=(17, 5))

                    for animal in AllData:
                        pos = AllData[animal]['informative_data'][posOrientation]
                        Nclust = AllData[animal]['SUA_data']['Nclust']
                        AllDepth = AllData[animal]['MUA_data']['AllDepth']
                        vMI = AllData[animal]['Statistics_data']['vMI']

                        plt.scatter(np.random.normal(loc=pos, scale=scale, size=(Nclust, 1)), AllDepth, c=vMI[condition][direction], cmap='coolwarm', s=s, alpha=alpha)
                        # plt.scatter(np.random.normal(loc=pos, scale=1, size=(Nclust, 1)), AllDepth, c=vMI[condition][direction], cmap='coolwarm', s=s, alpha=alpha)

                    plt.gca().invert_yaxis()
                    if posOrientation == 'ML_pos':
                        plt.gca().invert_xaxis()

                    plt.xlabel(f"{pos_title} position (mm)", fontsize=14)
                    plt.ylabel("Depth (µm)", fontsize=14)
                    plt.title(f"{direction} modulation of units in {pos_title} axis", fontsize=16)

                plt.colorbar(label=f"{direction} modulation index")

                if save:
                    direction_modulation_folder = os.path.join(analyse_path, 'Direction_modulation')
                    os.makedirs(direction_modulation_folder, exist_ok=True)
                    plt.savefig(os.path.join(direction_modulation_folder , f"{direction}_modulation_{pos_title}.png"))
                if show:
                    plt.show()
                else:
                    plt.close()
                print('\n')











# dirMI FUNCTION

def dirMI_function(AllData, analyse_path, save = False,s=35,alpha=0.6, show=True, scale=0.008):
    for pos_title, posOrientation in zip(['Mediolateral', 'Anteroposterior'], ['ML_pos', 'AP_pos']):
        for condition in ['second']:

            with load_theme("arctic_light"):
                plt.figure(figsize=(17, 5))

                for animal in AllData:
                    pos = AllData[animal]['informative_data'][posOrientation]
                    Nclust = AllData[animal]['SUA_data']['Nclust']
                    AllDepth = AllData[animal]['MUA_data']['AllDepth']
                    dirMI = AllData[animal]['Statistics_data']['dirMI']

                    plt.scatter(np.random.normal(loc=pos, scale=scale, size=(Nclust, 1)), AllDepth, c=dirMI[condition], cmap='coolwarm', s=s, alpha=alpha)
                    # plt.scatter((pos+random.uniform(-0.005,0.005))*np.ones(Nclust), AllDepth, c=dirMI[condition], cmap='coolwarm', s=s, alpha=alpha)

                plt.gca().invert_yaxis()
                if posOrientation == 'ML_pos':
                    plt.gca().invert_xaxis()

                plt.xlabel(f"{pos_title} position (mm)", fontsize=14)
                plt.ylabel("Depth (µm)", fontsize=14)
                plt.title(f"CW vs CCW preference in {pos_title} axis", fontsize=16)

            plt.colorbar(label=r"Direction Modulation Index : $\frac{CW - CCW}{CW + CCW}$")


                
            if save:
                direction_preference_folder = os.path.join(analyse_path, 'Direction_preference')
                os.makedirs(direction_preference_folder, exist_ok=True)
                plt.savefig(os.path.join(direction_preference_folder , f"Direction_preference_{pos_title}.png"))
            if show:
                plt.show()
            else:
                plt.close()
            print('\n')
















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










def scatter3D(x,y,z,colors,colorlabel,xlabel,ylabel,zlabel,title,filename,filepath,save=True,show=True,anim=True):

    theme = load_theme("arctic_light").set_overrides({
        "ytick.minor.visible": False,
        "xtick.minor.visible": False,
    })

    with theme:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        scatter = ax.scatter(x, y, z, c=colors, cmap='coolwarm', marker='o')

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

        if save & anim:
            os.makedirs(filepath, exist_ok=True)
            ani.save(os.path.join(filepath, filename), writer='pillow')
        elif save:
            os.makedirs(filepath, exist_ok=True)
            plt.savefig(os.path.join(filepath, filename), format='pdf')
        if show:
            plt.show()













        