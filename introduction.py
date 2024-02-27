import scipy.io
import os
import h5py
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
# from scipy.stats import sem
from scipy.stats import wilcoxon
import pandas as pd




def filepath(path, file):
    return os.path.join(path, file) 





def plots(row_number, col_number, *args, suptitle=None, **kwargs):
    if len(args) != row_number * col_number:
        raise ValueError("Le nombre d'arguments fournis ne correspond pas au nombre de sous-graphiques attendus.")
    
    fig, axs = plt.subplots(row_number, col_number)#, figsize=(width, height))
    
    # Flatten axes if there's only one row or one column
    if row_number == 1 or col_number == 1:
        axs = axs.reshape(-1)
    
    for i, plot_call in enumerate(args):
        if row_number > 1 and col_number > 1:  # If there's more than one row and more than one column
            col = i // row_number
            row = i % row_number
            ax = axs[row, col]
        else:  # If there's only one row or one column
            ax = axs[i]
        plot_call(ax)

    for ax in axs.flat:
        ax.set(**kwargs)

    plt.suptitle(suptitle)

    # plt.tight_layout()
    # fig.set_tight_layout(True)

    plt.show()





def plotVelocity(velocity, color='c', ax='', xlabel=True, ylabel=True):
    if ax:
        ax.plot(duration, velocity, color=color, linewidth='1')
        ax.set_ylim(-rotationSpeed, rotationSpeed)
        ax.set_xlim(-timeBef, timeAft)
        if ylabel:
            ax.set_ylabel(r'Mean Velocity $(^\circ/s)$')
        if xlabel:
            ax.set_xlabel(r'Time $(s)$')
        ax.set_yticks([-rotationSpeed, 0, rotationSpeed])
    else:
        plt.plot(duration, velocity, color=color, linewidth='1')
        plt.ylim(-rotationSpeed, rotationSpeed)
        plt.xlim(-timeBef, timeAft)
        if ylabel:
            plt.ylabel(r'Mean Velocity $(^\circ/s)$')
        if xlabel:
            plt.xlabel(r'Time $(s)$')
        plt.yticks([-rotationSpeed, 0, rotationSpeed])
        plt.show()






def plotRaster(spikeTimesObject, color='black', ax='', xlabel=True, ylabel=True, extra=None, title='', timeBef=timeBef, timeAft=timeAft):

    linelengths = 1

    if ax:
        # for trial_id, spike_times in enumerate(spikeTimesObject):
        #     ax.eventplot(spike_times, lineoffsets=trial_id+1, linelengths=linelengths, color=color)

        ax.eventplot(spikeTimesObject, linelengths=linelengths, colors=color)

        if ylabel:
            ax.set_ylabel('# Trial')
        if xlabel:
            ax.set_xlabel('Time (s)')
            
        ax.set_title('')
        ax.set_xlim(-timeBef,timeAft)
        
        ### annulate the offset due to python indexation
        def custom_formatter(x, pos):
            return f"{int(x) + 1}"
        ax.yaxis.set_major_formatter(plt.FuncFormatter(custom_formatter))

        ax.set_ylim(0-linelengths/2,len(spikeTimesObject)-1+linelengths/2)

        if extra is not None:
            extra()
        
        ax.set_title(title)
            
    else:
        # for trial_id, spike_times in enumerate(spikeTimesObject):
        #     plt.eventplot(spike_times, lineoffsets=trial_id+1, linelengths=linelengths, color=color)

        plt.eventplot(spikeTimesObject, linelengths=linelengths, colors=color)


        if extra is not None:
            extra()
        
        if ylabel:
            plt.ylabel('# Trial')
        if xlabel:
            plt.xlabel('Time (s)')
        
        plt.title('')
        plt.xlim(-timeBef,timeAft)

        ### annulate the offset due to python indexation
        def custom_formatter(x, pos):
            return f"{int(x) + 1}"
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(custom_formatter))

        plt.ylim(0-linelengths/2,len(spikeTimesObject)-1+linelengths/2)

        plt.title(title)

        plt.show()




def plotPSTH(StudiedSpikeTimes,color='k',shadedcolor='c',binResolution = 0.03,ax='',xlabel=True,ylabel=True):
    if ax:
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


        # plt.figure(figsize=(15,6))
        ax.plot(edges[:-1], Zscore, color=color)
        ax.fill_between(edges[:-1], Zscore-SEM, Zscore+SEM, alpha=0.1, color=shadedcolor)
        ax.set_xlim(-timeBef,timeAft)
        if ylabel:
            ax.set_ylabel('PSTH')
        if xlabel:
            ax.set_xlabel('Time (s)')
    else:
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


        # plt.figure(figsize=(15,6))
        plt.plot(edges[:-1], Zscore, color=color)
        plt.fill_between(edges[:-1], Zscore-SEM, Zscore+SEM, alpha=0.1, color=shadedcolor)
        plt.xlim(-timeBef,timeAft)
        if ylabel:
            plt.ylabel('PSTH')
        if xlabel:
            plt.xlabel('Time (s)')
        plt.show()








def figure1(unit,firstcolor='c',secondcolor='k',binResolution=0.03,suptitle='Figure 1 unit'):
      plots(3, 2,
            lambda ax: plotVelocity(MeanRotation['first']['CW'],color=firstcolor, ax=ax, xlabel=False),
            lambda ax: plotRaster(SpikeTimes['first']['CW'][unit], color=firstcolor, ax=ax, xlabel=False),
            lambda ax: plotPSTH(SpikeTimes['first']['CW'][unit],color=firstcolor, shadedcolor=firstcolor, ax=ax),
            lambda ax: plotVelocity(MeanRotation['second']['CCW'],color=secondcolor,ax=ax,ylabel=False,xlabel=False),
            lambda ax: plotRaster(SpikeTimes['second']['CCW'][unit], color=secondcolor, ax=ax,xlabel=False,ylabel=False),
            lambda ax: plotPSTH(SpikeTimes['second']['CCW'][unit],color=secondcolor, shadedcolor=secondcolor, ax=ax,ylabel=False),
            suptitle=suptitle
      )








def seeModulation(condition,modulated):
      position = np.where([(modulation[condition][neuron]['selectivity']==modulated.split()[0]) and (modulation[condition][neuron]['type']==modulated.split()[1]) for neuron in range(Nclust)])[0]
      print(f'The {modulated} modulated neurons are {position}')

      for unit in np.where([(modulation[condition][neuron]['selectivity']==modulated.split()[0]) and (modulation[condition][neuron]['type']==modulated.split()[1]) for neuron in range(Nclust)])[0]:
            plots(1,2,
                  lambda ax: plotRaster(SpikeTimes[condition]['CW'][unit], ax=ax),
                  lambda ax: plotRaster(SpikeTimes[condition]['CCW'][unit], ax=ax, ylabel=False)
                  )
            



