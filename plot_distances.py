import matplotlib.pyplot as plt 
import numpy as np


def plot_dists(metadata):
    fig = plt.figure()
    fig, axs = plt.subplots(2, 4)
    fig.suptitle("Distance Plots| System <PDB>")

#    for i in range(1, 3):
#        for j in range(1, 4):
    axs[0, 0].plot(metadata['distances']['stage1']['trapped water'])
    axs[0, 0].plot(metadata['distances']['stage1']['solvent water 1st'])
    axs[0, 0].plot(metadata['distances']['stage1']['solvent water 2nd'])

    axs[0, 1].plot(metadata['distances']['stage2']['trapped water'])
    axs[0, 1].plot(metadata['distances']['stage2']['solvent water 1st'])
    axs[0, 1].plot(metadata['distances']['stage2']['solvent water 2nd'])

    axs[0, 2].plot(metadata['distances']['stage3']['trapped water'])
    axs[0, 2].plot(metadata['distances']['stage3']['solvent water 1st'])
    axs[0, 2].plot(metadata['distances']['stage3']['solvent water 2nd'])

    axs[0, 3].plot(metadata['distances']['stage4']['trapped water'])
    axs[0, 3].plot(metadata['distances']['stage4']['solvent water 1st'])
    axs[0, 3].plot(metadata['distances']['stage4']['solvent water 2nd'])

    axs[1, 0].plot(metadata['distances']['stage5']['trapped water'])
    axs[1, 0].plot(metadata['distances']['stage5']['solvent water 1st'])
    axs[1, 0].plot(metadata['distances']['stage5']['solvent water 2nd'])

    axs[1, 1].plot(metadata['distances']['stage6']['trapped water'])
    axs[1, 1].plot(metadata['distances']['stage6']['solvent water 1st'])
    axs[1, 1].plot(metadata['distances']['stage6']['solvent water 2nd'])

    axs[1, 2].plot(metadata['distances']['stage7']['trapped water'])
    axs[1, 2].plot(metadata['distances']['stage7']['solvent water 1st'])
    axs[1, 2].plot(metadata['distances']['stage7']['solvent water 2nd'])

    
    for ax in axs.flat:
        ax.set(xlabel='time (ps)', ylabel='distance (nm)')
        ax.label_outer()

    plt.savefig('d1-7.png')
    
