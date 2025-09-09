import numpy as np
import pandas as pd

import numpy.linalg as LA
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, root
import sympy as sm
import os
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.gridspec as gridspec




all_params = pd.read_csv('../median_networks_for_clusters.csv')



simulation = pd.read_csv(
    '../output/MAP.csv',
    sep = ';')
simulation['time'] = simulation['Time']

import string

alphabet = list(string.ascii_uppercase)

label_dict = {
    'm': 'spiral sink',
    'b': 'attractor',
    'r': 'saddle',
    'g': 'repeller',
    'k': 'unlabelled',

}


ics = simulation[simulation['Time'] == 1]# .sample(frac = 0.1)

for _, row in all_params.iterrows():

    name = row['cluster']
    params=row[1:25]



    print(name)


    classified_unique_attr = pd.read_csv(
        f'../output/{name}_network_attractors.csv'
    )


    fig = plt.figure(figsize = (7, 2.5))
    # fig.set_layout_engine(layout='constrained')
    gs = gridspec.GridSpec(1, 4, width_ratios=[1,1,1, 0.2])


    ax0 = fig.add_subplot(gs[0], projection='3d')
    ax1 = fig.add_subplot(gs[1], projection='3d')
    ax2 = fig.add_subplot(gs[2], projection='3d')
    cbarax = fig.add_subplot(gs[3])

    lab = [i for i in classified_unique_attr['type']]

    g1 = list(classified_unique_attr['g1'])
    g2 = list(classified_unique_attr['g2'])
    g3 = list(classified_unique_attr['g3'])

    for l in set(lab):
        x = [g1[i] for i in range(len(lab)) if lab[i]==l]
        y = [g2[i] for i in range(len(lab)) if lab[i]==l]
        z = [g3[i] for i in range(len(lab)) if lab[i]==l]

        ax0.scatter3D(x, y, z,
                    c = l,
                    label = label_dict[l],
                        edgecolors='none',
                        s = 2)

    im1 = ax1.scatter3D(classified_unique_attr['g1'],
                    classified_unique_attr['g2'],
                    classified_unique_attr['g3'],
                    c = classified_unique_attr['Wnt'],
                    edgecolors='none',
                    s = 2)


    im2 = ax2.scatter3D(classified_unique_attr['g1'],
                    classified_unique_attr['g2'],
                    classified_unique_attr['g3'],
                    edgecolors='none',
                    c = classified_unique_attr['FGF'],
                    s=2)


    if name != 'MAP':
        ax0.set_title(f'{alphabet[int(name)]}.i: Attractor type')
        ax1.set_title(f'{alphabet[int(name)]}.ii: Wnt level')
        ax2.set_title(f'{alphabet[int(name)]}.iii: FGF level')
    else:
        ax0.set_title(f'A: Attractor type')
        ax1.set_title(f'B: Wnt level')
        ax2.set_title(f'C: FGF level')



    for axs in [ax0, ax1, ax2]:
        axs.set_xlim(0, 1.5)
        axs.set_ylim(0, 1.5)
        axs.set_zlim(0, 2)

        # Set labels closer to axes
        axs.set_xlabel('tbxta', fontsize=10, labelpad=-5)
        axs.set_ylabel('tbx16', fontsize=10, labelpad=-5)
        axs.set_zlabel('tbx6', fontsize=10, labelpad=-2)

        # Reduce tick padding
        axs.tick_params(pad=-2, labelsize=8)


    ax0.legend(
        loc='upper center', bbox_to_anchor=(0.5, -0.18),
        prop={'size': 8},
        ncol = len(set(lab))
        )

    cbarax.axis('off')

    cbar1 = plt.colorbar(im2, ax = cbarax, orientation="vertical",
                         fraction=0.2, # default 0.15 - why?!
                         # label = 'signal',
                         # pad=.9
                         )#
    cbar1.ax.tick_params(labelsize='xx-small')

    # Adjust spacing
    fig.subplots_adjust(left=0.05, right=0.9,
                        bottom=0.05, top=0.95,
                        wspace=0.5)
    fig.suptitle(f'Cluster {name}')

    # Add a border using a rectangle
    rect = plt.Rectangle(
        (0, 0), 1, 1,  # Position and size (normalized coordinates)
        transform=fig.transFigure,
        color="black",
        linewidth=2,
        fill=False
    )

    fig.patches.append(rect)

    plt.savefig(
        f'../graphics/phase_portrait/phase_port_with_sig_{name}.png',
        facecolor = 'white',
        dpi = 600
    )

    plt.clf()
    plt.close()



    # plt.clf()
    # plt.scatter(classified_ics['X'], classified_ics['Y'],
    #             c = classified_ics['attractor name'], s = 15)
    # ax = plt.gca()
    # ax.set_aspect('equal', adjustable='box')
    # plt.savefig(
    #     f'../graphics/phase_portrait/attractor_on_embryo_{name}.png'
    #     )

    # plt.clf()
    # plt.scatter(classified_ics['X'], classified_ics['Y'],
    #             c = classified_ics['final_attractor_g3'], s = 15,
    #             vmin = 0, vmax = 1.5)
    # plt.colorbar()
    # ax = plt.gca()
    # ax.set_aspect('equal', adjustable='box')
    # plt.savefig(
    #     f'../graphics/phase_portrait/attractor_on_embryo_tbx6_level_{name}.png'
    # )