import numpy as np
import pandas as pd

import numpy.linalg as LA
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, root
import sympy as sm


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy.ndimage.filters import uniform_filter1d
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import warnings

from phase_portrait_computation import *


def graph_pp(params, name):

    fig = plt.figure(
        figsize = (3*3, 4*3),
        dpi = 600)

    for i, row in ics.iterrows():
        wnt = row['Wnt']
        fgf = row['FGF']

        final_positions = ics_sim.sample(frac=1).apply(compute_attr,
                                                       args = (wnt, fgf, params),
                                                       axis = 1).tolist()
        final_positions = [x for x in final_positions if x is not None]

        attractors = np.unique(
            np.round(np.array(final_positions),
                    decimals = 5),
            axis =0)

        ss_types = []

        for ss in attractors:
            ss_type = classify_steady_state(ss, fgf, wnt, params)
            # print(ss, ss_type)
            ss_types.append(ss_type)


        traj_sim = ics_sim.sample(
            frac=.1
            ).apply(compute_trajectory,
                    args = (wnt, fgf, params),
                    axis = 1)

        trajectories = []
        for _, row in traj_sim.iterrows():
            traj = np.array([row['g1_traj'], row['g2_traj'], row['g3_traj']]).T
            trajectories.append(traj)

        # Create a Line3DCollection for all trajectories
        lines = Line3DCollection(trajectories, colors='k',
                                alpha=0.5, linewidths=0.1)


        axs = fig.add_subplot(4, 3,
                              i+1,
                              projection = '3d'
                              )
        # axs.set_title(alphabet[i], loc='left')
        title = f'{alphabet[i]}. AP={round(ics.iloc[i, 2])}um\n    Wnt={round(wnt, 2)}, FGF={round(fgf, 2)}'


        axs.set_title(title, loc='left')
        axs.scatter3D(
            attractors[:, 0],
            attractors[:, 1],
            attractors[:, 2],
            c = ss_types,
            s = 10
        )

        axs.add_collection(lines)
        axs.set_xlim(0, 1.5)
        axs.set_ylim(0, 1.5)
        axs.set_zlim(0, 2)

        axs.set_xlabel('tbxta')
        axs.set_ylabel('tbx16')
        axs.set_zlabel('tbx6')
    fig.subplots_adjust(left=0.05, right=0.90,
                        bottom=0.05, top=0.95,
                        wspace=0.4, hspace = 0.4)


    # fig.suptitle(f'Cluster {name} trajectories\n')
    plt.savefig(f'../graphics/phase_portrait/{name}_network_Phase_portrat_with_traj.png',
                dpi = 600, facecolor = 'white')




import string
alphabet = list(string.ascii_uppercase)


all_params = pd.read_csv('../median_networks_for_clusters.csv')

all_params = all_params.sort_values(by = 'll', ascending=False).reset_index(drop = True)

name = 'MAP'
params=all_params.iloc[0, 1:25]

print(all_params['ll'][0])

simulation = pd.read_csv(
    '../output/MAP.csv',
    sep = ';')
simulation['time'] = simulation['Time']

ics_sim = simulation[simulation['Time'] == 1]# .sample(frac = 0.1)

genes_groundtruth = []

sim_s = simulation[simulation['time'] == 1].reset_index(drop = True)

x = np.array(list(sim_s.sort_values("X")['X']))

ics = pd.DataFrame({'Wnt': [], 'FGF': [], })
ics['X'] = np.linspace(20, 170, 12)

for idx, gene in enumerate(['Wnt', 'FGF']):
    yhat_1 = uniform_filter1d(
        sim_s.sort_values("X")[gene],
        size = 100, mode = 'nearest')
    genes_groundtruth.append(yhat_1)

    gene_ics = []
    for xpos in ics['X']:
        idx = (np.abs(x - xpos)).argmin()
        gene_ics.append(list(yhat_1)[idx])
    ics[gene] = gene_ics


plt_w = 2.5

for idx, row in all_params.iterrows():
    name = row['cluster']
    params = row[1:25]
    print(name)
    graph_pp(params, name)

