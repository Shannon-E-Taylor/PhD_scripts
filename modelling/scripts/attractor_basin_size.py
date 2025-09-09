import numpy as np
import pandas as pd

import numpy.linalg as LA
from scipy.integrate import odeint
from scipy.optimize import fsolve, root
import sympy as sm

from scipy.ndimage.filters import uniform_filter1d
from scipy.interpolate import griddata

import warnings

from phase_portrait_computation import *


def compute_attr_basin(params, name):

    all_sims = []

    # For every combination of Wnt and Fgf, calculate and save flow and steatdy states
    for wnt in WNT_values:
        for fgf in FGF_values:

            traj_sim = ics.sample(
                frac=1
                ).apply(compute_tf,
                        args = (wnt, fgf, params),
                        axis = 1)


            traj_sim['Wnt_sim'] = wnt
            traj_sim['FGF_sim'] = fgf
            all_sims.append(traj_sim)

    all_sims_df = pd.concat(all_sims)
    all_sims_df.to_csv(
        f'../output/{name}_attr_basin.csv'
        )

    return(all_sims_df)


WNT_values = np.arange(0,1.3,0.1) #0.01
FGF_values = np.arange(0,1.3,0.1) #0.01

simulation = pd.read_csv(
    '../output/MAP.csv',
    sep = ';')
simulation['time'] = simulation['Time']


ics = simulation[simulation['Time'] == 1]#.sample(frac = 0.1)

all_params = pd.read_csv('../median_networks_for_clusters.csv')



import multiprocessing as mp

# Define a function that processes each row
def process_row(row):
    name = row['cluster']
    params = row[1:25]
    print(name)  # Debugging print statement
    compute_attr_basin(params, name)  # Call your function

# Create a pool of workers, using the number of available CPU cores
if __name__ == "__main__":
    num_workers = 8  # Eight rows
    with mp.Pool(num_workers) as pool:
        pool.map(process_row, [row for _, row in all_params.iterrows()])

