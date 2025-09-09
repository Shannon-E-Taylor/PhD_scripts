import numpy as np
import pandas as pd

import numpy.linalg as LA
from scipy.integrate import odeint
from scipy.optimize import fsolve, root
import sympy as sm

from scipy.ndimage.filters import uniform_filter1d
from scipy.interpolate import griddata

import warnings
import itertools

from phase_portrait_computation import *


def get_ss_from_grid(gene_expr_vals, wnt, fgf, t, params):

    gene_expr_vals_ = np.array(list(
        itertools.product(*[gene_expr_vals, gene_expr_vals, gene_expr_vals]))
        )

    psh = PSH(
        gene_expr_vals_.T, t, wnt, fgf, params)

    psh = np.array(psh)
    keep_idx = np.all(np.isclose(psh, 0, atol = 0.5), axis = 0)
    print(np.sum(keep_idx))

    expression_to_test = gene_expr_vals_[keep_idx]

    expression_to_test_df = pd.DataFrame(
    {'g1': expression_to_test[:, 0],
     'g2': expression_to_test[:, 1],
     'g3': expression_to_test[:, 2]
     })
    final_positions = expression_to_test_df.apply(compute_attr,
                                                args = (wnt, fgf, params),
                                                axis = 1)

    return(final_positions)



def save_attractors_by_sig(params, name):

    ss = np.array([]) # need to initialize and reset ss

    all_sims = []

    # For every combination of Wnt and Fgf, calculate and save flow and steatdy states
    for wnt in WNT_values:
        for fgf in FGF_values:

            final_positions = ics.apply(compute_attr,
                                               args = (wnt, fgf, params),
                                               axis = 1).tolist()
            final_positions = [x for x in final_positions if x is not None]

            fc_from_grid = get_ss_from_grid(gene_expr_vals,
                                            wnt, fgf, t,
                                            params
                                            ).tolist()

            fc = fc_from_grid + final_positions

            attractors = np.unique(
                np.round(np.array(fc),
                        decimals = 5),
                axis =0)

            ss_types = []

            for ss in attractors:
                ss_type = classify_steady_state(ss, fgf, wnt, params)
                ss_types.append(ss_type)

            attr_ = pd.DataFrame(
                {'g1': attractors[:, 0],
                'g2': attractors[:, 1],
                'g3': attractors[:, 2]}
            )

            attr_['Wnt'] = wnt
            attr_['FGF'] = fgf
            attr_['type'] = ss_types
            all_sims.append(attr_)

    all_sims_df = pd.concat(all_sims)
    all_sims_df.to_csv(
        f'../output/{name}_network_attractors.csv'
        )

    return(all_sims_df)


WNT_values = np.arange(0,1.3,0.1) #0.01
FGF_values = np.arange(0,1.3,0.1) #0.01

gene_expr_vals = np.arange(0, 5, 0.01)
t = np.arange(0, 10, 0.1)


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
    save_attractors_by_sig(params, name)  # Call your function

# Create a pool of workers, using the number of available CPU cores
if __name__ == "__main__":
    num_workers = 4  # Eight rows
    with mp.Pool(num_workers) as pool:
        pool.map(process_row, [row for _, row in all_params.iterrows()])