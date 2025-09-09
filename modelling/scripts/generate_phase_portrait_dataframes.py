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

import get_phase_portraits

### SCRIPTS

def PSH(s, t, wnt, fgf, params):

    B = np.array([wnt, fgf])

    W = np.array([[params[0] , params[1], params[2]],\
                 [ params[3],  params[4], params[5]],\
                 [ params[6],  params[7], params[8]]])

    E = np.array([[params[18],params[17]],[params[19],params[15]],[params[20],params[16]]])
    R = [params[9],  params[10],  params[11]]
    lmd = [params[12],  params[13] , params[14]]
    h = [params[21],  params[22] , params[23]]


    u0 = W[0][0]*s[0] + W[0][1]*s[1] + W[0][2]*s[2] + E[0].dot(B) + h[0]
    u1 = W[1][0]*s[0] + W[1][1]*s[1] + W[1][2]*s[2] + E[1].dot(B) + h[1]
    u2 = W[2][0]*s[0] + W[2][1]*s[1] + W[2][2]*s[2] + E[2].dot(B) + h[2]

    d_tbxta_dt = R[0] * g(u0)  - lmd[0]* s[0]
    d_tbx16_dt = R[1] * g(u1)  - lmd[1]* s[1]
    d_tbx24_dt = R[2] * g(u2)  - lmd[2]* s[2]

    dsdt = [d_tbxta_dt, d_tbx16_dt, d_tbx24_dt]

    return dsdt

# Function for all other genes
def fit_func2(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.0) / (2.0 * np.power(sig, 2.0)))

def fit_func2_symb(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.0) / (2.0 * np.power(sig, 2.0)))

# Function for tbxta
#def fit_func2_tbxta(x, a, b, c):
#    return np.power((1 - np.power(x, a)/(np.power(x, a) + np.power(b, a))), c)

def g(x):
    return (0.5* ((x / np.sqrt(x ** 2 + 1)) + 1))

def g_sym(x):
    return (0.5* ((x / sm.sqrt(x ** 2 + 1)) + 1))



def classify_ics(cell_ics, params, wnt, fgf, t = np.arange(0, 100, 0.01)):
    # I'm making ss a global variable
    # so that modifiying it within this code will hold
    # for the rest of the loop
    # this may be quite bad?
    global ss

    tbxta, tbx16, tbx24 = cell_ics['g1'], cell_ics['g2'], cell_ics['g3']

    # calculate long-term trajectory of cell
    traj = odeint(PSH, y0 = [tbxta, tbx16, tbx24], t = t,
              args = (wnt, fgf, params))

    flow = get_phase_portraits.PSH(traj[-1],
                                   t = 'd',
                                   wnt = wnt, fgf = fgf,
                                   params = params)
    # give this junk values
    final_attr = [200, 200, 200]
    # if we aren't close to a steady state
    # simulate again!
    if not np.isclose(flow, [0, 0, 0]).all():
        traj = odeint(PSH, y0 = traj[-1], t = np.arange(0, 1000, 0.01),
              args = (wnt, fgf, params))
        flow = get_phase_portraits.PSH(traj[-1],
                                   t = 'd',
                                   wnt = wnt, fgf = fgf,
                                   params = params)

    # break the code if we haven't hit a steady state yet
    # just for testing
    #assert np.isclose(flow, [0, 0, 0]).all()

    # see if attractor is one of the identified steady states
    found_ss = False
    for steady_state in ss:
        if (np.isclose(traj[-1], steady_state, atol = 0.01).all() and not found_ss):
            final_attr = steady_state
            found_ss = True
    # if the final timepoint of the trajectory is not an existing steady state,
    # figure out the nearby steady state
    # and add it to the steady state list
    # this should be more reliable
    if not found_ss:
        soln = root(PSH, traj[-1], args = (1, wnt, fgf, params))
        if soln.success:
            sol = soln.x
            flow = PSH(sol, 1, wnt, fgf, params)
            assert np.isclose(flow, [0, 0, 0], atol = 1e-7).all()
            final_attr = sol
            # add the solution to the steady states list
            # this will overwrite the global ss object
            # (I think...)
            if ss.shape[0]==0:
                ss = np.array([sol])
            else:
                ss = np.append(ss, [sol], axis = 0)
        else:
            print(soln)
            print(traj[-1])
            print(flow)
            print(fgf, wnt)
            print('Problem: solver failed when IDing new steady state')
            assert True == False #pfffft

    cell_ics['flow'] = flow
    # rounding as a crude way of keeping very close attractors the same
    cell_ics['final_attractor_g1'] = final_attr[0]
    cell_ics['final_attractor_g2'] = final_attr[1]
    cell_ics['final_attractor_g3'] = final_attr[2]
    cell_ics['Wnt'] = wnt
    cell_ics['FGF'] = fgf
    return cell_ics




# Define the space of your state variables
tol = 0.95
num_posi = 10
tbxta_ = np.linspace(0, 1, num_posi)
tbx16_ = np.linspace(0, 1, num_posi)
tbx24_ = np.linspace(0, 1, num_posi)

# Wnt and Fgf values to explore
WNT_values = np.arange(0,1.3,0.1) #0.01
FGF_values = np.arange(0,1.3,0.1) #0.01

tbxta, tbx16, tbx24 = np.meshgrid(tbxta_, tbx16_, tbx24_)

def attractors_by_ic_all_sig(params, name):

    ss = np.array([]) # need to initialize and reset ss

    all_sims = []

    FGF, Wnt = [], []

    steady_states = []
    steady_states_cells = []
    steady_states_root = []

    # For every combination of Wnt and Fgf, calculate and save flow and steatdy states
    for wnt in WNT_values:
        for fgf in FGF_values:
            # at the moment I'm not bothering to calculate the flow, ss etc
            # because we're going to do that again anyways

            ss = np.array([]) # need to initialize and reset ss

            #Flows.append(Flow)
            #steady_states.append(ss)
            #steady_states_cells.append(ss_cells)
            #steady_states_root.append(ss)

            df = ics.apply(classify_ics, args = (params, wnt, fgf), axis = 1)
            all_sims.append(df)

    all_sims_df = pd.concat(all_sims)
    all_sims_df.to_csv(f'../output/{name}_network_attractors.csv')

    return(all_sims_df)



def classify(ss, wnt, fgf):
    gene1 = sm.symbols("gene1", negative=False)
    gene2 = sm.symbols("gene2", negative=False)
    gene3 = sm.symbols("gene3", negative=False)

    eqMat = sm.Matrix(get_phase_portraits.PSH_sym([gene1, gene2, gene3], wnt, fgf, 1, params = params))
    Mat = sm.Matrix([gene1, gene2, gene3])
    jacMat = eqMat.jacobian(Mat)
    eigen = []

    for steady_state in ss:
        eqmat = jacMat.subs([(gene1, steady_state[0]), (gene2, steady_state[1]), (gene3, steady_state[2])])

        if len(eigen) == 0:
            eigen = list(eqmat.eigenvals().keys())
        else:
            eigen = np.vstack((eigen, list(eqmat.eigenvals().keys())))

    eigen = np.atleast_2d(eigen)

    color_ss = []

    for j in range(np.shape(ss)[0]):
        #if np.sum([abs(np.imag(complex(eigen[j][x].evalf()))) for x in range(0, 3)]) == 0:
        if (
                (np.imag(complex(eigen[j][0].evalf())) == 0)
                & (np.imag(complex(eigen[j][1].evalf())) == 0)
                & (np.imag(complex(eigen[j][2].evalf())) == 0)
            ):
            if (
                (np.real(complex(eigen[j][0].evalf())) > 0)
                & (np.real(complex(eigen[j][1].evalf())) > 0)
                & (np.real(complex(eigen[j][2].evalf())) > 0)
            ):
                # all positive real eigen values mean repeller (green)
                #ax.scatter(ss[j, 0], ss[j, 1], ss[j, 2], color="g")
                color_ss.append("g")
            elif (
                (np.real(complex(eigen[j][0].evalf())) < 0)
                & (np.real(complex(eigen[j][1].evalf())) < 0)
                & (np.real(complex(eigen[j][2].evalf())) < 0)
            ):
                # all negative real eigen values mean attracters (blue)
                #ax.scatter(ss[j, 0], ss[j, 1], ss[j, 2], color="b")
                color_ss.append("b")

            else:
                # any different signed real eigen values mean saddle point (red)
                #ax.scatter(ss[j, 0], ss[j, 1], ss[j, 2], color="r")
                color_ss.append("r")

        else:
            # any complex eigen value mean spiral (magenta)
            #ax.scatter(ss[j, 0], ss[j, 1], ss[j, 2], color="m")
            color_ss.append("m")

    return([wnt, fgf, color_ss.count('g'), color_ss.count('b'), color_ss.count('r'), color_ss.count('m')])


def classify_attr_by_row(row):
    wnt = row['Wnt']
    fgf = row['FGF']
    steady_state = [row['final_attractor_g1'], row['final_attractor_g2'], row['final_attractor_g3']]
    # double check this is actually an attractor!!
    #print(PSH(steady_state, t = 1, wnt = wnt, fgf = fgf, params = params))
    #assert np.isclose(PSH(steady_state, t = 1, wnt = wnt, fgf = fgf, params = params), [0, 0, 0]).all()

    gene1 = sm.symbols("gene1", negative=False)
    gene2 = sm.symbols("gene2", negative=False)
    gene3 = sm.symbols("gene3", negative=False)

    eqMat = sm.Matrix(get_phase_portraits.PSH_sym([gene1, gene2, gene3], wnt, fgf, 1, params = params))
    Mat = sm.Matrix([gene1, gene2, gene3])
    jacMat = eqMat.jacobian(Mat)
    eigen = []

    eqmat = jacMat.subs([(gene1, steady_state[0]), (gene2, steady_state[1]), (gene3, steady_state[2])])

    #if len(eigen) == 0:
    eigen = list(eqmat.eigenvals().keys())

    #eigen = np.atleast_2d(eigen)

    color_ss = []

    #if np.sum([abs(np.imag(complex(eigen[j][x].evalf()))) for x in range(0, 3)]) == 0:
    if (
            (np.imag(complex(eigen[0].evalf())) == 0)
            & (np.imag(complex(eigen[1].evalf())) == 0)
            & (np.imag(complex(eigen[2].evalf())) == 0)
        ):
        if (
            (np.real(complex(eigen[0].evalf())) > 0)
            & (np.real(complex(eigen[1].evalf())) > 0)
            & (np.real(complex(eigen[2].evalf())) > 0)
        ):
            # all positive real eigen values mean repeller (green)
            #ax.scatter(ss[j, 0], ss[j, 1], ss[j, 2], color="g")
            color_ss = "g"
        elif (
            (np.real(complex(eigen[0].evalf())) < 0)
            & (np.real(complex(eigen[1].evalf())) < 0)
            & (np.real(complex(eigen[2].evalf())) < 0)
        ):
            # all negative real eigen values mean attracters (blue)
            #ax.scatter(ss[j, 0], ss[j, 1], ss[j, 2], color="b")
            color_ss = "b"

        else:
            # any different signed real eigen values mean saddle point (red)
            #ax.scatter(ss[j, 0], ss[j, 1], ss[j, 2], color="r")
            color_ss = "r"

    else:
        # any complex eigen value mean spiral (magenta)
        #ax.scatter(ss[j, 0], ss[j, 1], ss[j, 2], color="m")
        color_ss = "m"

    row['attractor type'] = color_ss
    return(row)


all_params = pd.read_csv('../median_networks_for_clusters.csv')

ss = np.array([]) # need to initialize and reset ss



simulation = pd.read_csv(
    '../output/MAP.csv',
    sep = ';')
simulation['time'] = simulation['Time']


ics = simulation[simulation['Time'] == 1]#.sample(frac = 0.1)

for _, row in all_params.iterrows():

    name = row['cluster']
    params=row[1:25]
    print(name)

    all_sims_df = attractors_by_ic_all_sig(params, name)

    keep_cols = ['FGF', 'Wnt',
                 'final_attractor_g1',
                 'final_attractor_g2',
                 'final_attractor_g3'
                 ]

    unique_attr = all_sims_df[
        keep_cols
         ].reset_index(
             drop = True
             ).drop_duplicates(subset = keep_cols)

    classified_unique_attr = unique_attr.apply(classify_attr_by_row, axis = 1)
    classified_unique_attr.to_csv(f'../output/{name}_network_classified_unique_attractors.csv')


    # print(all_sims_df.columns)
    # print(classified_unique_attr.columns)

    classified_attractors = all_sims_df.merge(classified_unique_attr[[
        'final_attractor_g1', 'final_attractor_g2', 'final_attractor_g3', 'FGF', 'Wnt',
        'attractor type', # 'attractor name'
    ]], on = ['final_attractor_g1', 'final_attractor_g2', 'final_attractor_g3', 'FGF', 'Wnt'])

    classified_attractors.to_csv(f'../output/{name}_network_classified_attractors.csv')


    classified_ics = []

    for i, cell in ics.iterrows():
        cell = classify_ics(cell, params, cell['Wnt'], cell['FGF'])
        classified_ics.append(pd.DataFrame(cell).T)

    classified_ics = pd.concat(classified_ics).reset_index(drop = True)
    classified_ics = classified_ics.apply(classify_attr_by_row, axis = 1)

    classified_ics.to_csv(f'../output/{name}_network_classified_ics.csv')




