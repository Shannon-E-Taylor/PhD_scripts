import numpy as np
import pandas as pd

import numpy.linalg as LA
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, root
import sympy as sm


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import get_phase_portraits

from scipy.ndimage.filters import uniform_filter1d
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import warnings


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





def classify_steady_state(steady_state, fgf, wnt, params):

    gene1 = sm.symbols("gene1", negative=False)
    gene2 = sm.symbols("gene2", negative=False)
    gene3 = sm.symbols("gene3", negative=False)

    eqMat = sm.Matrix(get_phase_portraits.PSH_sym(
        [gene1, gene2, gene3],
        wnt, fgf, 1,
        params = params))
    Mat = sm.Matrix([gene1, gene2, gene3])
    jacMat = eqMat.jacobian(Mat)
    eigen = []

    eqmat = jacMat.subs([(gene1, steady_state[0]), (gene2, steady_state[1]), (gene3, steady_state[2])])
    eigen = list(eqmat.eigenvals().keys())

    # extract the imaginary component of the eigenvalue
    imag = [sm.im(i) for i in eigen]
    # extract real component
    real = [sm.re(i) for i in eigen]

    # if any of the imaginary components are non-zero,
    # we have a spiral attractor which we colour magenta
    if any(i!=0 for i in imag):
        return ('m')

    # if all real components are negative,
    # we have a normal attractor, in blue
    if all(i<0 for i in real):
        return('b')

    # if all real component are positive,
    # we have a repellor, in green
    if all(i>0 for i in real):
        return('g')

    #otherwise, its a saddle point, in red
    else:
        return ('r')


def compute_trajectory(cell_ics, wnt, fgf, params):

    tbxta, tbx16, tbx24 = cell_ics[['g1', 'g2', 'g3']]

    t = np.arange(0, 10, 0.001)

    # calculate long-term trajectory of cell
    traj = odeint(PSH, y0 = [tbxta, tbx16, tbx24], t = t,
            args = (wnt, fgf, params))

    cell_ics['g1_traj'] = traj[:, 0]
    cell_ics['g2_traj'] = traj[:, 1]
    cell_ics['g3_traj'] = traj[:, 2]
    return(cell_ics)


def compute_tf(cell_ics, wnt, fgf, params):

    tbxta, tbx16, tbx24 = cell_ics[['g1', 'g2', 'g3']]

    t = np.arange(0, 10, 0.001)

    # calculate long-term trajectory of cell
    traj = odeint(PSH, y0 = [tbxta, tbx16, tbx24], t = t,
            args = (wnt, fgf, params))

    cell_ics['g1_traj'] = traj[-1, 0]
    cell_ics['g2_traj'] = traj[-1, 1]
    cell_ics['g3_traj'] = traj[-1, 2]
    return(cell_ics)



def compute_attr(cell_ics, wnt, fgf, params):

    tbxta, tbx16, tbx24 = cell_ics[['g1', 'g2', 'g3']]

    # Try to find the steady state using ICs as inputs
    soln = root(PSH,
                np.array([tbxta, tbx16, tbx24]),
                args = (1, wnt, fgf, params))

    if soln.success:
        ss = soln.x

    # If this doesn't work, simulate it
    else:
        # # calculate long-term trajectory of cell
        traj = odeint(PSH, y0 = [tbxta, tbx16, tbx24], t = np.arange(0, 100, 0.1),
                args = (wnt, fgf, params))
        if np.isnan(traj).any():
            warnings.warn(f"NaN values detected in trajectory for initial conditions {cell_ics.values}. Skipping.", RuntimeWarning)
            return None
        soln = root(PSH, traj[-1], args = (1, wnt, fgf, params))

        if soln.success:
            ss = soln.x
        else:
            print(f'Error: no ss found for {tbxta, tbx16, tbx24}')
            print(wnt, fgf)
            print(soln)
            return None


    flow = PSH(ss, 't', wnt, fgf, params)
    if not np.allclose(flow, [0, 0, 0], atol=1e-6):
        # print(f"Steady-state verification failed for {ss}")
        warnings.warn(f"Steady-state verification failed for {[ss, wnt, fgf]}. Skipping.", RuntimeWarning)
        # assert True == False # stop it
        print(flow)
        return None

    return(ss)

