import get_phase_portraits
import numpy as np
import pandas as pd

import numpy.linalg as LA
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, root
import sympy as sm


def classify_ics(cell_ics, wnt, fgf, params, t = np.arange(0, 100, 0.01)): 
    # I'm making ss a global variable 
    # so that modifiying it within this code will hold 
    # for the rest of the loop
    # this may be quite bad?
    global ss
    
    tbxta, tbx16, tbx24 = cell_ics['g1'], cell_ics['g2'], cell_ics['g3']

    # calculate long-term trajectory of cell 
    traj = odeint(get_phase_portraits.PSH, y0 = [tbxta, tbx16, tbx24], t = t, 
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
        traj = odeint(get_phase_portraits.PSH, y0 = traj[-1], t = t, 
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
        soln = root(get_phase_portraits.PSH, traj[-1], args = (1, wnt, fgf, params))
        if soln.success: 
            sol = soln.x
            flow = get_phase_portraits.PSH(sol, 1, wnt, fgf, params) 
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


def classify_attr_by_row(row, params): 
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