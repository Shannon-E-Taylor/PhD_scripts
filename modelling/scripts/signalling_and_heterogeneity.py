import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation #animations
from IPython.display import HTML

import pickle

from scipy.ndimage.filters import uniform_filter1d
import matplotlib.ticker as ticker


import seaborn as sns


from numba import jit, cfunc, njit
from numbalsoda import lsoda, lsoda_sig


col_list = ['#ed392b', '#feb441', '#4090c5']
col_max = ['#af1117', '#ea6e13', '#0f59a3']
gene_list = ['tbxta', 'tbx16', 'tbx6']
emb_cmap = ['Reds', 'YlOrBr', 'Blues']


plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Lato', 'Tahoma', 'DejaVu Sans',
                               'Lucida Grande', 'Verdana']

plt.rcParams['font.size'] = 12

import string
alphabet = list(string.ascii_uppercase)


parameter_values = pd.read_csv('../median_networks_for_clusters.csv')

def get_dist_to_rolling_means(sim_s, genes):
    genes_groundtruth = []
    distances = []
    sim_s = sim_s.sort_values("X")
    for gene in genes:
            yhat_1 = uniform_filter1d(sim_s[gene], size = 75, mode = 'nearest')
            genes_groundtruth.append(yhat_1)
            distances.append((yhat_1 - sim_s[gene])**2)

    return(genes_groundtruth, distances)



fastmath = False

@njit(fastmath=fastmath)
def dot(x, y):
    s = 0
    for i in range(len(x)):
        s += x[i]*y[i]
    return s

@njit(fastmath=fastmath)
def g(x):
    return 0.5 * ((x / np.sqrt(x ** 2 + 1)) + 1)

@cfunc(lsoda_sig)
def PSH_lsoda(t, s, du, params):
    # unpack params
    W1 = [params[0], params[1], params[2]]
    W2 = [params[3], params[4], params[5]]
    W3 = [params[6], params[7], params[8]]

    E1 = [params[18], params[17]]
    E2 = [params[19], params[15]]
    E3 = [params[20], params[16]]

    R = [params[9], params[10], params[11]]
    lmd = [params[12], params[13], params[14]]
    h = [params[21], params[22], params[23]]

    B = [params[24], params[25]]

    # do computation
    u = [
        dot(W1, s) + dot(E1, B) + h[0],
        dot(W2, s) + dot(E2, B) + h[1],
        dot(W3, s) + dot(E3, B) + h[2]
    ]

    # we assign results to du rather than returning them
    du[0] = R[0] * g(u[0]) - lmd[0] * s[0]
    du[1] = R[1] * g(u[1]) - lmd[1] * s[1]
    du[2] = R[2] * g(u[2]) - lmd[2] * s[2]



funcptr = PSH_lsoda.address # address to ODE function

@njit(fastmath=fastmath)
def solve_lsoda(s, t_interval, params):
    # define jitted functions
    usol, success = lsoda(funcptr, s, t_interval, params,
        rtol = 1.0e-8, atol = 1.0e-8) # these are important for acceptable accuracy
    assert success
    return usol


class Simulate_on_tracks:
    '''
    The simulator takes as input a set of 24 parameters describing the GRN
    and uses ODEs with the initial and boundary conditions from the AGETs to
    simulate gene expression for the chosen 10 AGETs. The simulated gene expression
    is then added to the AGET, so that the AGET includes the original target as well
    as the simulated gene expression.
    '''

    def __init__(self, parameters, list_of_tracks):
        self.params = parameters
        self.list_of_tracks = list_of_tracks

    def simulation(self):
        params = self.params
        tracks_df = self.list_of_tracks

        # precompute t_intervals to speed up function
        base_t_intervals = []
        for i in range(1, 61): # number of timepoints
            t_interval = np.linspace(i*90/3600/3, (i+1)*90/3600/3, 10)
            base_t_intervals.append(t_interval)

        # precompute the Time_nPSM
        time_npsm = [i*90/3600/3 for i in range(1, 62)]


        # add the params to an empty array of shape (26) to pass to PSH_lsoda
        params_to_pass = np.empty((26, ))
        params_to_pass[0:24] = params

        df_out = []
        # Iterate through celltracks (AGETs)
        for df_celltrack in tracks_df:
            # add biological time for the ODEs
            # I precomputed this earlier to speed up the iterations
            df_celltrack = df_celltrack.reset_index(drop = True)
            if df_celltrack.shape[0] == 61:
                df_celltrack['Time_nPSM'] = time_npsm
            else:
                df_celltrack['Time_nPSM'] = df_celltrack['Time']*90/3600/3


            if (df_celltrack.shape[0]==0):
                continue

            # collect parameters / ICs
            simulated_expression = [df_celltrack.loc[0, ['g1', 'g2', 'g3']].values]
            B1 = df_celltrack['Wnt'].values
            B2 = df_celltrack['FGF'].values
            Time_nPSM = df_celltrack['Time_nPSM'].values
            time = df_celltrack['Time'].values

            # loop through positions in dataframe
            for index in df_celltrack.index[:-1]:
                # if the timestep = 1,
                # we can use the precomputed information above
                count_timesteps_between = time[index+1] - time[index]
                if count_timesteps_between != 1:
                    t_interval = np.linspace(Time_nPSM[index], Time_nPSM[index+1], 10*int(count_timesteps_between))
                else:
                    t_interval = base_t_intervals[int(time[index])-1]

                # define initial conditions
                # and signalling
                s0 = simulated_expression[-1]
                params_to_pass[24:26] =  [B1[index], B2[index]]

                simulated_expression.append(solve_lsoda(
                    s0,
                    t_interval,
                    params_to_pass
                    )[-1]
                )

            simulated_expression = np.array(simulated_expression)
            df_celltrack['g1_sim'] = simulated_expression[:, 0]
            df_celltrack['g2_sim'] = simulated_expression[:, 1]
            df_celltrack['g3_sim'] = simulated_expression[:, 2]
            df_out.append(df_celltrack) # append AGET with the simulated expression to the output list

        return df_out

def tuftean_axes(ax, minx, maxx, miny, maxy):
    # Remove axis lines.
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    ax.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    right = False,
    ) # labels along the bottom edge are off

    ax.spines['bottom'].set_bounds(minx, maxx)
    ax.spines['left'].set_bounds(miny, maxy)

simulation = pd.read_csv('../output/MAP.csv', sep = ';')
simulation['time'] = simulation['Time']

for _, row in parameter_values.iterrows():
    name = row['cluster']
    params=row[1:25]

    if name == 'MAP':
        alphabet_idx = 0
    elif name == '1':
        alphabet_idx = 1
    else:
        alphabet_idx = max(int(name) - 1, 0)
    print(alphabet_idx, name)


    genes_groundtruth = []
    sim_s = simulation[simulation['Time']==1]
    sim_s = sim_s.sort_values("X")
    for idx, gene in enumerate(['g1', 'g2', 'g3']):
            yhat_1 = uniform_filter1d(sim_s[gene], size = 75, mode = 'nearest')
            genes_groundtruth.append(yhat_1)
            sim_s[gene] = yhat_1


    tracks_t0 = sim_s['TrackID']
    simulation = simulation[simulation['TrackID'].isin(tracks_t0)]

    list_of_dfs = []

    df2=sim_s.set_index(['TrackID','Time'])
    simulation_smooth_IC = df2.combine_first(simulation.set_index(['TrackID','Time'])).reset_index()

    for track in tracks_t0:
        df = simulation_smooth_IC[simulation_smooth_IC['TrackID']==track].reset_index(drop = True)
        list_of_dfs.append(df)


    sim_smoothed = Simulate_on_tracks(params, list_of_dfs).simulation()
    sim_smoothed = pd.concat(sim_smoothed)
    list_of_dfs_nc_sig = []

    for track in tracks_t0:
        df = simulation_smooth_IC[simulation_smooth_IC['TrackID']==track].reset_index(drop = True)
        df['FGF'] = np.mean(df['FGF'])
        df['Wnt'] = np.mean(df['Wnt'])
        list_of_dfs_nc_sig.append(df)

    sim_smoothed_nc_sig = Simulate_on_tracks(params, list_of_dfs_nc_sig).simulation()
    sim_smoothed_nc_sig = pd.concat(sim_smoothed_nc_sig)


    list_of_dfs_original = []

    for track in tracks_t0:
        df = simulation[simulation['TrackID']==track].reset_index(drop = True)
        # df['FGF'] = np.mean(df['FGF'])
        # df['Wnt'] = np.mean(df['Wnt'])
        list_of_dfs_original.append(df)

    original_sim = Simulate_on_tracks(params,
                                      list_of_dfs_original).simulation()
    original_sim = pd.concat(original_sim)

    orig_mse_over_time = []
    smooth_ic_mse_over_time = []
    smooth_ic_flat_mse_over_time = []
    for time in range(1, 61):
        genes_smooth, distances_original_sim = get_dist_to_rolling_means(original_sim[original_sim['Time']==time], ['g1_sim', 'g2_sim', 'g3_sim'])
        genes_smooth, distances_smooth_ic = get_dist_to_rolling_means(sim_smoothed[sim_smoothed['Time']==time], ['g1_sim', 'g2_sim', 'g3_sim'])
        genes_smooth, distances_smooth_ic_nc_sig = get_dist_to_rolling_means(sim_smoothed_nc_sig[sim_smoothed_nc_sig['Time']==time], ['g1_sim', 'g2_sim', 'g3_sim'])


        orig_mse_over_time.append(np.mean(distances_original_sim, axis = 1))
        smooth_ic_mse_over_time.append(np.mean(distances_smooth_ic, axis = 1))
        smooth_ic_flat_mse_over_time.append(np.mean(distances_smooth_ic_nc_sig, axis = 1))


    fig, ax = plt.subplots(
        1, 3, figsize = (7, 2.5),
        sharex = True, sharey = True, tight_layout = True)

    for i in range(3):

        ax[i].plot(np.array(orig_mse_over_time)[:, i], label='original simulation', c = 'k')
        ax[i].plot(np.array(
            smooth_ic_flat_mse_over_time)[:, i],
            label = 'Changing sig.',
            c = '#315997')
        ax[i].plot(
            np.array(smooth_ic_mse_over_time)[:, i],
            label = 'Constant sig.',
            c = '#707f3b'
            )
        ax[i].set_xlabel('Time (AU)')
        tuftean_axes(ax[i], minx=0, maxx=61, miny=0, maxy=0.03)

    fig.suptitle(f'Network {name}: Mean squared error over simulation timecourse')

    ax[0].set_title(f'{alphabet[alphabet_idx]}i.      tbxta')
    ax[1].set_title(f'{alphabet[alphabet_idx]}ii.     tbx16')
    ax[2].set_title(f'{alphabet[alphabet_idx]}iii.    tbx6')

    # ax[0].legend(position)
    ax[0].set_ylabel('Mean squared error')

    # Add some spacing at the bottom of the figure
    fig.subplots_adjust(bottom=0.2)

    lgd = fig.legend(
        ax[0].get_legend_handles_labels()[1],
        bbox_to_anchor=(0.5, -0.07),
        loc='lower center',
        fontsize="10",
        ncol=3 # Number of columns in the legend
    )

    # add a bit more breathing room around the axes for the frames
    fig.subplots_adjust(top=0.85, bottom=0.3, left=0.2, hspace=0.8)

    fig.patch.set_linewidth(2)
    fig.patch.set_edgecolor('black')


    plt.savefig(f'../graphics/heterogeneity_over_time_{name}.png',
                dpi = 300,
                bbox_inches = 'tight')