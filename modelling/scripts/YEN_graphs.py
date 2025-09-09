from cmath import sin
import os
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd

from scipy.ndimage.filters import uniform_filter1d
import matplotlib.ticker as ticker

from numba import jit, cfunc, njit
from numbalsoda import lsoda, lsoda_sig

##########
# INPUTS #
##########

sim = pd.read_csv('../run_mcmc/many_aget_coordinated_analysis_2/Simulations/Network_9.csv')
netnum = '9'

params_all = np.array(pd.read_csv('../output/new_mcmc/in_progress_simulations_3.csv'))
params = params_all[8]

col_list = ['#ed392b', '#feb441', '#4090c5']

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Lato', 'Tahoma', 'DejaVu Sans',
                               'Lucida Grande', 'Verdana']

plt.rcParams['font.size'] = 12

#######
# scripts to run simulation
#######

#region

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



class Simulate_on_tracks_static:
    '''
    Static simulation.
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
            if df_celltrack.shape[0] == 61:
                df_celltrack['Time_nPSM'] = time_npsm
            else:
                df_celltrack['Time_nPSM'] = df_celltrack['Time']*90/3600/3

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
                # Set signalling to ICs!
                params_to_pass[24:26] =  [B1[0], B2[0]]

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
            if df_celltrack.shape[0] == 61:
                df_celltrack['Time_nPSM'] = time_npsm
            else:
                df_celltrack['Time_nPSM'] = df_celltrack['Time']*90/3600/3

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

#endregion




path = 'Input/List_of_all_cell_tracks.txt'

def run_sim_make_plot(netnum, params):
    # last entry appears garbled
    sim = Simulate_on_tracks(params, list_of_cell_tracks[:-1]).simulation()

    sim = pd.concat(sim)

    static = Simulate_on_tracks(parameters=params,
                                filename = path,
                                which_celltracks=[0, 1903]
                               )

    static.Moving = False
    static = static.simulation()

    static_df = pd.concat(static[0])#Kay's code also outputs other things we don't need
    moving_df = pd.concat(moving[0])


    print(str(netnum+1) + ": finished runing simulation")

    ####################
    # CONSTRUCT OUTPUT #
    ####################

    # now we just want to do some data analysis....
    merged_df = moving_df
    merged_df['g1_static'] = static_df['g1_sim']
    merged_df['g2_static'] = static_df['g2_sim']
    merged_df['g3_static'] = static_df['g3_sim']

    merged_df['g1_diff'] = merged_df['g1_sim'] - merged_df['g1_static']
    merged_df['g2_diff'] = merged_df['g2_sim'] - merged_df['g2_static']
    merged_df['g3_diff'] = merged_df['g3_sim'] - merged_df['g3_static']
    merged_df['total_diff'] = np.sqrt(
        merged_df['g1_diff']**2 +
        merged_df['g2_diff']**2 +
        merged_df['g3_diff']**2
    )
    #and save!
    merged_df.to_csv('Simulations/Network_' + str(netnum+1) + '.csv')

    make_plot(merged_df, str(netnum+1))

# run all simulations in parallel
# may as well use all the CPUs I've been given!
# Parallel(n_jobs=16)(delayed(run_sim_make_plot)(netnum, params) for netnum, params in enumerate(params_all))

##########################
# GRAPH 1
# Pattern output
##########################

#region

def tuftean_axes(ax, minx, maxx, miny, maxy, sim_s):
    # much of this code is cribbed from here
    # https://www.ajnisbet.com/blog/tufte-in-matplotlib#line

    # Remove axis lines.
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    right = False,
    ) # labels along the bottom edge are off

    ax.spines['bottom'].set_bounds(minx, maxx)
    ax.spines['left'].set_bounds(miny, maxy)

    # Adjust lower limits to let data breathe.
    ax.set_xlim([minx - 10, maxx + 10])
    ax.set_ylim([miny - 0.1, maxy + 0.1])

    get_percent_psm_ax(ax, sim_s, minx, maxx)

def get_percent_psm_ax(ax, sim_s, minx, maxx):
    percent_psm = ['0%', '25%', '50%', '75%', '100%']
    scale_factor = (max(sim_s['X'])- min(sim_s['X'])/min(sim_s['X']))
    minx = min(sim_s['X'])
    maxx = max(sim_s['X'])
    abs_position = [minx, (maxx - minx) / 4, (maxx - minx) / 2, (maxx - minx) * (3 / 4), maxx]
    ax.set_xticks(abs_position)
    ax.set_xticklabels(labels = percent_psm)

def plot_embryo(t, sim, cmax = None, s = 3):
    sim_s = sim[sim['Time']==t]

    if not cmax:
        cmax  = np.max([sim['g1_sim'], sim['g2_sim'], sim['g3_sim']])

    # intialization
    # fig.suptitle('Simulation output: Network ' + str(netnum), fontsize = 'x-large')
    for ax in axes.flatten():
        ax.clear()

    ###############
    # plot embryo #
    ###############

    embryo_pointsize = 10

    axes[0, 0].scatter(sim_s['X'], sim_s['Y'], c = sim_s['g1_sim'], cmap = emb_cmap[0], vmin = 0, vmax = cmax, s = embryo_pointsize)
    axes[0, 1].scatter(sim_s['X'], sim_s['Y'], c = sim_s['g2_sim'], cmap = emb_cmap[1], vmin = 0, vmax = cmax, s = embryo_pointsize)
    axes[0, 2].scatter(sim_s['X'], sim_s['Y'], c = sim_s['g3_sim'], cmap = emb_cmap[2], vmin = 0, vmax = cmax, s = embryo_pointsize)

    for ax in axes[0]:
        ax.axis("off")
        ax.axis('scaled')
        ax.set_ylim(min(sim['Y']) - 10, max(sim['Y']) - 10)
        ax.set_xlim(min(sim['X']) - 10, max(sim['X']) - 10)


    axes[0, 0].set_title("tbxta", style='italic')
    axes[0, 1].set_title("tbx16", style='italic')
    axes[0, 2].set_title("tbx6", style='italic')

    ##########################
    # Scatter plot of values #
    ##########################

    text_annotation = [150, 150, 10]
    annotation_from = [170, 170, 35]
    annotation_to = [100, 120, 100]

    genes_groundtruth = []

    for idx, gene in enumerate(['g1', 'g2', 'g3']):
            yhat_1 = uniform_filter1d(sim_s.sort_values("X")[gene], size = 75, mode = 'nearest')
            yhat_2 = uniform_filter1d(sim_s.sort_values("X")[gene + '_sim'], size = 75, mode = 'nearest')
            axes[1, idx].plot(sim_s.sort_values("X")['X'], yhat_1, color = 'black', label = 'ground truth\nrolling average')
            axes[1, idx].plot(sim_s.sort_values("X")['X'], yhat_2, color = col_max[idx], label = 'simulated data\nrolling average')
            axes[1, idx].scatter(sim_s['X'], sim_s[gene + '_sim'], color = col_list[idx],  s = s, label = 'simulated data')
            axes[1, idx].set_ylim(0, 1.7)
            tuftean_axes(axes[1, idx], min(sim['X']), max(sim['X']), 0, cmax, sim_s)
            axes[1, idx].set_title(gene_list[idx], style='italic')
            genes_groundtruth.append(yhat_1)
    axes[1, 0].legend(prop={'size': 9}, frameon=False)
    axes[1, 1].legend(prop={'size': 9}, frameon=False)

fig, axes = plt.subplots(
    2, 3, figsize = (25/2.5, 10/2.5),
    sharex = True, sharey = 'row',
    tight_layout = True,
    gridspec_kw={'height_ratios': [1.1, 1]}
    )

col_list = ['#ed392b', '#feb441', '#4090c5']
col_max = ['#af1117', '#ea6e13', '#0f59a3']
gene_list = ['tbxta', 'tbx16', 'tbx6']
emb_cmap = ['Reds', 'YlOrBr', 'Blues']

# col_list = ['darkgrey', 'darkgrey', 'darkgrey']
# # col_max = ['black', 'black', 'black']
# emb_cmap = ['Greys', 'Greys', 'Greys']

plot_embryo(61, sim)

plt.savefig('YEN_graphics/simulation_pattern.png', dpi = 300)

#endregion


##########################
# GRAPH 2
# autonomous differentation
##########################

#region

sim_t1 = sim[sim['Time'] == 1]
min_x = min(sim_t1['X'])
max_x = max(sim_t1['X'])

sim_t1['x_norm'] = (sim_t1['X'] - min_x) / (max_x - min_x)
data_subset = sim_t1[sim_t1['x_norm'] < 0.2]
data_subset = data_subset.reset_index(drop = True)

t_interval = np.linspace(0, 200, 10100)


list_of_df = []

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

    # # Adjust lower limits to let data breathe.
    # ax.set_xlim([minx-0.05, maxx*1.1])
    # ax.set_ylim([miny-0.05, maxy*1.1])



def run_sim(params):

    trajectories_tbxta = []
    trajectories_tbx16 = []
    trajectories_tbx24 = []

    #for every cell that made the cut
    for i in range(data_subset.shape[0]):
        # g1, g2, g3 are at these indices
        IC = np.array(data_subset.loc[i, ['g1', 'g2', 'g3']])

        params_to_pass = np.empty((26, ))
        params_to_pass[0:24] = params
        params_to_pass[24:26] =  data_subset.loc[i, ['Wnt', 'FGF']].values

        single_cell_trajectory = solve_lsoda(IC, t_interval, params_to_pass)
        trajectories_tbxta.append(single_cell_trajectory[:,0])
        trajectories_tbx16.append(single_cell_trajectory[:,1])
        trajectories_tbx24.append(single_cell_trajectory[:,2])

    return(
        [np.array(trajectories_tbxta),
        np.array(trajectories_tbx16),
        np.array(trajectories_tbx24)])




trajectories_tbxta, trajectories_tbx16, trajectories_tbx24 = run_sim(params)

print(trajectories_tbx24.shape)

# assert True == False
fig.clf()
fig, ax = plt.subplots(figsize = (4,4), tight_layout= True)
sns.stripplot(
    data=[trajectories_tbx24[:, 0], trajectories_tbx24[:, 30], trajectories_tbx24[:, 65], trajectories_tbx24[:, 100]],
    color = col_list[2], size = 2, ax = ax)
# sns.violinplot(data=[trajectories_tbx24[:, 0], trajectories_tbx24[:,30], trajectories_tbx24[:,65], trajectories_tbx24[:,100]], color = "cyan", ax = ax9, inner = None)
ax.set_xticklabels(['0','2','4','6'])


minx, maxx = -0.1, 3.1
miny, maxy = 0, 1.2

tuftean_axes(ax, minx, maxx, miny, maxy)

ax.set_ylabel('tbx6 expression')
ax.set_xlabel('Time since removal from PSM (hrs)')

# ax.set_title('Network 9')
plt.savefig('YEN_graphics/disass_experiment.png', dpi = 300)


fig = plt.figure(figsize = (3,3))
ax = fig.add_subplot(111, projection = '3d')
for cell in range(len(trajectories_tbxta)):
    ax.plot(
        trajectories_tbxta[cell], trajectories_tbx16[cell], trajectories_tbx24[cell],
        c = 'k', linewidth = 1, alpha = 0.1)
ax.scatter3D(
    trajectories_tbxta[:, -1], trajectories_tbx16[:, -1], trajectories_tbx24[:, -1],
    c = 'blue', alpha = 1)
# ax.set_title('Final relative expression, network ' + str(network_number))
plt.savefig('YEN_graphics/disass_experiment_trajectories.png')

fig = plt.figure(figsize = (5, 2.5))
ax1 = fig.add_subplot(121, projection = '3d')
ax2 = fig.add_subplot(122, projection = '3d')

for axs in [ax1, ax2]:
    axs.set_xlim3d(0, 1.5)
    axs.set_ylim3d(0, 1.5)
    axs.set_zlim3d(0, 1.5)
    axs.set_xlabel('tbxta')
    axs.set_ylabel('tbx16')
    axs.set_zlabel('tbx6')
    axs.set_xticks([0, 0.5, 1, 1.5])
    axs.set_yticks([0, 0.5, 1, 1.5])
    axs.set_zticks([0, 0.5, 1, 1.5])
    axs.set_facecolor('white')

for cell in range(len(trajectories_tbxta)):
    if trajectories_tbx24[cell][-1] > 0.4:
        ax1.plot(
            trajectories_tbxta[cell], trajectories_tbx16[cell], trajectories_tbx24[cell],
            c = 'k', linewidth = 1, alpha = 0.1)
        ax1.scatter3D(
            trajectories_tbxta[cell, -1], trajectories_tbx16[cell, -1], trajectories_tbx24[cell, -1],
            c = col_max[2], alpha = 1)
    else:
        ax2.plot(
            trajectories_tbxta[cell], trajectories_tbx16[cell], trajectories_tbx24[cell],
            c = 'k', linewidth = 1, alpha = 0.1)
        ax2.scatter3D(
            trajectories_tbxta[cell, -1], trajectories_tbx16[cell, -1], trajectories_tbx24[cell, -1],
            c = col_max[1], alpha = 1)
# ax.set_title('Final relative expression, network ' + str(network_number))
plt.savefig('YEN_graphics/disass_experiment_trajectories.svg')

#endregion

##########################
# GRAPH 4
# signalling perturbations
##########################

#sim_reversed = pd.read_csv('new_networks/network_9_reversed.csv')

sim_reversed = sim # last minute change...

size_x = 7

size_y = 1.8




genes = ['tbxta', 'tbx16', 'tbx6']

# df = df.groupby('TrackID').apply(get_age)
df = sim_reversed[sim_reversed['Time'] == 61]
x = np.array(list(df.sort_values("X")['X']))
#df = df[df['age'] > 50]

#ax[0, idx].set_title('Network ' + str(network_number))

minx, maxx, miny, maxy = min(df['X']), max(df['X']), 0, 1.2

s = 1
a = 0.5

fig, ax = plt.subplots(1, 3, figsize = (7, 2.5), sharex = True, sharey = True, tight_layout = True)

for idxx, gene in enumerate(['g1', 'g2', 'g3']):
    # tuftean_axes(ax[idxx], minx, maxx, miny, maxy)
    # ax[idxx].scatter(
    #     df['X'], df[gene + '_sim'],
    #     c = 'k',
    #     s = s, alpha = a, label = 'normal signalling')
    yhat_1 = uniform_filter1d(
        df.sort_values("X")[gene + '_sim'],
        size = 50, mode = 'nearest')
    ax[idxx].plot(x, yhat_1, c = 'k')
    ax[idxx].scatter(
        df['X'], df[gene + '_static'],
        c = col_list[idxx],
        s = s, alpha = a, label = 'reversed signalling')
    yhat_1 = uniform_filter1d(
        df.sort_values("X")[gene + '_static'],
        size = 50, mode = 'nearest')
    ax[idxx].plot(x, yhat_1,
        c = col_max[idxx], )
    # ax[idxx].legend(loc='upper center', bbox_to_anchor=(0.5, 0))
    ax[idxx].set_title(genes[idxx] + ' expression')
    ax[idxx].spines['top'].set_visible(False)
    ax[idxx].spines['right'].set_visible(False)
    ax[idxx].tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    right = False,
    ) # labels along the bottom edge are off

    ax[idxx].spines['bottom'].set_bounds(minx, maxx)
    ax[idxx].spines['left'].set_bounds(miny, maxy)
    ax[idxx].set_xlabel('AP position (um)')


ax[0].set_ylabel('Absolute expression (AU)')


plt.savefig('YEN_graphics/reversing_signalling_pattern.png', dpi = 300)

fig, ax = plt.subplots(
    1, 3,
    figsize = (7, 2.2),
    tight_layout = True,
    sharex = True, sharey = True
    )

# from matplotlib import colors
# divnorm=colors.TwoSlopeNorm(vmin=-0.5, vcenter=0., vmax=0.5)

minx, maxx, miny, maxy = min(df['X']), max(df['X']), -0.5, 0.5

for idxx, gene in enumerate(['g1', 'g2', 'g3']):
    ax[idxx].scatter(
        df['X'], #df['Y'],
        df[gene + '_diff'],
        s = s, alpha = a,
        c = col_list[idxx])
    ax[idxx].set_title(genes[idxx] + ' difference')
    # ax[idxx].legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
    ax[idxx].set_title(genes[idxx] + ' expression')
    ax[idxx].spines['top'].set_visible(False)
    ax[idxx].spines['right'].set_visible(False)
    ax[idxx].tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    right = False,
    ) # labels along the bottom edge are off

    ax[idxx].spines['bottom'].set_bounds(minx, maxx)
    ax[idxx].spines['left'].set_bounds(miny, maxy)
    ax[idxx].set_xlabel('AP position (um)')


ax[0].set_ylabel('Normal - perturbed\nsignalling')
ax[1].set_xlabel('AP position (um)')

plt.savefig('YEN_graphics/reversing_signalling_diff.png', dpi = 300)

###
# attr on embryo



fig, ax = plt.subplots(1,1, figsize = (4, 2), tight_layout = True, sharex = True, sharey = True)

df = pd.read_csv('new_networks/classified_ics9.csv')
df['attractor colour'] = col_max[1]

df.loc[df['attractor name'] == 'r', ['attractor colour']] = col_max[2]

ax.scatter(df['X'],
            df['Y'],
            c = df['attractor colour'],
            vmin = 0, vmax = 1.5,
            s = 3,
                #c = classified_unique_attr['attractor name']
            )

ax.set_aspect('equal')
ax.axis('off')

plt.savefig('YEN_graphics/Attractor_on_embryo.png', dpi = 300)



classified_attractors_2 = pd.read_csv('new_networks/classified_attr_9.csv')

WNT_values = np.unique(classified_attractors_2['Wnt'])
FGF_values = np.unique(classified_attractors_2['FGF'])


fig, ax = plt.subplots(1,2, figsize = (4.4, 2.2),
                      sharex = True, sharey = True,
                      tight_layout = True)

fargs = {'interpolation': 'none',
        'extent' : [WNT_values.min(), WNT_values.max(), FGF_values.min(), FGF_values.max()],
        'origin' : 'lower',
        #'cmap': plt.cm.get_cmap('Greys', ics.shape[0]),
        'alpha' : 1
        }

names = ['spirals', 'high_tbx6', 'other'] # I know that there are only 2 attractors for these networks


def extract_counts_of_attractors(df):
    df = df.reset_index(drop = True)
    fgf = df.loc[0, 'FGF']
    wnt = df.loc[0, 'Wnt']

    num_spirals = df[df['attractor name']=='m'].shape[0]
    num_high_tbxta = df[df['attractor name']=='k'].shape[0]
    num_high_tbx6 = df[df['attractor name']=='r'].shape[0]
    num_low_everything = df[df['attractor name']=='y'].shape[0]

    return pd.DataFrame({'FGF': [fgf], 'Wnt': [wnt],
                        'spirals': [num_spirals],
                        'high tbxta': [num_high_tbxta],
                        'high_tbx6': [num_high_tbx6],
                        'other': [num_low_everything]})

classified_attractors_2 = classified_attractors_2.groupby(['FGF', 'Wnt']).apply(extract_counts_of_attractors)
classified_attractors_2 = classified_attractors_2.reset_index(drop=True)

# the replace is to make all 0 objects white
# so it's obvious where a given attractor is never present
# might be better to modify the colormap tbh

ax[0].imshow(
    np.reshape(
        list(classified_attractors_2['other'].replace(0, np.nan)),
        (len(WNT_values), len(WNT_values))
                        ),
                        cmap = emb_cmap[1],
            vmin = 1, vmax = df.shape[0], **fargs)
im = ax[1].imshow(
    np.reshape(list(classified_attractors_2['high_tbx6'].replace(0, np.nan)),
                        (len(WNT_values), len(WNT_values))
                        ),
                        cmap = emb_cmap[2],
            vmin = 1, vmax = df.shape[0], **fargs)


ax[0].set_title('low tbx6')
ax[1].set_title('high tbx6')

plt.savefig('YEN_graphics/attr_counts.png', dpi = 300)