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


parameter_values = pd.read_csv('../median_networks_for_clusters.csv')


with open("../Input/List_of_all_cell_tracks.txt", "rb") as fp:   # Unpickling
    list_of_cell_tracks = pickle.load(fp)[0:1903]

# print(list_of_cell_tracks[0])

#############
# FUNCTIONS #
#############

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





class Simulate_on_tracks_static:
    '''
    static simulation.
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

#endregion


# Run MAP network.

all_params = pd.read_csv('../ALL_NETWORKS.csv')

all_params = all_params.sort_values(by = 'll', ascending=False).reset_index(drop = True)

name = 'MAP'
params=all_params.iloc[0, 1:25]

print(all_params['ll'][0])

print(name)
print(params.shape)


##################
# RUN SIMULATION #
##################

# last entry appears garbled
sim = Simulate_on_tracks(params, list_of_cell_tracks[:-1]).simulation()

sim = pd.concat(sim)

static_sim = pd.concat(
    Simulate_on_tracks_static(params,
                                list_of_cell_tracks[:-1]).simulation()
)

# now we just want to do some data analysis....
merged_df = sim
merged_df['g1_static'] = static_sim['g1_sim']
merged_df['g2_static'] = static_sim['g2_sim']
merged_df['g3_static'] = static_sim['g3_sim']

merged_df['g1_diff'] = merged_df['g1_sim'] - merged_df['g1_static']
merged_df['g2_diff'] = merged_df['g2_sim'] - merged_df['g2_static']
merged_df['g3_diff'] = merged_df['g3_sim'] - merged_df['g3_static']
merged_df['total_diff'] = np.sqrt(
    merged_df['g1_diff']**2 +
    merged_df['g2_diff']**2 +
    merged_df['g3_diff']**2
)
#and save!
merged_df.to_csv(f"../output/{name}.csv", sep = ";")



genes = ['tbxta', 'tbx16', 'tbx6']

# df = df.groupby('TrackID').apply(get_age)
df = merged_df[merged_df['Time'] == 61]
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


plt.savefig(f'../graphics/reversing_signalling_pattern_{name}.png', dpi = 300)

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

plt.savefig(f'../graphics/reversing_signalling_diff_{name}.png', dpi = 300)


import string
alphabet = list(string.ascii_uppercase)


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


    ##################
    # RUN SIMULATION #
    ##################

    # last entry appears garbled
    sim = Simulate_on_tracks(params, list_of_cell_tracks[:-1]).simulation()

    sim = pd.concat(sim)

    static_sim = pd.concat(
        Simulate_on_tracks_static(params,
                                  list_of_cell_tracks[:-1]).simulation()
    )

    # now we just want to do some data analysis....
    merged_df = sim
    merged_df['g1_static'] = static_sim['g1_sim']
    merged_df['g2_static'] = static_sim['g2_sim']
    merged_df['g3_static'] = static_sim['g3_sim']

    merged_df['g1_diff'] = merged_df['g1_sim'] - merged_df['g1_static']
    merged_df['g2_diff'] = merged_df['g2_sim'] - merged_df['g2_static']
    merged_df['g3_diff'] = merged_df['g3_sim'] - merged_df['g3_static']
    merged_df['total_diff'] = np.sqrt(
        merged_df['g1_diff']**2 +
        merged_df['g2_diff']**2 +
        merged_df['g3_diff']**2
    )
    #and save!
    merged_df.to_csv(f"../output/{name}.csv", sep = ";")



    genes = ['tbxta', 'tbx16', 'tbx6']

    # df = df.groupby('TrackID').apply(get_age)
    df = merged_df[merged_df['Time'] == 61]
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

    ax[0].set_title(f'{alphabet[alphabet_idx]}i.      tbxta')
    ax[1].set_title(f'{alphabet[alphabet_idx]}ii.     tbx16')
    ax[2].set_title(f'{alphabet[alphabet_idx]}iii.    tbx6')

    # add a bit more breathing room around the axes for the frames
    fig.subplots_adjust(top=0.85, bottom=0.3, left=0.2, hspace=0.8)
    fig.patch.set_linewidth(2)
    fig.patch.set_edgecolor('black')


    ax[0].set_ylabel('Absolute expression (AU)')
    fig.suptitle(f'Network {name}')


    plt.savefig(f'../graphics/reversing_signalling_pattern_{name}.png', dpi = 300)

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
        # ax[idxx].set_title(genes[idxx] + ' difference')
        # ax[idxx].legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
        # ax[idxx].set_title(genes[idxx] + ' expression')
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

    ax[0].set_title(f'{alphabet[alphabet_idx+6]}i.      tbxta')
    ax[1].set_title(f'{alphabet[alphabet_idx+6]}ii.     tbx16')
    ax[2].set_title(f'{alphabet[alphabet_idx+6]}iii.    tbx6')

    ax[0].set_ylabel('Normal - perturbed\nsignalling')
    ax[1].set_xlabel('AP position (um)')

    fig.suptitle(f'Network {name}')

    # add a bit more breathing room around the axes for the frames
    fig.subplots_adjust(top=0.85, bottom=0.3, left=0.2, hspace=0.8)

    fig.patch.set_linewidth(2)
    fig.patch.set_edgecolor('black')

    plt.savefig(f'../graphics/reversing_signalling_diff_{name}.png',
                dpi = 300,
                bbox_inches = 'tight')

