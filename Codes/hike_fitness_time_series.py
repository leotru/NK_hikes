# ------------------------------------------------------------------------
# Project NK-hikes/
# ------------------------------------------------------------------------
# hikes_fitness_time_series.py / version 0.1 
# Program to plot the fitness time series.
# ------------------------------------------------------------------------
# Created:  2021-01-19 
# by Leonardo Trujillo (leonardo.trujillo@inria.fr)
# You may use, share, or modify this file freely
# ------------------------------------------------------------------------
# GLOBAL NOMENCLATURE:
# Compile nk_walk.cpp and run the simulation as 
# ./nk_walk  -n $N -k $K -a $A -e $epi -snk $nk_seed -swlk $wlk_seed 
#            -t $T -m $M
# where, 
# N:........Genome length.
# K:........Number of epistatic interactions (0 <= K < N).
# A:........Alphabet size (default = 2).
# epi:......Type of epistatic interactions ('ADJ' (default) or 'RND').
# nk_seed:..PRNG seed value for the landscape (default = -1).
# wlk_seed:.PRNG seed value for the random walk (default = -1).
# T:........Number of steps in the random walk (default = 10000).
# M:........Mutation type (float) (0.0 <= M <= 1.0)," 
# ..........such that M = 0.0: point mutations only (default)," 
# ....................M = 1.0: inversions only."
# ------------------------------------------------------------------------
# IMPORTAN: The raw data comes from the simulations performed with nk_walk. 
# Before to run the present code, please verify the paths where 
# this program and your data are placed. For the data, for example, we 
# labeled the outcomes from the simulations as:
# N50_K20_T200000_epiADJ_M0.0_sample1.csv
# In this format the label sample is "synonymous" with hike
# ------------------------------------------------------------------------
import os
import glob
import random
import pandas as pd
import matplotlib.pyplot as plt
# ------------------------------------------------------------------------
# Parameters
Number_of_time_series_to_print = 4
N = 50
K = 2
T = 2000000
epi = 'ADJ'
M = 0.01
hike = '*'

random.seed(10)

# ------------------------------------------------------------------------
# Paths and file name
data_file = '../Data/BurstsCandidates/**/' # the /**/ is used to make glob search recursively
file_name = data_file + ("N%s_K%s_T%s_epi%s_M%s_sample%s.csv" % (N, K, T, epi, M, hike))

files_concerned = [name for name in  glob.iglob( file_name, recursive = True)]

# ------------------------------------------------------------------------
# Read and load the *.csv file with the raw data
def add_a_time_serie(ax, file_path):
    df = pd.read_csv(file_path, sep=';',header=None).values
    
    time = df[:,1]
    fitness = df[:,2]
    mut_type = df[:,5]
    total_time = df.shape[0]
    
    head, sample_name = os.path.split(file_path)
    for i in range(total_time):
        if int(mut_type[i]) == 1:
            col = 'red'
            mark = 'o'
        else:
            col = 'blue'
            mark = 's'

        ax.plot(time[i], fitness[i], marker = mark,color = col, markersize = 4.5, fillstyle='none')
    ax.plot(time, fitness, '-', label = sample_name)

#----------------------------------------------------------------------------
# Plot the fitness time series
fig, ax = plt.subplots()

k = 0 
while  (k < Number_of_time_series_to_print) and (files_concerned != []) :
    n = len(files_concerned)
    random_index = random.randint(0,n-1)
    target_file = files_concerned[random_index]
    files_concerned.remove(target_file)
    
    add_a_time_serie(ax, target_file)
    
    k+= 1
    
ax.set_xscale('log')
ax.grid(which = 'minor',axis='both',  alpha = 0.3)
ax.grid(which = 'major', axis='both', alpha = 0.7)        
ax.set_xlabel('generations',fontsize = 16)
ax.set_ylabel('fitness',fontsize = 16) 
ax.set_title = ("N%s_K%s_epi%s_M%s_hike%s" % (N, K, epi, M, hike))  

fig_file_name = ("figFitnessTimeSeries_N%s_K%s_epi%s_M%s_hike%s.png" % (N, K, epi, M, hike))
plt.legend()

plt.savefig(fig_file_name, dpi=300, bbox_inches='tight')

#plt.show()
        
# ------------------------------------------------------------------------
# Move the figure
os.system(f"mv {fig_file_name} ../Figures")

#-----------------------------------------------------------------------------
# Bye!  
#-----------------------------------------------------------------------------
