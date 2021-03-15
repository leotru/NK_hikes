# ------------------------------------------------------------------------
# Project NK-hikes/
# ------------------------------------------------------------------------
# hikes_final_time.py / version 0.1 
# This program plots the histogram of the final fitness per hike  
# ------------------------------------------------------------------------
# Created: 2021-03-10 
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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# ------------------------------------------------------------------------
# Parameters
N = 50
K = 20
T = 200000
epi = 'ADJ'
M = 0.0
realizations = 100
#title = '(e)'

# ------------------------------------------------------------------------
# Paths and file name
file_hike = 'Hikes1'

data_file_name = ("../Data/%s/N%s_K%s_T%s_epi%s_M%s/" % (file_hike,N, K, T, epi, M))
data_file = f"./{data_file_name}"

# ------------------------------------------------------------------------
# Lists initializations
fitness = []
time = []
max_time = []
mut_type = []
hike_list = []
last_genome = []
# ------------------------------------------------------------------------
# Read and load the *.csv file with the raw data
for i in range (realizations):
    sample = i    
    file_name = ("N%s_K%s_T%s_epi%s_M%s_sample%s.csv" % (N, K, T, epi, M, 
                                                         sample + 1))
    #df = pd.read_csv(file_name, sep=';',header=None).values
    df = pd.read_csv(f"{data_file}{file_name}", sep=';',header=None).values
    max_time.append([i,np.max(df[:,1])])
    
    for j in range(df.shape[0]):
        #---------------------------------------------------
        # Hike function data:
        # hike number | time | fitness | mut size | mut_type
        hike_list.append([i, df[j,1], df[j,2],abs(df[j,3]-df[j,4]), df[j,5]])
        time.append([i,df[j,1]])
    last_genome.append([i,df[j,0]])
    fitness.append([i,df[j,2]])

# ------------------------------------------------------------------------
# Visualise the histogram of final times by hike
max_t = np.array(max_time)

#plt.title(f"{title}", fontsize = 14)
plt.title(r"$N$" f" = {N}," r" $K$" f" = {K}," r" $epi$" f" = {epi}," r" $p$" f" = {M}," r" hikes" f" = {realizations} ", fontsize = 14)

plt.hist(max_t[:,1],bins=[2**k for k in range(20) ])

plt.xscale('log')
plt.xlim([1, 2**19])  
plt.xlabel('final-time by hike',fontsize = 14)
plt.grid(which = 'minor',axis='both',  alpha = 0.3)
plt.grid(which = 'major', axis='both', alpha = 0.7)  

fig_file_name = ("figFinalTimes_N%s_K%s_T%s_epi%s_M%s_hikes%s.png" % (N, K, T, epi, M, realizations))
plt.savefig(fig_file_name, dpi=300, bbox_inches='tight')

#plt.show()

# ------------------------------------------------------------------------
# Move the figure
os.system(f"mv {fig_file_name} ../Figures")

# ------------------------------------------------------------------------
# Bye!
# ------------------------------------------------------------------------
      