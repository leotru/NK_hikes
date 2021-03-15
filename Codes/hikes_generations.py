# ------------------------------------------------------------------------
# Project NK-hikes/
# ------------------------------------------------------------------------
# hike_generations.py / version 0.1 
# This program draws the scatter plot of hikes (and fitness) as a 
# function of hikes
# ------------------------------------------------------------------------
# Created: 2021-02-10 
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
#title = '(d)'

# ------------------------------------------------------------------------
# Paths and file name
file_hike = 'Hikes1'

data_file_name = ("../Data/%s/N%s_K%s_T%s_epi%s_M%s/" % (file_hike,N, K, T, epi, M))
data_file = f"./{data_file_name}"

# ------------------------------------------------------------------------
# Lists initializations
fitness = ([])
time = []
max_time = []
mut_type = []

# ------------------------------------------------------------------------
# Read and load the *.csv file with the raw data
for i in range (realizations):
    sample = i
    file_name = ("N%s_K%s_T%s_epi%s_M%s_sample%s.csv" % (N, K, T, epi, M, sample + 1))
    #df = pd.read_csv(file_name, sep = ';',header = None).values
    df = pd.read_csv(f"{data_file}{file_name}", sep=';',header=None).values
    max_time.append([i,df.shape[0]])    
    for j in range(df.shape[0]):
        fitness.append([i,df[j,2]])
        time.append([i,df[j,1]])
        mut_type.append([i,df[j,5]])

# ------------------------------------------------------------------------
# Plot generations vs hikes
t = np.array(time)
mut = np.array(mut_type)
max_t = np.array(max_time)

for k in range(int(t.shape[0])):
    if mut[k,1]==1:
        col = 'red'
        mark = 'o'
    else:
        col = 'blue'
        mark = 's'        
    plt.plot(t[k,0]+1,t[k,1],marker = mark,color = col, markersize = 2.5,fillstyle='none')

plt.grid(which = 'minor',axis='both',  alpha = 0.3)
plt.grid(which = 'major', axis='both', alpha = 0.7)        
plt.yscale('log')
plt.xlabel('hikes',fontsize = 16)
plt.ylabel('generations',fontsize = 16)

#plt.title(f"{title}", fontsize = 14)
plt.title(r"$N$" f" = {N}," r" $K$" f" = {K}," r" $epi$" f" = {epi}," r" $p$" f" = {M}," r" hikes" f" = {realizations} ", fontsize = 14)

fig_file_name = ("figGenerations_N%s_K%s_T%s_epi%s_M%s_hikes%s.png" % (N, K, T, epi, M, realizations))
plt.savefig(fig_file_name, dpi=300, bbox_inches='tight')

#plt.show()

# ------------------------------------------------------------------------
# Move the figure
os.system(f"mv {fig_file_name} ../Figures")

# ------------------------------------------------------------------------
# “Think of a ball of steel as large as the world, and a fly alighting 
# on it once every million years. When the ball of steel is rubbed away 
# by the friction, eternity will not even have begun.” 
# – The Picturegoers by David Lodge.
# ------------------------------------------------------------------------
# Bye!
# ------------------------------------------------------------------------
#    