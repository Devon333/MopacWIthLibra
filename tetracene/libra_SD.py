import os
import sys
import math
sys.path.append("../")
from libra_mopac import *

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
    

from libra_py import units as units
from libra_py import influence_spectrum as infsp
from libra_py import data_visualize, data_conv, data_read, data_stat, data_outs
from libra_py.workflows.nbra import mapping, step2_many_body, step3, step4, step3_many_body
import matplotlib.pyplot as plt
import numpy as np


plt.rc('axes', titlesize=24)      # fontsize of the axes title
plt.rc('axes', labelsize=20)      # fontsize of the x and y labels
plt.rc('legend', fontsize=20)     # legend fontsize
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels

plt.rc('figure.subplot', left=0.2)
plt.rc('figure.subplot', right=0.95)
plt.rc('figure.subplot', bottom=0.13)
plt.rc('figure.subplot', top=0.88)

colors = {}

colors.update({"11": "#8b1a0e"})  # red       
colors.update({"12": "#FF4500"})  # orangered 
colors.update({"13": "#B22222"})  # firebrick 
colors.update({"14": "#DC143C"})  # crimson   

colors.update({"21": "#5e9c36"})  # green
colors.update({"22": "#006400"})  # darkgreen  
colors.update({"23": "#228B22"})  # forestgreen
colors.update({"24": "#808000"})  # olive      

colors.update({"31": "#8A2BE2"})  # blueviolet
colors.update({"32": "#00008B"})  # darkblue  

colors.update({"41": "#2F4F4F"})  # darkslategray

clrs_index = ["11", "21", "31", "41", "12", "22", "32", "13","23", "14", "24"]

# setting number of files used to calculate NACs
start = 0
final = 1000

# setting time step
dt = 1*units.fs2au

# setting up params dict that is use by libra functions to read matrices
params = {}
if os.path.exists('*.png'):
  os.remove('*.png')
path = os.getcwd()
params["data_set_paths"] = [ "tetraceneTraj/tetracene_" ]  # where the step2 data is located
params["file_prefix"] = "tetracene_" # filename
params["data_dim"] = 6 # the number of rows/columsn in E_ks files
params["active_space"] = range(0,6) # alpha channel only here 
params["isnap"] = 0 # numbering of the first available file
params["fsnap"] = final # number of the final step

params["orbital_indices"]     = list( range(0,6) )    # orbital indices from waveplot
params["logfile_directory"]   = "tetraceneTraj/"
params["es_software"]         = "mopac"
params["isUKS"]               = 0
params["tolerance"]           = 0.0  # set this to 0.0 for cp2k
params["number_of_states"]    = 3 




# function that reads basis from mopac output files
sd_basis=read_mopac_SD_config(params)
print(sd_basis)
# reading matrices from tetraceneTraj/
params.update( {  "data_dim":6, "isnap":start,"fsnap":final,"active_space":range(0,6),"get_imag":0,"get_real":1,
                  "data_re_prefix": "tetraceneTraj/tetracene_S_", "data_re_suffix": "_re.txt",
                  "data_im_prefix": "tetraceneTraj/tetracene_S_", "data_im_suffix": "_im.txt"
               })
S = data_read.get_data(params)
print("GOT S")

params.update( {  "data_dim":6, "isnap":start,"fsnap":final,"active_space":range(0,6),"get_imag":0,"get_real":1,
                  "data_re_prefix": "tetraceneTraj/tetracene_St_", "data_re_suffix": "_re.txt",
                  "data_im_prefix": "tetraceneTraj/tetracene_St_", "data_im_suffix": "_im.txt"
               })
tim_ov = data_read.get_data(params)
print("GOT St")

params.update( {  "data_dim":10,'active_space':range(0,10),"isnap":start,"fsnap":final,
                  "data_re_prefix": "tetraceneTraj/tetracene_SD_mid_E_", "data_re_suffix": "_re.txt",
                  "data_im_prefix": "tetraceneTraj/tetracene_SD_mid_E_", "data_im_suffix": "_im.txt"
               })
Esd = data_read.get_data(params)
print("GOT SD_midpoints")
Esd[0].show_matrix()

raw_basis=get_raw_basis(params)
print(raw_basis)
basis= reindex_basis(raw_basis)

print(basis)

#exit()



#applying phase corrections to Kohn-Sham orbitals using libra functions
#step3.apply_orthonormalization_general(S, tim_ov)
#ks_phases=step3.apply_phase_correction_general(tim_ov)

Stsd=[]
Ssd=[]

#building lists of overlap and time overlap matrices using SD_ovlp
for time in range(final-start):
    Stsd.append(SD_ovlp(basis, tim_ov[time]))
    Ssd.append(SD_ovlp(basis, S[time]))
    
print("this is Stsd ", Stsd[0].show_matrix() )    
print("this is Ssd ", Ssd[0].show_matrix() )    

#step3.apply_orthonormalization_general( Ssd, Stsd )
#sd_phases=step3.apply_phase_correction_general( Stsd )
#print("the phases",sd_phases[0].get(1,0).real)



CIs=len(basis)
print(f"length of CI basis {CIs}")


dt = 1*units.fs2au
start_time = params["isnap"]
res_dir="tetraceneTraj"
nsteps = final-start
print("Outputting the CI data to the res directory..." )
# Printing matrices in results directory

for step in range(nsteps-1):
    Ssd[step].real().show_matrix("%s/S_sd_%d_re" % (res_dir, int(step)))
    Stsd[step].real().show_matrix("%s/St_sd_%d_re" % (res_dir, int(step)))

# Make the Hvib in the many-body basis
Hvib = []
sd_hvib = None
for step in range( len(Stsd) ):
    #Calculating CI NACs
    sd_nacs = (  0.5j / dt ) * CMATRIX ( ( Stsd[step] - Stsd[step].H() ).real() )    
    #Creating KS Vibrational Hamiltonian
    sd_hvib = Esd[step] - sd_nacs
    Hvib.append( sd_hvib )
    sd_hvib.real().show_matrix("%s/Hvib_sd_%d_re" % (res_dir, int( step )))
    sd_hvib.imag().show_matrix("%s/Hvib_sd_%d_im" % (res_dir, int( step )))
    sd_hvib.show_matrix()
    #exit()
Hvib.append( sd_hvib)    # appending the last element twice to make it nsteps
#Hvibsd.append( sd_hvib )




ntraj  = len(Hvib)
nsteps = len(Hvib)
nCIs   = Hvib[0].num_of_cols
nSDs   = Stsd[0].num_of_cols
print(Esd[0].num_of_cols)
print(Hvib[0].num_of_cols)
print(nSDs)
print(nCIs)
print(ntraj)
print(nsteps)
#exit()
# Make a list for the SD energies and populate it
SD_energy = []
md_time = list( range(nsteps) )

for sd_index in range( nCIs ):
    SD_energy.append( [] )
    for step in range( nsteps ):        
        En = Hvib[step].get( sd_index, sd_index ).real 
        E0 = Hvib[step].get( 0, 0 ).real
        SD_energy[ sd_index ].append( En - E0 )
        
SD_energy = np.array( SD_energy  )
md_time   = np.array( md_time )

# Functions to compute the time-averaged SD NACs and make a list of them
sd_res= Mat_avg(Hvib)
sd_tNACs = []


nstates = 10
orbitals = ["GS","42 -> 43","42 -> 44","42 -> 45","41 -> 43","41 -> 44","41 -> 45","40 -> 43","40 -> 44","40 -> 45"]
state_label = ["GS","2","3","4","5","6","7","8","9","10"]#['H-2','H-1','H','L','L+1','L+2'] 
state_label_x = ["GS","H->L","H->L+1","H->L+2","H-1->L","H-1->L+1","H-1->L+2","H-2->L","H-2->L+1","H-2->L+2"]


fig=plt.figure(num=None, figsize=(6,6), dpi=100, edgecolor='black', frameon=True)
plt.title('INDO SD NACs', fontsize=25)
ax = plt.subplot(111)
ax.set_xticks(np.arange(len(orbitals)))
ax.set_yticks(np.arange(len(orbitals)))
ax.set_xticklabels(orbitals,fontsize=15)
ax.set_yticklabels(orbitals,fontsize=15)
cb=ax.imshow(sd_res, origin='lower', cmap='plasma', interpolation='nearest')
plt.xticks(rotation=90)
plt.colorbar(cb,label="meV")

plt.tight_layout()
plt.savefig("SD_avgNAC.png")


# Figures - Time-Dependent SD Data
plt.figure(num=None, figsize=(6, 6), dpi=100, edgecolor='black', frameon=True)
ax=plt.subplot(1,1,1)
plt.title('INDO SD Energies', fontsize=25)
plt.xlabel('Time, fs')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
for state in range(1, nstates ):    
    ax.plot(md_time, SD_energy[state], label=orbitals[state], linewidth=1)
plt.tight_layout()
plt.savefig("SD_energytime.png")
figsize=(5,5)
fig_leg=plt.figure(figsize=figsize)
ax_leg = fig_leg.add_subplot(111)
ax_leg.legend(*ax.get_legend_handles_labels(),loc='center',frameon=False)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)
ax_leg.xaxis.set_visible(False)
ax_leg.yaxis.set_visible(False)
ax_leg.set_frame_on(False)
plt.savefig('SD_legend.png')
