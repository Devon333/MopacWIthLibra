import os
import sys
import math
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
final = 200
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



#sd_basis =read_mopac_num_of_SD_states(params)


# reading basis from mopac output file
sd_basis=read_mopac_SD_config(params)

# reading matrices from h2oTraj/
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
params.update( {  "data_dim":10, "isnap":start,"fsnap":final,"active_space":range(0,10),"get_imag":0,"get_real":1,
                  "data_re_prefix": "tetraceneTraj/tetracene_T_", "data_re_suffix": "_re.txt",
                  #"data_im_prefix": "Tetracene_indo_traj/h2_T_", "data_im_suffix": "_im.txt"
               })
T = data_read.get_data(params)
print("GOT CI_matrix")
params.update( {  "data_dim":10, "isnap":start,"fsnap":final,
                  "data_re_prefix": "tetraceneTraj/tetracene_SD_mid_E_", "data_re_suffix": "_re.txt",
                  "data_im_prefix": "tetraceneTraj/tetracene_SD_mid_E_", "data_im_suffix": "_im.txt"
               })
Esd = data_read.get_data(params)
print("GOT SD_midpoints")
params.update( {  "data_dim":10, "isnap":start,"fsnap":final,
                  "data_re_prefix": "tetraceneTraj/tetracene_CI_mid_E_", "data_re_suffix": "_re.txt",
                  "data_im_prefix": "tetraceneTraj/tetracene_CI_mid_E_", "data_im_suffix": "_im.txt"
               })
E = data_read.get_data(params)
print("GOT CI_midpoints")

raw_basis=get_raw_basis(params)
basis= reindex_basis(raw_basis)



#function that calculates overlap between Slater determinants using libra functions
def SD_ovlp(basis, time_ov_mat):#,phases):
  N, M = len(basis), len(basis)
  res = CMATRIX(N,M)
  for i in range(len(basis)):
      for j in range(len(basis)):
         test = mapping.ovlp_arb(basis[i], basis[j], time_ov_mat,False)
         #test *= phases.get(j,0).real
         res.set(i,j,test)
         #print(test) 
  return res

#function that populates matrices with Slater determinant energies using libra functions
def SD_energy(basis,E_mat):
  N, M = len(basis), len(basis)
  res = CMATRIX(N,M)
  for i in range(len(basis)):
     test = mapping.energy_arb(basis[i],E_mat)
     res.set(i,i,test)
     #print(test)
  #res.show_matrix() 
  return res

#function that populates matrix with CI energies energy usig libra functions
def CI_energy(basis,E_mat):
  N, M = len(basis), len(basis)
  res = CMATRIX(N,M)
  for i in range(len(basis)):
     test = mapping.energy_arb(basis[i],E_mat)
     res.set(i,i,test)
     #print(test)
  #res.show_matrix() 
  return res

def Mat_avg(cmatList):
    cols = cmatList[0].num_of_cols
    matAvg = np.zeros((cols,cols))
    steps = len(cmatList)
    for row in range(cols): 
        for col in range(cols):
            sumNacs = 0
            for time in range(steps):
                if row != col:
                    sumNacs += 27.2114*abs(cmatList[time].get(row,col))* 1000
            avg = sumNacs/steps
            print(avg)
            matAvg[row,col]=avg
    return matAvg




#applying phase corrections to Kohn-Sham orbitals using libra functions
step3.apply_orthonormalization_general(S, tim_ov)
ks_phases=step3.apply_phase_correction_general(tim_ov)


dt = 1*units.fs2au
start_time = params["isnap"]
res_dir="tetraceneTraj"
nsteps = final-start
print("Outputting the CI data to the res directory..." )
# Printing matrices in results directory

for step in range(nsteps-1):
    S[step].real().show_matrix("%s/S_ks_%d_re" % (res_dir, int(step)))
    tim_ov[step].real().show_matrix("%s/St_ks_%d_re" % (res_dir, int(step)))

# Make the Hvib in the many-body basis
Hvib = []
ks_hvib = None
for step in range( len(tim_ov) ):
    #Calculating CI NACs
    ks_nacs = (  0.5j / dt ) * CMATRIX ( ( tim_ov[step] - tim_ov[step].H() ).real() )    
    #Creating KS Vibrational Hamiltonian
    ks_hvib = E[step] - ks_nacs
    Hvib.append( ks_hvib )
    ks_hvib.real().show_matrix("%s/Hvib_ks_%d_re" % (res_dir, int( step )))
    ks_hvib.imag().show_matrix("%s/Hvib_ks_%d_im" % (res_dir, int( step )))

Hvib.append( ks_hvib)    # appending the last element twice to make it nsteps
#Hvibsd.append( sd_hvib )




ntraj  = len(Hvib)
nsteps = len(Hvib)
nCIs   = Hvib[0].num_of_cols

# Make a list for the SD energies and populate it
KS_energy = []
md_time = list( range(nsteps) )

for ks_index in range( nCIs ):
    KS_energy.append( [] )
    for step in range( nsteps ):        
        En = Hvib[step].get( ks_index, ks_index ).real 
        E0 = Hvib[step].get( 0, 0 ).real
        KS_energy[ ks_index ].append( En - E0 )
        
KS_energy = np.array( KS_energy  )
md_time   = np.array( md_time )

# Functions to compute the time-averaged SD NACs and make a list of them
ks_res= Mat_avg(Hvib)
sd_tNACs = []

#for i in range(nCIs):
#    ci_tNACs.append( [] )
#    sd_tNACs.append( [] )
#    for j in range(nCIs):
#        #ci_tNACs[i].append( ci_res.get(i,j).imag * 1000.0 / units.ev2Ha )        
#        sd_tNACs[i].append( sd_res.get(i,j).imag * 1000.0 /units.ev2Ha  )        
#        ci_tNACs[i].append( ci_res.get(i,j).imag * 1000.0 /units.ev2Ha  )        
#sd_tNACs = np.array(sd_tNACs)
#ci_tNACs = np.array(ci_tNACs)
#
#
nstates = 6
orbitals = ['H-2','H-1','H','L','L+1','L+2']#["GS","42 -> 43","42 -> 44","42 -> 45","41 -> 43","41 -> 44","41 -> 45","40 -> 43","40 -> 44","40 -> 45"]
state_label = ['H-2','H-1','H','L','L+1','L+2'] #["GS","2","3","4","5","6","7","8","9","10"]
state_label_x = ['H-2','H-1','H','L','L+1','L+2'] #["GS","H->L","H->L+1","H->L+2","H-1->L","H-1->L+1","H-1->L+2","H-2->L","H-2->L+1","H-2->L+2"]
## Figures - Time-Dependent CI Data
#plt.figure(num=None, figsize=(6,6), dpi=100, edgecolor='black', frameon=True)
#plt.title('INDO CI Energies', fontsize=25)
#plt.xlabel('Time, fs')
#plt.xticks(fontsize=15)
#plt.yticks(fontsize=15)
#for state in range(1,nstates ):    
#    plt.plot(md_time, CI_energy[state], label=state_label[state], linewidth=2)
#plt.tight_layout()
#plt.savefig("CI_energytime.png")

fig=plt.figure(num=None, figsize=(6,6), dpi=100, edgecolor='black', frameon=True)
plt.title('INDO KS NACs', fontsize=25)
ax = plt.subplot(111)
ax.set_xticks(np.arange(len(orbitals)))
ax.set_yticks(np.arange(len(orbitals)))
ax.set_xticklabels(state_label,fontsize=15)
ax.set_yticklabels(state_label,fontsize=15)
cb=ax.imshow(ks_res, origin='lower', cmap='plasma', interpolation='nearest')
plt.xticks(rotation=90)
plt.colorbar(cb,label="meV")

plt.tight_layout()
plt.savefig("KS_avgNAC.png")


# Figures - Time-Dependent SD Data
plt.figure(num=None, figsize=(6, 6), dpi=100, edgecolor='black', frameon=True)
ax=plt.subplot(1,1,1)
plt.title('INDO Energies', fontsize=25)
plt.xlabel('Time, fs')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
for state in range(1, nstates ):    
    ax.plot(md_time, KS_energy[state], label=orbitals[state], linewidth=1)
plt.tight_layout()
plt.savefig("KS_energytime.png")
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
plt.savefig('KS_legend.png')
