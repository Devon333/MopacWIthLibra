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
final = 98
# setting time step
dt = 1*units.fs2au
# setting up params dict that is use by libra functions to read matrices
params = {}
if os.path.exists('*.png'):
  os.remove('*.png')
path = os.getcwd()
params["data_set_paths"] = [ "h2oTraj/h2o_" ]  # where the step2 data is located
params["file_prefix"] = "h2o_" # filename
params["data_dim"] = 4 # the number of rows/columsn in E_ks files
params["active_space"] = range(0,4) # alpha channel only here 
params["isnap"] = 0 # numbering of the first available file
params["fsnap"] = final # number of the final step

params["orbital_indices"]     = list( range(0,4) )    # orbital indices from waveplot
params["logfile_directory"]   = "h2oTraj/"
params["es_software"]         = "mopac"
params["isUKS"]               = 0
params["tolerance"]           = 0.0  # set this to 0.0 for cp2k
params["number_of_states"]    = 3 



#sd_basis =read_mopac_num_of_SD_states(params)


# reading basis from mopac output file
sd_basis=read_mopac_SD_config(params)

# reading matrices from h2oTraj/
params.update( {  "data_dim":4, "isnap":start,"fsnap":final,"active_space":range(0,4),"get_imag":0,"get_real":1,
                  "data_re_prefix": "h2oTraj/h2o_S_", "data_re_suffix": "_re.txt",
                  "data_im_prefix": "h2oTraj/h2o_S_", "data_im_suffix": "_im.txt"
               })
S = data_read.get_data(params)
print("GOT S")

params.update( {  "data_dim":4, "isnap":start,"fsnap":final,"active_space":range(0,4),"get_imag":0,"get_real":1,
                  "data_re_prefix": "h2oTraj/h2o_St_", "data_re_suffix": "_re.txt",
                  "data_im_prefix": "h2oTraj/h2o_St_", "data_im_suffix": "_im.txt"
               })
tim_ov = data_read.get_data(params)
print("GOT St")
params.update( {  "data_dim":5, "isnap":start,"fsnap":final,"active_space":range(0,5),"get_imag":0,"get_real":1,
                  "data_re_prefix": "h2oTraj/h2o_T_", "data_re_suffix": "_re.txt",
                  #"data_im_prefix": "Tetracene_indo_traj/h2_T_", "data_im_suffix": "_im.txt"
               })
T = data_read.get_data(params)
print("GOT CI_matrix")
params.update( {  "data_dim":5, "isnap":start,"fsnap":final,
                  "data_re_prefix": "h2oTraj/h2o_SD_mid_E_", "data_re_suffix": "_re.txt",
                  "data_im_prefix": "h2oTraj/h2o_SD_mid_E_", "data_im_suffix": "_im.txt"
               })
Esd = data_read.get_data(params)
print("GOT SD_midpoints")
params.update( {  "data_dim":5, "isnap":start,"fsnap":final,
                  "data_re_prefix": "h2oTraj/h2o_CI_mid_E_", "data_re_suffix": "_re.txt",
                  "data_im_prefix": "h2oTraj/h2o_CI_mid_E_", "data_im_suffix": "_im.txt"
               })
E = data_read.get_data(params)
print("GOT CI_midpoints")

raw_basis=get_raw_basis(params)
basis= reindex_basis(raw_basis)



#function that calculates overlap between Slater determinants using libra functions
def SD_ovlp(basis, time_ov_mat):
  N, M = len(basis), len(basis)
  res = CMATRIX(N,M)
  for i in range(len(basis)):
      for j in range(len(basis)):
         test = mapping.ovlp_arb(basis[i], basis[j], time_ov_mat,False)
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



Stsd=[]
Ssd=[]
Sci=[]
Stci=[]

#applying phase corrections to Kohn-Sham orbitals using libra functions
step3.apply_orthonormalization_general(S, tim_ov)
step3.apply_phase_correction_general(tim_ov)

#building lists of overlap and time overlap matrices using function built above
for time in range(final-start):
    Stsd.append(SD_ovlp(basis, tim_ov[time]))
    Ssd.append(SD_ovlp(basis, S[time]))

print("this is Stsd ", Stsd[0].show_matrix() )    
print("this is Ssd ", Ssd[0].show_matrix() )    



T[0].show_matrix()
SD2CI=[]
T_co=[]
CIs=len(basis)
#print(f"length of CI basis {CIs}")
#exit()


# Building Transformation/CI matrix SD2CI
for time in range(final-start):
    SD2CI.append(CMATRIX(CIs,CIs))
    for col in range(CIs):
        for row in range(CIs):
            val= (T[time].get(row,col))#*(1.0+0.0j)
            SD2CI[time].set(row,col,val)
            print(f'this is a value put in the SD2CI matrix {val}')
    # normalizing rows of Transformation matrix
    for row in range(CIs):
        norm=0.0
        for col in range(CIs):
            norm += abs(SD2CI[time].get(row,col))**2
        norm = 1.0/math.sqrt(norm)
        SD2CI[time].scale(-1,col,norm*(1.0+0.0j))
#calculating CI overlap and time overlaps
for time in range(len(Ssd)-1):
    Stci.append( SD2CI[time].H() * Stsd[time] * SD2CI[time+1])    
    Stci[time].show_matrix()
    Sci.append( SD2CI[time].H() * Ssd[time] * SD2CI[time])    
    

step3.apply_orthonormalization_general( Ssd, Stsd )
#step3.apply_phase_correction_general( Stsd )

step3.apply_orthonormalization_general( Sci, Stci )
#step3.apply_phase_correction_general( Stci )



dt = 1*units.fs2au
start_time = params["isnap"]
res_dir="h2oTraj"
nsteps = final-start
print("Outputting the CI data to the res directory..." )
# Printing matrices in results directory
for step in range(nsteps-1):
    Ssd[step].real().show_matrix("%s/S_sd_%d_re" % (res_dir, int(step)))
    Sci[step].real().show_matrix("%s/S_ci_%d_re" % (res_dir, int(step)))
for step in range(nsteps-1):
    Stsd[step].real().show_matrix("%s/St_sd_%d_re" % (res_dir, int(step)))
    Stci[step].real().show_matrix("%s/St_ci_%d_re" % (res_dir, int(step)))

# Make the Hvib in the many-body basis
Hvib,Hvibsd = [[]], [[]]
ci_hvib = None
sd_hvib = None
for step in range( nsteps-1 ):
    #Calculating SD NACs
    sd_nacs = (  0.5j / dt ) * CMATRIX ( ( Stsd[step] - Stsd[step].T() ).real() )    
    #Calculating CI NACs
    ci_nacs = (  0.5j / dt ) * CMATRIX ( ( Stci[step] - Stci[step].T() ).real() )    
    #Creating CI Vibrational Hamiltonian
    ci_hvib = E[step] - ci_nacs
    #Creating SD Vibrational Hamiltonian 
    sd_hvib = Esd[step] - sd_nacs
    Hvibsd[0].append( sd_hvib )
    Hvib[0].append( ci_hvib )
    sd_hvib.real().show_matrix("%s/Hvib_sd_%d_re" % (res_dir, int( step )))
    sd_hvib.imag().show_matrix("%s/Hvib_sd_%d_im" % (res_dir, int( step )))
    ci_hvib.real().show_matrix("%s/Hvib_ci_%d_re" % (res_dir, int( step )))
    ci_hvib.imag().show_matrix("%s/Hvib_ci_%d_im" % (res_dir, int( step )))

Hvib[0].append( ci_hvib)    # appending the last element twice to make it nsteps
Hvibsd[0].append( sd_hvib )




ntraj  = len(Hvib)
nsteps = len(Hvib[0])
nCIs   = Hvib[0][0].num_of_cols

# Make a list for the SD energies and populate it
CI_energy = []
SD_energy = []
md_time = list( range(nsteps) )

for sd_index in range( nCIs ):
    CI_energy.append( [] )
    SD_energy.append( [] )
    for step in range( nsteps ):        
        Ensd = Hvibsd[0][step].get( sd_index, sd_index ).real 
        E0sd = Hvibsd[0][step].get( 0, 0 ).real
        En = Hvib[0][step].get( sd_index, sd_index ).real 
        E0 = Hvib[0][step].get( 0, 0 ).real
        CI_energy[ sd_index ].append( En - E0 )
        SD_energy[ sd_index ].append( Ensd - E0sd )        
        
CI_energy = np.array( CI_energy  )
SD_energy = np.array( SD_energy ) 
md_time   = np.array( md_time )

# Functions to compute the time-averaged SD NACs and make a list of them
sd_res = data_stat.cmat_stat2(Hvibsd[0],2)
sd_tNACs = []
# Functions to compute the time-averaged CI NACs and make a list of them
ci_res = data_stat.cmat_stat2(Hvib[0],2)
ci_tNACs = []

for i in range(nCIs):
    ci_tNACs.append( [] )
    sd_tNACs.append( [] )
    for j in range(nCIs):
        #ci_tNACs[i].append( ci_res.get(i,j).imag * 1000.0 / units.ev2Ha )        
        sd_tNACs[i].append( sd_res.get(i,j).imag * 1000.0 /units.ev2Ha  )        
        ci_tNACs[i].append( ci_res.get(i,j).imag * 1000.0 /units.ev2Ha  )        
sd_tNACs = np.array(sd_tNACs)
ci_tNACs = np.array(ci_tNACs)


nstates = 5
orbitals = ["GS","4 -> 5","4 -> 6","3 -> 5","3 -> 6","41 -> 44","41 -> 45","40 -> 43","40 -> 44","40 -> 45"]
state_label = ["GS","2","3","4","5","6","7","8","9","10"]
state_label_x =["GS","H->L","H->L+1","H->L+2","H-1->L","H-1->L+1","H-1->L+2","H-2->L","H-2->L+1","H-2->L+2"]
# Figures - Time-Dependent CI Data
plt.figure(num=None, figsize=(6,6), dpi=100, edgecolor='black', frameon=True)
plt.title('INDO CI Energies', fontsize=25)
plt.xlabel('Time, fs')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
for state in range(1,nstates ):    
    plt.plot(md_time, CI_energy[state], label=state_label[state], linewidth=2)
plt.tight_layout()
plt.savefig("CI_energytime.png")


fig=plt.figure(num=None, figsize=(6,6), dpi=100, edgecolor='black', frameon=True)
plt.title('INDO CINACs', fontsize=25)
ax = plt.subplot(111)
ax.set_xticks(np.arange(len(orbitals)))
ax.set_yticks(np.arange(len(orbitals)))
ax.set_xticklabels(state_label,fontsize=15)
ax.set_yticklabels(state_label,fontsize=15)
cb=ax.imshow(ci_tNACs, origin='lower', cmap='plasma', interpolation='nearest')
plt.xticks(rotation=90)
plt.colorbar(cb,label="meV")

plt.tight_layout()
plt.savefig("CI_avgNAC.png")


# Figures - Time-Dependent SD Data
plt.figure(num=None, figsize=(6, 6), dpi=100, edgecolor='black', frameon=True)
ax=plt.subplot(1,1,1)
plt.title('INDO Energies', fontsize=25)
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

fig=plt.figure(num=None, figsize=(6, 6), dpi=100, edgecolor='black', frameon=True)
plt.title('INDO SD NACs', fontsize=25)
ax = plt.subplot(111)
ax.set_xticks(np.arange(len(orbitals)))
ax.set_yticks(np.arange(len(orbitals)))
ax.set_xticklabels(labels=orbitals,fontsize=15)
ax.set_yticklabels(labels=orbitals,fontsize=15)
plt.xticks(rotation=90)
cb=ax.imshow(sd_tNACs,origin='lower', cmap='plasma', interpolation='nearest')
plt.colorbar(cb,label="meV")

plt.tight_layout()
plt.savefig("SD_avgNACs.png")


