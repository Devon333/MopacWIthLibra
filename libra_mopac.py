from libra_mopac import *
from liblibra_core import *
from libra_py.workflows.nbra import mapping, step2_many_body, step3, step4, step3_many_body
import numpy as np



def get_raw_basis(params):
    '''
    Main function that uses other functions to reads SD basis from 
    mopac output file returns a list of list containing SD representations
    '''
    output_file_name=params["logfile_directory"]+params["file_prefix"]+"0"+".out"
    output = open(output_file_name)
    min_orb=1000000
    max_orb=0
    sd_basis=read_mopac_SD_config(params)
    HOMO_INDEX=0
    line = output.readline()
    while line:
        line = output.readline()
        if "SHELL" in line:
            line = output.readline()
            split_line = line.split()
            #print(line)
            HOMO_INDEX=int(split_line[1])
            break
    #    line = output.readline()
    output.close()
    
    for i in range(1,len(sd_basis)):
        temp_bas=sd_basis[i][0]
        #print(temp_bas)
        if temp_bas[0] < min_orb:
            min_orb=temp_bas[0]
            print(min_orb)
        if temp_bas[0] > max_orb:
            max_orb=temp_bas[1]
            print(max_orb)
            #print(sd_basis[i][0])


    orbital_range=list(range(min_orb,max_orb))
    gs_basis = []
    for i in orbital_range:
        gs_basis.append(i)
        gs_basis.append(i*-1)
    print(gs_basis)
    
    basis=[]
    basis.append(gs_basis) 
    for i in range(1,len(sd_basis)):
       exc = sd_basis[i][0]
       j=0
       #temp_lis=[]
       basis.append([])
       while j < len(gs_basis):
           if gs_basis[j] == exc[0]:
               temp_ele =exc[1]
               basis[i].append(temp_ele)
           else:
               basis[i].append(gs_basis[j])
           j+=1
       #basis.append(temp_lis)
    print(basis[0])
    print(basis[0][1])
    print("this it the basis",basis)

    if len(basis) > 0:
        return basis
    else:
        print("did not create sd basis")


def reindex_basis(basis):
    HOMO = max(basis[0])
    lis=[]
    for i in basis[0]:
        if i >0:
            lis.append(i)
    
    HOMO_N=min(lis)
    print("HOMO index ",HOMO)
    homo_index=HOMO-HOMO_N
    print("homo-N index ",HOMO_N)
    print("new HOMO index ",homo_index)
    reindexed_basis=[]
    for sd in range(len(basis)):
        temp_li=[]
        reindexed_basis.append([])
        for orb in range(len(basis[0])):
            if basis[sd][orb] >0:
                reindexed_basis[sd].append(basis[sd][orb]-(HOMO_N-1))
                #temp_li.append(orb-HOMO_N)
            if basis[sd][orb] <0:
                reindexed_basis[sd].append((abs(basis[sd][orb])-(HOMO_N-1))*-1)
                #temp_li.append((abs(orb)-HOMO_N)*-1)
        #reindexed_basis.append(temp_li)

    print("re-indexed basis",reindexed_basis)
    return reindexed_basis




def read_num_of_MOs_output(params):
    '''
    function that reads the number of MOs for output file
    '''
    output_file_name=params["logfile_directory"]+params["file_prefix"]+"0"+".out"
    output = open(output_file_name)
    num_of_mos=0
    line = output.readline()
    while line:
        line = output.readline()
        if "TOTAL" in line:
            split_line = line.split()
            #print(line)
            num_of_mos=int(split_line[1])
            break
    #    line = output.readline()
    output.close()
    if num_of_mos > 0:
        return num_of_mos
    else:
        print("number of molecular orbitals not found")



def read_mopac_MO_eigenvalues(params):
    '''
    function reads MO energies/eigenvalues and returna a list
    '''
    output_file_name=params["logfile_directory"]+params["file_prefix"]+"0"+".out"
    output = open(output_file_name)
    num_of_mos = read_num_of_MOs_output(params)
    energies=[""]
    mo_coeff=[[] for i in range(num_of_mos+1)]
    line = output.readline()
    counter=0
    orb_count=0
    while line:
        #print(line)
        line=output.readline()
        if "Root No." in line:
            line = output.readline()
            line = output.readline()
            line = output.readline()
            line = output.readline()
            for energy in line.split():
                energies.append(float(energy)) 
    output.close()
    if len(energies) > 0:
        return energies
    else:
        print("molecular orbital energies not found")
        exit()



def read_mopac_MO_eigenvectors(params):
    '''
    function that reads molecular orbital vectors from mopac output file and returns 
    a list of list containing molecular orbital vectors
    '''
    output_file_name=params["logfile_directory"]+params["file_prefix"]+"0"+".out"
    output = open(output_file_name)
    num_of_mos = read_num_of_MOs_output(params)
    energies=[""]
    mo_coeff=[[] for i in range(num_of_mos+1)]
    line = output.readline()
    counter=1
    orb_count=0
    while line:
        #print(line)
        line=output.readline()
        if "Root No." in line:
            line = output.readline()
            line = output.readline()
            line = output.readline()
            line = output.readline()
            line = output.readline()
            line = output.readline()
            line_count = line.split()
            orb_count = len(line_count[3:])
            #print(mo_coeff)
            for atomic_orb in range(num_of_mos):
                split_line = line.split()
                split_line = split_line[3:]
                #print(split_line)
                for ato in range(len(split_line)):
                    mo_coeff[counter+ato].append(float(split_line[ato]))
                line = output.readline()
            #print("counter: ",counter)
            counter += orb_count
    output.close()
    #print(mo_coeff) 
    if len(mo_coeff) > 0:
        return mo_coeff
    else:
        print("molecular orbital coefficients not found")
        exit()



def read_mopac_num_of_SD_states(params):
    '''
    function reads the number of SDs from output file
    '''
    output_file_name=params["logfile_directory"]+params["file_prefix"]+"0"+".out"
    output = open(output_file_name)
    line = output.readline()
    num_of_SD_states=0
    while line:
        line=output.readline()
        if "CI excitations=" in line:
            split_line = line.split()
            #print(split_line)
            num_of_SD_states = int(split_line[4])
            line = output.readline()
            
    output.close()
    #print(num_of_SD_states)
    if num_of_SD_states > 0:
        return num_of_SD_states
    else:
        print("number of molecular orbitals not found")
        exit()


def read_mopac_SD_energies(params):
    '''
    function returns a list of SD energies from mopac output file
    '''
    output_file_name=params["logfile_directory"]+params["file_prefix"]+"0"+".out"
    output = open(output_file_name)
    line = output.readline()
    num_of_SD_states=read_mopac_num_of_SD_states(params)
    sd_energies=[]
    while line:
        #print(line)
        line=output.readline()
        if "CI excitations=" in line:
            line=output.readline()
            line=output.readline()
            line=output.readline()
            line=output.readline()
            line=output.readline()
            split_line = line.split()
            for state in range(num_of_SD_states):
                sd_energies.append(float(split_line[1]))
                line = output.readline()
                split_line = line.split()

    #print(sd_energies) 
    output.close()
    if len(sd_energies) > 0:
        return sd_energies
    else:
        print("SD energies not found")
        exit()


def read_mopac_SD_config(params):
    '''
    Function that reads in SD configurations from output fle
    '''
    output_file_name=params["logfile_directory"]+params["file_prefix"]+"0"+".out"
    output = open(output_file_name)
    line = output.readline()
    num_of_SD_states=read_mopac_num_of_SD_states(params)
    sd_configs=[]
    while line:
        #print(line)
        line=output.readline()
        if "CI excitations=" in line:
            line=output.readline()
            line=output.readline()
            line=output.readline()
            line=output.readline()
            line=output.readline()
            split_line = line.split()
            for state in range(num_of_SD_states):
               if float(split_line[1]) ==0.00:
                   sd_configs.append([[0,0],"gs"])
               else:
                   if len(split_line) == 13:
                       sd_configs.append([[int(split_line[9]),int(split_line[11])],"alp"])
                   if len(split_line) == 14:
                       sd_configs.append([[int(split_line[10]),int(split_line[12])],"alp"])
                       
               line = output.readline()
               split_line = line.split()
    #print(sd_configs) 
    output.close()
    if len(sd_configs) > 0:
        return sd_configs
    else:
        print("SD configs not found")
        exit()



def read_mopac_CI_matrix(params):
    '''
    Function that reads CI matrix from mopac output
    '''
    output_file_name=params["logfile_directory"]+params["file_prefix"]+"0"+".out"
    output = open(output_file_name)
    num_of_sd = read_mopac_num_of_SD_states(params)
    line = output.readline()
    counter=1
    orb_count=0
    ci_mat=[]
    while line:
        #print(line)
        if "CI Matrix" in line:
            line = output.readline()
            #print(mo_coeff)
            for state in range(num_of_sd):
                split_line = line.split()
                temp_li=[]
                for i in split_line:
                    temp_li.append(float(i))
                ci_mat.append((temp_li))
                line= output.readline()
                split_line = line.split() 
        line = output.readline()
    output.close()
    if len(ci_mat) > 0:
        return ci_mat
    else:
        print("CI matrix not found")
        exit()
    



def read_mopac_CI_energies(params):
    '''
    Function that reads CI energies from mopac output file
    returns a list of CI energies
    '''
    output_file_name=params["logfile_directory"]+params["file_prefix"]+"0"+".out"
    output = open(output_file_name)
    line = output.readline()
    num_of_SD_states=read_mopac_num_of_SD_states(params)
    CI_energies=[0.00]
    count =1
    while line:
        #print(line)
        line=output.readline()
        if "CI trans." in line:
            line=output.readline()
            line=output.readline()
            line=output.readline()
            split_line = line.split()
            print(split_line)
            print(num_of_SD_states)
            for state in range(1,num_of_SD_states):
               print(split_line[0])
               CI_energies.append(float(split_line[1]))
               #else:
               #    if len(split_line) == 13:
               #        CI_energies.append(float(split_line[1]))
               #    if len(split_line) == 14:
               #        CI_energies.append(float(split_line[1]))
               #        
               line = output.readline()
               split_line = line.split()
    print(CI_energies) 
    output.close()
    N, M = len(CI_energies), len(CI_energies)
    res = CMATRIX(N,M)
    for i in range(len(CI_energies)):
       test = CI_energies[i]
       res.set(i,i,test)
       print(test)
    res.show_matrix()
    if len(CI_energies) > 0:
        return res
    else:
        print("CI Energies not found")
        exit()


#function that calculates overlap between Slater determinants using libra functions
def SD_ovlp(basis, time_ov_mat):
  """basis -> list[list]
     time_ov_mat -> list[CMATRIX]
     takes in list of SD basis lists (list with indexes of KS orbitals used to describe 
     a single SD) and computes overlaps between SD bases by pulling the appropriate indices
     from the time doverlap matrix (time_ov_mat)
     returns a matrix with overlaps of SD overlaps
  """
  N, M = len(basis), len(basis)
  res = CMATRIX(N,M)
  for i in range(len(basis)):
      for j in range(len(basis)):
         test = mapping.ovlp_arb(basis[i], basis[j], time_ov_mat,False)
         #test *= phases.get(j,0).real
         res.set(i,j,test)
         #print(test) 
  return res


def Mat_avg(cmatList):
    """cmatList -> list[CMATRIX]
    takes in a list of cmatrices and returns the average value of each matrix element 
    over all of the in a numpy matrix
    """
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

