import matplotlib.pyplot as plt
import numpy as np
import sys

start=0
final=2
nsteps=final-start
names = ['tetraceneTraj/tetracene_T_'+str(count)+'_re.txt' for count in range(start,final)]
#names = ['tetracene/h2o_T_'+str(count)+'_re.txt' for count in range(start,final)]
#state_label = ["GS","2","3","4","5","6","7","8","9","10"]
#orbitals=["GS","H->L","H->L+1","H->L+2","H-1->L","H-1->L+1","H-1->L+2","H-2->L","H-2->L+1","H-2->L+2"]
state_label = ["2","3","4","5","6","7","8","9","10"]
orbitals=["H->L","H->L+1","H->L+2","H-1->L","H-1->L+1","H-1->L+2","H-2->L","H-2->L+1","H-2->L+2"]
#state_label = ["1","2","3","4","5",'6','7',8]
#orbitals=["H->L","H->L+1","H->L+2","H-1->L","H-1->L+1"]
#image = np.matrix(np.zeros((10,10)))

avg_matrix= np.matrix(np.zeros((10,10)))
for name in names:
  with open(name,'r') as fi:
    lines = fi.readlines()
    image = np.matrix([[float(ele) for ele in line.split()] for line in lines])
    #print(image)
    for i in range(len(image)):
      for j in range(len(image)):
        avg_matrix[i,j] += abs(image[i,j])
    fi.close()
  #image = lines

for i in range(len(image)):
  for j in range(len(image)):
    avg_matrix[i,j] /= nsteps

fig=plt.figure(num=None, figsize=(6, 6), dpi=100, edgecolor='black', frameon=True)
plt.title('INDO CI Coefficients', fontsize=25)
#ax = fig.add_subplot(111)

li = []
for i in range(len(avg_matrix)):
    co_sum=0
    for j in range(len(avg_matrix)):
      co_sum += avg_matrix[i,j]
    print(co_sum)
    avg_matrix[i]=avg_matrix[i]*1/(co_sum)**(1/2)
    li.append((co_sum))  
plt.scatter(range(10),li)
plt.ylim(0,max(li))
plt.show()
print(avg_matrix)
li2=[]
for i in range(len(avg_matrix)):
    for j in range(len(avg_matrix)):
      li2.append(avg_matrix[i,j])


plt.scatter(range(10),li2[:10],color='k')
plt.scatter(range(10,20),li2[10:20],color='r')
plt.scatter(range(20,30),li2[20:30],color='b')
plt.scatter(range(30,40),li2[30:40],color='g')
plt.scatter(range(40,50),li2[40:50],color='k')
plt.scatter(range(50,60),li2[50:60],color='r')
plt.scatter(range(60,70),li2[60:70],color='b')
plt.scatter(range(70,80),li2[70:80],color='g')
plt.scatter(range(80,90),li2[80:90],color='k')
plt.scatter(range(90,100),li2[90:100],color='r')
plt.ylim(0,max(li2))
plt.show()

#ax.set_xticks(np.arange(5))
#ax.set_yticks(np.arange(5))
#ax.set_xticklabels(labels=orbitals,fontsize=15)
#ax.set_yticklabels(labels=state_label,fontsize=15)
#plt.xticks(rotation=90)
##plt.xticks([0,1,2,3,4,5,6,7,8,9])
##plt.yticks([0,1,2,3,4,5,6,7,8,9])
#cb=ax.imshow(avg_matrix[1:,1:],origin='lower', cmap='plasma', interpolation='nearest')
#plt.colorbar(cb,label="arbitrary units")
#
##cb.ax.tick_params(labelsize=30)
#
#plt.tight_layout()
##plt.show()
#plt.savefig("INDO_avgCoeffs.png")
##plt.imshow(avg_matrix,cmap='jet',origin='lower')
##plt.show()
##print(image)
#
