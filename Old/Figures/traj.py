import matplotlib
import scipy
import pylab 
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import glob
from parameters import *


pvelr=scipy.genfromtxt(path+'trajectories.dat',comments='#')
t,npar,coor = pvelr[:,0],pvelr[:,1],pvelr[:,2:]

fig = plt.figure()

ax = fig.add_subplot(111)
for j in range(0,len(coor[0,:])/3):
 r=np.sqrt(coor[:,j*3]*coor[:,j*3]+coor[:,j*3+1]*coor[:,j*3+1]+coor[:,j*3+2]*coor[:,j*3+2])
 print(j)
 print(r)
 ax.plot(t, r, c='k')
#ax.set_ylim([np.min(r),np.max(r)])
ax.set_ylim([0,1])
ax.set_ylabel('Radius')
ax.set_xlabel('Time')
ax.set_xlim([np.min(t),np.max(t)])
#ax.set_aspect((np.max(t)-np.min(t))/(np.max(r)-np.min(r)))
ax.set_aspect((np.max(t)-np.min(t)))
plt.savefig(path+'traj.jpg')
plt.show()
