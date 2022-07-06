import matplotlib
import scipy
import pylab 
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import glob
from parameters import *

#path='../DATA/galxoc/galxoc3/'

#fig = plt.figure()

#ax = fig.add_subplot(111)

import os
import glob

files = glob.glob(path+'*.jpg')
for f in files:
    os.remove(f)

for fil in sorted(glob.glob(path+'Snap*')):
 print(fil)
 pvelr=scipy.genfromtxt(fil,comments='#')
 x,y,z,vx,vy,vz,e1,e2,e3,m = pvelr[:,0],pvelr[:,1],pvelr[:,2],pvelr[:,3],pvelr[:,4],pvelr[:,5],pvelr[:,6],pvelr[:,7],pvelr[:,8],pvelr[:,9]

 fig = plt.figure()

 ax = fig.add_subplot(221)
 ax.scatter(x, y, c='r',s=0.1)
 ax1 = fig.add_subplot(222)
 ax1.scatter(x, z, c='r',s=0.1) 
 ax2 = fig.add_subplot(223)
 ax2.scatter(y, z, c='r',s=0.1)
 ax.set_ylim([-15, 15])
 ax.set_xlabel('X')
 ax.set_ylabel('Y')
 ax.set_xlim([-15,15])
 ax.set_aspect(1.)
 ax1.set_ylim([-15, 15])
 ax1.set_xlim([-15,15])
 ax1.set_xlabel('X')
 ax1.set_ylabel('Z')
 ax1.set_aspect(1.)
 ax2.set_ylim([-15,15])
 ax2.set_xlim([-15,15])
 ax2.set_xlabel('Y')
 ax2.set_ylabel('Z')
 ax2.set_aspect(1.)

 pylab.savefig(fil+'.jpg')
# plt.show()
 plt.clf()
 plt.cla()
 plt.close()

