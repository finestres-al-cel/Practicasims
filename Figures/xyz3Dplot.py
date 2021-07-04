import matplotlib
import scipy
import pylab 
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import glob
from parameters import *
import os
import glob

files = glob.glob(path+'*.jpg')
for f in files:
    os.remove(f)

for fil in sorted(glob.glob(path+'Snap*')):
 print fil
 pvelr=scipy.genfromtxt(fil,comments='#')
 x,y,z,vx,vy,vz,e1,e2,e3,m = pvelr[:,0],pvelr[:,1],pvelr[:,2],pvelr[:,3],pvelr[:,4],pvelr[:,5],pvelr[:,6],pvelr[:,7],pvelr[:,8],pvelr[:,9]

 fig = plt.figure()

 ax=fig.add_subplot(111,projection='3d')
 ax.scatter(x,y,z,c='r',s=1)
# ax.set_xlabel('x [kpc]')
# ax.set_ylabel('y [kpc]')
# ax.set_zlabel('z [kpc]')
 ax.set_ylim([-15,15])
 ax.set_xlim([-15,15])
 ax.set_zlim([-15,15])
 ax.grid(False)
 ax.set_aspect(1.)
# plt.show()
 pylab.savefig(fil+'3D.jpg')
# plt.show()
 plt.clf()
 plt.cla()
 plt.close()
# plt.show()
