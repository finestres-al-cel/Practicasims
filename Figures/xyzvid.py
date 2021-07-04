import matplotlib
matplotlib.use('TkAgg') # do this before importing pylab
from PIL import Image
import scipy
import pylab 
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import glob
from parameters import *

fig = plt.figure()
ax = fig.add_subplot(111)

#path='../DATA/galxoc/galxoc3/'

def animate():
  for i in range(0,10):
    filenames=sorted(glob.glob(path+'*.jpg'))
    im=plt.imshow(Image.open(filenames[0]))
    for filename in filenames[:]:
        image=Image.open(filename)
        im.set_data(image)
        fig.canvas.manager.window.after(100)        
        fig.canvas.draw() 
                

win = fig.canvas.manager.window
fig.canvas.manager.window.after(100, animate)
plt.show()
