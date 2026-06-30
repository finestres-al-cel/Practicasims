import matplotlib
matplotlib.use('TkAgg') # do this before importing pylab
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import glob
from parameters import path

fig = plt.figure()
ax = fig.add_subplot(111)

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
