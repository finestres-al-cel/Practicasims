import matplotlib
matplotlib.use("TkAgg")

from PIL import Image
import matplotlib.pyplot as plt
import glob
from parameters import path

filenames = sorted(glob.glob(path + "*.jpg"))

fig, ax = plt.subplots()

im = ax.imshow(Image.open(filenames[0]))
ax.axis("off")

plt.show(block=False)

for filename in filenames:
    im.set_data(Image.open(filename))
    fig.canvas.draw_idle()
    plt.pause(0.1)
