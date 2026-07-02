import matplotlib
matplotlib.use("TkAgg")

from matplotlib.animation import FuncAnimation
from PIL import Image
import matplotlib.pyplot as plt
import glob
from parameters import path

filenames = sorted(glob.glob(path + "*3D.jpg"))

fig, ax = plt.subplots()

im = ax.imshow(Image.open(filenames[0]))
ax.axis("off")

def update(frame):
    im.set_data(Image.open(filenames[frame]))
    return [im]

ani = FuncAnimation(
    fig,
    update,
    frames=len(filenames),
    interval=100,
    blit=True,
    repeat=True
)

plt.show()
