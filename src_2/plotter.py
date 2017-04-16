import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()

def transform3d(x, y, z, N):
	return N*N*x + N*y + z;

def get_matrix(file_name):
	data_flat = np.genfromtxt(file_name)
	N = int(data_flat.shape[0]**(1/2))
	data = data_flat.reshape((N, N))
	return data[:, :]

files = glob("*.dat")
number_files = len(files)

data = [get_matrix("%d_energy.dat"%i) for i in range(number_files)]
data[0][0,0] = 1e20
im = ax.imshow(np.log10(data[0]),  cmap = 'jet')
fig.colorbar(im)

def animate(i):
	global ax, data, im
	im.set_data(np.log10(data[i]))
	ax.set_xlabel("%d"%i)

ani = FuncAnimation(fig, animate, frames = number_files, interval = 100)
# ani.save("explosion.gif", writer = 'imagemagick')
plt.show()
