import numpy as np
import matplotlib.pyplot as plt

data_flat = np.genfromtxt("energy.dat")
N = int(data_flat.shape[0]**(1/3))
data = np.zeros((N, N, N))

def transform3d(x, y, z):
	return N*N*x + N*y + z;

for i in range(N):
	for j in range(N):
		for k in range(N):
			index = transform3d(i, j, k)
			data[i, j, k] = data_flat[index]

N_half = int(N/2)
temp = data[:, :, N_half]
plt.imshow(temp, cmap = 'jet')
plt.colorbar()
plt.show()
