import sys
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from matplotlib.animation import FuncAnimation


fig, axes = plt.subplots(2, 2, sharex=True)
axes = axes.reshape(4)

labels = ["Energy", "Density", "Speed", "Pressure"]
N = len(labels)

if len(sys.argv) == 1:
    data = np.genfromtxt("tube.dat")
    for i in range(N):
        axes[i].plot(data[:, 0], data[:, i+1], lw=1)
        axes[i].set_ylabel(labels[i])
        if i > 1:
            axes[i].set_xlabel("$x$")
    fig.tight_layout()
    fig.savefig("tube.png")

elif sys.argv[1] == "animate":
    files = glob("*_data.dat")
    number = len(files)
    whole = []
    plots = []

    frames = 50
    step = int(number/frames)
    for i in range(0, number, step):
        data = np.genfromtxt("%d_data.dat"%i)
        whole.append(data)

    for i in range(N):
        plot = axes[i].plot([], [], lw=1)[0]
        plots.append(plot)
        axes[i].set_ylabel(labels[i])
        axes[i].set_xlim(0, 1)
        axes[i].set_ylim(int(min(data[:, i+1])), np.ceil(1.2*max(data[:, i+1])))
        if i > 1:
            axes[i].set_xlabel("$x$")

    def update(j):
        data = whole[j]
        for i in range(N):
            plots[i].set_data(data[:, 0], data[:, i+1])

    fig.tight_layout()
    ani = FuncAnimation(fig, update, frames = frames)
    ani.save("animated.gif", writer="imagemagick", dpi=70, fps=10)
