import sys
import numpy as np
from glob import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def calculateU4(Gamma, beta, p3):
	num = 1-Gamma
	denom = rho5*(p3 + Gamma*p5)
	return (p3 - p5)*np.sqrt(num/denom)

def calculateU2(Gamma, beta, p3):
	Gamma2 = Gamma**2
	num = (1-Gamma2)*p1**(1/gamma)
	denom = Gamma2*rho1
	return (p1**beta - p3**beta)*np.sqrt(num/denom)

def calculateU3(Gamma, beta, p3):
	denom = 0.5*rho5*((gamma+1)*p3 + (gamma-1)*p5)
	return u5 + (p3 - p5)/np.sqrt(denom)

def get_speeds(Gamma, beta, p3):
	u2 = calculateU2(Gamma, beta, p3)
	u4 = calculateU4(Gamma, beta, p3)
	return u4-u2

def find_p3(Gamma, beta, thresh=1e-4):
	a = p5
	b = p1
	ans = 1
	p3 = (a+b)/2
	while abs(ans) > thresh:
		ans = get_speeds(Gamma, beta, a)
		fa = get_speeds(Gamma, beta, p3)
		if ans == 0:
			break
		elif fa*ans < 0:
			b = p3
		else:
			a = p3
		p3 = (a+b)/2
	return p3

def analytic_sod():
	n = 1000
	x = np.linspace(0, 1.0, n)
	Gamma = (gamma-1)/(gamma+1)
	beta = 0.5*(gamma-1)/gamma
	x_final = 0.9
	c1 = np.sqrt(gamma*p1/rho1)

	p3 = find_p3(Gamma, beta)
	p4 = p3
	rho3 = rho1*(p3/p1)**(1/gamma)
	rho4 = rho5*(p4+Gamma*p5)/(p5+Gamma*p4)
	u3 = calculateU3(Gamma, beta, p3)
	u4 = u3

	density_relation = (rho3/rho5 - Gamma)/(1+Gamma*rho3/rho5)
	vpost = 2*(np.sqrt(gamma)/(gamma-1))*(1-p3**((gamma-1)/(2*gamma)))
	vshock = vpost*(density_relation/(density_relation -1))
	x1 = x0 - t*c1
	x2 = x0 + (vpost - c1 + ((gamma - 1)/2)*vpost)*t
	x3 = x0 + vpost*t
	x4 = x0 + vshock*t
	pressure = np.zeros(n)
	speed = np.zeros(n)
	density = np.zeros(n)
	for i in range(n):
		if x[i] < x1:
			pressure[i] = p1
			speed[i] = u1
			density[i] = rho1
		elif x[i] < x2:
			c = Gamma*((x0-x[i])/t) + (1-Gamma)*c1
			density[i] = rho1*(c/c1)**(2/(gamma-1))
			pressure[i] = p1*(density[i]/rho1)**gamma
			speed[i] = (1-Gamma)*(c1-(x0-x[i])/t)
			x4 = x_final

		elif x[i] < x3:
			pressure[i] = p3
			density[i] = rho3
			speed[i] = u3

		elif x[i] < x4:
			pressure[i] = p4
			density[i] = rho4
			speed[i] = u4
		else:
			pressure[i] = p5
			density[i] = rho5
		label = "Exact"
	return [x, density, speed, pressure]

if sys.argv[1] == '1':
	files = glob("*_sedov.dat")
	number_files = len(files)
	fig, axes = plt.subplots(number_files, sharex=True)
	for i in range(number_files):
		y = np.genfromtxt(files[i])
		x = np.linspace(0, 1, y.shape[0])*128
		# print(len(y), len(x))
		if number_files == 1:
			axes.plot(x, y)
			axes.set_ylabel("Density (kg/m$^3$)")
		else:
			axes[i].plot(x, y)
			axes[i].set_ylabel("Density (kg/m$^3$)")
	if number_files > 1:
		axes[-1].set_xlabel("Position (m)")
	fig.tight_layout()
	fig.savefig("sedov.pdf")

elif sys.argv[1] == '0':
	fig, axes = plt.subplots(3, sharex=True)
	axes = axes.reshape(3)

	labels = ["Density", "Speed", "Pressure"]
	N = len(labels)
	with open("tube.dat", "r") as f:
		for i in range(9):
			text = f.readline()
			exec(text)
	analytic = analytic_sod()
	data = np.genfromtxt("tube.dat", skip_header=9)
	for i in range(N):
		axes[i].plot(data[:, 0], data[:, i+2], "-o", label="Numerical", lw=0.5, ms=2, markevery=5)
		axes[i].plot(analytic[0], analytic[i+1], label="Exact")
		axes[i].set_ylabel(labels[i])
		axes[i].legend()
		if i > 1:
			axes[i].set_xlabel("$x$")
	fig.tight_layout()
	fig.savefig("shock.pdf")
