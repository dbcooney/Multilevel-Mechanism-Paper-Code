"""
Script used to generate Figure 3.2, which illustrates steady state densities for the
multilevel dynamics under a range of within-group assortment probabilities.
"""


import matplotlib.pyplot as plt
import numpy as np

import scipy.integrate as spi
import os

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


"""
Parameters describing game payoffs and the initial density. 
"""

alpha = -1.0
beta = -1.0
gamma = 1.5
theta = 2.0

r_any_coop = (-beta) / ( gamma + alpha - beta)
r_all_coop = (-alpha - beta) / (gamma - beta) 

def G_r(x, r, gamma, alpha, beta):
	return (gamma + alpha * r) * x + (1.0 - r) * alpha * (x ** 2.0)
	


"""
Function for plotting unnormalized version of steady state density for various values 
of the assortment probability r. 
"""

def steady_density(x,lamb, r, beta, alpha, gamma, theta):
	a = (1.0 - r) * alpha
	b = beta + r * (gamma + alpha - beta)
	gamma = gamma + r * alpha
	
	if r > r_any_coop and x <= -(b/a):
		return 0.0
		
	elif r < r_any_coop:
		x_power = (1.0 / np.abs(b)) * (lamb * (gamma + a) - (np.abs(b) - a) * theta) - 1.0
		bminusax_power = -1.0 * (1.0 / np.abs(b)) * (lamb * (gamma + np.abs(b) + a) + a * theta) - 1.0
		return (x ** x_power) * ((1.0 - x) ** (theta - 1.0)) * ((np.abs(b) - a * x) ** bminusax_power)
		
	elif r_any_coop < r < r_all_coop:
		x_power = (1.0 / b) * (lamb*(np.abs(a) - gamma) + (np.abs(a) - b) * theta) - 1.0
		bminusax_power = (1.0 / b) * (lamb * (gamma - np.abs(a) - b) - np.abs(a) * theta) - 1.0
		return (x ** x_power) * ((1.0 - x) ** (theta - 1.0)) * ((np.abs(a) * x - b) ** bminusax_power)	

steady_state_vec = np.vectorize(steady_density)
lamb = 10.


"""
Plotting the normalized steady state densities for various assortment probabilities r
and fixed game parameters and Hoelder exponent near x=1.
"""

r = 0.0

x_step = 0.001
x_range = np.arange(x_step,1.0+x_step,x_step)

holderr0 = steady_state_vec(x_range,lamb,r,beta,alpha,gamma,theta)
holderr0int = spi.simps(holderr0,x_range)
holderr0 = holderr0 / holderr0int
#print  spi.simps(holderr0,x_range)
#print holderr0
plt.plot(x_range,holderr0, lw = 4., color = plt.cm.YlOrRd(0.), label = r"$r = 0.0$")

r = 0.15

x_range = np.arange(0.000,1.0+x_step,x_step)

holderr2 = steady_state_vec(x_range,lamb,r,beta,alpha,gamma,theta)
holderr2int = spi.simps(holderr2,x_range)
holderr2 = holderr2/ holderr2int
#print  spi.simps(holderr2,x_range)
#print holderr2
plt.plot(x_range,holderr2, lw = 4., color = plt.cm.YlOrRd(0.2), label = r"$r = 0.15$")


r = 0.3

holderr4 = steady_state_vec(x_range,lamb,r,beta,alpha,gamma,theta)
holderr4int = spi.simps(holderr4,x_range)
holderr4 = holderr4/ holderr4int
#print  spi.simps(holderr4,x_range)
#print holderr4
plt.plot(x_range,holderr4, lw = 4., color = plt.cm.YlOrRd(0.4), label = r"$r = 0.3$")

r = 0.45

holderr6 = steady_state_vec(x_range,lamb,r,beta,alpha,gamma,theta)
holderr6int = spi.simps(holderr6,x_range)
holderr6 = holderr6/ holderr6int
#print  spi.simps(holderr4,x_range)
#print holderr4
plt.plot(x_range,holderr6, lw = 4., color = plt.cm.YlOrRd(0.6), label = r"$r = 0.45$")



r = 0.6

holderr8 = steady_state_vec(x_range,lamb,r,beta,alpha,gamma,theta)
holderr8int = spi.simps(holderr8,x_range)
holderr8 = holderr8/ holderr8int
#print  spi.simps(holderr4,x_range)
print holderr8
plt.plot(x_range,holderr8, lw = 4., color = plt.cm.YlOrRd(0.8), label = r"$r = 0.6$")


r = 0.7

holderr10 = steady_state_vec(x_range,lamb,r,beta,alpha,gamma,theta)
holderr10int = spi.simps(holderr10,x_range)
holderr10 = holderr10/ holderr10int
#print  spi.simps(holderr4,x_range)
#print holderr10
plt.plot(x_range,holderr10, lw = 4., color = plt.cm.YlOrRd(1.0), label = r"$r = 0.7$")
plt.axvline(x= -beta / alpha + (r / (1.0 - r)) * ((gamma + alpha) / (-alpha)), lw = 4., color = plt.cm.YlOrRd(1.0), ls = '--')



plt.axis([0.0,1.0,0.0,15.0])
plt.legend(loc = "upper center")


plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.xlabel(r"Fraction of Cooperators ($x$)", fontsize = 20., labelpad = 20.)
plt.ylabel(r"Probability Density", fontsize = 20.)

plt.tight_layout()




script_folder = os.getcwd()
mechanism_folder = os.path.dirname(script_folder)

plt.savefig(mechanism_folder + "/Figures/rdensities.png")

plt.show()