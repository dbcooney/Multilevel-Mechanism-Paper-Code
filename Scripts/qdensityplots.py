"""
Script used to generate Figure 3.2, which illustrates steady state densities for the
multilevel dynamics under a range of defector defection probabilities $q$.
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

q_any_coop = (-beta) / ( gamma + alpha - beta)
q_all_coop = (-alpha - beta) / (gamma - beta) 


def G_q(x, q, gamma, alpha, beta):
	return gamma * (1 - q) * x + (alpha + q * gamma) * (x ** 2.0)

	
	
"""
Function for plotting unnormalized version of steady state density for various values 
of the defector detection probability $q$.  
"""
	
	
def steady_density(x,lamb, r, beta, alpha, gamma, theta):
	a = alpha + q * gamma
	b = (1.0 - q) * beta
	
	x_eq = -(1.0 - q) * beta * (1.0 / (alpha + q * gamma))
	
	E = 1.0 / (x_eq * (1.0 - x_eq))
	Delta_pi = -(b + a)
	G_1 = G_q(1.0,r,gamma,alpha,beta)
	G_eq = G_q(x_eq,r,gamma,alpha,beta)
	Delta_G = G_1 - G_eq
	
	x_power = (1.0 / np.abs(b)) * (lamb * G_1 - Delta_pi * theta) - 1.0
	one_minus_x_power = theta - 1.0
	bminusax_power = (lamb * Delta_G - Delta_pi * theta) * (E / np.abs(a)) - 1.0
	
	if q > q_all_coop:
		return 0.0
	else:
		return (x ** x_power) * ((1.0 - x) ** (theta - 1.0)) * ((np.abs(b + a * x)) ** bminusax_power)

steady_state_vec = np.vectorize(steady_density)


"""
Plotting the normalized steady state densities for various defector detection
probabilities $q$ and fixed game parameters and Hoelder exponent near x=1.
"""


lamb = 10.
q = 0.0

x_step = 0.001
x_range = np.arange(x_step,1.0+x_step,x_step)

holderr0 = steady_state_vec(x_range,lamb,q,beta,alpha,gamma,theta)
holderr0int = spi.simps(holderr0,x_range)
holderr0 = holderr0 / holderr0int
#print  spi.simps(holderr0,x_range)
#print holderr0
plt.plot(x_range,holderr0, lw = 4., color = plt.cm.YlOrRd(0.), label = r"$q = 0.0$")

q = 0.125

x_range = np.arange(0.000,1.0+x_step,x_step)

holderr2 = steady_state_vec(x_range,lamb,q,beta,alpha,gamma,theta)
holderr2int = spi.simps(holderr2,x_range)
holderr2 = holderr2/ holderr2int
#print  spi.simps(holderr2,x_range)
#print holderr2
plt.plot(x_range,holderr2, lw = 4., color = plt.cm.YlOrRd(0.2), label = r"$q = 0.125$")


q = 0.25

holderr4 = steady_state_vec(x_range,lamb,q,beta,alpha,gamma,theta)
holderr4int = spi.simps(holderr4,x_range)
holderr4 = holderr4/ holderr4int
#print  spi.simps(holderr4,x_range)
#print holderr4
plt.plot(x_range,holderr4, lw = 4., color = plt.cm.YlOrRd(0.4), label = r"$q = 0.25$")

q = 0.375

holderr6 = steady_state_vec(x_range,lamb,q,beta,alpha,gamma,theta)
holderr6int = spi.simps(holderr6,x_range)
holderr6 = holderr6/ holderr6int
#print  spi.simps(holderr4,x_range)
#print holderr4
plt.plot(x_range,holderr6, lw = 4., color = plt.cm.YlOrRd(0.6), label = r"$q = 0.375$")



q = 0.5

holderr8 = steady_state_vec(x_range,lamb,q,beta,alpha,gamma,theta)
holderr8int = spi.simps(holderr8,x_range)
holderr8 = holderr8/ holderr8int
#print  spi.simps(holderr4,x_range)
print holderr8
plt.plot(x_range,holderr8, lw = 4., color = plt.cm.YlOrRd(0.8), label = r"$q = 0.5$")




q = 0.625

holderr10 = steady_state_vec(x_range,lamb,q,beta,alpha,gamma,theta)
holderr10int = spi.simps(holderr10,x_range)
holderr10 = holderr10/ holderr10int
#print  spi.simps(holderr4,x_range)
print holderr10
plt.plot(x_range,holderr10, lw = 4., color = plt.cm.YlOrRd(1.0), label = r"$q = 0.625$")



plt.axis([0.0,1.0,0.0,12.0])
plt.legend(loc = "upper center")


plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.xlabel(r"Fraction of Cooperators ($x$)", fontsize = 20., labelpad = 20.)
plt.ylabel(r"Probability Density", fontsize = 20.)

plt.tight_layout()



print q_all_coop

script_folder = os.getcwd()
mechanism_folder = os.path.dirname(script_folder)

plt.savefig(mechanism_folder + "/Figures/qdensities.png")

plt.show()