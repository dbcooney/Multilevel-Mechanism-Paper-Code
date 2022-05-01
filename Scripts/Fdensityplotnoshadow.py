"""
This sript is used to generate Figure 4.2(left), which illustrates steady state densities
 for the multilevel dynamics under various levels of the other-regarding preference
 parameter $F$. For this panel, we consider an example in which average payoff in a group
 is maximized by the all-cooperator composition. 
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
gamma = 2.5
theta = 2.0

F_any_coop = (-beta) / ( gamma - beta)
F_all_coop = (-alpha - beta) / (gamma + alpha - beta) 

def G_F(x, F, gamma, alpha, beta):
	return gamma * x + alpha * (x ** 2.0)
	
def threshold_lambda(x_eq, F, alpha, beta, gamma):
	numerator = (-alpha - beta - F * (gamma - beta)) * theta
	if r > x_any_coop and r < x_all_coop:
		denominator = G_F(1,F,gamma,alpha,beta) - G_F(x_eq,F,gamma,alpha,beta)
	elif r < x_any_coop:
		denominator = G_F(1,F,gamma,alpha,beta) - G_F(0.0,F,gamma,alpha,beta)
	else:
		return 0.0
	return numerator / denominator


"""
Function for plotting unnormalized version of steady state density for various values 
of the weight parameter $F$ corresponding to other-regarding preference. 
"""
	
def steady_density(x,lamb, F, beta, alpha, gamma, theta):
	a = alpha
	b = beta + (F / (1.0 + F)) * (gamma - 2.0 * beta)
	gamma = gamma

	
	if F > F_any_coop and x <= -(b/a):
		return 0.0
	elif F < F_any_coop:
		x_power = (1.0 / np.abs(b)) * (lamb * (gamma + a) - (np.abs(b) - a) * theta) - 1.0
		bminusax_power = -1.0 * (1.0 / np.abs(b)) * (lamb * (gamma + np.abs(b) + a) + a * theta) - 1.0
		return (x ** x_power) * ((1.0 - x) ** (theta - 1.0)) * ((np.abs(b) - a * x) ** bminusax_power)
	elif F_any_coop < F < F_all_coop:
		x_power = (1.0 / b) * (lamb*(np.abs(a) - gamma) + (np.abs(a) - b) * theta) - 1.0
		bminusax_power = (1.0 / b) * (lamb * (gamma - np.abs(a) - b) - np.abs(a) * theta) - 1.0
		return (x ** x_power) * ((1.0 - x) ** (theta - 1.0)) * ((np.abs(a) * x - b) ** bminusax_power)
		
		
steady_state_vec = np.vectorize(steady_density)



"""
Plotting the normalized steady state densities for various weights other-regarding
preference $F$ and fixed game parameters and Hoelder exponent near $x=1$.
"""

lamb = 5.
F = 0.0

x_step = 0.001
x_range = np.arange(x_step,1.0+x_step,x_step)

holderr0 = steady_state_vec(x_range,lamb,F,beta,alpha,gamma,theta)
holderr0int = spi.simps(holderr0,x_range)
holderr0 = holderr0 / holderr0int
#print  spi.simps(holderr0,x_range)
#print holderr0
plt.plot(x_range,holderr0, lw = 4., color = plt.cm.YlOrRd(0.), label = r"$F = 0.0$")

F = 0.12

x_range = np.arange(0.000,1.0+x_step,x_step)

holderr2 = steady_state_vec(x_range,lamb,F,beta,alpha,gamma,theta)
holderr2int = spi.simps(holderr2,x_range)
holderr2 = holderr2/ holderr2int
#print  spi.simps(holderr2,x_range)
#print holderr2
plt.plot(x_range,holderr2, lw = 4., color = plt.cm.YlOrRd(0.2), label = r"$F = 0.12$")


F = 0.24

holderr4 = steady_state_vec(x_range,lamb,F,beta,alpha,gamma,theta)
holderr4int = spi.simps(holderr4,x_range)
holderr4 = holderr4/ holderr4int
#print  spi.simps(holderr4,x_range)
#print holderr4
plt.plot(x_range,holderr4, lw = 4., color = plt.cm.YlOrRd(0.4), label = r"$F = 0.24$")

F = 0.32

holderr6 = steady_state_vec(x_range,lamb,F,beta,alpha,gamma,theta)
holderr6int = spi.simps(holderr6,x_range)
holderr6 = holderr6/ holderr6int
#print  spi.simps(holderr4,x_range)
#print holderr4
plt.plot(x_range,holderr6, lw = 4., color = plt.cm.YlOrRd(0.6), label = r"$F = 0.32$")




F = 0.48

holderr8 = steady_state_vec(x_range,lamb,F,beta,alpha,gamma,theta)
holderr8int = spi.simps(holderr8,x_range)
holderr8 = holderr8/ holderr8int
#print  spi.simps(holderr4,x_range)
print holderr8
plt.plot(x_range,holderr8, lw = 4., color = plt.cm.YlOrRd(0.8), label = r"$F = 0.48$")



F = 0.6

holderr10 = steady_state_vec(x_range,lamb,F,beta,alpha,gamma,theta)
holderr10int = spi.simps(holderr10,x_range)
holderr10 = holderr10/ holderr10int
#print  spi.simps(holderr4,x_range)
#print holderr10
plt.plot(x_range,holderr10, lw = 4., color = plt.cm.YlOrRd(1.0), label = r"$F = 0.6$")


plt.axis([0.0,1.0,0.0,16.0])
plt.legend(loc = "upper left")

plt.axvline(x=gamma / (-2.0 * alpha), color = 'Gray', lw = 4., ls = '--')

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.xlabel(r"Fraction of Cooperators ($x$)", fontsize = 20., labelpad = 20.)
plt.ylabel(r"Probability Density", fontsize = 20.)



plt.tight_layout()

script_folder = os.getcwd()
mechanism_folder = os.path.dirname(script_folder)

plt.savefig(mechanism_folder + "/Figures/Fdensityplotnoshadow.png")

plt.show()