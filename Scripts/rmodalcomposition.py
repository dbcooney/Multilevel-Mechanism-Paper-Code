"""
This script is used to generate Figure 3.3, which provides a comparison between the 
group composition maximizing average payoff and the modal fraction of cooperations in the
steady state density, plotted as a function of the assortment probability $r$. We consider
both the case of a finite stength of between-group selection (left panel) and in limit
of infinitely strong between-group competition (right panel). 
"""

import matplotlib.pyplot as plt
import numpy as np
import os

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

alpha = -1.0
beta = -1.0
gamma = 1.5
theta = 2.0

r_any_coop = (-beta) / ( gamma + alpha - beta)
r_all_coop = (-alpha - beta) / (gamma - beta) 

def G_r(x, r, gamma, alpha, beta):
	return (gamma + alpha * r) * x + (1.0 - r) * alpha * (x ** 2.0)
	
"""
Defining threshold selection strengths to achieve cooperation at steady state, and the 
corresponding threshold to have a steady state with modal composition exceeding the
level of cooperation achieved at the within-group equilibrium. 
"""	
	
def threshold_lambda(x_eq, r, alpha, beta, gamma):
	numerator = (-alpha - beta - r * (gamma - beta)) * theta
	if r > r_any_coop and r < r_all_coop:
		denominator = G_r(1,r,gamma,alpha,beta) - G_r(x_eq,r,gamma,alpha,beta)
	elif r < r_any_coop:
		denominator = G_r(1,r,gamma,alpha,beta) - G_r(0.0,r,gamma,alpha,beta)
	else:
		return 0.0
	return numerator / denominator
	
def peak_threshold_lambda(x_eq, r, alpha, beta, gamma):
	numerator = (-alpha - beta - r * (gamma - beta)) * theta
	b = beta + r * (gamma + alpha - beta)
	if r > r_any_coop and r < r_all_coop:
		denominator = G_r(1,r,gamma,alpha,beta) - G_r(x_eq,r,gamma,alpha,beta)
		peak_cushion = b * x_eq * (1.0 - x_eq)
	elif r < r_any_coop:
		denominator = G_r(1,r,gamma,alpha,beta) - G_r(0.0,r,gamma,alpha,beta)
		peak_cushion = b
	else:
		return 0.0
	return (numerator + peak_cushion) / denominator
	
def fitness(lamb,r,theta,alpha,beta,gamma):
	return gamma + alpha + (theta / lamb) * (alpha + beta + r* (gamma - beta))
	
lamb_holder = []


r_list =  np.arange(0.0,1.01,0.01)
r_list_list = [r for r in r_list]


plt.figure(1)


"""
Computing modal composition of cooperators at density steady state as a function of
the assortment probability $r$ and for fixed finite between-group selection strength 
$\lambda$. 
"""


lamb = 8.0


def best_group(r,alpha,gamma):
	if r < 2.0 - (gamma / np.abs(alpha)):
		return (gamma - np.abs(alpha) * r) / (2.0 * (1.0 - r) * np.abs(alpha))
	else:
		return 1.0
		
def density_peak(lamb, r, beta, alpha, gamma, theta):
	a = (1.0 - r) * alpha
	b = beta + r * (gamma + alpha - beta)
	
	x_eq = -beta / alpha + (r / (1.0 - r)) * ((gamma + alpha) / (-alpha))
	
	E = 1.0 / ( x_eq * (1.0 - x_eq))
	Delta_pi = -(b + a)
	G_1 = G_r(1.0,r,gamma,alpha,beta)
	G_eq = G_r(x_eq,r,gamma,alpha,beta)
	Delta_G = G_1 - G_eq
	
	constant = lamb * G_1 - Delta_pi * theta + b
	linear = (a/b - 1.0) * lamb * G_1 + lamb * Delta_G * E + (1.0 - (a/b)) * Delta_pi * theta - (2.0 - theta) * b  + 2.0 * a - Delta_pi * theta * E
	quadratic = -lamb * G_1 * (a / b) - lamb * Delta_G * E - (3.0 - theta) * a  + Delta_pi * theta  * (a/b +  E)
	
	root_plus = (-linear + np.sqrt(linear**(2.0) - 4.0 * constant * quadratic)) / (2.0 * quadratic)
	root_minus = (-linear - np.sqrt(linear**(2.0) - 4.0 * constant * quadratic)) / (2.0 * quadratic)
	
	if r < r_any_coop and lamb < peak_threshold_lambda(x_eq, r, alpha, beta, gamma) + np.abs(b):
		return 0.0	
	elif r > r_any_coop and r < r_all_coop and lamb < peak_threshold_lambda(x_eq, r, alpha, beta, gamma):
		return x_eq
	elif r >= r_all_coop:
		return 1.0
	elif quadratic != 0.0:
		return root_minus
	else:
		return constant / (-linear)
	
		
max_list = [best_group(r,alpha,gamma) for r in r_list]
peak_list = [density_peak(lamb,r,beta,alpha,gamma,theta) for r in r_list]
#	x_eq = -beta / alpha + (r / (1.0 - r)) * ((gamma + alpha) / (-alpha))
threshold_list = [threshold_lambda(-beta / alpha + (r / (1.0 - r)) * ((gamma + alpha) / (-alpha)),r,alpha,beta,gamma) for r in r_list]
x_eq_list = [-beta / alpha + (r / (1.0 - r)) * ((gamma + alpha) / (-alpha)) for r in r_list]

"""
Plotting modal composition at steady state and composition that maximizes collective
payoff. 
"""

plt.plot(r_list,max_list,lw = 3., label = r"$x^*$")
plt.plot(r_list,peak_list,lw = 3., label = r"$\hat{x}^{\lambda}$")
plt.plot(r_list,x_eq_list, lw = 3., ls = '--', color = 'k', label = r"$x_{eq}^r$")


plt.axis([0.0,1.0,0.0,1.01])

plt.xlabel(r"Assortment Probability ($r$)", fontsize = 20.,labelpad =20)
plt.ylabel(r"Fraction of Cooperators ($x$)", fontsize = 20.)

plt.axvline(r_any_coop, lw = 3.0, color = "gray", ls = "--")
plt.axvline(r_all_coop, lw = 3.0, color = "gray", ls = "--")


plt.legend(loc = "lower right", prop={'size': 16})

plt.tight_layout()



script_folder = os.getcwd()
mechanism_folder = os.path.dirname(script_folder)

plt.savefig(mechanism_folder + "/Figures/rpeakfixedlambdaunghostfirst.png")




plt.figure(2)

"""
Creating similar plots in the limit of infinite strength of between-group competition 
($\lambda \to \infty$). 
"""

def infinite_density_peak(lamb, r, beta, alpha, gamma, theta):
	a = (1.0 - r) * alpha
	b = beta + r * (gamma + alpha - beta)
	
	x_eq = -beta / alpha + (r / (1.0 - r)) * ((gamma + alpha) / (-alpha))
	
	E = 1.0 / (x_eq * (1.0 - x_eq))
	Delta_pi = -(b + a)
	G_1 = G_r(1.0,r,gamma,alpha,beta)
	G_eq = G_r(x_eq,r,gamma,alpha,beta)
	Delta_G = G_1 - G_eq
	
	possible_peak = -G_1 / ((a/b) * G_1 + Delta_G * E)
	if r >= r_all_coop:
		return 1.0
	if possible_peak < 1.0:
		return possible_peak
	else:
		return 1.0
		
infinity_peak_list = [infinite_density_peak(lamb,r,beta,alpha,gamma,theta) for r in r_list]
	
plt.plot(r_list,infinity_peak_list,lw = 4., color ='g', label = r"Modal Group Type as $\lambda \to \infty$")
plt.plot(r_list,max_list,lw = 4., color = 'b', label = r'Group with Maximal Average Payoff $x^*$')
	
plt.axis([0.0,1.0,0.0,1.01])

plt.xlabel(r"Assortment Probability ($r$)", fontsize = 20.,labelpad =20)
plt.ylabel(r"Fraction of Cooperators ($x$)", fontsize = 20.)

plt.axvline(r_any_coop, lw = 3.0, color = "gray", ls = "--")
plt.axvline(r_all_coop, lw = 3.0, color = "gray", ls = "--")

plt.legend(loc = "lower right")

plt.tight_layout()

plt.savefig(mechanism_folder + "/Figures/rpeakinfinitelambdaunghostfirst.png")
		

plt.show()

\