"""
This script is used to generate Figure 4.3, which provides a comparison between the 
group composition maximizing average payoff and the modal fraction of cooperations in the
steady state density, plotted as a function of the weight parameter $F$ characterizing
other-regarding preference. 
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
Parameters characterizing game payoffs and initial distribution.
"""
alpha = -1.0
beta = -1.0
gamma = 1.5
theta = 2.0

lamb = 15.


F_step = 0.01
F_range = np.arange(0.,1.0+F_step,F_step)



def F_any_coop(alpha,beta,gamma):
	return  (-beta) / ( gamma - beta)
	
def F_all_coop(alpha,beta,gamma):
	return (-alpha - beta) / (gamma + alpha - beta)
	 
def F_AC(alpha,beta,gamma):
	return 1. / (F_all_coop(alpha,beta,gamma))

def x_eqF(alpha,beta,gamma,F):
	return (beta + (F / (1. + F)) * (gamma - 2. * beta) ) / (-alpha)
	
	
"""
Defining threshold to have a steady state with modal composition exceeding the
level of cooperation achieved at the within-group equilibrium. 
"""	
		
	
def peak_threshold(alpha,beta,gamma,F,theta):
	betaF = beta + (F / (1. + F)) * (gamma - 2. * beta)
	
	if F < F_any_coop(alpha,beta,gamma):
		num = -(betaF + alpha) * theta
		denom = gamma + alpha
		return (num - betaF) / denom
	elif F < F_AC(alpha,beta,gamma) and F < F_all_coop(alpha,beta,gamma):
		num = -alpha * theta
		denom = gamma + alpha - betaF
		return (num + betaF) / denom
	else:
		return 0.
		
peak_threshold_vec = np.vectorize(peak_threshold)


"""
Formulas for modal composition at steady state when modified multilevel dynamics are in
the PD and HD regimes under our model of other-regarding preference. 
"""


def PD_peak(alpha,beta,gamma,lamb,theta,F):
	betaF = beta  + (F / (1. + F)) * (gamma - 2. * beta)
	if lamb < peak_threshold_vec(alpha,beta,gamma,F,theta):
		return 0.
	elif lamb > peak_threshold_vec(alpha,beta,gamma,F,theta) and F < min(F_AC(alpha,beta,gamma),F_all_coop(alpha,beta,gamma)):
		betaF = beta  + (F / (1. + F)) * (gamma - 2. * beta)
		a = -(3. + lamb) * alpha
		b = 2. * (alpha - betaF) - lamb * gamma
		c = lamb * (gamma + alpha) + (betaF + alpha) * theta + betaF
		print a, b, c
		return (-b - np.sqrt((b**2) - 4. * a * c)) / (2. * a)
	
def HD_peak(alpha,beta,gamma,lamb,theta,F):
	betaF = beta  + (F / (1. + F)) * (gamma - 2. * beta)
	if lamb < peak_threshold_vec(alpha,beta,gamma,F,theta):
		return x_eqF(alpha,beta,gamma,F)
	elif lamb >= peak_threshold_vec(alpha,beta,gamma,F,theta) and F > F_any_coop(alpha,beta,gamma):
		a = -(3. + lamb) * alpha
		b = 2. * (alpha - betaF) - lamb * gamma
		c = lamb * (gamma + alpha) + (betaF + alpha) * theta + betaF
		return (-b - np.sqrt((b**2) - 4. * a * c)) / (2. * a)
		
		
"""
Calculating the modal level of cooperation at steady state, incorporating the different
cases for the other-regarding preference parameter $F$ in which the dynamics resemble
PD games, HD games, and anti-coordination games.
"""
		
def all_peak(alpha,beta,gamma,lamb,theta,F):
	if F <= F_any_coop(alpha,beta,gamma):
		return PD_peak(alpha,beta,gamma,lamb,theta,F)
	elif F >= F_AC(alpha,beta,gamma):
		return x_eqF(alpha,beta,gamma,F)
	elif F <= F_all_coop(alpha,beta,gamma):
		return HD_peak(alpha,beta,gamma,lamb,theta,F)
	else:
		return 1.
		
all_peak_vec = np.vectorize(all_peak)
		
	
PD_peak_vec = np.vectorize(PD_peak)
print PD_peak_vec(alpha,beta,gamma,lamb,theta,0.1)


"""
Plotting modal composition at steady state and composition that maximizes collective
payoff. 
"""


plt.plot(F_range,all_peak_vec(alpha,beta,gamma,lamb,theta,F_range), lw = 5., color = 'g')
plt.plot(F_range,gamma / (-2. * alpha) * np.ones(len(F_range)), lw = 5., color = 'b')
plt.plot(F_range,x_eqF(alpha,beta,gamma,F_range), lw = 5., ls = '--', color = 'k')

plt.axvline(x = F_any_coop(alpha,beta,gamma), lw = 5., ls = '--', color = 'gray', alpha = 0.8)
plt.axvline(x = F_AC(alpha,beta,gamma), lw = 5., ls = '--', color = 'gray', alpha = 0.8)


plt.annotate(r"$x^*_F$",xy = (0.2, 0.77), fontsize = 20.)
plt.annotate(r"$\hat{x}^{\lambda}_F$",xy = (0.2, 0.25), fontsize = 20.)
plt.annotate(r"$F_W^s$", xy = (0.32,0.6), fontsize = 20.)
plt.annotate(r"$F_W^{AC}$", xy = (0.65,0.6), fontsize = 20.)
plt.annotate(r"$x_{eq}^{F}$", xy = (0.52,0.1), fontsize = 20.)
plt.axis([0.,1.,0.,1.])

plt.xlabel(r"Fraternity Parameter ($F$)",fontsize = 20.,labelpad = 10.)
plt.ylabel(r"Fraction of Cooperators ($x$)", fontsize = 20.)

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.tight_layout()

plt.savefig("Fmodelpeakplot.png")

plt.show()



script_folder = os.getcwd()
mechanism_folder = os.path.dirname(script_folder)

plt.savefig(mechanism_folder + "/Figures/Fmodelpeakplot.png")

		

