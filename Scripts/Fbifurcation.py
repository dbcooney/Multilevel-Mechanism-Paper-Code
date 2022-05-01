"""
Script used to generate Figure 4.1, presenting a bifurcation diagram for the
within-group dynamics and the group composition maximizing average payoff under our
model of other-regarding preferences.
"""

import matplotlib.pyplot as plt
import numpy as np
import os

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)



"""
Game-theoretic parameters for our model. 
"""

alpha = -1.0
beta = -1.0
gamma = 1.5


step_size = 0.010


"""
Calculating the interior within-group equilibrium and the composition maximizing
average payoff of group members under the model for other-regarding preferences.
"""

def interior_equilibrium(F,alpha,beta,gamma):
	numerator = beta + (F / (1.0 + F)) * (gamma - 2.0 * beta)
	denominator = -alpha
	return numerator / denominator

def group_maximizer(alpha,gamma):
	interior = gamma / (-2.0 * alpha)
	maximizer = np.minimum(interior,1.0)
	return maximizer


F_any_coop = (-beta) / ( gamma - beta)
F_all_coop = (-alpha - beta) / (gamma + alpha - beta) 


F_between_group = gamma / (-2.0 * alpha)

step_size = 0.010

"""
Calculating the stable and unstable curves for the bifurcation diagram and the curve
describing the composition maximizing collective payoff in groups, described as
functions of the parameter $F$ characterizing the weight placed on opponent payoff. 
"""

r_holder = np.arange(0.0,1.0 + step_size,step_size)
group = [group_maximizer(alpha,gamma)] * len(r_holder)
plt.plot(r_holder,group, lw = 4., color = 'g', ls = '-')

zero_stable_holder = np.arange(0.0,F_any_coop+step_size,step_size)
zero_stable = np.zeros(len(zero_stable_holder))

zero_unstable_holder = np.arange(F_any_coop,1.0+step_size,step_size)
zero_unstable = np.zeros(len(zero_unstable_holder))

plt.plot(zero_stable_holder,zero_stable,lw=4.,ls = '-', color = 'b')
plt.plot(zero_unstable_holder,zero_unstable,lw=4.,ls = '--', color = 'b')

one_unstable_holder = np.arange(0.0,F_all_coop+step_size,step_size)
one_unstable = np.ones(len(one_unstable_holder))

one_stable_holder = np.arange(F_all_coop,1.0+step_size,step_size)
one_stable = np.ones(len(one_stable_holder))

plt.plot(zero_stable_holder,zero_stable,lw=4.,ls = '-', color = 'b')
plt.plot(zero_unstable_holder,zero_unstable,lw=4.,ls = '--', color = 'b')

plt.plot(one_stable_holder,one_stable,lw=4.,ls = '-', color = 'b')
plt.plot(one_unstable_holder,one_unstable,lw=4.,ls = '--', color = 'b')

print zero_unstable_holder


interior_eq_holder = np.arange(F_any_coop,F_all_coop+step_size,step_size)
interior_eq = interior_equilibrium(interior_eq_holder,alpha,beta,gamma)

interior_eq_holder[-1] = F_all_coop
interior_eq[-1] = interior_equilibrium(F_all_coop,alpha,beta,gamma)

plt.plot(interior_eq_holder,interior_eq, lw = 4., ls = '-', color = 'b')

print interior_eq_holder
print interior_equilibrium(F_all_coop,alpha,beta,gamma)


"""
Plotting the curves calculated above. 
"""


plt.axvline(F_any_coop, lw = 3., ls = '-.',color = 'gray')
if gamma + 2.0 * alpha < 0:
	plt.axvline(1.0 / F_all_coop, lw = 3., ls = '-.',color = 'gray')
else:
	plt.axvline(F_all_coop, lw = 3., ls = '-.',color = 'gray')



print r_holder
print group


plt.annotate(r"$x^{*}_{F}$", xy = (0.13,group_maximizer(alpha,gamma) - 0.09), fontsize = 24.)




plt.annotate(r"$F_{W}^{s}$", xy = (F_any_coop - 0.075,0.25), fontsize = 24., rotation = 0)
if gamma + 2.0 * alpha < 0:
	plt.annotate(r"$F_{W}^{AC}$", xy = (1/(F_all_coop) + 0.02,0.25), fontsize = 24., rotation = 0)
	plt.annotate(r"$x^{eq}_{F}$", xy = (0.56,0.335), fontsize = 24.)
	
else:
	plt.annotate(r"$F_{W}^{a}$", xy = (F_all_coop + 0.02,0.25), fontsize = 24., rotation = 0)
	plt.annotate(r"$x^{eq}_{F}$", xy = (0.39,0.41), fontsize = 24.)
	

plt.axis([0.0,1.0,-0.02,1.02])

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.xlabel(r"Opponent Payoff Weight ($F$)", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Fraction of Cooperators ($x$)", fontsize = 20.)

plt.tight_layout()

script_folder = os.getcwd()
mechanism_folder = os.path.dirname(script_folder)

if gamma == 1.5:
	plt.savefig(mechanism_folder + "/Figures/Fbifurcationshadow.png")
if gamma == 2.5:
	plt.savefig(mechanism_folder + "/Figures/Fbifurcationnoshadow.png")




plt.show()