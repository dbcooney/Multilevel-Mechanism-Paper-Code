"""
This script is used to generate Figure 3.1, presenting a bifurcation diagram for the
within-group dynamics and the group composition maximizing average payoff under the
assortment model.
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
Game-theoretic parameters for our model. We use beta = -1.0 for left panel, and 
beta = -0.25 in right panel.
"""

alpha = -1.0
beta = -1.0
gamma = 1.5


"""
Calculating the interior within-group equilibrium and the composition maximizing
average payoff of group members under the assortment process.
"""

def interior_equilibrium(r,alpha,beta,gamma):
	numerator = beta + r * (gamma + alpha - beta)
	denominator = -alpha * (1 - r)
	return numerator / denominator

def group_maximizer(r,alpha,gamma):
	interior = ((1 - r) * gamma + r * (gamma + alpha))/ (- 2 * alpha * (1-r))
	maximizer = np.minimum(interior,1.0)
	return maximizer
	
def Group_payoff(x,alpha,gamma,r):
	return r * (gamma + alpha) * x + (1.0 - r) * (gamma * x + alpha * x * x)
	


r_any_coop = (-beta) / ( gamma + alpha - beta)
r_all_coop = (-alpha - beta) / (gamma - beta) 

r_between_group = (gamma + 2.0 * alpha) / (alpha)

step_size = 0.010


"""
Calculating the stable and unstable curves for the bifurcation diagram and the curve
describing the composition maximizing collective payoff in groups, described as
functions of the assortment probability. 
"""

r_short_holder = np.arange(0.0, (-alpha - beta) / (gamma - beta),step_size)
r_holder = np.arange(0.0,1.0 + step_size,step_size)
group = group_maximizer(r_holder,alpha,gamma)
plt.plot(r_holder,group, lw = 4., color = 'g', ls = '-')

zero_stable_holder = np.arange(0.0,r_any_coop+step_size,step_size)
zero_stable = np.zeros(len(zero_stable_holder))

zero_unstable_holder = np.arange(r_any_coop,1.0+step_size,step_size)
zero_unstable = np.zeros(len(zero_unstable_holder))

plt.plot(zero_stable_holder,zero_stable,lw=4.,ls = '-', color = 'b')
plt.plot(zero_unstable_holder,zero_unstable,lw=4.,ls = '--', color = 'b')

one_unstable_holder = np.arange(0.0,r_all_coop+step_size,step_size)
one_unstable = np.ones(len(one_unstable_holder))

one_stable_holder = np.arange(r_all_coop,1.0+step_size,step_size)
one_stable = np.ones(len(one_stable_holder))

plt.plot(zero_stable_holder,zero_stable,lw=4.,ls = '-', color = 'b')
plt.plot(zero_unstable_holder,zero_unstable,lw=4.,ls = '--', color = 'b')

plt.plot(one_stable_holder,one_stable,lw=4.,ls = '-', color = 'b')
plt.plot(one_unstable_holder,one_unstable,lw=4.,ls = '--', color = 'b')

print zero_unstable_holder

interior_eq_holder = np.arange(r_any_coop,r_all_coop+step_size,step_size)
interior_eq = interior_equilibrium(interior_eq_holder,alpha,beta,gamma)

interior_eq_holder[-1] = r_all_coop
interior_eq[-1] = interior_equilibrium(r_all_coop,alpha,beta,gamma)


"""
Plotting the curves calculated above. 
"""

plt.plot(interior_eq_holder,interior_eq, lw = 4., ls = '-', color = 'b')

print interior_eq_holder
print interior_equilibrium(r_all_coop,alpha,beta,gamma)

plt.axvline(r_any_coop, lw = 3., ls = '-.',color = 'gray')
plt.axvline(r_all_coop, lw = 3., ls = '-.',color = 'gray')
plt.axvline(r_between_group, lw = 3., ls = '-.',color = 'gray')


print r_holder
print group

plt.annotate(r"$x^{eq}_{r}$", xy = (0.71,0.7), fontsize = 20.)
plt.annotate(r"$x^{*}_{r}$", xy = (0.13,0.81), fontsize = 20.)
#plt.annotate(r"$x^{eq}_{r}$", xy = (0.63,0.81), fontsize = 20.)


plt.annotate(r"$r_{B}$", xy = (r_between_group - 0.05,0.25), fontsize = 20., rotation = 90 )
plt.annotate(r"$r_{W}^{s}$", xy = (r_any_coop - 0.05,0.25), fontsize = 20., rotation = 90 )
plt.annotate(r"$r_{W}^{a}$", xy = (r_all_coop + 0.02,0.25), fontsize = 20., rotation = 90 )

plt.axis([0.0,1.0,-0.02,1.02])

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.xlabel(r"Assortment Probability ($r$)", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Fraction of Cooperators ($x$)", fontsize = 20.)

plt.tight_layout()



script_folder = os.getcwd()
mechanism_folder = os.path.dirname(script_folder)

if beta == -1.0:
	plt.savefig(mechanism_folder + "/Figures/rbifurcationtype1.png")
elif beta == -0.25:
	plt.savefig(mechanism_folder + "/Figures/rbifurcationtype2.png")



plt.show()