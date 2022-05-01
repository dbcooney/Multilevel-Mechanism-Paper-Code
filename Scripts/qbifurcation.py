"""
This script is used to generate Figure 5.1, presenting a bifurcation diagram for the
within-group dynamics and the group composition maximizing average payoff under our
simplified model for indirect reciprocity.
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
alpha = -1.1
beta = -1.0
gamma = 1.5



"""
Calculating the interior within-group equilibrium and the composition maximizing
average payoff of group members under the indirect reciprocity model.
"""

def interior_equilibrium(q,alpha,beta,gamma):
	numerator = -beta * (1 - q)
	denominator = q * gamma + alpha
	return numerator / denominator

def group_maximizer(q,alpha,gamma):
	interior = ((1-q) * gamma)/ (- 2 * (alpha + q * gamma))
	#maximizer = np.minimum(interior,1.0)
	if q < (gamma + 2 * alpha) / (-gamma):
		return interior
	else:
		return 1.0


q_all_coop = (-alpha - beta) / (gamma - beta) 
q_between_group = (gamma + 2.0 * alpha) / (-gamma)

step_size = 0.010

"""
Calculating the stable and unstable curves for the bifurcation diagram and the curve
describing the composition maximizing collective payoff in groups, described as
functions of the defector detection probability $q$. 
"""

q_holder = np.arange(0.0,1.0 + step_size,step_size)
group_maximizer_vec = np.vectorize(group_maximizer)
group = group_maximizer_vec(q_holder,alpha,gamma)
plt.plot(q_holder,group, lw = 4., color = 'g', ls = '-')

zero_stable_holder = np.arange(0.0,1.0+step_size,step_size)
zero_stable = np.zeros(len(zero_stable_holder))


plt.plot(zero_stable_holder,zero_stable,lw=4.,ls = '-', color = 'b')
#plt.plot(zero_unstable_holder,zero_unstable,lw=4.,ls = '--', color = 'b')

one_unstable_holder = np.arange(0.0,q_all_coop+step_size,step_size)
one_unstable = np.ones(len(one_unstable_holder))

one_stable_holder = np.arange(q_all_coop,1.0+step_size,step_size)
one_stable = np.ones(len(one_stable_holder))

plt.plot(zero_stable_holder,zero_stable,lw=4.,ls = '-', color = 'b')
#plt.plot(zero_unstable_holder,zero_unstable,lw=4.,ls = '--', color = 'b')

plt.plot(one_stable_holder,one_stable,lw=4.,ls = '-', color = 'b')
plt.plot(one_unstable_holder,one_unstable,lw=4.,ls = '--', color = 'b')

#print zero_unstable_holder

interior_eq_holder = np.arange(q_all_coop,1.0+step_size,step_size)
interior_eq = interior_equilibrium(interior_eq_holder,alpha,beta,gamma)

#interior_eq_holder[-1] = 1.0
#interior_eq[-1] = interior_equilibrium(q_all_coop,alpha,beta,gamma)

"""
Plotting the curves calculated above. 
"""


plt.plot(interior_eq_holder,interior_eq, lw = 4., ls = '--', color = 'b')

print interior_eq_holder
print interior_equilibrium(q_all_coop,alpha,beta,gamma)

#plt.axvline(x__coop, lw = 3., ls = '-.',color = 'gray')
plt.axvline(q_all_coop, lw = 3., ls = '-.',color = 'gray')
plt.axvline(q_between_group, lw = 3., ls = '-.',color = 'gray')


print q_holder
print group

plt.annotate(r"$x^{eq}_{q}$", xy = (0.875,0.7), fontsize = 20.)
plt.annotate(r"$x^{*}_{q}$", xy = (0.13,0.75), fontsize = 20.)



plt.annotate(r"$q_{B}$", xy = (q_between_group - 0.055,0.25), fontsize = 20., rotation = 90 )
plt.annotate(r"$q_{W}^{a}$", xy = (q_all_coop - 0.055,0.25), fontsize = 20., rotation = 90 )

plt.axis([0.0,1.0,-0.02,1.02])

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.xlabel(r"Detection Probability ($q$)", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Fraction of Cooperators ($x$)", fontsize = 20.)

plt.tight_layout()


script_folder = os.getcwd()
mechanism_folder = os.path.dirname(script_folder)

plt.savefig(mechanism_folder + "/Figures/qbifurcation.png")


plt.show()