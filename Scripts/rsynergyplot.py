"""
Script used for generating Figure 6.1, which illustrates the synergy between within-group
assortment and multilevel selection for promoting cooperation. We plot the regions of
parameter space in which cooperation can be achieved via multilevel selection for a given
level of between-group selection $\lambda$ and assortment probability $r$. 
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

alpha = -1.0
beta = -1.0
gamma = 1.5
theta = 2.0

r_any_coop = -beta / (gamma + alpha - beta)
r_all_coop = (-alpha - beta) / (gamma - beta) 
lamb_PD = (-(beta + alpha) * theta) / (gamma + alpha)

r_step = 0.001
r_range = np.arange(0.,1. + r_step, r_step)
r_short_range = np.arange(0.,r_any_coop+r_step,r_step)
r_medium_range = np.arange(0.,r_all_coop,r_step)



	
def G_r(alpha,beta,gamma,r,x):
	return r * (gamma + alpha)* x + (1. - r) * (gamma * x + alpha * x * x)
	
def x_eq(alpha,beta,gamma,r):
	return -beta / alpha + (r / (1.0 - r)) * ((gamma + alpha) / (-alpha))
	
	
"""
Calculating the threshold between-group selection strength $\lambda^*_r$ required to
achieve cooperation via multilevel selection for a given within-group assortment
probability $r$. 
"""	
	
def lamb_r(alpha,beta,gamma,r,theta):
	num = (-(beta + alpha) - r * (gamma - beta)) * theta
	if  r < -beta / (gamma + alpha - beta):
		#num = (-(beta + alpha) - r * (gamma - beta)) * theta
		denom = gamma + alpha
		print num/denom
		return num / denom
	elif r < (-alpha - beta) / (gamma - beta):
		denom = gamma + alpha - G_r(alpha,beta,gamma,r,x_eq(alpha,beta,gamma,r))
		print num/denom
		print r
		return num / denom
	else:
		return -10.
		
lamb_r_vec = np.vectorize(lamb_r)


plt.axvline(x = r_any_coop, lw = 5., ls = '--', color = 'k')
plt.axvline(x = r_all_coop, lw = 5., ls = '--', color = 'k')
plt.plot(r_short_range,lamb_PD * np.ones(len(r_short_range)),lw = 5., ls ='--', color = 'k')
print r_any_coop


lamb_r_0 = lamb_r(alpha,beta,gamma,0.,theta) * np.ones(len(r_short_range))


"""
Plotting the region of parameter space in which cooperation can be supported by a
combination of within-group assortment and multilevel selection, even when neither 
mechanism can support cooperation when acting on their own. 
"""

plt.fill_between(r_short_range,lamb_r_vec(alpha,beta,gamma,r_short_range,theta),lamb_r_0, color = 'green', alpha = 0.5)


plt.plot(r_short_range,lamb_r_vec(alpha,beta,gamma,r_short_range,theta), lw = 5., ls = '--', color = 'g')

plt.axis([0.,0.78,0.,12.])

plt.xlabel(r"Assortment Probability ($r$)", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Between-Group Selection Strength ($\lambda$)", fontsize = 20.)

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.annotate(r"$\lambda^*_{PD}$", xy = (0.15,8.3), fontsize = 20.)
plt.annotate(r"$\lambda^*_{r}$", xy = (0.17,5.), fontsize = 20.)
plt.annotate(r"$r_W^s$", xy = (0.61,10), fontsize = 20.)

plt.annotate(r"I", xy = (0.4,10), fontsize = 24)
plt.annotate(r"II", xy = (0.4,6), fontsize = 24)
plt.annotate(r"III", xy = (0.4,2), fontsize = 24)
plt.annotate(r"IV", xy = (0.72,6), fontsize =24)

plt.tight_layout()


script_folder = os.getcwd()
mechanism_folder = os.path.dirname(script_folder)

plt.savefig(mechanism_folder + "/Figures/rsynergyplot.png")


plt.show()