"""
This script is used to generate Figure 5.3, which presents the average payoff achieved at
steady state for the multilevel dynamics with direct reciprocity, plotted as a function
of the continuation probability / discount rate $\delta$ and the relative strength of 
between-group competition $\lambda$. 
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
Parameters describing game payoffs and the initial density. 
"""

alpha = -1.0
beta = -1.0
gamma = 1.5
theta = 2.0


"""
Calculation ofquantities for the multilevel dynamics with repeated interactions, 
including the average payoff at steady state.
"""

delta_all_coop = (-alpha - beta) / (gamma - beta)

def G_delta(x, delta, gamma, alpha, beta):
	return gamma *  x + (alpha + (delta / (1 - delta)) * (gamma + alpha)) * (x ** 2.0)
	
def threshold_lambda(x_eq, delta, alpha, beta, gamma):
	numerator = (-alpha - beta - (delta / (1 - delta)) * (gamma + alpha)) * theta
	denominator = G_delta(1,delta,gamma,alpha,beta)
	return numerator / denominator
	
def fitness(lamb,delta,theta,alpha,beta,gamma):
	return (1.0 / (1.0 - delta)) * (gamma + alpha) + (theta / lamb) * (beta + alpha + (delta / (1 - delta)) * (gamma + alpha))
	
lamb_holder = []


delta_list =  np.arange(0.0,0.825,0.01)
delta_list_list = [delta for delta in delta_list]
for delta in delta_list:
	if delta < delta_all_coop:
		x_eq = 0.0
	else:
		x_eq = 1.0
	print x_eq


print lamb_holder
print G_delta(delta_all_coop,delta,gamma,alpha,beta)

plt.figure(2)


threshold_holder = []


lamb_list = np.arange(0.0,40.0,0.5)
lamb_list_list = [lamb for lamb in lamb_list]
steady_fitness = np.zeros([len(lamb_list),len(delta_list)])

"""
Constructing array of average payoffs at steady states for various relative levels of 
between-group competition $\lambda$ and discount rate / continuating probability
parameter $\delta$. 
"""


for delta in delta_list:
	if delta < delta_all_coop:
		x_eq = 0.0
	else:
		x_eq = 1.0
	for lamb in lamb_list:
		if lamb >= threshold_lambda(x_eq,delta,alpha,beta,gamma) and delta < delta_all_coop:
			steady_fitness[lamb_list_list.index(lamb),delta_list_list.index(delta)] = fitness(lamb,delta,theta,alpha,beta,gamma)
			#print steady_fitness
		elif delta >= delta_all_coop:
			steady_fitness[lamb_list_list.index(lamb),delta_list_list.index(delta)] = G_delta(1.0,delta,gamma,alpha,beta)
		else:
			steady_fitness[lamb_list_list.index(lamb),delta_list_list.index(delta)] = G_delta(x_eq, delta, gamma, alpha, beta)
			
	threshold_holder.append(threshold_lambda(x_eq,delta,alpha,beta,gamma))
			

"""
Plotting heatmap for average payoff at steady state. 
"""

fig = plt.imshow(np.flipud(steady_fitness), "jet")
plt.colorbar(fig,fraction=0.04, pad=0.04)


plt.axvline(100.0*delta_all_coop, lw = 6.0, color = "gray", ls = "--")
x_labels = delta_list[::10]
x_labels = np.around(x_labels, decimals = 1)
x_label_loc = [100.0 * delta for delta in delta_list]
x_label_loc = x_label_loc[::10]
print x_labels
print x_label_loc

#y_labels = lamb_list
y_labels = np.flip(lamb_list[::4],0)
print y_labels
print lamb_list[::4]
y_labels_loc = np.arange(0.0,42.0,2.0)
y_lables = np.arange(40.0,-2.0,-2.0)
plt.xticks(x_label_loc,x_labels)
plt.yticks(2.0*y_labels_loc,y_lables)


plt.xlabel(r"Continuation Probability / Discount Rate ($\delta$)", fontsize = 20.,labelpad =20)
plt.ylabel(r"Relative Selection Strength ($\lambda$)", fontsize = 20.)



plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

threshold_holder = [80 - 2.5 * y if 80 - 2.5 * y <= 79.5 else 79.5 for y in threshold_holder]


plt.plot(threshold_holder, lw = 6., ls = '--', color = 'gray')


plt.tight_layout()


script_folder = os.getcwd()
mechanism_folder = os.path.dirname(script_folder)

plt.savefig(mechanism_folder + "/Figures/deltafitnessheatmapjet.png")

plt.show()

		

#plt.show()	