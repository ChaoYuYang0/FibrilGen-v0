import numpy as np
from func import *

def generate_a_fibril_sample(r,p):
	coord,i = [],0
	for angle in range(0,360,10):
		theta = np.pi*angle/180
		c,s = np.cos(theta),np.sin(theta)
		coord += [[r*c,r*s,i*p]]
		i += 1
	return np.array(coord)

def generate_fibrils(lor,lop):
	acc = []
	for (r,p) in zip(lor,lop):
		acc += [generate_a_fibril_sample(r,p)]
	return np.stack(acc,axis=0)

fibril_traj = generate_fibrils([20,18,16],[3,4,5])
fibril = fibril_conformation([[0,35]])
morph = fibril.get_morph(fibril_traj)
print ('Radii: %s and Pitches: %s'%(morph['radius'],morph['pitch']))