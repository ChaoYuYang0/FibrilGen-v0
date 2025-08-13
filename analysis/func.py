
import numpy as np

def flatten_one_order(x):
	acc = []
	for xi in x:
		acc += xi
	return acc

def cal_cos(v1,v2):
	v_dot = np.sum(v1*v2,axis=1)
	norm_v1 = np.linalg.norm(v1,axis=1)
	norm_v2 = np.linalg.norm(v2,axis=1)
	return v_dot/(norm_v1*norm_v2)

class fibril_conformation():
	def __init__(self,s):
		self.s = s

	def get_neighbor_k(self,lo_r,k):
		acc = []
		for r in lo_r:
			acc += [[i+j for j in range(k)] for i in range(r[0],r[1]+2-k)]
		return acc

	def fit_axis(self,coord_all,idx):
		# Get COM
		idx_all = flatten_one_order(idx)
		com = np.mean(coord_all[idx_all],axis=0).reshape((1,3))
		# Get t
		acc_t = []
		for idx_i in idx:
			idx_c = np.argmin(np.linalg.norm(coord_all[idx_i]-com,axis=1))
			t = np.array([i for i in range(len(idx_i))])-idx_c
			acc_t += t.tolist()
		# Get b
		X = np.stack([np.ones(len(acc_t)),acc_t],axis=0).T
		Y = coord_all[idx_all] 
		b = np.dot(np.linalg.inv(np.dot(X.T,X)),np.dot(X.T,Y))
		return b[1].reshape(1,3),b[0].reshape(1,3)

	def cal_radius(self,coord_all,p,axis):
		cos = cal_cos(coord_all-p,axis)
		sin = (1-cos**2)**0.5
		dist = np.linalg.norm(coord_all-p,axis=1)*sin
		return dist

	def cal_pitch(self,coord_all,n,axis):
		p1 = coord_all[n[:,0].tolist()]
		p2 = coord_all[n[:,1].tolist()]
		v12 = p2-p1
		p = cal_cos(v12,axis)*np.linalg.norm(v12,axis=1)
		return p

	def get_morph(self,traj):
		idx = [[i for i in range(si[0],si[1]+1)] for si in self.s]
		idx_all = flatten_one_order(idx)
		n = np.array(self.get_neighbor_k(self.s,2))
		lo_radius,lo_pitch = [],[]
		for i in range(traj.shape[0]):
			coord_all = traj[i]
			axis,p = self.fit_axis(coord_all,idx)
			print (axis,p)
			# Radius
			radius = self.cal_radius(coord_all[idx_all],p,axis)
			lo_radius += [np.mean(radius,axis=0)]
			# Pitch
			lo_pitch += [np.linalg.norm(axis)]
		return {'radius':np.array(lo_radius),'pitch':np.array(lo_pitch)}



