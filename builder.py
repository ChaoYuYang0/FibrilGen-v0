# run builder.py

import numpy as np
import pymol
import math

def get_ca(name):
	pymol.cmd.select('lo_ca','name ca and '+name)
	pos_ca = pymol.cmd.get_coords('lo_ca',1)
	return np.array(pos_ca)

def get_com(name):
	return np.mean(np.array(pymol.cmd.get_coords(name,1)),0)

def get_boundary(lo_name):
	lop = []
	for name in lo_name:
		lop += pymol.cmd.get_coords(name,1).tolist()
	lop = np.array(lop)
	return [max(lop[:,0]),min(lop[:,0]),max(lop[:,1]),min(lop[:,1]),max(lop[:,2]),min(lop[:,2])]

def rotate_coordinate(name,coord):
	mf = coord[0].tolist()+[0]+coord[1].tolist()+[0]+coord[2].tolist()+[0]+[0,0,0,1]
	pymol.cmd.transform_selection(name,mf)

def check_clash(lon,dist_tolorence):
	for name in lon:
		pymol.cmd.select('contacts',name+' around '+str(dist_tolorence))
		print (name,pymol.cmd.get_coords('contacts',1))


def example(morphology):
	## Reset
	pymol.cmd.delete('all')
	pymol.cmd.reset()
	## Load PDB
	pymol.cmd.load('examples/input/capF8.pdb')
	## Select a reference coordinate
	pymol.cmd.select('p1','resi 21-30')
	pymol.cmd.select('p2','resi 31-40')
	pymol.cmd.select('p3','resi 111-120')
	pymol.cmd.select('p4','resi 121-130')
	pymol.cmd.select('po1','resi 22 and name ca')
	pymol.cmd.select('po2','resi 29 and name ca')
	pymol.cmd.select('po3','resi 62 and name ca')

	## Create a periodic unit
	unit = create_sheet_unit('p1','p2','p3','p4','po1','po2','po3')

	## Create a fibril object
	fibril = create_fibril(unit)
	## Change default geo parameters
	# fibril.tilt_s1, fibril.tilt_s2 = 0, 0
	# fibril.dist_tolorence = 1.0

	## --- Examples of structures ---
	## Build a plain sheet (num of units per sheet)
	if morphology == 'a_sheet':
		fibril.build_a_flat_sheet(5)
	## Build a plain sheet (stacking pattern, num of units per sheet)
	elif morphology == 's_sheet':
		fibril.build_a_stacked_sheet([[0,1],[1,1]],5)								
	## Build a rod (tilt angle, num of units per sheet, twist sign)
	elif morphology == 'a_rod':
		fibril.build_a_rod(20,25,1) 	
	## Build a stacked rod (tilt angle, stacking pattern, num of units per sheet, twist sign)								
	elif morphology == 's_rod':
		fibril.build_a_stacked_rod(10,[[0,1],[1,1]],15,1)	
	## Build a ribbon (tilt angle, radius, num of units per sheet, twist sign)
	elif morphology == 'a_ribbon':
		fibril.build_a_ribbon(30,30,20,-1)
	## Build a stacked ribbon (tilt angle, radius, stacking angle, stacking number, num of units per sheet, twist sign)
	elif morphology == 's_ribbon':
		fibril.build_a_stacked_ribbon(10,30,70,2,5,1)
	else:
		print ('Unknown morphology!')
	## --- End ---
	
	fibril.get_dimension() # Print the refined fibril dimension
	pymol.cmd.zoom()


class create_sheet_unit():
	# INPUT (sheet 1 peptide 1, sheet 1 peptide 2, sheet 2 peptide 1, sheet 2 peptide 2, point for x head, point for x tail, point for y)
	def __init__(self,pep1_s1,pep2_s1,pep1_s2,pep2_s2,po1,po2,po3):
		# Create a unit
		pymol.cmd.create('s1_pep1',pep1_s1,0,0,1)
		pymol.cmd.create('s1_pep2',pep2_s1,0,0,1)
		pymol.cmd.create('s2_pep1',pep1_s2,0,0,1)
		pymol.cmd.create('s2_pep2',pep2_s2,0,0,1)
		# Create a cooridinate
		coord = self.get_coordinate_by_xy(po1,po2,po3)
		# Center the unit
		pymol.cmd.group('unit','s1_pep1 s1_pep2 s2_pep1 s2_pep2')
		com = np.mean(get_ca('unit'),0)
		pymol.cmd.translate((-com).tolist(),'s1_pep1')
		pymol.cmd.translate((-com).tolist(),'s1_pep2')
		pymol.cmd.translate((-com).tolist(),'s2_pep1')
		pymol.cmd.translate((-com).tolist(),'s2_pep2')
		# Align the unit to the coordinate
		rotate_coordinate('s1_pep1',coord)
		rotate_coordinate('s1_pep2',coord)
		rotate_coordinate('s2_pep1',coord)
		rotate_coordinate('s2_pep2',coord)
		# Set dimensions
		ca_s1 = np.mean(get_ca('s1_pep1').tolist()+get_ca('s1_pep2').tolist(),0)
		ca_s2 = np.mean(get_ca('s2_pep1').tolist()+get_ca('s2_pep2').tolist(),0)
		self.d = (ca_s2-ca_s1)[2]
		self.b1 = np.mean(get_ca('s1_pep2'),0)[1]-np.mean(get_ca('s1_pep1'),0)[1]
		self.b2 = np.mean(get_ca('s2_pep2'),0)[1]-np.mean(get_ca('s2_pep1'),0)[1]
		self.b = min(self.b1,self.b2)
		boundary = get_boundary(['s1_pep1','s1_pep2'])
		self.l = boundary[0]-boundary[1]
		box_boundaries = get_boundary(['s1_pep1','s1_pep2','s2_pep1','s2_pep2'])
		self.box_w = box_boundaries[4]-box_boundaries[5]
		self.box_l = box_boundaries[0]-box_boundaries[1]

	def get_coordinate_by_xy(self,po1,po2,po3):
		pos_po1 = np.array(pymol.cmd.get_coords(po1,1)).reshape(3)
		pos_po2 = np.array(pymol.cmd.get_coords(po2,1)).reshape(3)
		pos_po3 = np.array(pymol.cmd.get_coords(po3,1)).reshape(3)
		x = pos_po2 - pos_po1
		y = pos_po3 - pos_po1
		z = np.cross(x,y)
		x = np.cross(y,z) 
		unit_x,unit_y,unit_z = x/np.linalg.norm(x),y/np.linalg.norm(y),z/np.linalg.norm(z)
		return [unit_x,unit_y,unit_z]


class create_fibril():
	def __init__(self,unit): 
		self.unit = unit
		# Default tilt angle
		self.tilt_s1 = 0
		self.tilt_s2 = 0
		# Default tolorence
		self.dist_tolorence = 0.6

	def set_max_twist(self,angle_z,z_sign,y_sign):
		self.max_twist = 0
		for i in range(21):
			angle = 0 + i*0.4
			if self.twist_y_wo_intersection(angle_z,angle,z_sign,y_sign):
				self.max_twist = angle
			else:
				return
		return

	def twist_y_wo_intersection(self,angle_z,angle_y,z_sign,y_sign):
		def idx_collection(x):
			acc = ''
			for xi in x:
				acc += (str(xi)+'+')
			return acc[:-1]

		def get_sidechain_coord(name):
			pymol.cmd.select('sele_ca',name+' and name ca')
			pymol.cmd.select('sele_cb',name+' and name cb')
			pymol.cmd.select('sele_sc',name+' and sidechain and not (name ca+cb+ha)')
			idx_ca = [item[1] for item in pymol.cmd.index('sele_ca')]
			idx_cb = [item[1] for item in pymol.cmd.index('sele_cb')]
			idx_sc = [item[1] for item in pymol.cmd.index('sele_sc')]
			dic_coord,acc = {},0
			for (idx1,idx2) in zip(idx_ca,idx_ca[1:]+[idx_ca[-1]+50]):
				c_cb = [idx for idx in idx_cb if idx1<idx<idx2]
				c_sc = [idx for idx in idx_sc if idx1<idx<idx2]
				if (c_cb != []) and (c_sc != []):
					coord_ca = pymol.cmd.get_coords(name+' and index '+str(idx1),1)
					coord_cb = pymol.cmd.get_coords(name+' and index '+str(c_cb[0]),1)
					coord_sc = pymol.cmd.get_coords(name+' and index '+idx_collection(c_sc),1)
					dic_coord[acc] = [coord_ca,coord_cb,coord_sc,idx_sc]
				acc += 1
			return dic_coord


		def is_clash(name):
			pymol.cmd.select('contacts',name+' around '+str(self.dist_tolorence))
			if pymol.cmd.get_coords('contacts',1) is not None:
				return 1
			else:
				return 0

		# Testing by a far away unit
		i = 0
		z_sign1 = z_sign[0]
		name_pep1,name_pep2,name_pep3 = 'test_s1_pep1','test_s1_pep2','test_s1_pep3'
		pymol.cmd.create(name_pep1,'s1_pep1')
		self.affine_transformation(name_pep1,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*2*i,[0,0,0],[0,self.unit.b1*2*i+200,0])
		pymol.cmd.create(name_pep2,'s1_pep2')
		self.affine_transformation(name_pep2,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*(2*i+1),[0,-self.unit.b1,0],[0,self.unit.b1*(2*i+1)+200,0])
		pymol.cmd.create(name_pep3,'s1_pep1')
		self.affine_transformation(name_pep3,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*2*(i+1),[0,0,0],[0,self.unit.b1*2*(i+1)+200,0])

		z_sign2 = z_sign[1]
		name_pep1,name_pep2,name_pep3 = 'test_s2_pep1','test_s2_pep2','test_s2_pep3'
		pymol.cmd.create(name_pep1,'s2_pep1')
		self.affine_transformation(name_pep1,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*2*i,[0,0,0],[0,self.unit.b2*2*i+200,0])
		pymol.cmd.create(name_pep2,'s2_pep2')
		self.affine_transformation(name_pep2,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*(2*i+1),[0,-self.unit.b2,0],[0,self.unit.b2*(2*i+1)+200,0])
		pymol.cmd.create(name_pep3,'s2_pep1')
		self.affine_transformation(name_pep3,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*2*(i+1),[0,0,0],[0,self.unit.b2*2*(i+1)+200,0])

		# Check whether in close contacts
		for name in ['test_s1_pep1','test_s1_pep2','test_s1_pep3','test_s2_pep1','test_s2_pep2','test_s2_pep3']:
			if is_clash(name):
				pymol.cmd.delete('test*')
				return 0
				
		pymol.cmd.delete('test*')
		return 1

	def build_a_flat_sheet(self,num_half):
		for i in range(num_half):
			name_pep1,name_pep2 = 'p_s1_pep1_'+str(i),'p_s1_pep2_'+str(i)
			pymol.cmd.create(name_pep1,'s1_pep1')
			pymol.cmd.translate([0,self.unit.b1*2*i,0],name_pep1)
			pymol.cmd.create(name_pep2,'s1_pep2')
			pymol.cmd.translate([0,self.unit.b1*2*i,0],name_pep2)
		for i in range(num_half):
			name_pep1,name_pep2 = 'p_s2_pep1_'+str(i),'p_s2_pep2_'+str(i)
			pymol.cmd.create(name_pep1,'s2_pep1')
			pymol.cmd.translate([0,self.unit.b2*2*i,0],name_pep1)
			pymol.cmd.create(name_pep2,'s2_pep2')
			pymol.cmd.translate([0,self.unit.b2*2*i,0],name_pep2)
		pymol.cmd.color('green','p_s1_*')
		pymol.cmd.color('orange','p_s2_*')
		pymol.cmd.group('plain_sheet','p_*')
		self.set_dimension(0,0,0,self.unit.b)
		return

	def build_a_stacked_sheet(self,stacking,num_half):
		# Build a position matrix
		stacking = np.array(stacking)
		stack_z,stack_x = stacking.shape
		pos_matrix = np.zeros((stack_z,stack_x,2))
		started_z,started_x = (stack_z-1)/2.0,(stack_x-1)/2.0
		dist_z = self.unit.box_w
		dist_x = self.unit.box_l
		for j in range(stack_z):
			for k in range(stack_x):
				pos_matrix[j,k,0] = (started_z-j)*dist_z
				pos_matrix[j,k,1] = (-started_x+k)*dist_x
		idx_unit = 0
		for (pos,stack_this) in zip(pos_matrix.reshape((-1,2)).tolist(),stacking.flatten()):
			if stack_this:
				# Get param
				pos_z,pos_x = pos
				# Build the structure
				for i in range(num_half):
					name = str(i)+'_'+str(idx_unit)
					name_pep1,name_pep2 = 'sp_s1_pep1_'+name,'sp_s1_pep2_'+name
					pymol.cmd.create(name_pep1,'s1_pep1')
					pymol.cmd.translate([pos_x,self.unit.b1*2*i,pos_z],name_pep1)
					pymol.cmd.create(name_pep2,'s1_pep2')
					pymol.cmd.translate([pos_x,self.unit.b1*2*i,pos_z],name_pep2)
				for i in range(num_half):
					name = str(i)+'_'+str(idx_unit)
					name_pep1,name_pep2 = 'sp_s2_pep1_'+name,'sp_s2_pep2_'+name
					pymol.cmd.create(name_pep1,'s2_pep1')
					pymol.cmd.translate([pos_x,self.unit.b2*2*i,pos_z],name_pep1)
					pymol.cmd.create(name_pep2,'s2_pep2')
					pymol.cmd.translate([pos_x,self.unit.b2*2*i,pos_z],name_pep2)
				idx_unit += 1
		pymol.cmd.color('green','sp_s1_*')
		pymol.cmd.color('orange','sp_s2_*')
		pymol.cmd.group('plain_sheet','sp_*')
		self.set_dimension(0,0,0,self.unit.b)
		return


	def build_a_rod(self,angle_z,num_half,sign):
		radius = self.unit.d/2.0
		# Check geometry
		param = self.refine_theta(angle_z*np.pi/180,radius,sign)
		if param == None:
			print ('Please decrease tilt angle!')
			print ('Stop to update ... ')
			return
		else:
			theta_z,theta_y,y = param
			angle_z,angle_y = theta_z*180/np.pi,theta_y*180/np.pi
			z_sign1,y_sign = -sign,-sign
			# Build the structure
			for i in range(num_half):
				name_pep1,name_pep2 = 'nr_s1_pep1_'+str(i),'nr_s1_pep2_'+str(i)
				pymol.cmd.create(name_pep1,'s1_pep1')
				self.affine_transformation(name_pep1,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*2*i,[0,0,0],[0,y*2*i,0])
				pymol.cmd.create(name_pep2,'s1_pep2')
				self.affine_transformation(name_pep2,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*(2*i+1),[0,-self.unit.b1,0],[0,y*(2*i+1),0])
			z_sign2,y_sign = sign,-sign
			for i in range(num_half):
				name_pep1,name_pep2 = 'nr_s2_pep1_'+str(i),'nr_s2_pep2_'+str(i)
				pymol.cmd.create(name_pep1,'s2_pep1')
				self.affine_transformation(name_pep1,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*2*i,[0,0,0],[0,y*2*i,0])
				pymol.cmd.create(name_pep2,'s2_pep2')
				self.affine_transformation(name_pep2,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*(2*i+1),[0,-self.unit.b2,0],[0,y*(2*i+1),0])
			pymol.cmd.color('green','nr_s1_*')
			pymol.cmd.color('orange','nr_s2_*')
			pymol.cmd.group('a_rod','nr_*')
			self.set_dimension(radius,theta_z,theta_y,y)
			return

	def build_a_stacked_rod(self,angle_z,stacking,num_half,sign):
		# Build a position matrix
		stacking = np.array(stacking)
		stack_z,stack_x = stacking.shape
		pos_matrix = np.zeros((stack_z,stack_x,2))
		started_z,started_x = (stack_z-1)/2.0,(stack_x-1)/2.0
		theta_z = angle_z*np.pi/180
		dist_z = self.unit.box_w
		dist_x = self.unit.box_l*np.cos(theta_z)+self.unit.b*np.sin(theta_z)
		for j in range(stack_z):
			for k in range(stack_x):
				pos_matrix[j,k,0] = (started_z-j)*dist_z
				pos_matrix[j,k,1] = (-started_x+k)*dist_x
		# Check geometry
		radius_matrix = (pos_matrix[:,:,0]**2+pos_matrix[:,:,1]**2)**0.5
		param = self.refine_stack_rod(angle_z*np.pi/180,radius_matrix.reshape(-1),sign)
		if (param == None):
			print ('Please decrease tilt angle!')
			print ('Stop to update ... ')
			return
		else:
			# Stacking
			theta_z,theta_y,y = param
			angle_z,angle_y = theta_z*180/np.pi,theta_y*180/np.pi
			idx_unit = 0
			for (pos,stack_this) in zip(pos_matrix.reshape((-1,2)).tolist(),stacking.flatten()):
				if stack_this:
					# Get param
					pos_z,pos_x = pos
					y_sign = -sign
					z_sign1 = np.sign(pos_z-self.unit.d/2.0)*sign
					# Build the structure
					for i in range(num_half):
						name = str(i)+'_'+str(idx_unit)
						name_pep1,name_pep2 = 'snr_s1_pep1_'+name,'snr_s1_pep2_'+name
						pymol.cmd.create(name_pep1,'s1_pep1')
						pymol.cmd.translate([pos_x,0,pos_z],name_pep1)
						self.affine_transformation(name_pep1,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*2*i,[0,0,0],[0,y*2*i,0])
						pymol.cmd.create(name_pep2,'s1_pep2')
						pymol.cmd.translate([pos_x,0,pos_z],name_pep2)
						self.affine_transformation(name_pep2,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*(2*i+1),[0,-self.unit.b1,0],[0,y*(2*i+1),0])
					z_sign2 = np.sign(pos_z+self.unit.d/2.0)*sign
					# print (angle_y,z_sign1,z_sign2)
					for i in range(num_half):
						name = str(i)+'_'+str(idx_unit)
						name_pep1,name_pep2 = 'snr_s2_pep1_'+name,'snr_s2_pep2_'+name
						pymol.cmd.create(name_pep1,'s2_pep1')
						pymol.cmd.translate([pos_x,0,pos_z],name_pep1)
						self.affine_transformation(name_pep1,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*2*i,[0,0,0],[0,y*2*i,0])
						pymol.cmd.create(name_pep2,'s2_pep2')
						pymol.cmd.translate([pos_x,0,pos_z],name_pep2)
						self.affine_transformation(name_pep2,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*(2*i+1),[0,-self.unit.b2,0],[0,y*(2*i+1),0])
					idx_unit += 1
			pymol.cmd.color('green','snr_s1_*')
			pymol.cmd.color('orange','snr_s2_*')
			pymol.cmd.group('s_rod','snr_*')
			self.set_dimension(max(radius_matrix.reshape(-1)),theta_z,theta_y,y)
			return

	def build_a_ribbon(self,angle_z,radius,num_half,sign):
		radius_offset = self.unit.d/2.0
		# Check geometry
		param = self.refine_theta_radius(angle_z*np.pi/180,radius,sign)
		if param == None:
			print ('Please decrease tilt angle or try another radius!')
			print ('Stop to update ... ')
			return
		# Get param
		else:
			theta_z,theta_y,radius,y = param
			angle_z,angle_y = theta_z*180/np.pi,theta_y*180/np.pi
			y_sign = -sign
			# Build the structure
			z_sign1 = np.sign(radius-self.unit.d/2.0)*sign
			for i in range(num_half):
				name_pep1,name_pep2 = 'r_s1_pep1_'+str(i),'r_s1_pep2_'+str(i)
				pymol.cmd.create(name_pep1,'s1_pep1')
				pymol.cmd.translate([0,0,radius],name_pep1)
				self.affine_transformation(name_pep1,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*2*i,[0,0,0],[0,y*2*i,0])
				pymol.cmd.create(name_pep2,'s1_pep2')
				pymol.cmd.translate([0,0,radius],name_pep2)
				self.affine_transformation(name_pep2,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*(2*i+1),[0,-self.unit.b1,0],[0,y*(2*i+1),0])
			z_sign2 = np.sign(radius+self.unit.d/2.0)*sign
			for i in range(num_half):
				name_pep1,name_pep2 = 'r_s2_pep1_'+str(i),'r_s2_pep2_'+str(i)
				pymol.cmd.create(name_pep1,'s2_pep1')
				pymol.cmd.translate([0,0,radius],name_pep1)
				self.affine_transformation(name_pep1,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*2*i,[0,0,0],[0,y*2*i,0])
				pymol.cmd.create(name_pep2,'s2_pep2')
				pymol.cmd.translate([0,0,radius],name_pep2)
				self.affine_transformation(name_pep2,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*(2*i+1),[0,-self.unit.b2,0],[0,y*(2*i+1),0])
				
			pymol.cmd.color('green','r_s1_*')
			pymol.cmd.color('orange','r_s2_*')
			pymol.cmd.group('a_ribbon','r_*')
			self.set_dimension(radius,theta_z,theta_y,y)
			return

	def build_a_stacked_ribbon(self,angle_z,radius,angle_stack,num_stack,num_half,sign):
		radius_offset = self.unit.d/2.0
		# Check inner sheet geometry and edge contact
		param = self.refine_stack_ribbon(angle_z*np.pi/180,radius,sign,angle_stack*np.pi/180,num_stack)
		if param == None:
			print ('Please decrease tilt angle, try another radius, or try another stack angle!')
			print ('Stop to update ... ')
			return		
		# Get param
		else:
			theta_z,theta_y,radius,y = param
			angle_z,angle_y = theta_z*180/np.pi,theta_y*180/np.pi
			y_sign = -sign
			# Stacking
			for j in range(num_stack):
				# Build the structure
				z_sign1 = np.sign(radius-self.unit.d/2.0)*sign
				for i in range(num_half):
					name = str(i)+'_'+str(j)
					name_pep1,name_pep2 = 'sr_s1_pep1_'+name,'sr_s1_pep2_'+name
					pymol.cmd.create(name_pep1,'s1_pep1')
					pymol.cmd.translate([0,0,radius],name_pep1)
					self.affine_transformation(name_pep1,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*2*i-angle_stack*j,[0,0,0],[0,y*2*i,0])
					pymol.cmd.create(name_pep2,'s1_pep2')
					pymol.cmd.translate([0,0,radius],name_pep2)
					self.affine_transformation(name_pep2,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*(2*i+1)-angle_stack*j,[0,-self.unit.b1,0],[0,y*(2*i+1),0])
				z_sign2 = np.sign(radius+self.unit.d/2.0)*sign
				for i in range(num_half):
					name = str(i)+'_'+str(j)
					name_pep1,name_pep2 = 'sr_s2_pep1_'+name,'sr_s2_pep2_'+name
					pymol.cmd.create(name_pep1,'s2_pep1')
					pymol.cmd.translate([0,0,radius],name_pep1)
					self.affine_transformation(name_pep1,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*2*i-angle_stack*j,[0,0,0],[0,y*2*i,0])
					pymol.cmd.create(name_pep2,'s2_pep2')
					pymol.cmd.translate([0,0,radius],name_pep2)
					self.affine_transformation(name_pep2,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*(2*i+1)-angle_stack*j,[0,-self.unit.b2,0],[0,y*(2*i+1),0])
			pymol.cmd.color('green','sr_s1_*')
			pymol.cmd.color('orange','sr_s2_*')
			pymol.cmd.group('s_ribbon','sr_*')
			self.set_dimension(radius,theta_z,theta_y,y)
			return


	def set_dimension(self,radius,theta_z,theta_y,y):
		self.radius = radius
		self.angle_z = theta_z*180/np.pi
		self.angle_y = theta_y*180/np.pi
		if self.angle_y < 10**-3:
			self.period = float('inf')
			self.pitch = float('inf')
		else:
			self.period = 2*np.pi/theta_y
			self.pitch = y*self.period
		return

	def get_dimension(self):
		print ('Radius: '+str(self.radius/10.0)+' nm')
		if self.period != float('inf'):
			print ('Pitch length: '+str(self.pitch/10)+' nm')
		print ('Tilt angle: '+str(self.angle_z)+' degree')
		print ('Twist angle: '+str(self.angle_y)+' degree')
		print ('Period: '+str(self.period)+' peptides')
		return


	def refine_theta(self,theta_z,radius,sign):
		# Refine the structure by 40 iterations
		for i in range(40):
			self.set_max_twist(theta_z*180/np.pi,[np.sign(-self.unit.d/2.0)*sign,np.sign(self.unit.d/2.0)*sign],sign)
			# Estimate the twist angle
			k = self.unit.b*np.sin(theta_z)
			theta_y = np.arccos(1-0.5*(k/radius)**2)
			y = self.unit.b*np.cos(theta_z)
			small_enough_twist = (theta_y < self.max_twist*np.pi/180)
			if small_enough_twist:
				return [theta_z,theta_y,y]
			else:
				if theta_z > 0.02:
					theta_z -= 0.02

		self.set_dimension(radius,theta_z,theta_y,y)
		return None

	def refine_theta_radius(self,theta_z,radius,sign):
		# Refine the structure by 40 iterations
		for i in range(40):
			self.set_max_twist(theta_z*180/np.pi,[sign,sign],sign)
			# Estimate the twist angle
			k = self.unit.b*np.sin(theta_z)
			theta_y = np.arccos(1-0.5*(k/radius)**2)
			y = self.unit.b*np.cos(theta_z)
			small_enough_twist = (theta_y < self.max_twist*np.pi/180)
			if small_enough_twist:
				return [theta_z,theta_y,radius,y]
			else:
				radius += 1
				if theta_z > 0.02:
					theta_z -= 0.02
		self.set_dimension(radius,theta_z,theta_y,y)
		return None


	def refine_stack_rod(self,theta_z,lo_radius,sign):
		# Refine the structure by 40 iterations
		radius = min(lo_radius)	# check from the outermost sheet
		for i in range(40):
			self.set_max_twist(theta_z*180/np.pi,[np.sign(-self.unit.d/2.0)*sign,np.sign(self.unit.d/2.0)*sign],sign)
			k = self.unit.b*np.sin(theta_z)
			theta_y = np.arccos(1-0.5*(k/radius)**2)
			y = self.unit.b*np.cos(theta_z)
			small_enough_twist = (theta_y < self.max_twist*np.pi/180)
			if small_enough_twist:
				return [theta_z,theta_y,y]
			else:
				if theta_z > 0.02:
					theta_z -= 0.02
		self.set_dimension(radius,theta_z,theta_y,y)
		return None

	def refine_stack_ribbon(self,theta_z,radius,sign,theta_stack,num_stack):
		def edge_contact(radius,theta_z,theta_y,y,theta_stack,num_stack):
			edge = (radius**2+(self.unit.l*np.cos(theta_z)/2.0)**2)**0.5
			theta_shift = 2*np.arcsin(self.unit.l*np.cos(theta_z)/(2*edge))
			# Check intersection on the cross-section
			if (theta_stack < theta_shift) or (num_stack*theta_stack > 2*np.pi):
				return 'too_close'
			else:
				# Check whether peptides from neighbor sheets intersect
				y_shift = self.unit.l*np.sin(theta_z)/2.0
				pos_edge1 = np.array([np.cos(theta_stack)*edge,y_shift,np.sin(theta_stack)*edge])
				lo_contact = []
				for i in range(int(2*np.pi/theta_y)):
					pos_edge2 = np.array([np.cos(theta_shift+theta_y*i)*edge,-y_shift+y*i,np.sin(theta_shift+theta_y*i)*edge])
					contact = np.linalg.norm(pos_edge1-pos_edge2)
					if (1 < contact):
						lo_contact += [contact]
						continue
					else:
						return 'too_close'
				# Check whether existing peptides from neighbor sheets are in close contact
				if min(lo_contact) < 4:
					# print ('contact',lo_contact)
					return 'good_dist'
				else:
					return 'too_far'

		# Refine the structure by 40 iterations
		for i in range(40):
			self.set_max_twist(theta_z*180/np.pi,[sign,sign],sign)
			# print (self.max_twist)
			# Estimate the twist angle
			k = self.unit.b*np.sin(theta_z)
			theta_y = np.arccos(1-0.5*(k/radius)**2)
			y = self.unit.b*np.cos(theta_z)
			small_enough_twist = (theta_y < self.max_twist*np.pi/180)
			contact = edge_contact(radius-self.unit.d/2.0,theta_z,theta_y,y,theta_stack,num_stack)
			if (small_enough_twist and contact=='good_dist'):
				return [theta_z,theta_y,radius,y]
			elif (not small_enough_twist and contact=='too_close'):
				radius += 1
				if theta_z > 0.02:
					theta_z -= 0.02
			elif (small_enough_twist and contact=='too_far'):
				radius -= 1
				theta_z += 0.02
			else:
				radius += (np.random.randint(3)-1)
		self.set_dimension(radius,theta_z,theta_y,y)
		return None



	def affine_transformation(self,name,angle_z,angle_y,translation1,translation2):
		# pymol.cmd.translate(translation1,name)
		y = np.mean(get_ca(name),0)[1]
		pymol.cmd.translate([0,-y,0],name)
		pymol.cmd.rotate('z',angle_z,name,-1,1,None,[0,0,0])
		pymol.cmd.rotate('y',angle_y,name,-1,1,None,[0,0,0])
		pymol.cmd.translate(translation2,name)









