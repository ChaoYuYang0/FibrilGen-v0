# run builder.py

import numpy as np
import pymol
import math
from pymol.cgo import *

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

def get_bounding_vertices(box_boundaries):
	minX,minY,minZ,maxX,maxY,maxZ = box_boundaries[1],box_boundaries[3],box_boundaries[5],box_boundaries[0],box_boundaries[2],box_boundaries[4]
	loPoints = [[minX, minY, minZ],[minX, minY, maxZ],\
				[minX, maxY, minZ], [minX, maxY, maxZ],\
				[maxX, minY, minZ], [maxX, minY, maxZ],\
				[maxX, maxY, minZ], [maxX, maxY, maxZ]]
	return loPoints

def draw_box(loPoints,boxName):
	linewidth, r, g, b = 2.0, 1.0, 0.0, 0.0
	boundingBox = [
		LINEWIDTH, float(linewidth),

		BEGIN, LINES,
		COLOR, float(r), float(g), float(b),

		VERTEX, loPoints[0][0], loPoints[0][1], loPoints[0][2],       #1
		VERTEX, loPoints[1][0], loPoints[1][1], loPoints[1][2],       #2

		VERTEX, loPoints[2][0], loPoints[2][1], loPoints[2][2],       #3
		VERTEX, loPoints[3][0], loPoints[3][1], loPoints[3][2],       #4

		VERTEX, loPoints[4][0], loPoints[4][1], loPoints[4][2],       #5
		VERTEX, loPoints[5][0], loPoints[5][1], loPoints[5][2],       #6

		VERTEX, loPoints[6][0], loPoints[6][1], loPoints[6][2],       #7
		VERTEX, loPoints[7][0], loPoints[7][1], loPoints[7][2],       #8


		VERTEX, loPoints[0][0], loPoints[0][1], loPoints[0][2],       #1
		VERTEX, loPoints[4][0], loPoints[4][1], loPoints[4][2],       #5

		VERTEX, loPoints[2][0], loPoints[2][1], loPoints[2][2],       #3
		VERTEX, loPoints[6][0], loPoints[6][1], loPoints[6][2],       #7

		VERTEX, loPoints[3][0], loPoints[3][1], loPoints[3][2],       #4
		VERTEX, loPoints[7][0], loPoints[7][1], loPoints[7][2],       #8

		VERTEX, loPoints[1][0], loPoints[1][1], loPoints[1][2],       #2
		VERTEX, loPoints[5][0], loPoints[5][1], loPoints[5][2],       #6


		VERTEX, loPoints[0][0], loPoints[0][1], loPoints[0][2],       #1
		VERTEX, loPoints[2][0], loPoints[2][1], loPoints[2][2],       #3

		VERTEX, loPoints[4][0], loPoints[4][1], loPoints[4][2],       #5
		VERTEX, loPoints[6][0], loPoints[6][1], loPoints[6][2],       #7

		VERTEX, loPoints[1][0], loPoints[1][1], loPoints[1][2],       #2
		VERTEX, loPoints[3][0], loPoints[3][1], loPoints[3][2],       #4

		VERTEX, loPoints[5][0], loPoints[5][1], loPoints[5][2],       #6
		VERTEX, loPoints[7][0], loPoints[7][1], loPoints[7][2],       #8

		END
		]
	cmd.load_cgo(boundingBox,boxName)
	return

def example(morphology):
	## Reset
	pymol.cmd.delete('all')
	pymol.cmd.reset()
	## Load PDB
	pymol.cmd.load('structures/input/capF8_bilayer.pdb')
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
		fibril.build_a_ribbon(10,30,60,-1)
	## Build a stacked ribbon (tilt angle, radius, stacking angle, stacking number, num of units per sheet, twist sign)
	elif morphology == 's_ribbon':
		fibril.build_a_stacked_ribbon(10,30,70,2,20,1)
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
		self.box_boundaries = get_boundary(['s1_pep1','s1_pep2','s2_pep1','s2_pep2'])
		self.box_w = self.box_boundaries[4]-self.box_boundaries[5]
		self.box_l = self.box_boundaries[0]-self.box_boundaries[1]
		# Draw bounding box
		draw_box(get_bounding_vertices(self.box_boundaries),'UnitBox')


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

	def check_unit(self,angle_z,angle_y,z_sign,y_sign,radius):
		max_twist = 0
		for i in range(21):
			angle = 0 + i*0.4
			if self.unit_is_not_clashed(angle_z,angle,z_sign,y_sign,radius):
				max_twist = angle
		if max_twist > angle_y:
			return 1
		else:
			return 0

	def unit_is_not_clashed(self,angle_z,angle_y,z_sign,y_sign,radius):
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
		pymol.cmd.translate([0,0,radius],name_pep1)
		self.affine_transformation(name_pep1,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*2*i,[0,0,0],[0,self.unit.b1*2*i+200,0])
		pymol.cmd.create(name_pep2,'s1_pep2')
		pymol.cmd.translate([0,0,radius],name_pep2)
		self.affine_transformation(name_pep2,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*(2*i+1),[0,-self.unit.b1,0],[0,self.unit.b1*(2*i+1)+200,0])
		pymol.cmd.create(name_pep3,'s1_pep1')
		pymol.cmd.translate([0,0,radius],name_pep3)
		self.affine_transformation(name_pep3,z_sign1*(angle_z+self.tilt_s1),y_sign*angle_y*2*(i+1),[0,0,0],[0,self.unit.b1*2*(i+1)+200,0])

		z_sign2 = z_sign[1]
		name_pep1,name_pep2,name_pep3 = 'test_s2_pep1','test_s2_pep2','test_s2_pep3'
		pymol.cmd.create(name_pep1,'s2_pep1')
		pymol.cmd.translate([0,0,radius],name_pep1)
		self.affine_transformation(name_pep1,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*2*i,[0,0,0],[0,self.unit.b2*2*i+200,0])
		pymol.cmd.create(name_pep2,'s2_pep2')
		pymol.cmd.translate([0,0,radius],name_pep2)
		self.affine_transformation(name_pep2,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*(2*i+1),[0,-self.unit.b2,0],[0,self.unit.b2*(2*i+1)+200,0])
		pymol.cmd.create(name_pep3,'s2_pep1')
		pymol.cmd.translate([0,0,radius],name_pep3)
		self.affine_transformation(name_pep3,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*2*(i+1),[0,0,0],[0,self.unit.b2*2*(i+1)+200,0])

		# Check whether in close contacts
		for name in ['test_s1_pep1','test_s1_pep2','test_s1_pep3','test_s2_pep1','test_s2_pep2','test_s2_pep3']:
			if is_clash(name):
				pymol.cmd.delete('test*')
				return 0
				
		pymol.cmd.delete('test*')
		return 1		


	def build_a_flat_sheet(self,num_half):
		# Build the first sheet of the bilayer
		for i in range(num_half):
			name_pep1,name_pep2 = 'p_s1_pep1_'+str(i),'p_s1_pep2_'+str(i)
			pymol.cmd.create(name_pep1,'s1_pep1')
			pymol.cmd.translate([0,self.unit.b1*2*i,0],name_pep1)
			pymol.cmd.create(name_pep2,'s1_pep2')
			pymol.cmd.translate([0,self.unit.b1*2*i,0],name_pep2)
		# Build the second sheet of the bilayer
		for i in range(num_half):
			name_pep1,name_pep2 = 'p_s2_pep1_'+str(i),'p_s2_pep2_'+str(i)
			pymol.cmd.create(name_pep1,'s2_pep1')
			pymol.cmd.translate([0,self.unit.b2*2*i,0],name_pep1)
			pymol.cmd.create(name_pep2,'s2_pep2')
			pymol.cmd.translate([0,self.unit.b2*2*i,0],name_pep2)
		# Build box representation
		for i in range(num_half):
			name_box = 'vis_'+str(i)
			self.affine_transformation_a_box(name_box,0,0,[0,0,0],[0,self.unit.b1*2*i,0])
		pymol.cmd.color('green','p_s1_*')
		pymol.cmd.color('orange','p_s2_*')
		pymol.cmd.group('plain_sheet','p_*')
		pymol.cmd.group('visBox','vis_*')
		self.set_dimension(0,0,0,self.unit.b)
		return

	def build_a_stacked_sheet(self,stacking,num_half):
		# Get a position matrix
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
				for i in range(num_half):
					name_box = 'vis_'+str(i)+'_'+str(idx_unit)
					self.affine_transformation_a_box(name_box,0,0,[0,0,0],[pos_x,self.unit.b1*2*i,pos_z])
				idx_unit += 1
		pymol.cmd.color('green','sp_s1_*')
		pymol.cmd.color('orange','sp_s2_*')
		pymol.cmd.group('plain_sheet','sp_*')
		pymol.cmd.group('visBox','vis_*')
		self.set_dimension(0,0,0,self.unit.b)
		return


	def build_a_rod(self,angle_z,num_half,sign):
		radius = self.unit.d/2.0
		# Refine the input geometry
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
			for i in range(num_half):
				name_box = 'vis_'+str(i)
				self.affine_transformation_a_box(name_box,0,y_sign*theta_y*(2*i+0.5),[0,0,0],[0,y*2*i,0])

			pymol.cmd.color('green','nr_s1_*')
			pymol.cmd.color('orange','nr_s2_*')
			pymol.cmd.group('a_rod','nr_*')
			pymol.cmd.group('visBox','vis_*')
			self.set_dimension(radius,theta_z,theta_y,y)
			return

	def build_a_stacked_rod(self,angle_z,stacking,num_half,sign):
		# Get a position matrix
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
		# Refine the input geometry
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
					for i in range(num_half):
						name = str(i)+'_'+str(idx_unit)
						name_pep1,name_pep2 = 'snr_s2_pep1_'+name,'snr_s2_pep2_'+name
						pymol.cmd.create(name_pep1,'s2_pep1')
						pymol.cmd.translate([pos_x,0,pos_z],name_pep1)
						self.affine_transformation(name_pep1,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*2*i,[0,0,0],[0,y*2*i,0])
						pymol.cmd.create(name_pep2,'s2_pep2')
						pymol.cmd.translate([pos_x,0,pos_z],name_pep2)
						self.affine_transformation(name_pep2,z_sign2*(angle_z+self.tilt_s2),y_sign*angle_y*(2*i+1),[0,-self.unit.b2,0],[0,y*(2*i+1),0])
					for i in range(num_half):
						name_box = 'vis_'+str(i)+'_'+str(idx_unit)
						self.affine_transformation_a_box(name_box,0,y_sign*theta_y*(2*i+0.5),[pos_x,0,pos_z],[0,y*2*i,0])
					idx_unit += 1
			pymol.cmd.color('green','snr_s1_*')
			pymol.cmd.color('orange','snr_s2_*')
			pymol.cmd.group('s_rod','snr_*')
			pymol.cmd.group('visBox','vis_*')
			self.set_dimension(max(radius_matrix.reshape(-1)),theta_z,theta_y,y)
			return

	def build_a_ribbon(self,angle_z,radius,num_half,sign):
		radius_offset = self.unit.d/2.0
		# Refine the input geometry
		param = self.refine_theta_radius(angle_z*np.pi/180,radius,sign)
		if param == None:
			print ('Please decrease tilt angle or try another radius!')
			print ('Stop to update ... ')
			return
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
			tilt_s1,tilt_s2 = self.tilt_s1*np.pi/180,self.tilt_s2*np.pi/180
			for i in range(num_half):
				name_box = 'vis_'+str(i)
				self.affine_transformation_a_box(name_box,(z_sign1/2.0+z_sign2/2.0)*(theta_z+tilt_s1/2.0+tilt_s2/2.0),y_sign*theta_y*(2*i+0.5),[0,0,radius],[0,y*2*i,0])

			pymol.cmd.color('green','r_s1_*')
			pymol.cmd.color('orange','r_s2_*')
			pymol.cmd.group('a_ribbon','r_*')
			pymol.cmd.group('visBox','vis_*')
			self.set_dimension(radius,theta_z,theta_y,y)
			return

	def build_a_stacked_ribbon(self,angle_z,radius,angle_stack,num_stack,num_half,sign):
		radius_offset = self.unit.d/2.0
		# Refine the input geometry
		param = self.refine_stack_ribbon(angle_z*np.pi/180,radius,sign,angle_stack*np.pi/180,num_stack)
		if param == None:
			print ('Please decrease tilt angle, try another radius, or try another stack angle!')
			print ('Stop to update ... ')
			return		
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
				tilt_s1,tilt_s2,theta_stack = self.tilt_s1*np.pi/180,self.tilt_s2*np.pi/180,angle_stack*np.pi/180
				for i in range(num_half):
					name_box = 'vis_'+str(i)+'_'+str(j)
					self.affine_transformation_a_box(name_box,(z_sign1/2.0+z_sign2/2.0)*(theta_z+tilt_s1/2.0+tilt_s2/2.0),y_sign*theta_y*(2*i+0.5)-theta_stack*j,[0,0,radius],[0,y*2*i,0])
			pymol.cmd.color('green','sr_s1_*')
			pymol.cmd.color('orange','sr_s2_*')
			pymol.cmd.group('s_ribbon','sr_*')
			pymol.cmd.group('visBox','vis_*')
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
			# Calculate the twist angle
			k = self.unit.b*np.sin(theta_z)
			theta_y = np.arccos(1-0.5*(k/radius)**2)
			y = self.unit.b*np.cos(theta_z)
			# Check clashes
			good_unit = self.check_unit(theta_z*180/np.pi,theta_y*180/np.pi,[np.sign(-self.unit.d/2.0)*sign,np.sign(self.unit.d/2.0)*sign],sign,0)
			if good_unit:
				return [theta_z,theta_y,y]
			else:
				if theta_z > 0.02:
					theta_z -= 0.02

		self.set_dimension(radius,theta_z,theta_y,y)
		return None

	def refine_theta_radius(self,theta_z,radius,sign):
		# Refine the structure by 40 iterations
		for i in range(40):
			# Calculate the twist angle
			k = self.unit.b*np.sin(theta_z)
			theta_y = np.arccos(1-0.5*(k/radius)**2)
			y = self.unit.b*np.cos(theta_z)
			# Check clashes
			good_unit = self.check_unit(theta_z*180/np.pi,theta_y*180/np.pi,[sign,sign],sign,radius)
			if good_unit:
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
			k = self.unit.b*np.sin(theta_z)
			theta_y = np.arccos(1-0.5*(k/radius)**2)
			y = self.unit.b*np.cos(theta_z)
			# Check clashes
			good_unit = 1
			for a_radius in lo_radius:
				good_unit *= self.check_unit(theta_z*180/np.pi,theta_y*180/np.pi,[np.sign(-self.unit.d/2.0)*sign,np.sign(self.unit.d/2.0)*sign],sign,a_radius)
			if good_unit:
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
			# Calculate the twist angle
			k = self.unit.b*np.sin(theta_z)
			theta_y = np.arccos(1-0.5*(k/radius)**2)
			y = self.unit.b*np.cos(theta_z)
			# Check clashes
			good_unit = self.check_unit(theta_z*180/np.pi,theta_y*180/np.pi,[sign,sign],sign,radius)
			contact = edge_contact(radius-self.unit.d/2.0,theta_z,theta_y,y,theta_stack,num_stack)
			if (good_unit and contact=='good_dist'):
				return [theta_z,theta_y,radius,y]
			elif (not good_unit and contact=='too_close'):
				radius += 1
				if theta_z > 0.02:
					theta_z -= 0.02
			elif (good_unit and contact=='too_far'):
				radius -= 1
				theta_z += 0.02
			else:
				radius += (np.random.randint(3)-1)
		self.set_dimension(radius,theta_z,theta_y,y)
		return None

	def affine_transformation_a_box(self,name,theta_z,theta_y,translation1,translation2):
		vertices = np.array(get_bounding_vertices(self.unit.box_boundaries))
		coord = (vertices+np.array(translation1).reshape((-1,3))).T
		c,s = np.cos(theta_z),np.sin(theta_z)
		matrix_z = np.array([[c,-s,0],[s,c,0],[0,0,1]])
		coord = np.dot(matrix_z,coord)
		c,s = np.cos(theta_y),np.sin(theta_y)
		matrix_y = np.array([[c,0,s],[0,1,0],[-s,0,c]])
		coord = np.dot(matrix_y,coord)
		coord = coord.T+np.array(translation2).reshape((-1,3))
		draw_box(coord.tolist(),name)

	def affine_transformation(self,name,angle_z,angle_y,translation1,translation2):
		# pymol.cmd.translate(translation1,name)
		y = np.mean(get_ca(name),0)[1]
		pymol.cmd.translate([0,-y,0],name)
		pymol.cmd.rotate('z',angle_z,name,-1,1,None,[0,0,0])
		pymol.cmd.rotate('y',angle_y,name,-1,1,None,[0,0,0])
		pymol.cmd.translate(translation2,name)









