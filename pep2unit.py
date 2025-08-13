# run pep2unit.py

import numpy as np
import pymol
import math

def get_ca(name):
	pymol.cmd.select('lo_ca','name ca and '+name)
	pos_ca = pymol.cmd.get_coords('lo_ca',1)
	return np.array(pos_ca)

def get_boundary(lo_name):
	lop = []
	for name in lo_name:
		lop += pymol.cmd.get_coords(name,1).tolist()
	lop = np.array(lop)
	return [max(lop[:,0]),min(lop[:,0]),max(lop[:,1]),min(lop[:,1]),max(lop[:,2]),min(lop[:,2])]

def rotate_coordinate(name,coord):
	mf = coord[0].tolist()+[0]+coord[1].tolist()+[0]+coord[2].tolist()+[0]+[0,0,0,1]
	pymol.cmd.transform_selection(name,mf)


def example():
	## Reset
	pymol.cmd.delete('all')
	pymol.cmd.reset()
	## Load PDB
	pymol.cmd.load('examples/input/AL1.pdb')
	## Select a reference coordinate
	pymol.cmd.select('po1','resi 10 and name ca')
	pymol.cmd.select('po2','resi 2 and name ca')
	pymol.cmd.select('po3','resi 7 and name O')
	pymol.cmd.select('po4','resi 8 and name N')

	## Create a periodic unit
	unit = create_pep_unit('AL1','po1','po2','po3','po4')
	# unit.rotate_sidechain() ## list of angles

	## --- Examples of sheet structures ---
	sheet = create_sheet(unit,[0,11],[0,11])
	# INPUT: ([aaa/apa/aap/app/paa/ppa/pap/ppp,sidechain flip],num of units per sheet)
	sheet.build_a_plain_sheet(['pap','d'],5)
	## --- End ---
	
	sheet.get_dimension()
	pymol.cmd.zoom()


class create_pep_unit():
	# INPUT (point for the head of x, point for the tail of x, point for the head of y, point for the tail of y)
	def __init__(self,name,po1,po2,po3,po4):
		# Create a cooridinate
		coord = self.get_coordinate_by_xy(po1,po2,po3,po4)
		# Align the unit to the coordinate
		rotate_coordinate(name,coord)
		self.name = name
		self.set_boundary()

	def get_coordinate_by_xy(self,po1,po2,po3,po4):
		pos_po1 = np.array(pymol.cmd.get_coords(po1,1)).reshape(3)
		pos_po2 = np.array(pymol.cmd.get_coords(po2,1)).reshape(3)
		pos_po3 = np.array(pymol.cmd.get_coords(po3,1)).reshape(3)
		pos_po4 = np.array(pymol.cmd.get_coords(po4,1)).reshape(3)
		x = pos_po2 - pos_po1
		y = pos_po4 - pos_po3
		z = np.cross(x,y)
		y = np.cross(z,x)
		unit_x,unit_y,unit_z = x/np.linalg.norm(x),y/np.linalg.norm(y),z/np.linalg.norm(z)
		return [unit_x,unit_y,unit_z]

	def rotate_sidechain(self,ax,rot,idx):
		def rotate_atom(ax,idx,theta,ori):
			init_pos = np.array(pymol.cmd.get_coords('index '+idx,1))-ori_i
			c,s = np.cos(theta),np.sin(theta)
			if ax == 'x':
				m = np.array([[1,0,0],[0,c,-s],[0,s,c]])
			elif ax == 'y':
				m = np.array([[c,0,s],[0,1,0],[-s,0,c]])
			else:
				m = np.array([[c,-s,0],[s,c,0],[0,0,1]])
			final_pos = np.dot(m,init_pos.T).T
			pymol.cmd.translate((final_pos-init_pos)[0].tolist(),'index '+idx)

		# rotate around axis ax
		for i in idx:
			resi_i,rot_i = str(i+1),rot[i]
			# get origin
			ori_i = get_ca('resi '+resi_i)	
			# get side chain idx
			lo_idx = pymol.cmd.index('resi '+resi_i+' and not name CA+N+H+C+O+H1+H2+H3+OXT')
			# rotate side chain
			for idx in lo_idx:
				rotate_atom(ax,str(idx[1]),rot_i*np.pi/180,ori_i)
			self.set_boundary()

	def set_boundary(self):
		# Set dimensions
		boundary = get_boundary([self.name])
		self.y = 4.8 # Set peptide to peptide dist as 0.48 nm 
		self.z = boundary[4]-boundary[5]


class create_sheet():
	def __init__(self,unit,resi,resi_a): 
		self.unit = unit
		self.resi_start = resi[0]
		self.resi_end = resi[1]
		self.resi_a_start = resi_a[0]
		self.resi_a_end = resi_a[1]

	def build_a_plain_sheet(self,alignment,num_half):
		b_alignment,s_alignment = alignment
		b_flip_angle = self.get_b_flip_angle(b_alignment)
		s_flip_angle = self.get_s_flip_angle(s_alignment)
		z_offset = self.unit.z/2.0
		for i in range(num_half):
			name_pep1,name_pep2 = 'p_s1_pep1_'+str(i),'p_s1_pep2_'+str(i)
			pymol.cmd.create(name_pep1,self.unit.name)
			self.affine_transformation_backbone(name_pep1,[0,0])
			pymol.cmd.translate([0,self.unit.y*2*i,-z_offset],name_pep1)
			pymol.cmd.create(name_pep2,self.unit.name)
			self.affine_transformation_backbone(name_pep2,b_flip_angle[0])
			pymol.cmd.translate([0,self.unit.y*(2*i+1),-z_offset],name_pep2)
		for i in range(num_half):
			name_pep1,name_pep2 = 'p_s2_pep1_'+str(i),'p_s2_pep2_'+str(i)
			pymol.cmd.create(name_pep1,self.unit.name)
			self.affine_transformation_sidechain(name_pep1,s_flip_angle)
			self.affine_transformation_backbone(name_pep1,b_flip_angle[1])
			pymol.cmd.translate([0,self.unit.y*2*i,z_offset],name_pep1)
			pymol.cmd.create(name_pep2,self.unit.name)
			self.affine_transformation_sidechain(name_pep2,s_flip_angle)
			self.affine_transformation_backbone(name_pep2,b_flip_angle[2])
			pymol.cmd.translate([0,self.unit.y*(2*i+1),z_offset],name_pep2)
		pymol.cmd.color('green','p_s1_*')
		pymol.cmd.color('orange','p_s2_*')
		pymol.cmd.group('plain_sheet','p_*')
		lon = ['p_s1_pep1_'+str(i) for i in range(num_half)]+['p_s1_pep2_'+str(i) for i in range(num_half)]+\
				['p_s2_pep1_'+str(i) for i in range(num_half)]+['p_s2_pep2_'+str(i) for i in range(num_half)]
		self.set_dimension(lon)
		return

	def set_dimension(self,lon):
		boundary = get_boundary(lon)
		self.x = boundary[0]-boundary[1]
		self.y = boundary[2]-boundary[3]
		self.z = boundary[4]-boundary[5]
		return

	def get_dimension(self):
		print ('Dimension x: '+str(self.x)+' nm')
		print ('Dimension y: '+str(self.y)+' nm')
		print ('Dimension z: '+str(self.z)+' nm')
		return

	def get_b_flip_angle(self,b_alignment):
		# assign flip angle for [s1_pep2,s2_pep1,s2_pep2]
		if b_alignment=='aaa':
			return [[0,180],[180,0],[180,180]]
		elif b_alignment=='apa':
			return [[0,180],[180,180],[180,0]]
		elif b_alignment=='aap':
			return [[0,180],[180,0],[180,0]]
		elif b_alignment=='app':
			return [[0,180],[180,180],[180,180]]
		elif b_alignment=='paa':
			return [[0,0],[180,0],[180,180]]
		elif b_alignment=='ppa':
			return [[0,0],[180,180],[180,0]]
		elif b_alignment=='pap':
			return [[0,0],[180,0],[180,0]]
		elif b_alignment=='ppp':
			return [[0,0],[180,180],[180,180]]
		else:
			return None

	def get_s_flip_angle(self,s_alignment):
		if s_alignment == 's':
			return [0,0]
		elif s_alignment == 'd':
			return [180,180]
		else:
			return None

	def affine_transformation_sidechain(self,name,angle): # rotate z-axis than x-axis
		angle_y,angle_z = angle
		pymol.cmd.rotate('z',angle_z,name)
		pymol.cmd.rotate('y',angle_y,name)
		# Correct pos
		pos_com = np.mean(get_ca(name)[self.resi_start:self.resi_end],0)
		pymol.cmd.translate([0,-pos_com[1],0],name)	


	def affine_transformation_backbone(self,name,angle): # rotate y-axis than z-axis
		angle_y,angle_z = angle
		pymol.cmd.rotate('y',angle_y,name)
		pymol.cmd.rotate('z',angle_z,name)
		# Correct pos
		pos_com = np.mean(get_ca(name)[self.resi_start:self.resi_end],0)
		pymol.cmd.translate([0,-pos_com[1],0],name)
		if  (angle_y+angle_z)%360 != 0:
			pos_x_init = get_ca(name)[self.resi_a_end][0]
		else:
			pos_x_init = get_ca(name)[self.resi_a_start][0]
		pymol.cmd.translate([-pos_x_init,0,0],name)			

	def show_unit(self):
		pymol.cmd.color('cyan','p_s1_pep2_*')
		pymol.cmd.color('purple','p_s2_pep2_*')
		pymol.cmd.hide('all')
		pymol.cmd.show('sticks','backbone or resn ace+nhe')
		return





