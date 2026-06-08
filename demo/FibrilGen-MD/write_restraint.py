
def make_a_list(lo_range,lo_value,idx):
	acc = []
	if not idx:
		for (r,v) in zip (lo_range,lo_value):
			choices = len(v)
			acc += [v[(i-r[0])%choices] for i in range(r[0],r[1])]
	else:
		for r in lo_range:
			acc += [i for i in range(r[0],r[1])]		
	return acc

# INPUT structure information
lo_range = [[i*10,(i+1)*10] for i in range(8)]				# The range of beta-stand index in every beta-sheet 
space_ = [[163],[165],[163],[165],[163],[165],[163],[165]]	# The number of atoms in a beta-stand for each b-sheet
offset_head_ = [[66],[66],[66],[66],[66],[66],[66],[66]]	# The first central atom in a beta-stand for each b-sheet
offset_tail_ = [[83],[84],[83],[84],[83],[84],[83],[84]]	# The second central atom in a beta-stand for each b-sheet

# INPUT flat-bottomed harmonic potential
b1,b2,b3,b4 = 0.0,2.0,9.0,11.0								# 0-2 Å for the lower bounds and 9-11 Å for the upper bounds of the NMR restraint
k_init = 10.												# The initial force constant
k_step0 = [k_init]
s1,s2,s3 = [0,2500,2501,5000],[0,50000,50001,100000]\
			,[0,20000,20001,500000]							# The number of steps for a constant force constant
k_step3 = [9.,8.,7.,6.,5.,4.,3.,2.,1.,0.]					# The decreasing force constants
k_step_split3 = [[20001,22000],[22001,24000],\
					[24001,26000],[26001,28000],\
					[28001,30000],[30001,32000],\
					[32001,34000],[34001,36000],\
					[36001,38000],[38001,500000]] 			# The number of steps for each decreasing force constant

# OUTPUT writing the restraints iteratively
lo_space = make_a_list(lo_range,space_,False)
lo_head = make_a_list(lo_range,offset_head_,False)
lo_tail = make_a_list(lo_range,offset_tail_,False)
lo_end = [r[1] for r in lo_range]
for (root_out,s,k_step) in zip(['dist1.RST','dist2.RST','dist3.RST'],[s1,s2,s3],[k_step0,k_step0,k_step3]):
	c_step0,c_step1,c_step2,c_step3 = s
	if len(k_step) == 1:
		k_step_split = [[c_step2,c_step3]] 
	else:
		k_step_split = k_step_split3

	with open(root_out,'w') as f:
		print (root_out)
		acc_atom_idx,acc_r_idx = 0,0
		for (offset_head,offset_tail,n_offset_head,n_offset_tail,space) in zip(lo_head[:-1],lo_tail[:-1],lo_head[1:],lo_tail[1:],lo_space[:-1]):			
			head_atom = acc_atom_idx+offset_head+1
			tail_atom = acc_atom_idx+offset_tail+1
			next_head_atom = acc_atom_idx + space + n_offset_head + 1
			next_tail_atom = acc_atom_idx + space + n_offset_tail + 1
			acc_atom_idx += space
			acc_r_idx += 1
			if acc_r_idx not in lo_end:
				f.write(f'&rst iat={str(head_atom)},{str(next_head_atom)}, r1={str(b1)}, r2={str(b2)}, r3={str(b3)}, r4={str(b4)}, rk2={str(k_init)}, rk3={str(k_init)}, nstep1={str(c_step0)}, nstep2={str(c_step1)}, &end/ \n')
				f.write(f'&rst iat={str(tail_atom)},{str(next_tail_atom)}, r1={str(b1)}, r2={str(b2)}, r3={str(b3)}, r4={str(b4)}, rk2={str(k_init)}, rk3={str(k_init)}, nstep1={str(c_step0)}, nstep2={str(c_step1)}, &end/ \n')
				for (k,ks) in zip(k_step,k_step_split):
					f.write(f'&rst iat={str(head_atom)},{str(next_head_atom)}, nstep1={str(ks[0])}, nstep2={str(ks[1])}, ifvari=1, r1a={str(b1)}, r2a={str(b2)}, r3a={str(b3)}, r4a={str(b4)}, rk2a={str(k)}, rk3a={str(k)}, &end/ \n')
					f.write(f'&rst iat={str(tail_atom)},{str(next_tail_atom)}, nstep1={str(ks[0])}, nstep2={str(ks[1])}, ifvari=1, r1a={str(b1)}, r2a={str(b2)}, r3a={str(b3)}, r4a={str(b4)}, rk2a={str(k)}, rk3a={str(k)}, &end/ \n')

