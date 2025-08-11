# FibrilGen
FibrilGen is a Python library that builds various cross-beta structures with a set of controllable parameters. Such cross-beta structures can be purely hypothetical or modelled based on experimental data (e.g., cryo-EM or ssNMR data). Combining FibrilGen with MD simulations constitutes a modelling pipeline for screening geometrically feasible and energetically favorable cross-beta structures, which enables a systematic assessment of peptide packing and peptide-based nanomaterials.

## PyMol versions
v2.3.5 (commercial), v1.7.4.5 (educational), v2.3.0 (open-source)

## Make example structures 
### 1. Import FibrilGen functions
```bash
# Change the current directory to FibrilGen
cd [directory of FibrilGen]
# Import FibrilGen functions 
run builder.py
```
### 2. Build example cross-beta structures
```bash
# Demonstrate building a beta-sheet structure
example('a_sheet')
# Demonstrate building a stacked beta-sheet structure
example('s_sheet')
# Demonstrate building a rod structure
example('a_rod')
# Demonstrate building a stacked rod structure
example('s_rod')
# Demonstrate building a ribbon structure
example('a_ribbon')
#  Demonstrate building a stacked ribbon (/tube) structure
example('s_ribbon')
```

## Make your structures
### 1. Modify bilayer_init.py to initialize a bilayer beta-sheet structure
```bash
# Import pep2unit library
pymol.cmd.run('pep2unit.py')
# Load a beta-strand
pymol.cmd.load('examples/input/AL1.pdb')
# Select a reference coordinate system from the beta-strand
pymol.cmd.select('po1','resi 10 and name ca')
pymol.cmd.select('po2','resi 2 and name ca')
pymol.cmd.select('po3','resi 7 and name O')
pymol.cmd.select('po4','resi 8 and name N')
# Assign a reference coordinate system
unit = create_pep_unit('AL1','po1','po2','po3','po4')
# Create a sheet object
sheet = create_sheet(unit,[0,11],[0,11])
# Build a bilayer structure 
sheet.build_a_plain_sheet(['pap','d'],10) # Parameters of ([backbone alignment (e.g., aaa,apa,aap,app,paa,ppa,pap,ppp), beta-sheets arranged face-to-face or face-to-back] and the number of peptides in each sheet
```
### 2. Execute bilayer_init.py in the PyMOL command line
```bash
run bilayer_init.py
```
### 3. Modify fibril_init.py to initialize a cross-beta nanostructure
```bash
# Import FibrilGen library
pymol.cmd.run('builder.py')
# Load a bilayer structure
pymol.cmd.load('examples/input/capF8.pdb')
# Select four peptides p1,p2,p3,p4 to construct a 2x2 assembly unit
pymol.cmd.select('p1','resi 21-30')
pymol.cmd.select('p2','resi 31-40')
pymol.cmd.select('p3','resi 111-120')
pymol.cmd.select('p4','resi 121-130')
# Select three atoms po1,po2,po3 to construct a reference coordinate system
pymol.cmd.select('po1','resi 22 and name ca')
pymol.cmd.select('po2','resi 29 and name ca')
pymol.cmd.select('po3','resi 62 and name ca')
# Create a periodic unit
unit = create_sheet_unit('p1','p2','p3','p4','po1','po2','po3')
# Create a fibril object
fibril = create_fibril(unit)
# Build a cross-beta structure (e.g., a ribbon)
fibril.build_a_stacked_ribbon(20,30,90,3,20,1) # This is an example of calling the function build_a_stacked_ribbon with the initial geometrical parameters (tilt angle, radius, angle for rotational stacking, number of rotational stacking, number of units in the bilayer, and the direction of twist)
```

### 4. Execute fibril_init.py in the PyMOL command line
```bash
run fibril_init.py
```

## Controllable parameters
FibrilGen offers a fibril assembly space that relates to a set of input parameters. Each FibrilGen fibril model represents a stacking pattern and an automatic adjustment of helical twist for a compact assembly. 
1. fibril.build_a_flat_sheet (N). The function “build_a_flat_sheet” takes only one parameter, N, to specify the repeat of the 2x2 unit along the beta-sheet axis (Figure 2a).
2. fibril.build_a_stacked_sheet(K, N). The function “build_a_stacked_sheet” takes input parameters N and K to stack the 2x2 unit along the fibril long axis (Figure 2a) and on the fibril cross-section (Figure 2b), respectively.
3. fibril.build_a_rod(theta_z, N, the sign of theta_y). The function “build_a_rod” takes an input parameter N to stack the 2x2 unit along the beta-sheet axis (Figure 2a). An initial helical twist is assigned with a tilt angle theta_z and the direction (assigned to 1 or -1) of the twist angle theta_y (Figure 2c).
4. fibril.build_a_stacked_rod(theta_z, K, N, the sign of theta_y). The function “build_a_stacked_rod” takes input parameters N and K for a linear stacking of the 2x2 unit along the fibril long axis (Figure 2a) and on the fibril cross-section (Figure 2b). An initial helical twist is assigned with a tilt angle theta_z and the direction (assigned as 1 or -1) of the twist angle theta_y (Figure 2c).
5. fibril.build_a_ribbon(theta_z, r_y, N, the sign of theta_y). The function “build_a_ribbon” takes an input parameter N to stack the 2x2 unit along the fibril long axis (Figure 2a). An initial helical twist is assigned with a tilt angle theta_z, a radius r_y and the direction (assigned as 1 or -1) of the twist angle theta_y (Figure 2c).
6. fibril.build_a_stacked_ribbon(theta_z, r_y, theta_s, M, the sign of theta_y). The function “build_a_stacked_ribbon” takes input parameter N to stack the 2x2 unit along the fibril long axis (Figure 2a). The rotational stacking on the fibril cross-section with an incremental rotation angle theta_s repeated for M times is assigned (Figure 2b). An initial helical twist is assigned with a tilt angle theta_z, a radius r_y and the direction (assigned as 1 or -1) of the twist angle theta_y (Figure 2c). Here tube is a special case that theta_∙M=360°.


## Contact
Chao-Yu Yang email address: chao-yu.yang@postgrad.manchester.ac.uk

## Citation
Chao-Yu Yang, Aline F. Miller, Alberto Saiani, Richard A. Bryce. FibrilGen: a Python package for atomistic modelling of peptide b-sheet nanostructures. Submitted.
