# FibrilGen
FibrilGen is a Python library that builds various cross-beta structures, where a set of parameters can control the assembled morphology. Such cross-beta structures can be compared with experimental data (e.g., reconstructed electron density from cryo-EM data) or be analyzed with MD simulations. This tool enables combinatorial peptide assembly into cross-beta nanostructures for systematically investigating or designing novel fibril structures.

## PyMol versions
v2.3.5 (commercial), v1.7.4.5 (educational), v2.3.0 (open-source)

## Generate cross-beta nanostructures from PyMOL command line
```bash
# Change the current directory to FibrilGen
cd [directory of FibrilGen]
# Import FibrilGen functions 
run builder.py
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
# Demonstrate building a stacked ribbon structure
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

### 2. Modify fibril_init.py to initialize a cross-beta nanostructure
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
fibril.build_a_stack_ribbon(20,30,90,3,20,1) # This is an example of calling the function build_a_stack_ribbon with the initial geometrical parameters (tilt angle, radius, angle for rotational stacking, number of rotational stacking, number of units in the bilayer, and the direction of twist)
```

### 3. Execute scripts in the PyMOL command line
```bash
run fibril_init.py
```

## Citation

