pymol.cmd.run('pep2unit.py')

## Reset
pymol.cmd.delete('all')
pymol.cmd.reset()
## Load PDB
pymol.cmd.load('structures/input/AL1.pdb')
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


sheet.get_dimension() # Print the refined fibril dimension
pymol.cmd.zoom()