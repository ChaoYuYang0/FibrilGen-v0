

pymol.cmd.run('builder.py')

## Reset
pymol.cmd.delete('all')
pymol.cmd.reset()

## Load PDB
pymol.cmd.load('structures/input/capF8_bilayer.pdb')

## Select a reference unit
pymol.cmd.select('p1','resi 21-30')
pymol.cmd.select('p2','resi 31-40')
pymol.cmd.select('p3','resi 111-120')
pymol.cmd.select('p4','resi 121-130')
## Select a reference coordinate
pymol.cmd.select('po1','resi 22 and name ca')
pymol.cmd.select('po2','resi 29 and name ca')
pymol.cmd.select('po3','resi 62 and name ca')

## Create a periodic unit
unit = create_sheet_unit('p1','p2','p3','p4','po1','po2','po3')
## Create a fibril object
fibril = create_fibril(unit)

## --- Examples of structures ---
# fibril.build_a_flat_sheet(10)
## Build a plain sheet (stacking pattern, num of units per sheet)
# fibril.build_a_stacked_sheet([[0,1],[1,1]],10)								
## Build a rod (tilt angle, num of units per sheet, twist sign)
# fibril.build_a_rod(20,50,1) 	
## Build a stacked rod (tilt angle, stacking pattern, num of units per sheet, twist sign)								
# fibril.build_a_stacked_rod(10,[[0,1],[1,1]],10,1)	
## Build a ribbon (tilt angle, radius, num of units per sheet, twist sign)
# fibril.build_a_ribbon(30,60,10,-1)
## Build a stacked ribbon (tilt angle, radius, stacking angle, stacking number, num of units per sheet, twist sign)
fibril.build_a_stacked_ribbon(10,30,90,3,20,1)
## --- End ---

fibril.get_dimension() # Print the refined fibril dimension
pymol.cmd.zoom()