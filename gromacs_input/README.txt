ACEMD to GROMACS Conversion - README
============================================================

OVERVIEW
------------------------------------------------------------
This directory contains GROMACS input files converted from
ACEMD format. The conversion includes:
  - Coordinate file (GRO format)
  - Topology template
  - MD parameter files (MDP)
  - Helper scripts for full topology conversion

IMPORTANT: TOPOLOGY CONVERSION
------------------------------------------------------------
The PSF/PRM files use CHARMM36 force field. You need to:

METHOD 1: Using ParmEd (Recommended for command-line)
  1. Install ParmEd:
     pip install parmed

  2. Run the converter:
     python3 parmed_converter.py

  This will generate proper GROMACS topology files.

METHOD 2: Using CHARMM-GUI (Web-based, most reliable)
  1. Visit: https://charmm-gui.org/?doc=input/pdbreader
  2. Upload your PDB file
  3. Select 'GROMACS' as output format
  4. Download and extract the topology files
  5. Copy the topology to this directory

METHOD 3: Manual setup with GROMACS CHARMM port
  1. Download CHARMM36 for GROMACS:
     http://mackerell.umaryland.edu/charmm_ff.shtml
  2. Use pdb2gmx with CHARMM36 force field

RUNNING GROMACS SIMULATION
------------------------------------------------------------
After obtaining proper topology files:

1. Energy Minimization:
   gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr
   gmx mdrun -v -deffnm em

2. NVT Equilibration:
   gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr
   gmx mdrun -v -deffnm nvt

3. NPT Equilibration:
   gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -t nvt.cpt
   gmx mdrun -v -deffnm npt

4. Production MD:
   gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -t npt.cpt
   gmx mdrun -v -deffnm md

NOTES
------------------------------------------------------------
- Original ACEMD run time: 28ns
- Box size: [62.23, 62.23, 62.23] Å
- Thermostat: True
- All GROMACS commands assume you have GROMACS installed
- Adjust MDP parameters as needed for your system
- Monitor equilibration carefully before production

TROUBLESHOOTING
------------------------------------------------------------
If you encounter errors:
- Ensure CHARMM36 force field is properly installed
- Check that all residues are recognized
- Verify box dimensions and PBC settings
- Consult GROMACS manual: manual.gromacs.org

For questions about conversion, check:
- GROMACS documentation
- CHARMM-GUI tutorials
- ParmEd documentation: parmed.github.io
