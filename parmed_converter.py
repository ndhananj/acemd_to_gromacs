#!/usr/bin/env python3
"""
Convert CHARMM PSF/PDB to GROMACS using ParmEd
Install: pip install parmed
"""

import parmed as pmd
import sys

try:
    # Load CHARMM structure
    print('Loading PSF file: dhfr.psf')
    psf = pmd.load_file('dhfr.psf')
    print('Loading PDB coordinates: dhfr.pdb')
    pdb = pmd.load_file('dhfr.pdb')

    # Transfer coordinates
    psf.coordinates = pdb.coordinates
    psf.box = pdb.box

    # Save as GROMACS
    print('Saving GROMACS topology: system.top')
    psf.save('system.top', format='gromacs')
    print('Saving GROMACS coordinates: system.gro')
    psf.save('system.gro', format='gro')

    print('\nConversion successful!')
    print('Generated files: system.top, system.gro')

except ImportError:
    print('ERROR: ParmEd not installed')
    print('Install with: pip install parmed')
    sys.exit(1)
except Exception as e:
    print(f'ERROR: {e}')
    sys.exit(1)
