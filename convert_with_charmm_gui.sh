#!/bin/bash
# Helper script for CHARMM to GROMACS conversion
# This script provides guidance for using CHARMM-GUI

echo 'CHARMM to GROMACS Conversion Helper'
echo '===================================='
echo ''
echo 'OPTION 1: Using CHARMM-GUI (Recommended)'
echo '-----------------------------------------'
echo '1. Go to https://charmm-gui.org/?doc=input/pdbreader'
echo '2. Upload your PDB file: dhfr.pdb'
echo '3. Select GROMACS as the output format'
echo '4. Download the generated topology files'
echo ''
echo 'OPTION 2: Using ParmEd (Command-line)'
echo '--------------------------------------'
echo 'Install ParmEd: pip install parmed'
echo 'Then run the parmed_converter.py script'
echo ''
echo 'OPTION 3: Manual conversion'
echo '---------------------------'
echo 'Use the generated topology template and adjust manually'
