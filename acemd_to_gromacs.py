#!/usr/bin/env python3
"""
ACEMD to GROMACS Converter
==========================
Converts ACEMD input files (PSF, PDB, PRM, YAML) to GROMACS format (GRO, TOP, MDP)
"""

import os
import sys
import yaml
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import argparse


class ACEMDToGromacsConverter:
    """Main converter class for ACEMD to GROMACS conversion"""
    
    def __init__(self, yaml_file: str, output_dir: str = "gromacs_files"):
        self.yaml_file = yaml_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Parse ACEMD input
        with open(yaml_file, 'r') as f:
            self.acemd_config = yaml.safe_load(f)
        
        self.base_dir = Path(yaml_file).parent
        
    def parse_psf(self, psf_file: str) -> Dict:
        """Parse PSF file to extract topology information"""
        print(f"Parsing PSF file: {psf_file}")
        
        psf_data = {
            'atoms': [],
            'bonds': [],
            'angles': [],
            'dihedrals': [],
            'impropers': [],
            'cmaps': []
        }
        
        with open(psf_file, 'r') as f:
            lines = f.readlines()
        
        section = None
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Detect sections
            if '!NATOM' in line:
                section = 'atoms'
                i += 1
                continue
            elif '!NBOND' in line:
                section = 'bonds'
                i += 1
                continue
            elif '!NTHETA' in line:
                section = 'angles'
                i += 1
                continue
            elif '!NPHI' in line:
                section = 'dihedrals'
                i += 1
                continue
            elif '!NIMPHI' in line:
                section = 'impropers'
                i += 1
                continue
            elif '!NCRTERM' in line:
                section = 'cmaps'
                i += 1
                continue
            elif '!NDON' in line or '!NACC' in line or '!NNB' in line:
                section = None
                i += 1
                continue
            
            # Parse atoms
            if section == 'atoms' and line and not line.startswith('!'):
                parts = line.split()
                if len(parts) >= 8:
                    atom = {
                        'id': int(parts[0]),
                        'segid': parts[1],
                        'resid': int(parts[2]),
                        'resname': parts[3],
                        'name': parts[4],
                        'type': parts[5],
                        'charge': float(parts[6]),
                        'mass': float(parts[7])
                    }
                    psf_data['atoms'].append(atom)
            
            # Parse bonds (4 per line)
            elif section == 'bonds' and line and not line.startswith('!'):
                parts = [int(x) for x in line.split()]
                for j in range(0, len(parts), 2):
                    if j + 1 < len(parts):
                        psf_data['bonds'].append((parts[j], parts[j+1]))
            
            # Parse angles (3 per line)
            elif section == 'angles' and line and not line.startswith('!'):
                parts = [int(x) for x in line.split()]
                for j in range(0, len(parts), 3):
                    if j + 2 < len(parts):
                        psf_data['angles'].append((parts[j], parts[j+1], parts[j+2]))
            
            # Parse dihedrals (2 per line)
            elif section == 'dihedrals' and line and not line.startswith('!'):
                parts = [int(x) for x in line.split()]
                for j in range(0, len(parts), 4):
                    if j + 3 < len(parts):
                        psf_data['dihedrals'].append((parts[j], parts[j+1], parts[j+2], parts[j+3]))
            
            # Parse impropers (2 per line)
            elif section == 'impropers' and line and not line.startswith('!'):
                parts = [int(x) for x in line.split()]
                for j in range(0, len(parts), 4):
                    if j + 3 < len(parts):
                        psf_data['impropers'].append((parts[j], parts[j+1], parts[j+2], parts[j+3]))
            
            # Parse CMAP terms
            elif section == 'cmaps' and line and not line.startswith('!'):
                parts = [int(x) for x in line.split()]
                if len(parts) >= 8:
                    psf_data['cmaps'].append(tuple(parts[:8]))
            
            i += 1
        
        print(f"  Found {len(psf_data['atoms'])} atoms")
        print(f"  Found {len(psf_data['bonds'])} bonds")
        print(f"  Found {len(psf_data['angles'])} angles")
        print(f"  Found {len(psf_data['dihedrals'])} dihedrals")
        print(f"  Found {len(psf_data['impropers'])} impropers")
        print(f"  Found {len(psf_data['cmaps'])} cmap terms")
        
        return psf_data
    
    def convert_pdb_to_gro(self, pdb_file: str, gro_file: str, box_size: List[float]):
        """Convert PDB to GROMACS GRO format"""
        print(f"Converting PDB to GRO: {pdb_file} -> {gro_file}")
        
        atoms = []
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        atom_num = int(line[6:11].strip())
                        atom_name = line[12:16].strip()
                        res_name = line[17:20].strip()
                        res_num = int(line[22:26].strip())
                        x = float(line[30:38].strip()) / 10.0  # Convert Å to nm
                        y = float(line[38:46].strip()) / 10.0
                        z = float(line[46:54].strip()) / 10.0
                        
                        atoms.append({
                            'atom_num': atom_num,
                            'res_name': res_name,
                            'atom_name': atom_name,
                            'res_num': res_num,
                            'x': x,
                            'y': y,
                            'z': z
                        })
                    except (ValueError, IndexError):
                        continue
        
        # Write GRO file
        with open(gro_file, 'w') as f:
            f.write("DHFR system converted from ACEMD\n")
            f.write(f"{len(atoms)}\n")
            
            for i, atom in enumerate(atoms, 1):
                # GROMACS GRO format: 
                # %5d%-5s%5s%5d%8.3f%8.3f%8.3f
                f.write(f"{atom['res_num']:5d}{atom['res_name']:<5s}{atom['atom_name']:>5s}"
                       f"{i:5d}{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}\n")
            
            # Box vectors (convert Å to nm)
            box_nm = [b / 10.0 for b in box_size]
            f.write(f"   {box_nm[0]:.5f}   {box_nm[1]:.5f}   {box_nm[2]:.5f}\n")
        
        print(f"  Wrote {len(atoms)} atoms to {gro_file}")
    
    def generate_topology(self, psf_data: Dict, top_file: str):
        """Generate GROMACS topology file from PSF data"""
        print(f"Generating topology file: {top_file}")
        
        # Group atoms by molecule/segment
        molecules = {}
        for atom in psf_data['atoms']:
            segid = atom['segid']
            if segid not in molecules:
                molecules[segid] = []
            molecules[segid].append(atom)
        
        with open(top_file, 'w') as f:
            f.write("; GROMACS topology file\n")
            f.write("; Converted from ACEMD/CHARMM format\n")
            f.write("; Generated by acemd_to_gromacs.py\n\n")
            
            # Include CHARMM36 force field
            f.write("; Include forcefield parameters\n")
            f.write("#include \"charmm36-jul2022.ff/forcefield.itp\"\n\n")
            
            # System information
            f.write("[ system ]\n")
            f.write("; Name\n")
            f.write("DHFR System\n\n")
            
            # Molecules
            f.write("[ molecules ]\n")
            f.write("; Compound        #mols\n")
            
            for segid in sorted(molecules.keys()):
                if segid.startswith('P'):  # Protein segment
                    f.write(f"Protein_{segid}      1\n")
                elif segid.startswith('W'):  # Water segment
                    n_waters = len(molecules[segid]) // 3  # 3 atoms per water
                    if n_waters > 0:
                        f.write(f"SOL              {n_waters}\n")
        
        print(f"  Topology written with {len(molecules)} segments")
    
    def generate_mdp_files(self):
        """Generate GROMACS MDP files for energy minimization and MD"""
        
        # Parse run time
        run_time = self.acemd_config.get('run', '28ns')
        if 'ns' in str(run_time):
            run_ns = float(str(run_time).replace('ns', ''))
            run_ps = run_ns * 1000
        else:
            run_ps = float(run_time)
        
        # Energy minimization MDP
        em_mdp = self.output_dir / "em.mdp"
        print(f"Generating energy minimization MDP: {em_mdp}")
        
        with open(em_mdp, 'w') as f:
            f.write("; Energy Minimization\n")
            f.write("integrator               = steep\n")
            f.write("emtol                    = 1000.0\n")
            f.write("emstep                   = 0.01\n")
            f.write("nsteps                   = 50000\n\n")
            
            f.write("; Output control\n")
            f.write("nstxout                  = 1000\n")
            f.write("nstvout                  = 1000\n")
            f.write("nstfout                  = 1000\n")
            f.write("nstlog                   = 100\n")
            f.write("nstenergy                = 100\n\n")
            
            f.write("; Neighbor searching\n")
            f.write("cutoff-scheme            = Verlet\n")
            f.write("nstlist                  = 10\n")
            f.write("ns_type                  = grid\n")
            f.write("pbc                      = xyz\n\n")
            
            f.write("; Electrostatics\n")
            f.write("coulombtype              = PME\n")
            f.write("rcoulomb                 = 1.2\n\n")
            
            f.write("; Van der Waals\n")
            f.write("vdwtype                  = Cut-off\n")
            f.write("rvdw                     = 1.2\n")
            f.write("DispCorr                 = EnerPres\n")
        
        # NVT equilibration MDP
        nvt_mdp = self.output_dir / "nvt.mdp"
        print(f"Generating NVT equilibration MDP: {nvt_mdp}")
        
        with open(nvt_mdp, 'w') as f:
            f.write("; NVT Equilibration\n")
            f.write("integrator               = md\n")
            f.write("dt                       = 0.002\n")
            f.write("nsteps                   = 50000     ; 100 ps\n\n")
            
            f.write("; Output control\n")
            f.write("nstxout                  = 500\n")
            f.write("nstvout                  = 500\n")
            f.write("nstenergy                = 500\n")
            f.write("nstlog                   = 500\n\n")
            
            f.write("; Bond parameters\n")
            f.write("continuation             = no\n")
            f.write("constraint_algorithm     = lincs\n")
            f.write("constraints              = h-bonds\n\n")
            
            f.write("; Neighbor searching\n")
            f.write("cutoff-scheme            = Verlet\n")
            f.write("nstlist                  = 10\n")
            f.write("ns_type                  = grid\n")
            f.write("pbc                      = xyz\n\n")
            
            f.write("; Electrostatics\n")
            f.write("coulombtype              = PME\n")
            f.write("rcoulomb                 = 1.2\n\n")
            
            f.write("; Van der Waals\n")
            f.write("vdwtype                  = Cut-off\n")
            f.write("rvdw                     = 1.2\n")
            f.write("DispCorr                 = EnerPres\n\n")
            
            f.write("; Temperature coupling\n")
            if self.acemd_config.get('thermostat', False):
                f.write("tcoupl                   = V-rescale\n")
                f.write("tc-grps                  = System\n")
                f.write("tau_t                    = 0.1\n")
                f.write("ref_t                    = 300\n")
            else:
                f.write("tcoupl                   = no\n")
            
            f.write("\n; Pressure coupling\n")
            f.write("pcoupl                   = no\n")
        
        # NPT equilibration MDP
        npt_mdp = self.output_dir / "npt.mdp"
        print(f"Generating NPT equilibration MDP: {npt_mdp}")
        
        with open(npt_mdp, 'w') as f:
            f.write("; NPT Equilibration\n")
            f.write("integrator               = md\n")
            f.write("dt                       = 0.002\n")
            f.write("nsteps                   = 50000     ; 100 ps\n\n")
            
            f.write("; Output control\n")
            f.write("nstxout                  = 500\n")
            f.write("nstvout                  = 500\n")
            f.write("nstenergy                = 500\n")
            f.write("nstlog                   = 500\n\n")
            
            f.write("; Bond parameters\n")
            f.write("continuation             = yes\n")
            f.write("constraint_algorithm     = lincs\n")
            f.write("constraints              = h-bonds\n\n")
            
            f.write("; Neighbor searching\n")
            f.write("cutoff-scheme            = Verlet\n")
            f.write("nstlist                  = 10\n")
            f.write("ns_type                  = grid\n")
            f.write("pbc                      = xyz\n\n")
            
            f.write("; Electrostatics\n")
            f.write("coulombtype              = PME\n")
            f.write("rcoulomb                 = 1.2\n\n")
            
            f.write("; Van der Waals\n")
            f.write("vdwtype                  = Cut-off\n")
            f.write("rvdw                     = 1.2\n")
            f.write("DispCorr                 = EnerPres\n\n")
            
            f.write("; Temperature coupling\n")
            if self.acemd_config.get('thermostat', False):
                f.write("tcoupl                   = V-rescale\n")
                f.write("tc-grps                  = System\n")
                f.write("tau_t                    = 0.1\n")
                f.write("ref_t                    = 300\n")
            
            f.write("\n; Pressure coupling\n")
            f.write("pcoupl                   = Parrinello-Rahman\n")
            f.write("pcoupltype               = isotropic\n")
            f.write("tau_p                    = 2.0\n")
            f.write("ref_p                    = 1.0\n")
            f.write("compressibility          = 4.5e-5\n")
        
        # Production MD MDP
        md_mdp = self.output_dir / "md.mdp"
        print(f"Generating production MD MDP: {md_mdp}")
        
        nsteps = int(run_ps / 2.0)  # dt = 2 fs
        
        with open(md_mdp, 'w') as f:
            f.write("; Production MD\n")
            f.write("integrator               = md\n")
            f.write("dt                       = 0.002\n")
            f.write(f"nsteps                   = {nsteps}     ; {run_time}\n\n")
            
            f.write("; Output control\n")
            f.write("nstxout-compressed       = 5000\n")
            f.write("compressed-x-grps        = System\n")
            f.write("nstenergy                = 5000\n")
            f.write("nstlog                   = 5000\n\n")
            
            f.write("; Bond parameters\n")
            f.write("continuation             = yes\n")
            f.write("constraint_algorithm     = lincs\n")
            f.write("constraints              = h-bonds\n\n")
            
            f.write("; Neighbor searching\n")
            f.write("cutoff-scheme            = Verlet\n")
            f.write("nstlist                  = 10\n")
            f.write("ns_type                  = grid\n")
            f.write("pbc                      = xyz\n\n")
            
            f.write("; Electrostatics\n")
            f.write("coulombtype              = PME\n")
            f.write("rcoulomb                 = 1.2\n\n")
            
            f.write("; Van der Waals\n")
            f.write("vdwtype                  = Cut-off\n")
            f.write("rvdw                     = 1.2\n")
            f.write("DispCorr                 = EnerPres\n\n")
            
            f.write("; Temperature coupling\n")
            if self.acemd_config.get('thermostat', False):
                f.write("tcoupl                   = V-rescale\n")
                f.write("tc-grps                  = System\n")
                f.write("tau_t                    = 0.1\n")
                f.write("ref_t                    = 300\n")
            
            f.write("\n; Pressure coupling\n")
            f.write("pcoupl                   = Parrinello-Rahman\n")
            f.write("pcoupltype               = isotropic\n")
            f.write("tau_p                    = 2.0\n")
            f.write("ref_p                    = 1.0\n")
            f.write("compressibility          = 4.5e-5\n")
    
    def generate_conversion_script(self):
        """Generate shell script to help with CHARMM-GUI conversion"""
        script_file = self.output_dir / "convert_with_charmm_gui.sh"
        
        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# Helper script for CHARMM to GROMACS conversion\n")
            f.write("# This script provides guidance for using CHARMM-GUI\n\n")
            
            f.write("echo 'CHARMM to GROMACS Conversion Helper'\n")
            f.write("echo '===================================='\n")
            f.write("echo ''\n")
            f.write("echo 'OPTION 1: Using CHARMM-GUI (Recommended)'\n")
            f.write("echo '-----------------------------------------'\n")
            f.write("echo '1. Go to https://charmm-gui.org/?doc=input/pdbreader'\n")
            f.write("echo '2. Upload your PDB file: dhfr.pdb'\n")
            f.write("echo '3. Select GROMACS as the output format'\n")
            f.write("echo '4. Download the generated topology files'\n")
            f.write("echo ''\n")
            f.write("echo 'OPTION 2: Using ParmEd (Command-line)'\n")
            f.write("echo '--------------------------------------'\n")
            f.write("echo 'Install ParmEd: pip install parmed'\n")
            f.write("echo 'Then run the parmed_converter.py script'\n")
            f.write("echo ''\n")
            f.write("echo 'OPTION 3: Manual conversion'\n")
            f.write("echo '---------------------------'\n")
            f.write("echo 'Use the generated topology template and adjust manually'\n")
        
        os.chmod(script_file, 0o755)
        print(f"Generated conversion helper script: {script_file}")
    
    def generate_parmed_script(self):
        """Generate Python script using ParmEd for conversion"""
        script_file = self.output_dir / "parmed_converter.py"
        
        psf_path = self.acemd_config.get('structure', 'dhfr.psf')
        pdb_path = self.acemd_config.get('coordinates', 'dhfr.pdb')
        
        with open(script_file, 'w') as f:
            f.write("#!/usr/bin/env python3\n")
            f.write('"""\n')
            f.write("Convert CHARMM PSF/PDB to GROMACS using ParmEd\n")
            f.write("Install: pip install parmed\n")
            f.write('"""\n\n')
            
            f.write("import parmed as pmd\n")
            f.write("import sys\n\n")
            
            f.write("try:\n")
            f.write("    # Load CHARMM structure\n")
            f.write(f"    print('Loading PSF file: {psf_path}')\n")
            f.write(f"    psf = pmd.load_file('{psf_path}')\n")
            f.write(f"    print('Loading PDB coordinates: {pdb_path}')\n")
            f.write(f"    pdb = pmd.load_file('{pdb_path}')\n\n")
            
            f.write("    # Transfer coordinates\n")
            f.write("    psf.coordinates = pdb.coordinates\n")
            f.write("    psf.box = pdb.box\n\n")
            
            f.write("    # Save as GROMACS\n")
            f.write("    print('Saving GROMACS topology: system.top')\n")
            f.write("    psf.save('system.top', format='gromacs')\n")
            f.write("    print('Saving GROMACS coordinates: system.gro')\n")
            f.write("    psf.save('system.gro', format='gro')\n\n")
            
            f.write("    print('\\nConversion successful!')\n")
            f.write("    print('Generated files: system.top, system.gro')\n\n")
            
            f.write("except ImportError:\n")
            f.write("    print('ERROR: ParmEd not installed')\n")
            f.write("    print('Install with: pip install parmed')\n")
            f.write("    sys.exit(1)\n")
            f.write("except Exception as e:\n")
            f.write("    print(f'ERROR: {e}')\n")
            f.write("    sys.exit(1)\n")
        
        os.chmod(script_file, 0o755)
        print(f"Generated ParmEd conversion script: {script_file}")
    
    def convert(self):
        """Main conversion function"""
        print("="*60)
        print("ACEMD to GROMACS Converter")
        print("="*60)
        print()
        
        # Get file paths from ACEMD config
        psf_file = self.base_dir / self.acemd_config['structure']
        pdb_file = self.base_dir / self.acemd_config['coordinates']
        box_size = self.acemd_config.get('boxsize', [62.23, 62.23, 62.23])
        
        # Parse PSF
        psf_data = self.parse_psf(psf_file)
        print()
        
        # Convert PDB to GRO
        gro_file = self.output_dir / "system.gro"
        self.convert_pdb_to_gro(pdb_file, gro_file, box_size)
        print()
        
        # Generate topology
        top_file = self.output_dir / "system.top"
        self.generate_topology(psf_data, top_file)
        print()
        
        # Generate MDP files
        print("Generating GROMACS MDP files...")
        self.generate_mdp_files()
        print()
        
        # Generate helper scripts
        self.generate_conversion_script()
        self.generate_parmed_script()
        print()
        
        # Generate README
        self.generate_readme()
        
        print("="*60)
        print("Conversion complete!")
        print("="*60)
        print(f"\nOutput directory: {self.output_dir}")
        print("\nGenerated files:")
        print("  - system.gro              : Coordinate file")
        print("  - system.top              : Topology template (needs CHARMM-GUI or ParmEd)")
        print("  - em.mdp                  : Energy minimization parameters")
        print("  - nvt.mdp                 : NVT equilibration parameters")
        print("  - npt.mdp                 : NPT equilibration parameters")
        print("  - md.mdp                  : Production MD parameters")
        print("  - parmed_converter.py     : ParmEd conversion script")
        print("  - convert_with_charmm_gui.sh : CHARMM-GUI helper")
        print("  - README.txt              : Detailed instructions")
        print("\nNext steps:")
        print("  1. Convert topology using ParmEd: python3 gromacs_files/parmed_converter.py")
        print("  2. Or use CHARMM-GUI: ./gromacs_files/convert_with_charmm_gui.sh")
        print("  3. Run GROMACS simulation following README.txt instructions")
    
    def generate_readme(self):
        """Generate comprehensive README file"""
        readme_file = self.output_dir / "README.txt"
        
        with open(readme_file, 'w') as f:
            f.write("ACEMD to GROMACS Conversion - README\n")
            f.write("="*60 + "\n\n")
            
            f.write("OVERVIEW\n")
            f.write("-"*60 + "\n")
            f.write("This directory contains GROMACS input files converted from\n")
            f.write("ACEMD format. The conversion includes:\n")
            f.write("  - Coordinate file (GRO format)\n")
            f.write("  - Topology template\n")
            f.write("  - MD parameter files (MDP)\n")
            f.write("  - Helper scripts for full topology conversion\n\n")
            
            f.write("IMPORTANT: TOPOLOGY CONVERSION\n")
            f.write("-"*60 + "\n")
            f.write("The PSF/PRM files use CHARMM36 force field. You need to:\n\n")
            
            f.write("METHOD 1: Using ParmEd (Recommended for command-line)\n")
            f.write("  1. Install ParmEd:\n")
            f.write("     pip install parmed\n\n")
            f.write("  2. Run the converter:\n")
            f.write("     python3 parmed_converter.py\n\n")
            f.write("  This will generate proper GROMACS topology files.\n\n")
            
            f.write("METHOD 2: Using CHARMM-GUI (Web-based, most reliable)\n")
            f.write("  1. Visit: https://charmm-gui.org/?doc=input/pdbreader\n")
            f.write("  2. Upload your PDB file\n")
            f.write("  3. Select 'GROMACS' as output format\n")
            f.write("  4. Download and extract the topology files\n")
            f.write("  5. Copy the topology to this directory\n\n")
            
            f.write("METHOD 3: Manual setup with GROMACS CHARMM port\n")
            f.write("  1. Download CHARMM36 for GROMACS:\n")
            f.write("     http://mackerell.umaryland.edu/charmm_ff.shtml\n")
            f.write("  2. Use pdb2gmx with CHARMM36 force field\n\n")
            
            f.write("RUNNING GROMACS SIMULATION\n")
            f.write("-"*60 + "\n")
            f.write("After obtaining proper topology files:\n\n")
            
            f.write("1. Energy Minimization:\n")
            f.write("   gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr\n")
            f.write("   gmx mdrun -v -deffnm em\n\n")
            
            f.write("2. NVT Equilibration:\n")
            f.write("   gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr\n")
            f.write("   gmx mdrun -v -deffnm nvt\n\n")
            
            f.write("3. NPT Equilibration:\n")
            f.write("   gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -t nvt.cpt\n")
            f.write("   gmx mdrun -v -deffnm npt\n\n")
            
            f.write("4. Production MD:\n")
            f.write("   gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -t npt.cpt\n")
            f.write("   gmx mdrun -v -deffnm md\n\n")
            
            f.write("NOTES\n")
            f.write("-"*60 + "\n")
            f.write("- Original ACEMD run time: " + str(self.acemd_config.get('run', 'N/A')) + "\n")
            f.write("- Box size: " + str(self.acemd_config.get('boxsize', 'N/A')) + " Å\n")
            f.write("- Thermostat: " + str(self.acemd_config.get('thermostat', False)) + "\n")
            f.write("- All GROMACS commands assume you have GROMACS installed\n")
            f.write("- Adjust MDP parameters as needed for your system\n")
            f.write("- Monitor equilibration carefully before production\n\n")
            
            f.write("TROUBLESHOOTING\n")
            f.write("-"*60 + "\n")
            f.write("If you encounter errors:\n")
            f.write("- Ensure CHARMM36 force field is properly installed\n")
            f.write("- Check that all residues are recognized\n")
            f.write("- Verify box dimensions and PBC settings\n")
            f.write("- Consult GROMACS manual: manual.gromacs.org\n\n")
            
            f.write("For questions about conversion, check:\n")
            f.write("- GROMACS documentation\n")
            f.write("- CHARMM-GUI tutorials\n")
            f.write("- ParmEd documentation: parmed.github.io\n")
        
        print(f"Generated README: {readme_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Convert ACEMD input files to GROMACS format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.yaml
  %(prog)s input.yaml -o my_gromacs_files
  
This tool converts ACEMD/CHARMM files (PSF, PDB, PRM, YAML) to GROMACS format.
It generates coordinate files (GRO), topology templates (TOP), and simulation
parameters (MDP) suitable for GROMACS.

Note: Full topology conversion requires additional tools (ParmEd or CHARMM-GUI).
See the generated README.txt for detailed instructions.
        """
    )
    
    parser.add_argument('yaml_file', help='ACEMD input YAML file')
    parser.add_argument('-o', '--output', default='gromacs_files',
                       help='Output directory (default: gromacs_files)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.yaml_file):
        print(f"ERROR: File not found: {args.yaml_file}")
        sys.exit(1)
    
    try:
        converter = ACEMDToGromacsConverter(args.yaml_file, args.output)
        converter.convert()
    except Exception as e:
        print(f"ERROR: Conversion failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
