# ACEMD to GROMACS Conversion Utility

## Overview

This Python utility (`acemd_to_gromacs.py`) converts ACEMD molecular dynamics input files to GROMACS format. It handles:

- **Input formats**: PSF, PDB, PRM (CHARMM force fields), YAML configuration
- **Output formats**: GRO (coordinates), TOP (topology), MDP (simulation parameters)

## What Was Converted

Your DHFR system was successfully converted with:

- **23,558 atoms** (protein + water)
- **23,592 bonds**
- **11,584 angles**
- **6,701 dihedrals**
- **436 improper dihedrals**
- **157 CMAP terms** (backbone corrections)
- **Box size**: 62.23 × 62.23 × 62.23 Å
- **Simulation time**: 28 ns (from ACEMD config)

## Generated Files

### Coordinate and Topology Files
1. **system.gro** - Coordinates in GROMACS format (23,558 atoms)
2. **system.top** - Topology template (requires additional processing)

### Simulation Parameter Files (MDP)
3. **em.mdp** - Energy minimization (50,000 steps)
4. **nvt.mdp** - NVT equilibration (100 ps, constant volume)
5. **npt.mdp** - NPT equilibration (100 ps, constant pressure)
6. **md.mdp** - Production MD (28 ns run)

### Helper Scripts
7. **parmed_converter.py** - Automated topology conversion using ParmEd
8. **convert_with_charmm_gui.sh** - Instructions for web-based conversion
9. **README.txt** - Detailed usage instructions

## Key Features of the Utility

### 1. PSF Parser
- Extracts complete topology: atoms, bonds, angles, dihedrals, impropers, CMAP terms
- Reads CHARMM atom types and charges
- Preserves segmentation information

### 2. Coordinate Conversion
- Converts PDB to GRO format
- Handles Ångström to nanometer conversion
- Preserves residue numbering and names
- Includes box dimensions

### 3. Parameter Mapping
- Translates ACEMD simulation settings to GROMACS MDP format
- Maps thermostat settings (V-rescale for temperature coupling)
- Converts run time specifications
- Sets up PME electrostatics and cutoffs

### 4. Multi-Stage Simulation Setup
Creates a complete workflow:
- Energy minimization (steep descent)
- NVT equilibration (temperature stabilization)
- NPT equilibration (pressure/density stabilization)
- Production MD (data collection)

## Usage

### Basic Usage
```bash
python3 acemd_to_gromacs.py input.yaml
```

### Custom Output Directory
```bash
python3 acemd_to_gromacs.py input.yaml -o my_output_folder
```

### Command-Line Help
```bash
python3 acemd_to_gromacs.py --help
```

## Next Steps: Topology Conversion

**IMPORTANT**: The generated `system.top` is a template. You need to complete the topology conversion using one of these methods:

### Method 1: ParmEd (Recommended)
```bash
# Install ParmEd
pip install parmed

# Run the generated converter script
cd gromacs_files
python3 parmed_converter.py
```

This will create proper GROMACS topology files that include all force field parameters.

### Method 2: CHARMM-GUI (Web-based)
1. Visit: https://charmm-gui.org/?doc=input/pdbreader
2. Upload your `dhfr.pdb` file
3. Select "GROMACS" as the output format
4. Download the generated topology files
5. Replace the template `system.top` with the downloaded file

### Method 3: GROMACS pdb2gmx
```bash
# Download CHARMM36 force field for GROMACS
# From: http://mackerell.umaryland.edu/charmm_ff.shtml

# Use pdb2gmx with CHARMM36
gmx pdb2gmx -f dhfr.pdb -ff charmm36-jul2022 -water tip3p
```

## Running GROMACS Simulation

Once you have proper topology files:

```bash
# 1. Energy Minimization
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr
gmx mdrun -v -deffnm em

# 2. NVT Equilibration (100 ps)
gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

# 3. NPT Equilibration (100 ps)
gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -t nvt.cpt
gmx mdrun -v -deffnm npt

# 4. Production MD (28 ns)
gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -t npt.cpt
gmx mdrun -v -deffnm md
```

## Technical Details

### Force Field
- Original: CHARMM36 (all-atom protein force field with CMAP)
- Target: CHARMM36 for GROMACS
- Includes: Protein parameters, lipids, water (TIP3P), ions

### Simulation Parameters
- **Integrator**: Leap-frog (md)
- **Time step**: 2 fs
- **Constraints**: H-bonds (LINCS algorithm)
- **Electrostatics**: PME (Particle Mesh Ewald)
- **Cutoffs**: 1.2 nm (Coulomb and VdW)
- **Temperature**: 300 K (V-rescale thermostat, τ = 0.1 ps)
- **Pressure**: 1 bar (Parrinello-Rahman barostat, τ = 2.0 ps)

### Coordinate Conversion
- **Units**: Ångström → nanometer (÷10)
- **Format**: PDB (fixed-width) → GRO (GROMACS format)
- **Precision**: 3 decimal places in GRO (0.001 nm = 0.01 Å resolution)

## Customization

You can modify the generated MDP files to:
- Change simulation length (adjust `nsteps`)
- Modify output frequency (`nstxout-compressed`, `nstenergy`)
- Adjust temperature or pressure
- Change thermostat/barostat coupling constants
- Modify cutoff distances
- Enable/disable dispersion correction

## Limitations and Notes

1. **Topology Conversion**: The utility generates a topology template but requires external tools (ParmEd or CHARMM-GUI) for complete parameter assignment

2. **Force Field Files**: GROMACS needs CHARMM36 force field files installed or included

3. **System-Specific Adjustments**: You may need to adjust:
   - Temperature coupling groups (currently set to "System")
   - Pressure coupling for non-cubic boxes
   - Ion parameters if present
   - Custom residues or ligands

4. **CMAP Terms**: CHARMM36's CMAP cross-terms are preserved but require proper force field files

## Troubleshooting

### Common Issues

**Problem**: `Fatal error: Residue 'XXX' not found in residue topology database`
- **Solution**: Ensure CHARMM36 force field is properly installed in GROMACS
- Use CHARMM-GUI to generate complete topology

**Problem**: `Cannot find position restraint file`
- **Solution**: Remove or comment out position restraint includes if not needed

**Problem**: `Box size mismatch`
- **Solution**: Verify box dimensions in GRO file match topology expectations

**Problem**: `Atom type not recognized`
- **Solution**: Ensure you're using CHARMM36-compatible force field files

## References

- **GROMACS Manual**: https://manual.gromacs.org
- **CHARMM-GUI**: https://charmm-gui.org
- **ParmEd**: https://parmed.github.io
- **CHARMM Force Field**: http://mackerell.umaryland.edu/charmm_ff.shtml

## License and Citation

This utility is provided as-is for academic and research use. If you use it in published work, please cite:

- GROMACS: Abraham et al., SoftwareX (2015)
- CHARMM36: Best et al., J. Chem. Theory Comput. (2012)
- ParmEd (if used): github.com/ParmEd/ParmEd

## Support

For issues with:
- **This utility**: Check the generated README.txt and error messages
- **GROMACS**: Consult GROMACS manual and user forums
- **Force field**: CHARMM website and forums
- **Topology conversion**: CHARMM-GUI tutorials

---

**Version**: 1.0  
**Last Updated**: October 2025  
**Tested With**: GROMACS 2020+, CHARMM36-jul2022
