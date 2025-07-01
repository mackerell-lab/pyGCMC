#!/usr/bin/env python3
"""
Drude Force Field Analysis and Testing

This module provides analysis tools for understanding Drude polarizable force field
PDB files and implements comprehensive testing for the 4wp7_drude.pdb system.

Drude Theory Background:
- Drude oscillator model for electronic polarization
- Each polarizable atom has an associated Drude particle (virtual electron cloud)
- Drude particles are connected to parent atoms by harmonic springs
- Lone pairs (LP) represent electron density in specific directions
- Enables accurate modeling of polarization effects in molecular simulations

Key Features of Drude PDB Files:
1. Parent atoms (N, CA, CB, etc.) - the nuclear centers
2. Drude particles (DN, DCA, DCB, etc.) - polarizable electron clouds
3. Lone pairs (LPOA, LPOB, etc.) - directional electron density
4. Springs connect parent-Drude pairs with force constants
5. Anisotropic polarizabilities for different atom types
"""

import os
import sys
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class DrudeAtom:
    """Represents a Drude atom system (parent + Drude particle + lone pairs)"""
    parent_name: str
    parent_coords: Tuple[float, float, float] 
    drude_name: str = None
    drude_coords: Tuple[float, float, float] = None
    lone_pairs: List[Tuple[str, Tuple[float, float, float]]] = None
    residue_name: str = ""
    residue_number: int = 0
    
    def __post_init__(self):
        if self.lone_pairs is None:
            self.lone_pairs = []
    
    @property
    def has_drude_particle(self) -> bool:
        """Check if this atom has an associated Drude particle"""
        return self.drude_name is not None
    
    @property 
    def polarization_displacement(self) -> float:
        """Calculate displacement between parent and Drude particle"""
        if not self.has_drude_particle:
            return 0.0
        
        px, py, pz = self.parent_coords
        dx, dy, dz = self.drude_coords
        return ((px-dx)**2 + (py-dy)**2 + (pz-dz)**2)**0.5


class DrudePDBAnalyzer:
    """Analyzer for Drude force field PDB files"""
    
    # Standard amino acid atom types that commonly have Drude particles
    POLARIZABLE_ATOMS = {
        'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'CZ', 'OH', 'OG', 'OD1', 'OD2', 
        'OE1', 'OE2', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'SD', 'SG'
    }
    
    def __init__(self, pdb_file: str):
        self.pdb_file = pdb_file
        self.drude_atoms: Dict[int, DrudeAtom] = {}
        self.residue_stats: Dict[str, Dict] = defaultdict(lambda: {
            'total_atoms': 0, 'parent_atoms': 0, 'drude_particles': 0, 'lone_pairs': 0
        })
        self.total_stats = {
            'total_atoms': 0, 'parent_atoms': 0, 'drude_particles': 0, 
            'lone_pairs': 0, 'polarizable_systems': 0
        }
        
    def analyze(self) -> Dict:
        """Perform comprehensive analysis of Drude PDB file"""
        atoms_by_residue = defaultdict(list)
        
        # Read and categorize all atoms
        with open(self.pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    atom_data = self._parse_atom_line(line)
                    if atom_data:
                        atoms_by_residue[atom_data['res_key']].append(atom_data)
                        self.total_stats['total_atoms'] += 1
        
        # Process each residue to identify Drude systems
        for res_key, atoms in atoms_by_residue.items():
            self._process_residue(res_key, atoms)
            
        return self._generate_summary()
    
    def _parse_atom_line(self, line: str) -> Dict:
        """Parse ATOM line from PDB file"""
        try:
            return {
                'serial': int(line[6:11].strip()),
                'name': line[12:16].strip(), 
                'res_name': line[17:21].strip(),
                'res_num': int(line[22:26].strip()),
                'x': float(line[30:38].strip()),
                'y': float(line[38:46].strip()),
                'z': float(line[46:54].strip()),
                'res_key': f"{line[17:21].strip()}_{int(line[22:26].strip())}"
            }
        except (ValueError, IndexError):
            return None
    
    def _process_residue(self, res_key: str, atoms: List[Dict]):
        """Process a residue to identify Drude atom systems"""
        res_name = atoms[0]['res_name']
        
        # Group atoms by type
        parent_atoms = {}
        drude_particles = {}
        lone_pairs = {}
        
        for atom in atoms:
            name = atom['name']
            coords = (atom['x'], atom['y'], atom['z'])
            
            if name.startswith('D') and len(name) > 1:
                # Drude particle (D + parent atom name)
                parent_name = name[1:]  # Remove 'D' prefix
                drude_particles[parent_name] = (name, coords)
                self.residue_stats[res_name]['drude_particles'] += 1
                self.total_stats['drude_particles'] += 1
                
            elif name.startswith('LP'):
                # Lone pair
                lone_pairs[name] = coords
                self.residue_stats[res_name]['lone_pairs'] += 1
                self.total_stats['lone_pairs'] += 1
                
            else:
                # Parent atom
                parent_atoms[name] = coords
                self.residue_stats[res_name]['parent_atoms'] += 1
                self.total_stats['parent_atoms'] += 1
        
        # Create DrudeAtom objects for polarizable systems
        for parent_name, parent_coords in parent_atoms.items():
            drude_info = drude_particles.get(parent_name)
            
            drude_atom = DrudeAtom(
                parent_name=parent_name,
                parent_coords=parent_coords,
                drude_name=drude_info[0] if drude_info else None,
                drude_coords=drude_info[1] if drude_info else None,
                residue_name=res_name,
                residue_number=atoms[0]['res_num']
            )
            
            # Associate lone pairs (simple heuristic - same residue)
            for lp_name, lp_coords in lone_pairs.items():
                drude_atom.lone_pairs.append((lp_name, lp_coords))
            
            if drude_atom.has_drude_particle:
                self.total_stats['polarizable_systems'] += 1
                
            self.drude_atoms[atoms[0]['serial']] = drude_atom
        
        self.residue_stats[res_name]['total_atoms'] = len(atoms)
    
    def _generate_summary(self) -> Dict:
        """Generate comprehensive analysis summary"""
        return {
            'file_info': {
                'filename': os.path.basename(self.pdb_file),
                'total_atoms': self.total_stats['total_atoms'],
                'total_residues': len(self.residue_stats)
            },
            'drude_composition': {
                'parent_atoms': self.total_stats['parent_atoms'],
                'drude_particles': self.total_stats['drude_particles'], 
                'lone_pairs': self.total_stats['lone_pairs'],
                'polarizable_systems': self.total_stats['polarizable_systems']
            },
            'residue_breakdown': dict(self.residue_stats),
            'polarization_analysis': self._analyze_polarization(),
            'force_field_validation': self._validate_drude_structure()
        }
    
    def _analyze_polarization(self) -> Dict:
        """Analyze polarization characteristics"""
        displacements = []
        polarizable_residues = set()
        
        for drude_atom in self.drude_atoms.values():
            if drude_atom.has_drude_particle:
                displacement = drude_atom.polarization_displacement
                displacements.append(displacement)
                polarizable_residues.add(drude_atom.residue_name)
        
        return {
            'avg_displacement': sum(displacements) / len(displacements) if displacements else 0,
            'max_displacement': max(displacements) if displacements else 0,
            'polarizable_residue_types': len(polarizable_residues),
            'displacement_distribution': {
                'small (<0.1A)': sum(1 for d in displacements if d < 0.1),
                'medium (0.1-0.5A)': sum(1 for d in displacements if 0.1 <= d < 0.5),
                'large (>0.5A)': sum(1 for d in displacements if d >= 0.5)
            }
        }
    
    def _validate_drude_structure(self) -> Dict:
        """Validate Drude force field structure consistency"""
        issues = []
        
        # Check for orphaned Drude particles
        orphaned_drude = self.total_stats['drude_particles'] - self.total_stats['polarizable_systems'] 
        if orphaned_drude > 0:
            issues.append(f"Found {orphaned_drude} orphaned Drude particles")
        
        # Check polarization coverage
        coverage = (self.total_stats['polarizable_systems'] / 
                   self.total_stats['parent_atoms'] * 100) if self.total_stats['parent_atoms'] > 0 else 0
        
        return {
            'validation_issues': issues,
            'polarization_coverage': f"{coverage:.1f}%",
            'structure_integrity': 'PASS' if not issues else 'ISSUES_FOUND'
        }


def print_analysis_report(analysis: Dict):
    """Print formatted analysis report"""
    print("="*60)
    print("DRUDE FORCE FIELD PDB ANALYSIS REPORT")
    print("="*60)
    
    # File info
    file_info = analysis['file_info']
    print(f"\nFile: {file_info['filename']}")
    print(f"Total Atoms: {file_info['total_atoms']:,}")
    print(f"Total Residues: {file_info['total_residues']}")
    
    # Drude composition
    comp = analysis['drude_composition'] 
    print(f"\nDrude Force Field Composition:")
    print(f"  Parent Atoms: {comp['parent_atoms']:,}")
    print(f"  Drude Particles: {comp['drude_particles']:,}")
    print(f"  Lone Pairs: {comp['lone_pairs']:,}")
    print(f"  Polarizable Systems: {comp['polarizable_systems']:,}")
    
    # Polarization analysis
    pol = analysis['polarization_analysis']
    print(f"\nPolarization Analysis:")
    print(f"  Average Displacement: {pol['avg_displacement']:.3f} Å")
    print(f"  Maximum Displacement: {pol['max_displacement']:.3f} Å")
    print(f"  Polarizable Residue Types: {pol['polarizable_residue_types']}")
    
    dist = pol['displacement_distribution']
    print(f"  Displacement Distribution:")
    print(f"    Small (<0.1Å): {dist['small (<0.1A)']}")
    print(f"    Medium (0.1-0.5Å): {dist['medium (0.1-0.5A)']}")
    print(f"    Large (>0.5Å): {dist['large (>0.5A)']}")
    
    # Validation
    val = analysis['force_field_validation']
    print(f"\nForce Field Validation:")
    print(f"  Polarization Coverage: {val['polarization_coverage']}")
    print(f"  Structure Integrity: {val['structure_integrity']}")
    if val['validation_issues']:
        print("  Issues Found:")
        for issue in val['validation_issues']:
            print(f"    - {issue}")
    
    # Top residue types by polarization
    residues = analysis['residue_breakdown']
    sorted_residues = sorted(residues.items(), 
                           key=lambda x: x[1]['drude_particles'], reverse=True)
    
    print(f"\nTop Polarizable Residue Types:")
    for i, (res_name, stats) in enumerate(sorted_residues[:10]):
        if stats['drude_particles'] > 0:
            print(f"  {i+1}. {res_name}: {stats['drude_particles']} Drude particles, "
                  f"{stats['lone_pairs']} lone pairs")


if __name__ == "__main__":
    # Example usage
    pdb_file = "/home/zhaomt/gcmc/test107/pygcmc_dev/tests/data/4wp7/4wp7_drude.pdb"
    
    if os.path.exists(pdb_file):
        print("Analyzing Drude PDB file...")
        analyzer = DrudePDBAnalyzer(pdb_file)
        analysis = analyzer.analyze()
        print_analysis_report(analysis)
    else:
        print(f"PDB file not found: {pdb_file}")