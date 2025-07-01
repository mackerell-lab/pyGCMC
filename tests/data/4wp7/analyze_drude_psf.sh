#!/bin/bash
# Script to analyze Drude PSF file structure and properties
# This script provides detailed analysis of the 4wp7_drude.psf file
# including topology, charges, polarizability, and force field parameters

DRUDE_PSF="4wp7_drude.psf"
DRUDE_PDB="4wp7_drude.pdb"

# Function to analyze Drude PSF file
analyze_drude_psf() {
    echo "=========================================="
    echo "=== 4wp7 DRUDE PSF File Analysis ==="
    echo "File: $DRUDE_PSF"
    echo "=========================================="
    
    if [ ! -f "$DRUDE_PSF" ]; then
        echo "ERROR: Drude PSF file not found!"
        return 1
    fi

    echo "Total lines in file:"
    wc -l $DRUDE_PSF

    echo
    echo "PSF File Header:"
    head -1 $DRUDE_PSF
    echo "Features: $(head -1 $DRUDE_PSF | tr ' ' '\n' | grep -v PSF | tr '\n' ' ')"

    echo
    echo "=== PSF SECTION ANALYSIS ==="
    
    # Extract section counts from PSF file
    NATOM=$(grep "!NATOM" $DRUDE_PSF | awk '{print $1}')
    NBOND=$(grep "!NBOND" $DRUDE_PSF | awk '{print $1}')
    NTHETA=$(grep "!NTHETA" $DRUDE_PSF | awk '{print $1}')
    NPHI=$(grep "!NPHI" $DRUDE_PSF | awk '{print $1}')
    NIMPHI=$(grep "!NIMPHI" $DRUDE_PSF | awk '{print $1}')
    NDON=$(grep "!NDON" $DRUDE_PSF | awk '{print $1}')
    NACC=$(grep "!NACC" $DRUDE_PSF | awk '{print $1}')
    NUMLP=$(grep "!NUMLP" $DRUDE_PSF | awk '{print $1}')
    NUMLPH=$(grep "!NUMLP" $DRUDE_PSF | awk '{print $2}')
    NUMANISO=$(grep "!NUMANISO" $DRUDE_PSF | awk '{print $1}')
    NCRTERM=$(grep "!NCRTERM" $DRUDE_PSF | awk '{print $1}')
    
    echo "Atoms (NATOM): $NATOM"
    echo "Bonds (NBOND): $NBOND"  
    echo "Angles (NTHETA): $NTHETA"
    echo "Dihedrals (NPHI): $NPHI"
    echo "Impropers (NIMPHI): $NIMPHI"
    echo "Donors (NDON): $NDON"
    echo "Acceptors (NACC): $NACC"
    echo "Lone Pairs (NUMLP): $NUMLP"
    echo "Lone Pair Hosts (NUMLPH): $NUMLPH"
    echo "Anisotropic (NUMANISO): $NUMANISO"
    echo "Cross-terms (NCRTERM): $NCRTERM"

    echo
    echo "=== DRUDE-SPECIFIC PSF ANALYSIS ==="
    
    # Analyze Drude-specific atom types from PSF
    ATOM_START=$(grep -n "!NATOM" $DRUDE_PSF | cut -d: -f1)
    BOND_START=$(grep -n "!NBOND" $DRUDE_PSF | cut -d: -f1)
    ATOM_END=$((BOND_START - 1))
    
    echo "Analyzing atoms from line $((ATOM_START + 1)) to $ATOM_END..."
    
    # Extract atom section and analyze (exclude empty lines)
    sed -n "$((ATOM_START + 1)),$((ATOM_END))p" $DRUDE_PSF | grep -v "^$" > /tmp/drude_atoms.tmp
    
    # Count different atom types using robust identification
    # Drude particles: name starts with D AND mass ≈ 0.4 (more robust than type check)
    PSF_DRUDE_PARTICLES=$(awk '{if($5 ~ /^D/ && $8 > 0.3 && $8 < 0.5) print $5}' /tmp/drude_atoms.tmp | wc -l)
    PSF_LONE_PAIRS=$(awk '{if($5 ~ /^LP/) print $5}' /tmp/drude_atoms.tmp | wc -l)
    PSF_PARENT_ATOMS=$(awk '{if($5 !~ /^D/ && $5 !~ /^LP/) print $5}' /tmp/drude_atoms.tmp | wc -l)
    
    echo "PSF Parent atoms: $PSF_PARENT_ATOMS"
    echo "PSF Drude particles (DRUD type): $PSF_DRUDE_PARTICLES"
    echo "PSF Lone pairs (LP* types): $PSF_LONE_PAIRS"
    echo "PSF Total: $((PSF_PARENT_ATOMS + PSF_DRUDE_PARTICLES + PSF_LONE_PAIRS))"
    
    # Analyze atom type distribution
    echo
    echo "=== ATOM TYPE DISTRIBUTION ==="
    echo "Top 15 atom types:"
    awk '{print $5}' /tmp/drude_atoms.tmp | sort | uniq -c | sort -nr | head -15
    
    # Analyze charge distribution  
    echo
    echo "=== CHARGE ANALYSIS ==="
    TOTAL_CHARGE=$(awk '{sum += $7} END {printf "%.3f", sum}' /tmp/drude_atoms.tmp)
    DRUDE_CHARGE=$(awk '{if($5 ~ /^D/ && $8 > 0.3 && $8 < 0.5) sum += $7} END {printf "%.3f", sum}' /tmp/drude_atoms.tmp)
    PARENT_CHARGE=$(awk '{if($5 !~ /^D/ && $5 !~ /^LP/) sum += $7} END {printf "%.3f", sum}' /tmp/drude_atoms.tmp)
    LP_CHARGE=$(awk '{if($5 ~ /^LP/) sum += $7} END {printf "%.3f", sum}' /tmp/drude_atoms.tmp)
    
    echo "Total system charge: $TOTAL_CHARGE"
    echo "Parent atoms charge: $PARENT_CHARGE"
    echo "Drude particles charge: $DRUDE_CHARGE"
    echo "Lone pairs charge: $LP_CHARGE"
    
    # Check charge neutrality
    if (( $(echo "$TOTAL_CHARGE < 0.1 && $TOTAL_CHARGE > -0.1" | bc -l) )); then
        echo "Charge neutrality: ✓ PASS (|q| < 0.1)"
    else
        echo "Charge neutrality: ✗ WARNING (|q| = $TOTAL_CHARGE)"
    fi
    
    # Analyze polarizability (alpha values in column 11 for parent atoms)
    echo
    echo "=== POLARIZABILITY ANALYSIS ==="
    # Alpha values are stored with parent atoms, not Drude particles
    PARENT_ALPHAS=$(awk '{if($6 != "DRUD" && $5 !~ /^LP/ && $11 > 0) print $11}' /tmp/drude_atoms.tmp)
    ALPHA_COUNT=$(echo "$PARENT_ALPHAS" | grep -v "^$" | wc -l)
    if [ $ALPHA_COUNT -gt 0 ]; then
        AVG_ALPHA=$(echo "$PARENT_ALPHAS" | awk '{sum += $1; count++} END {if(count > 0) printf "%.3f", sum/count; else print "0"}')
        MAX_ALPHA=$(echo "$PARENT_ALPHAS" | sort -n | tail -1)
        MIN_ALPHA=$(echo "$PARENT_ALPHAS" | sort -n | head -1)
        echo "Parent atoms with polarizability: $ALPHA_COUNT"
        echo "Average polarizability (α): $AVG_ALPHA Å³"
        echo "Max polarizability: $MAX_ALPHA Å³"
        echo "Min polarizability: $MIN_ALPHA Å³"
        
        # Analyze polarizability distribution
        echo "Polarizability distribution:"
        echo "$PARENT_ALPHAS" | awk '{
            if($1 < 0.5) small++
            else if($1 < 1.0) medium++
            else if($1 < 2.0) large++
            else xlarge++
        } END {
            printf "  Small (α < 0.5): %d\n", small+0
            printf "  Medium (0.5 ≤ α < 1.0): %d\n", medium+0
            printf "  Large (1.0 ≤ α < 2.0): %d\n", large+0
            printf "  X-Large (α ≥ 2.0): %d\n", xlarge+0
        }'
        
        # Show sample alpha values by atom type
        echo "Sample alpha values by atom type:"
        awk '{if($6 != "DRUD" && $5 !~ /^LP/ && $11 > 0) print $5, $11}' /tmp/drude_atoms.tmp | sort | uniq | head -10
    else
        echo "No polarizability values found in parent atoms"
    fi
    
    # Analyze Thole screening parameters (column 10 for parent atoms)
    echo
    echo "=== THOLE SCREENING ANALYSIS ==="
    THOLE_VALUES=$(awk '{if($6 != "DRUD" && $5 !~ /^LP/ && $10 != 0) print $10}' /tmp/drude_atoms.tmp)
    THOLE_COUNT=$(echo "$THOLE_VALUES" | grep -v "^$" | wc -l)
    if [ $THOLE_COUNT -gt 0 ]; then
        AVG_THOLE=$(echo "$THOLE_VALUES" | awk '{sum += $1; count++} END {if(count > 0) printf "%.3f", sum/count; else print "0"}')
        echo "Parent atoms with Thole parameters: $THOLE_COUNT"
        echo "Average Thole parameter: $AVG_THOLE"
        
        # Note: Thole values are often negative in PSF files
        echo "Thole parameter range:"
        echo "$THOLE_VALUES" | sort -n | awk 'NR==1{min=$1} END{print "  Min: " min "  Max: " $1}'
    else
        echo "No Thole screening parameters found"
    fi
    
    # Analyze Drude particle masses (should be ~0.4 amu)
    echo
    echo "=== DRUDE PARTICLE MASS ANALYSIS ==="
    DRUDE_MASSES=$(awk '{if($5 ~ /^D/ && $8 > 0.3 && $8 < 0.5) print $8}' /tmp/drude_atoms.tmp)
    DRUDE_MASS_COUNT=$(echo "$DRUDE_MASSES" | grep -v "^$" | wc -l)
    if [ $DRUDE_MASS_COUNT -gt 0 ]; then
        AVG_DRUDE_MASS=$(echo "$DRUDE_MASSES" | awk '{sum += $1; count++} END {if(count > 0) printf "%.3f", sum/count; else print "0"}')
        echo "Drude particles: $DRUDE_MASS_COUNT"
        echo "Average Drude mass: $AVG_DRUDE_MASS amu (should be ~0.4)"
        
        # Check for mass consistency
        MASS_CONSISTENCY=$(echo "$DRUDE_MASSES" | awk '{if($1 < 0.3 || $1 > 0.5) bad++} END {if(bad) print "INCONSISTENT"; else print "CONSISTENT"}')
        echo "Mass consistency: $MASS_CONSISTENCY"
    fi
    
    # Note about spring constants
    echo
    echo "=== DRUDE SPRING CONSTANT INFO ==="
    echo "Note: Drude spring constants (k) are not stored per-atom in PSF files."
    echo "They are typically defined globally in parameter files as:"
    echo "  BOND DRUD X   500.000 0.000  (k=500 kcal/mol/Å²)"
    echo "  or KDRUDE = 500.0 in CHARMM"

    echo
    echo "=== CONNECTIVITY ANALYSIS ==="
    
    # Analyze bonds
    echo "Total bonds in system: $NBOND"
    echo "Bond density: $(echo "scale=2; $NBOND / $NATOM" | bc -l) bonds per atom"
    
    # Calculate average coordination
    echo "Average coordination: $(echo "scale=2; 2 * $NBOND / $NATOM" | bc -l)"
    
    echo
    echo "=== LONE PAIR ANALYSIS ==="
    echo "Lone pairs (NUMLP): $NUMLP"
    echo "Lone pair hosts (NUMLPH): $NUMLPH"
    if [ $NUMLP -gt 0 ] && [ $NUMLPH -gt 0 ]; then
        echo "Average lone pairs per host: $(echo "scale=2; $NUMLP / $NUMLPH" | bc -l)"
    fi
    
    # Analyze lone pair types
    echo "Lone pair types in PSF:"
    awk '{if($5 ~ /^LP/) print $5}' /tmp/drude_atoms.tmp | sort | uniq -c | sort -nr
    
    echo
    echo "=== ANISOTROPY ANALYSIS ==="
    echo "Anisotropic sites: $NUMANISO"
    if [ $NUMANISO -gt 0 ] && [ $PSF_DRUDE_PARTICLES -gt 0 ]; then
        echo "Anisotropic fraction: $(echo "scale=2; $NUMANISO * 100 / $PSF_DRUDE_PARTICLES" | bc -l)% of Drude particles"
    fi
    
    echo
    echo "=== CROSS-TERM (CMAP) ANALYSIS ==="
    echo "Cross-terms: $NCRTERM"
    if [ $NCRTERM -gt 0 ]; then
        echo "CMAP density: $(echo "scale=2; $NCRTERM * 100 / $NATOM" | bc -l) cross-terms per 100 atoms"
    fi
    
    echo
    echo "=== HYDROGEN BONDING ANALYSIS ==="
    echo "Hydrogen bond donors: $NDON"
    echo "Hydrogen bond acceptors: $NACC"
    echo "Total H-bonding sites: $((NDON + NACC))"
    if [ $NATOM -gt 0 ]; then
        echo "H-bonding density: $(echo "scale=2; ($NDON + $NACC) * 100 / $NATOM" | bc -l)% of atoms"
    fi

    echo
    echo "=== PSF VALIDATION SUMMARY ==="
    if [ -f "$DRUDE_PDB" ]; then
        echo "PSF vs PDB consistency check:"
        PDB_ATOMS=$(grep "^ATOM" $DRUDE_PDB | wc -l)
        PDB_DRUDE=$(grep "^ATOM" $DRUDE_PDB | awk '{print $3}' | grep "^D" | grep -v "^D$" | wc -l)
        PDB_LP=$(grep "^ATOM" $DRUDE_PDB | awk '{print $3}' | grep "^LP" | wc -l)
        
        echo "- PSF atoms: $NATOM vs PDB atoms: $PDB_ATOMS"
        echo "- PSF Drude particles: $PSF_DRUDE_PARTICLES vs PDB Drude particles: $PDB_DRUDE"
        echo "- PSF lone pairs: $PSF_LONE_PAIRS vs PDB lone pairs: $PDB_LP"
        
        # Consistency checks
        if [ $NATOM -eq $PDB_ATOMS ]; then
            echo "Atom count consistency: ✓ PASS"
        else
            echo "Atom count consistency: ✗ MISMATCH"
        fi
        
        if [ $PSF_DRUDE_PARTICLES -eq $PDB_DRUDE ]; then
            echo "Drude particle consistency: ✓ PASS"
        else
            echo "Drude particle consistency: ✗ MISMATCH"
        fi
        
        if [ $PSF_LONE_PAIRS -eq $PDB_LP ]; then
            echo "Lone pair consistency: ✓ PASS"
        else
            echo "Lone pair consistency: ✗ MISMATCH"
        fi
    else
        echo "PDB file not found - skipping consistency check"
    fi
    
    # PSF internal consistency
    PSF_TOTAL_CHECK=$((PSF_PARENT_ATOMS + PSF_DRUDE_PARTICLES + PSF_LONE_PAIRS))
    if [ $PSF_TOTAL_CHECK -eq $NATOM ]; then
        echo "PSF internal consistency: ✓ PASS"
    else
        echo "PSF internal consistency: ✗ MISMATCH ($PSF_TOTAL_CHECK vs $NATOM)"
    fi
    
    echo
    echo "=== FORCE FIELD QUALITY METRICS ==="
    
    # Calculate various quality metrics
    if [ $PSF_DRUDE_PARTICLES -gt 0 ] && [ $PSF_PARENT_ATOMS -gt 0 ]; then
        POLARIZATION_RATIO=$(echo "scale=3; $PSF_DRUDE_PARTICLES / $PSF_PARENT_ATOMS" | bc -l)
        echo "Polarization ratio (Drude/Parent): $POLARIZATION_RATIO"
        
        if (( $(echo "$POLARIZATION_RATIO > 0.3 && $POLARIZATION_RATIO < 0.8" | bc -l) )); then
            echo "Polarization coverage: ✓ GOOD (30-80% coverage)"
        elif (( $(echo "$POLARIZATION_RATIO >= 0.2" | bc -l) )); then
            echo "Polarization coverage: ○ MODERATE (≥20% coverage)"
        else
            echo "Polarization coverage: ✗ LOW (<20% coverage)"
        fi
    fi
    
    # Topology richness
    if [ $NATOM -gt 0 ]; then
        TOPOLOGY_RICHNESS=$(echo "scale=3; ($NBOND + $NTHETA + $NPHI) / $NATOM" | bc -l)
        echo "Topology richness (bonds+angles+dihedrals per atom): $TOPOLOGY_RICHNESS"
    fi
    
    # Clean up temp files
    rm -f /tmp/drude_atoms.tmp
}

# Main execution
echo "Starting Drude PSF analysis..."
echo "Date: $(date)"
echo

analyze_drude_psf

echo
echo "=========================================="
echo "Analysis complete. This data can be used to validate"
echo "PSF parsers and understand Drude force field structure."
echo "=========================================="