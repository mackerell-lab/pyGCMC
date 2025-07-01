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
    
    # Count different atom types by atom name (column 5) and type (column 6)
    PSF_PARENT_ATOMS=$(awk '{if($5 !~ /^D/ && $5 !~ /^LP/) print $5}' /tmp/drude_atoms.tmp | wc -l)
    PSF_DRUDE_PARTICLES=$(awk '{if($6 == "DRUD") print $5}' /tmp/drude_atoms.tmp | wc -l)
    PSF_LONE_PAIRS=$(awk '{if($5 ~ /^LP/) print $5}' /tmp/drude_atoms.tmp | wc -l)
    
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
    DRUDE_CHARGE=$(awk '{if($6 == "DRUD") sum += $7} END {printf "%.3f", sum}' /tmp/drude_atoms.tmp)
    PARENT_CHARGE=$(awk '{if($6 != "DRUD" && $5 !~ /^LP/) sum += $7} END {printf "%.3f", sum}' /tmp/drude_atoms.tmp)
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
    
    # Analyze polarizability (alpha values in column 7)
    echo
    echo "=== POLARIZABILITY ANALYSIS ==="
    DRUDE_ALPHAS=$(awk '{if($6 == "DRUD" && $8 > 0) print $8}' /tmp/drude_atoms.tmp)
    ALPHA_COUNT=$(echo "$DRUDE_ALPHAS" | grep -v "^$" | wc -l)
    if [ $ALPHA_COUNT -gt 0 ]; then
        AVG_ALPHA=$(echo "$DRUDE_ALPHAS" | awk '{sum += $1; count++} END {if(count > 0) printf "%.3f", sum/count; else print "0"}')
        MAX_ALPHA=$(echo "$DRUDE_ALPHAS" | sort -n | tail -1)
        MIN_ALPHA=$(echo "$DRUDE_ALPHAS" | sort -n | head -1)
        echo "Drude particles with polarizability: $ALPHA_COUNT"
        echo "Average polarizability (α): $AVG_ALPHA"
        echo "Max polarizability: $MAX_ALPHA"
        echo "Min polarizability: $MIN_ALPHA"
        
        # Analyze polarizability distribution
        echo "Polarizability distribution:"
        echo "$DRUDE_ALPHAS" | awk '{
            if($1 < 0.5) small++
            else if($1 < 1.0) medium++
            else large++
        } END {
            printf "  Small (α < 0.5): %d\n", small
            printf "  Medium (0.5 ≤ α < 1.0): %d\n", medium  
            printf "  Large (α ≥ 1.0): %d\n", large
        }'
    else
        echo "No explicit polarizability values found in DRUD atoms"
    fi
    
    # Analyze Drude spring constants (if present in column 8)
    echo
    echo "=== DRUDE SPRING CONSTANT ANALYSIS ==="
    SPRING_CONSTANTS=$(awk '{if($6 == "DRUD" && $9 > 0) print $9}' /tmp/drude_atoms.tmp)
    SPRING_COUNT=$(echo "$SPRING_CONSTANTS" | grep -v "^$" | wc -l)
    if [ $SPRING_COUNT -gt 0 ]; then
        AVG_SPRING=$(echo "$SPRING_CONSTANTS" | awk '{sum += $1; count++} END {if(count > 0) printf "%.3f", sum/count; else print "0"}')
        echo "Drude particles with spring constants: $SPRING_COUNT"
        echo "Average spring constant (k): $AVG_SPRING kcal/mol/Å²"
    else
        echo "No explicit spring constants found in DRUD atoms"
    fi

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