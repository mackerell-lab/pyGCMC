#!/bin/bash
# Script to verify atom counts in 4wp7 PDB files
# This provides the ground truth data used in the PDB parser tests
# Run this script to get exact counts directly from the raw PDB files

REGULAR_PDB="4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb"
DRUDE_PDB="4wp7_drude.pdb"

# Function to analyze regular PDB file
analyze_regular_pdb() {
    echo "=== 4wp7 Regular PDB File Analysis ==="
    echo "File: $REGULAR_PDB"
    echo

    echo "Total lines in file:"
    wc -l $REGULAR_PDB

    echo
    echo "Total ATOM records:"
    grep "^ATOM" $REGULAR_PDB | wc -l

    echo
    echo "Atom counts by residue type:"
    grep "^ATOM" $REGULAR_PDB | awk '{print $4}' | sort | uniq -c | sort -nr

    echo
    echo "=== GCMC Molecule Counts (used in tests) ==="
    echo "ACEY: $(grep "^ATOM" $REGULAR_PDB | awk '{print $4}' | grep -c "^ACEY$")"
    echo "BENX: $(grep "^ATOM" $REGULAR_PDB | awk '{print $4}' | grep -c "^BENX$")"
    echo "DMEE: $(grep "^ATOM" $REGULAR_PDB | awk '{print $4}' | grep -c "^DMEE$")"
    echo "FORM: $(grep "^ATOM" $REGULAR_PDB | awk '{print $4}' | grep -c "^FORM$")"

    echo
    echo "=== Protein Residue Counts ==="
    PROTEIN_RESIDUES="ALA ARG ASN ASP CYS GLN GLU GLY HSD ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"
    PROTEIN_TOTAL=0

    for residue in $PROTEIN_RESIDUES; do
        count=$(grep "^ATOM" $REGULAR_PDB | awk '{print $4}' | grep -c "^${residue}$")
        echo "$residue: $count"
        PROTEIN_TOTAL=$((PROTEIN_TOTAL + count))
    done

    echo
    echo "Total protein atoms: $PROTEIN_TOTAL"
    echo "Total atoms: $(grep "^ATOM" $REGULAR_PDB | wc -l)"
    echo "Non-protein atoms: $(($(grep "^ATOM" $REGULAR_PDB | wc -l) - PROTEIN_TOTAL))"

    echo
    echo "=== Crystal Cell Information ==="
    grep "^CRYST1" $REGULAR_PDB
}

# Function to analyze Drude PDB file
analyze_drude_pdb() {
    echo
    echo "=========================================="
    echo "=== 4wp7 DRUDE PDB File Analysis ==="
    echo "File: $DRUDE_PDB"
    echo
    
    if [ ! -f "$DRUDE_PDB" ]; then
        echo "ERROR: Drude PDB file not found!"
        return 1
    fi

    echo "Total lines in file:"
    wc -l $DRUDE_PDB

    echo
    echo "Total ATOM records:"
    TOTAL_ATOMS=$(grep "^ATOM" $DRUDE_PDB | wc -l)
    echo $TOTAL_ATOMS

    echo
    echo "=== DRUDE FORCE FIELD ANALYSIS ==="
    
    # Count parent atoms (standard atom names)
    PARENT_ATOMS=$(grep "^ATOM" $DRUDE_PDB | awk '{print $3}' | grep -v "^D" | grep -v "^LP" | wc -l)
    echo "Parent atoms (nuclear centers): $PARENT_ATOMS"
    
    # Count Drude particles (start with 'D', length > 1)
    DRUDE_PARTICLES=$(grep "^ATOM" $DRUDE_PDB | awk '{print $3}' | grep "^D" | grep -v "^D$" | wc -l)
    echo "Drude particles (electron clouds): $DRUDE_PARTICLES"
    
    # Count lone pairs (start with 'LP')
    LONE_PAIRS=$(grep "^ATOM" $DRUDE_PDB | awk '{print $3}' | grep "^LP" | wc -l)
    echo "Lone pairs (directional e- density): $LONE_PAIRS"
    
    # Verify total
    CALCULATED_TOTAL=$((PARENT_ATOMS + DRUDE_PARTICLES + LONE_PAIRS))
    echo "Calculated total: $CALCULATED_TOTAL"
    echo "Verification: $([ $TOTAL_ATOMS -eq $CALCULATED_TOTAL ] && echo "✓ PASS" || echo "✗ MISMATCH")"

    echo
    echo "=== DRUDE PARTICLE TYPES ==="
    echo "Unique Drude particle names:"
    grep "^ATOM" $DRUDE_PDB | awk '{print $3}' | grep "^D" | grep -v "^D$" | sort | uniq -c | sort -nr
    
    echo
    echo "=== LONE PAIR TYPES ==="
    echo "Unique lone pair names:"
    grep "^ATOM" $DRUDE_PDB | awk '{print $3}' | grep "^LP" | sort | uniq -c | sort -nr
    
    echo
    echo "=== POLARIZABLE RESIDUE TYPES ==="
    echo "Residues with Drude particles:"
    grep "^ATOM" $DRUDE_PDB | awk '{if($3 ~ /^D/ && $3 != "D") print $4}' | sort | uniq -c | sort -nr
    
    POLARIZABLE_RESIDUE_COUNT=$(grep "^ATOM" $DRUDE_PDB | awk '{if($3 ~ /^D/ && $3 != "D") print $4}' | sort | uniq | wc -l)
    echo "Total polarizable residue types: $POLARIZABLE_RESIDUE_COUNT"

    echo
    echo "=== PARENT-DRUDE PAIRING ANALYSIS ==="
    
    # Get sample of parent-Drude pairs for validation
    echo "Sample parent-Drude pairs (first 10):"
    grep "^ATOM" $DRUDE_PDB | awk '{
        if($3 ~ /^D/ && $3 != "D") {
            parent_name = substr($3, 2)  # Remove D prefix
            printf "Drude: %s -> Parent: %s (Residue: %s %s)\n", $3, parent_name, $4, $5
        }
    }' | head -10
    
    # Count orphaned Drude particles (Drude particles without corresponding parent)
    echo
    echo "Checking for orphaned Drude particles..."
    ORPHANED_COUNT=0
    grep "^ATOM" $DRUDE_PDB | awk '{
        res_key = $4 "_" $5
        if($3 ~ /^D/ && $3 != "D") {
            parent_name = substr($3, 2)
            drude_list[res_key, parent_name] = $3
        } else if($3 !~ /^LP/) {
            parent_list[res_key, $3] = 1
        }
    } END {
        orphaned = 0
        for(key in drude_list) {
            if(!(key in parent_list)) {
                orphaned++
            }
        }
        print "Orphaned Drude particles: " orphaned
    }'

    echo
    echo "=== FORCE FIELD COMPLETENESS ==="
    
    # Define potentially polarizable atom types (matching test definition)
    POLARIZABLE_TYPES="N CA C O CB CG CD CE CZ OH OG OD1 OD2 OE1 OE2"
    
    # Count potentially polarizable atoms
    POTENTIALLY_POLARIZABLE=0
    for atom_type in $POLARIZABLE_TYPES; do
        count=$(grep "^ATOM" $DRUDE_PDB | awk -v type="$atom_type" '{if($3 == type) print}' | wc -l)
        POTENTIALLY_POLARIZABLE=$((POTENTIALLY_POLARIZABLE + count))
    done
    
    echo "Potentially polarizable atoms: $POTENTIALLY_POLARIZABLE"
    echo "Actually polarized (Drude particles): $DRUDE_PARTICLES"
    
    # Calculate coverage percentage
    if [ $POTENTIALLY_POLARIZABLE -gt 0 ]; then
        COVERAGE=$(echo "scale=1; $DRUDE_PARTICLES * 100 / $POTENTIALLY_POLARIZABLE" | bc -l)
        echo "Polarization coverage: $COVERAGE%"
        
        # Check if coverage is in expected range
        COVERAGE_INT=$(echo "$COVERAGE" | cut -d. -f1)
        if [ $COVERAGE_INT -ge 20 ] && [ $COVERAGE_INT -le 200 ]; then
            echo "Coverage status: ✓ PASS (within 20-200% range)"
        else
            echo "Coverage status: ✗ WARNING (outside 20-200% range)"
        fi
    fi

    echo
    echo "=== TEST VALIDATION SUMMARY ==="
    echo "Values for test comparison:"
    echo "- total_atoms: $TOTAL_ATOMS"
    echo "- parent_atoms: $PARENT_ATOMS" 
    echo "- drude_particles: $DRUDE_PARTICLES"
    echo "- lone_pairs: $LONE_PAIRS"
    echo "- polarizable_residue_types: $POLARIZABLE_RESIDUE_COUNT"
    echo "- polarization_coverage: $COVERAGE%"
}

# Main execution
if [ -f "$REGULAR_PDB" ]; then
    analyze_regular_pdb
else
    echo "Warning: Regular PDB file $REGULAR_PDB not found"
fi

if [ -f "$DRUDE_PDB" ]; then
    analyze_drude_pdb
else
    echo "Warning: Drude PDB file $DRUDE_PDB not found"
fi

echo
echo "=========================================="
echo "This script provides the ground truth data used in pygcmc PDB parser tests"
echo "to avoid circular validation (testing pygcmc with data derived from pygcmc)."