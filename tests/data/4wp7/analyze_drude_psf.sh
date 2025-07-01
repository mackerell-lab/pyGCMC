#!/bin/bash
# Enhanced Script to analyze Drude PSF file structure and properties
# This script provides comprehensive analysis of the 4wp7_drude.psf file
# including topology, charges, polarizability, force field parameters,
# molecular integrity validation, and detailed statistical analysis

DRUDE_PSF="4wp7_drude.psf"
DRUDE_PDB="4wp7_drude.pdb"
ENHANCED_VERSION="2.0"

# Function to analyze Drude PSF file
analyze_drude_psf() {
    echo "=========================================="
    echo "=== 4wp7 DRUDE PSF File Analysis v$ENHANCED_VERSION ==="
    echo "File: $DRUDE_PSF"
    echo "Date: $(date)"
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
    echo "=== ENHANCED FILE INTEGRITY CHECK ==="
    
    # Check file format and encoding
    FILE_ENCODING=$(file -bi "$DRUDE_PSF" | cut -d';' -f2 | cut -d'=' -f2)
    echo "File encoding: $FILE_ENCODING"
    
    # Verify PSF file structure
    SECTION_ORDER=(NATOM NBOND NTHETA NPHI NIMPHI NDON NACC NUMLP NUMANISO NCRTERM)
    PREV_LINE=0
    SECTIONS_VALID=true
    
    for section in "${SECTION_ORDER[@]}"; do
        SECTION_LINE=$(grep -n "!$section" "$DRUDE_PSF" | cut -d: -f1)
        if [ -n "$SECTION_LINE" ]; then
            if [ "$SECTION_LINE" -le "$PREV_LINE" ]; then
                echo "✗ ERROR: Section !$section out of order (line $SECTION_LINE)"
                SECTIONS_VALID=false
            fi
            PREV_LINE=$SECTION_LINE
        fi
    done
    
    if [ "$SECTIONS_VALID" = true ]; then
        echo "✓ PSF section ordering: PASS"
    else
        echo "✗ PSF section ordering: FAIL"
    fi
    
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
    
    # Enhanced Drude particle identification and validation
    echo
    echo "=== ENHANCED DRUDE PARTICLE VALIDATION ==="
    
    # More robust Drude particle identification with multiple criteria
    awk '{
        # Criteria: name starts with D, mass ~0.4, AND type is DRUD
        if($5 ~ /^D/ && $8 > 0.3 && $8 < 0.5 && $6 == "DRUD") {
            drude_count++
            mass_sum += $8
            charge_sum += $7
            # Store details for optional verbose output
            if(ENVIRON["VERBOSE"] == "1" && drude_count <= 10) {
                print $1, $3, $4, $5, $7, $8  # atom_id, res_id, res_name, atom_name, charge, mass
            }
        }
    } END {
        if(drude_count > 0) {
            printf "Drude particles found: %d\n", drude_count
            printf "Average Drude mass: %.4f amu\n", mass_sum/drude_count
            printf "Total Drude charge: %.4f e\n", charge_sum
            if(ENVIRON["VERBOSE"] == "1" && drude_count > 10) {
                printf "(Showing first 10 of %d particles in verbose mode)\n", drude_count
            }
        }
    }' /tmp/drude_atoms.tmp
    
    # Verify Drude-parent pairing
    echo
    echo "Analyzing Drude-parent atom pairing..."
    DRUDE_PARENT_PAIRS=$(awk '{
        if($5 ~ /^D/ && $8 > 0.3 && $8 < 0.5) {
            drude_res[$3]++
        } else if($5 !~ /^LP/) {
            parent_res[$3]++
        }
    } END {
        paired = 0
        unpaired = 0
        for(res in drude_res) {
            if(res in parent_res) paired++
            else unpaired++
        }
        printf "Residues with proper Drude pairing: %d\n", paired
        if(unpaired > 0) printf "WARNING: Unpaired Drude residues: %d\n", unpaired
    }' /tmp/drude_atoms.tmp)
    
    echo "$DRUDE_PARENT_PAIRS"
    
    # Analyze atom type distribution with categories
    echo
    echo "=== ENHANCED ATOM TYPE DISTRIBUTION ==="
    echo "Categorized atom types:"
    awk '{
        type = $5
        if(type ~ /^D/) category = "Drude"
        else if(type ~ /^LP/) category = "LonePair"
        else if(type ~ /^H/) category = "Hydrogen"
        else if(type ~ /^C/) category = "Carbon"
        else if(type ~ /^N/) category = "Nitrogen"
        else if(type ~ /^O/) category = "Oxygen"
        else if(type ~ /^S/) category = "Sulfur"
        else category = "Other"
        
        types[type]++
        categories[category]++
    } END {
        printf "\n--- By Category ---\n"
        for(cat in categories) {
            printf "%-12s: %6d atoms\n", cat, categories[cat]
        }
        printf "\n--- Top 20 Atom Types ---\n"
        PROCINFO["sorted_in"] = "@val_num_desc"
        count = 0
        for(type in types) {
            if(++count <= 20) printf "%-8s: %6d\n", type, types[type]
        }
    }' /tmp/drude_atoms.tmp
    
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
    
    # Enhanced polarizability analysis with detailed statistics
    echo
    echo "=== ENHANCED POLARIZABILITY ANALYSIS ==="
    
    # Extract and analyze polarizability data with residue information
    awk '{
        if($6 != "DRUD" && $5 !~ /^LP/ && $11 > 0) {
            alpha = $11
            atom_type = $5
            res_name = $4
            res_id = $3
            
            alphas[NR] = alpha
            alpha_sum += alpha
            alpha_count++
            
            # Track by atom type
            type_alpha_sum[atom_type] += alpha
            type_alpha_count[atom_type]++
            
            # Track by residue type
            res_alpha_sum[res_name] += alpha
            res_alpha_count[res_name]++
            
            # Store for detailed output
            printf "%s %s %s %.4f\n", res_id, res_name, atom_type, alpha > "/tmp/alpha_details.tmp"
        }
    } END {
        if(alpha_count > 0) {
            # Calculate statistics
            avg = alpha_sum / alpha_count
            
            # Sort alphas for median calculation
            asort(alphas)
            if(alpha_count % 2 == 0) {
                median = (alphas[alpha_count/2] + alphas[alpha_count/2 + 1]) / 2
            } else {
                median = alphas[(alpha_count+1)/2]
            }
            
            # Calculate standard deviation
            for(i in alphas) {
                diff = alphas[i] - avg
                variance += diff * diff
            }
            stdev = sqrt(variance / alpha_count)
            
            printf "=== Statistical Summary ===\n"
            printf "Total polarizable atoms: %d\n", alpha_count
            printf "Average polarizability: %.4f Å³\n", avg
            printf "Median polarizability: %.4f Å³\n", median
            printf "Std deviation: %.4f Å³\n", stdev
            printf "Min polarizability: %.4f Å³\n", alphas[1]
            printf "Max polarizability: %.4f Å³\n", alphas[alpha_count]
            
            # Distribution analysis
            printf "\n=== Distribution Analysis ===\n"
            for(i in alphas) {
                alpha = alphas[i]
                if(alpha < 0.5) dist["Small (< 0.5)"]++
                else if(alpha < 1.0) dist["Medium (0.5-1.0)"]++
                else if(alpha < 1.5) dist["Large (1.0-1.5)"]++
                else if(alpha < 2.0) dist["X-Large (1.5-2.0)"]++
                else dist["Huge (> 2.0)"]++
            }
            
            for(range in dist) {
                printf "%-20s: %4d (%.1f%%)\n", range, dist[range], 100*dist[range]/alpha_count
            }
            
            # Top atom types by polarizability
            printf "\n=== Top 10 Atom Types by Average Polarizability ===\n"
            for(type in type_alpha_sum) {
                type_avg[type] = type_alpha_sum[type] / type_alpha_count[type]
            }
            
            # Sort by average alpha
            PROCINFO["sorted_in"] = "@val_num_desc"
            count = 0
            for(type in type_avg) {
                if(++count <= 10) {
                    printf "%-8s: %.4f Å³ (n=%d)\n", type, type_avg[type], type_alpha_count[type]
                }
            }
            
            # Top residue types by total polarizability
            printf "\n=== Top 10 Residue Types by Total Polarizability ===\n"
            PROCINFO["sorted_in"] = "@val_num_desc"
            count = 0
            for(res in res_alpha_sum) {
                if(++count <= 10) {
                    printf "%-8s: %.2f Å³ total (n=%d atoms, avg=%.3f)\n", 
                           res, res_alpha_sum[res], res_alpha_count[res], 
                           res_alpha_sum[res]/res_alpha_count[res]
                }
            }
        } else {
            printf "No polarizability values found in parent atoms\n"
        }
    }' /tmp/drude_atoms.tmp
    
    # Analyze polarizability patterns
    if [ -f "/tmp/alpha_details.tmp" ]; then
        echo
        echo "=== Polarizability Pattern Analysis ==="
        
        # Check for residue-specific patterns
        echo "Checking for residue-specific polarizability patterns..."
        sort -k2,2 -k4,4n /tmp/alpha_details.tmp | awk '{
            res = $2
            alpha = $4
            if(res != prev_res && NR > 1) {
                if(count > 2) {
                    avg = sum/count
                    printf "  %s: %.3f Å³ average (n=%d)\n", prev_res, avg, count
                }
                sum = 0
                count = 0
            }
            sum += alpha
            count++
            prev_res = res
        } END {
            if(count > 2) {
                avg = sum/count
                printf "  %s: %.3f Å³ average (n=%d)\n", prev_res, avg, count
            }
        }' | head -10
        
        rm -f /tmp/alpha_details.tmp
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
    echo "=== ENHANCED CONNECTIVITY ANALYSIS ==="
    
    # Analyze bonds with more detail
    echo "Total bonds in system: $NBOND"
    echo "Bond density: $(echo "scale=3; $NBOND / $NATOM" | bc -l) bonds per atom"
    
    # Calculate average coordination
    AVG_COORD=$(echo "scale=3; 2 * $NBOND / $NATOM" | bc -l)
    echo "Average coordination: $AVG_COORD"
    
    # Extract and analyze bond section
    if [ $NBOND -gt 0 ]; then
        BOND_END=$(grep -n "!NTHETA" $DRUDE_PSF | cut -d: -f1)
        BOND_END=$((BOND_END - 1))
        
        # Sample bond analysis (first 1000 bonds)
        sed -n "$((BOND_START + 1)),$((BOND_START + 1000))p" $DRUDE_PSF | awk '{
            for(i=1; i<=NF; i+=2) {
                if(i+1 <= NF) {
                    atom1 = $i
                    atom2 = $(i+1)
                    if(atom1 != "" && atom2 != "") bond_count++
                }
            }
        } END {
            printf "Sample bond count (first 1000): %d\n", bond_count
        }'
        
        echo "Analyzing bond types..."
        # This would require atom type information from bonds
    fi
    
    # Analyze angles and dihedrals
    echo
    echo "Angle statistics:"
    echo "- Total angles: $NTHETA"
    echo "- Angles per atom: $(echo "scale=3; $NTHETA / $NATOM" | bc -l)"
    
    echo
    echo "Dihedral statistics:"
    echo "- Total dihedrals: $NPHI"
    echo "- Dihedrals per atom: $(echo "scale=3; $NPHI / $NATOM" | bc -l)"
    echo "- Total impropers: $NIMPHI"
    echo "- Impropers per atom: $(echo "scale=3; $NIMPHI / $NATOM" | bc -l)"
    
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
    
    echo
    echo "=========================================="
    echo "=== ENHANCED ANALYSIS SUMMARY ==="
    echo "=========================================="
    
    # Create comprehensive summary
    echo
    echo "1. SYSTEM COMPOSITION:"
    echo "   - Total atoms: $NATOM"
    echo "   - Parent atoms: $PSF_PARENT_ATOMS ($(echo "scale=1; $PSF_PARENT_ATOMS * 100 / $NATOM" | bc -l)%)"
    echo "   - Drude particles: $PSF_DRUDE_PARTICLES ($(echo "scale=1; $PSF_DRUDE_PARTICLES * 100 / $NATOM" | bc -l)%)"
    echo "   - Lone pairs: $PSF_LONE_PAIRS ($(echo "scale=1; $PSF_LONE_PAIRS * 100 / $NATOM" | bc -l)%)"
    
    echo
    echo "2. DRUDE POLARIZATION:"
    if [ $PSF_PARENT_ATOMS -gt 0 ] && [ $PSF_DRUDE_PARTICLES -gt 0 ]; then
        POLARIZATION_COVERAGE=$(echo "scale=1; $PSF_DRUDE_PARTICLES * 100 / $PSF_PARENT_ATOMS" | bc -l)
        echo "   - Polarization coverage: $POLARIZATION_COVERAGE% of parent atoms"
        echo "   - Polarizable atoms: $ALPHA_COUNT"
        echo "   - Average polarizability: Check detailed analysis above"
    fi
    
    echo
    echo "3. CONNECTIVITY METRICS:"
    echo "   - Bond density: $(echo "scale=2; $NBOND / $NATOM" | bc -l) bonds/atom"
    echo "   - Angle density: $(echo "scale=2; $NTHETA / $NATOM" | bc -l) angles/atom"
    echo "   - Dihedral density: $(echo "scale=2; $NPHI / $NATOM" | bc -l) dihedrals/atom"
    
    echo
    echo "4. SPECIAL FEATURES:"
    echo "   - Anisotropic sites: $NUMANISO"
    echo "   - CMAP cross-terms: $NCRTERM"
    echo "   - H-bond donors: $NDON"
    echo "   - H-bond acceptors: $NACC"
    
    echo
    echo "5. VALIDATION STATUS:"
    echo -n "   - Charge neutrality: "
    if (( $(echo "$TOTAL_CHARGE < 0.1 && $TOTAL_CHARGE > -0.1" | bc -l) )); then
        echo "✓ PASS"
    else
        echo "✗ FAIL (q=$TOTAL_CHARGE)"
    fi
    
    echo -n "   - PSF internal consistency: "
    if [ $PSF_TOTAL_CHECK -eq $NATOM ]; then
        echo "✓ PASS"
    else
        echo "✗ FAIL"
    fi
    
    if [ -f "$DRUDE_PDB" ]; then
        echo -n "   - PSF-PDB consistency: "
        if [ $NATOM -eq $PDB_ATOMS ] && [ $PSF_DRUDE_PARTICLES -eq $PDB_DRUDE ]; then
            echo "✓ PASS"
        else
            echo "✗ FAIL"
        fi
    fi
}

# Function to generate machine-readable output
generate_json_summary() {
    echo
    echo "=== JSON SUMMARY (for automated parsing) ==="
    cat << EOF
{
  "file": "$DRUDE_PSF",
  "version": "$ENHANCED_VERSION",
  "date": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "atoms": {
    "total": $NATOM,
    "parent": $PSF_PARENT_ATOMS,
    "drude": $PSF_DRUDE_PARTICLES,
    "lonepair": $PSF_LONE_PAIRS
  },
  "connectivity": {
    "bonds": $NBOND,
    "angles": $NTHETA,
    "dihedrals": $NPHI,
    "impropers": $NIMPHI
  },
  "special": {
    "anisotropic": $NUMANISO,
    "cmap": $NCRTERM,
    "donors": $NDON,
    "acceptors": $NACC
  },
  "charge": {
    "total": $TOTAL_CHARGE,
    "parent": $PARENT_CHARGE,
    "drude": $DRUDE_CHARGE
  }
}
EOF
}

# Parse command line options
VERBOSE=0
JSON_OUTPUT=0
HELP=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --verbose|-v)
            VERBOSE=1
            export VERBOSE
            shift
            ;;
        --json)
            JSON_OUTPUT=1
            shift
            ;;
        --help|-h)
            HELP=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            HELP=1
            shift
            ;;
    esac
done

# Show help if requested
if [ $HELP -eq 1 ]; then
    echo "Usage: $0 [OPTIONS]"
    echo
    echo "Analyze Drude PSF file structure and properties"
    echo
    echo "Options:"
    echo "  --verbose, -v    Show detailed output including atom listings"
    echo "  --json          Generate machine-readable JSON summary"
    echo "  --help, -h      Show this help message"
    echo
    echo "Examples:"
    echo "  $0                  # Standard analysis"
    echo "  $0 --verbose        # Include detailed atom listings"
    echo "  $0 --json           # Add JSON summary at end"
    echo "  $0 --verbose --json # Both verbose and JSON output"
    exit 0
fi

# Main execution
echo "Starting Enhanced Drude PSF analysis v$ENHANCED_VERSION..."
echo "Date: $(date)"
if [ $VERBOSE -eq 1 ]; then
    echo "Mode: VERBOSE"
fi
echo

analyze_drude_psf

# Generate JSON summary if requested
if [ $JSON_OUTPUT -eq 1 ]; then
    generate_json_summary
fi

echo
echo "=========================================="
echo "Enhanced analysis complete (v$ENHANCED_VERSION)"
echo "This data provides comprehensive validation for"
echo "PSF parsers and Drude force field understanding."
echo "Use --verbose flag for detailed atom listings."
echo "Use --json flag for machine-readable output."
echo "=========================================="