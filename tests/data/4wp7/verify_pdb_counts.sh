#!/bin/bash
# Script to verify atom counts in 4wp7 PDB file
# This provides the ground truth data used in the PDB parser tests
# Run this script to get exact counts directly from the raw PDB file

PDB_FILE="4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb"

echo "=== 4wp7 PDB File Analysis ==="
echo "File: $PDB_FILE"
echo

echo "Total lines in file:"
wc -l $PDB_FILE

echo
echo "Total ATOM records:"
grep "^ATOM" $PDB_FILE | wc -l

echo
echo "Atom counts by residue type:"
grep "^ATOM" $PDB_FILE | awk '{print $4}' | sort | uniq -c | sort -nr

echo
echo "=== GCMC Molecule Counts (used in tests) ==="
echo "ACEY: $(grep "^ATOM" $PDB_FILE | awk '{print $4}' | grep -c "^ACEY$")"
echo "BENX: $(grep "^ATOM" $PDB_FILE | awk '{print $4}' | grep -c "^BENX$")"
echo "DMEE: $(grep "^ATOM" $PDB_FILE | awk '{print $4}' | grep -c "^DMEE$")"
echo "FORM: $(grep "^ATOM" $PDB_FILE | awk '{print $4}' | grep -c "^FORM$")"

echo
echo "=== Protein Residue Counts ==="
PROTEIN_RESIDUES="ALA ARG ASN ASP CYS GLN GLU GLY HSD ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"
PROTEIN_TOTAL=0

for residue in $PROTEIN_RESIDUES; do
    count=$(grep "^ATOM" $PDB_FILE | awk '{print $4}' | grep -c "^${residue}$")
    echo "$residue: $count"
    PROTEIN_TOTAL=$((PROTEIN_TOTAL + count))
done

echo
echo "Total protein atoms: $PROTEIN_TOTAL"
echo "Total atoms: $(grep "^ATOM" $PDB_FILE | wc -l)"
echo "Non-protein atoms: $(($(grep "^ATOM" $PDB_FILE | wc -l) - PROTEIN_TOTAL))"

echo
echo "=== Crystal Cell Information ==="
grep "^CRYST1" $PDB_FILE

echo
echo "This script provides the ground truth data used in pygcmc PDB parser tests"
echo "to avoid circular validation (testing pygcmc with data derived from pygcmc)."