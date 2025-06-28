#!/usr/bin/env python3
"""
Compare all OpenMM test functions between original and modular versions
"""

import re
import os

def extract_function_definitions(file_path):
    """Extract all test function definitions from a Python file."""
    functions = {}
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Find all test function definitions
    pattern = r'def (test_\w+)\(.*?\):(.*?)(?=\ndef |\nclass |\Z)'
    matches = re.findall(pattern, content, re.DOTALL)
    
    for func_name, func_body in matches:
        # Clean up the function body
        functions[func_name] = func_body.strip()
    
    return functions

def compare_functions(original_functions, modular_functions, file_description):
    """Compare functions between original and modular versions."""
    print(f"\n{'='*60}")
    print(f"Comparing {file_description}")
    print(f"{'='*60}")
    
    all_functions = set(original_functions.keys()) | set(modular_functions.keys())
    differences = []
    
    for func_name in sorted(all_functions):
        if func_name not in original_functions:
            print(f"❌ MISSING in original: {func_name}")
            differences.append(f"Missing in original: {func_name}")
        elif func_name not in modular_functions:
            print(f"❌ MISSING in modular: {func_name}")
            differences.append(f"Missing in modular: {func_name}")
        else:
            original_body = original_functions[func_name]
            modular_body = modular_functions[func_name]
            
            if original_body == modular_body:
                print(f"✅ IDENTICAL: {func_name}")
            else:
                print(f"❌ DIFFERENT: {func_name}")
                differences.append(f"Different implementation: {func_name}")
                
                # Show the differences
                print(f"   Original length: {len(original_body)} chars")
                print(f"   Modular length: {len(modular_body)} chars")
                
                # Find specific differences
                original_lines = original_body.split('\n')
                modular_lines = modular_body.split('\n')
                
                for i, (orig_line, mod_line) in enumerate(zip(original_lines, modular_lines)):
                    if orig_line.strip() != mod_line.strip():
                        print(f"   Line {i+1} differs:")
                        print(f"     Original: {orig_line.strip()}")
                        print(f"     Modular:  {mod_line.strip()}")
                        break
    
    return differences

def main():
    base_dir = "/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation"
    
    # Define file mappings
    comparisons = [
        {
            "original": f"{base_dir}/test_openmm_naive_nonbonded.py",
            "modular": f"{base_dir}/vsOpenmm/naive_comparison.py",
            "description": "Naive Nonbonded Tests"
        },
        {
            "original": f"{base_dir}/test_openmm_nonbonded_file.py", 
            "modular": f"{base_dir}/vsOpenmm/file_based.py",
            "description": "File-based Tests"
        }
    ]
    
    all_differences = []
    
    for comparison in comparisons:
        if os.path.exists(comparison["original"]) and os.path.exists(comparison["modular"]):
            original_functions = extract_function_definitions(comparison["original"])
            modular_functions = extract_function_definitions(comparison["modular"])
            
            differences = compare_functions(original_functions, modular_functions, comparison["description"])
            all_differences.extend(differences)
        else:
            print(f"❌ Missing file(s) for {comparison['description']}")
            if not os.path.exists(comparison["original"]):
                print(f"   Original not found: {comparison['original']}")
            if not os.path.exists(comparison["modular"]):
                print(f"   Modular not found: {comparison['modular']}")
    
    # For test_openmm_nonbonded.py, we need to check across multiple modular files
    large_original = f"{base_dir}/test_openmm_nonbonded.py"
    if os.path.exists(large_original):
        print(f"\n{'='*60}")
        print("Comparing Large Original File (test_openmm_nonbonded.py)")
        print(f"{'='*60}")
        
        original_functions = extract_function_definitions(large_original)
        
        # Combine all modular functions
        modular_files = [
            f"{base_dir}/vsOpenmm/simple_interactions.py",
            f"{base_dir}/vsOpenmm/intermediate_tests.py", 
            f"{base_dir}/vsOpenmm/method_comparisons.py",
            f"{base_dir}/vsOpenmm/analysis_tests.py",
            f"{base_dir}/vsOpenmm/periodic_tests.py"
        ]
        
        all_modular_functions = {}
        for modular_file in modular_files:
            if os.path.exists(modular_file):
                modular_functions = extract_function_definitions(modular_file)
                all_modular_functions.update(modular_functions)
                print(f"   Loaded {len(modular_functions)} functions from {os.path.basename(modular_file)}")
        
        differences = compare_functions(original_functions, all_modular_functions, "Large Original vs All Modular")
        all_differences.extend(differences)
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    
    if all_differences:
        print(f"❌ Found {len(all_differences)} differences:")
        for diff in all_differences:
            print(f"   - {diff}")
    else:
        print("✅ All functions are IDENTICAL between original and modular versions!")
    
    return len(all_differences) == 0

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)