#!/usr/bin/env python
"""
Test cavity bias + thermal wavelength (Λ³) combination for detailed balance
Tests the critical case that was previously untested and prone to sign errors
"""
import pytest
import numpy as np
import math
import pygcmc


def test_cavity_bias_with_lambda_detailed_balance():
    """Test detailed balance with both cavity bias and thermal wavelength enabled.
    
    This is the critical combination that was previously untested and prone to 
    sign errors between the two cavity bias implementations.
    """
    # Setup ideal gas system for precise theoretical comparison
    state = pygcmc.MCState()
    L = 3.0  # nm
    state.info.box = np.array([L, L, L])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]      # Ideal gas - no interactions
    ff.ljSigma = [0.1]    # Small but non-zero to avoid numerical issues
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -5.0   # Higher value for better insertion rates with Λ³
    params.useCavityBias = True
    params.cavityGridSpacing = 0.2
    params.probeRadius = 0.15
    params.thermalLambdaNm = 0.3      # Non-default value to trigger WithLambda path
    params.seed = 12345               # Reproducible results
    params.updateDerivedParameters()
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Equilibrate to get some particles in the system
    np.random.seed(12345)
    for _ in range(500):
        if np.random.random() < 0.6:
            mover.attemptInsertion(state)
        else:
            mover.attemptDeletion(state)
    
    # Test microstate pairing detailed balance
    beta = 1.0 / (8.314e-3 * params.temperature)
    V = L**3
    ratios = []
    
    print(f"\\nCavity Bias + Λ³ Detailed Balance Test:")
    print(f"Temperature: {params.temperature} K")
    print(f"Chemical potential: {params.chemicalPotential} kJ/mol") 
    print(f"Thermal wavelength: {params.thermalLambdaNm} nm")
    print(f"Box volume: {V} nm³")
    print("=" * 60)
    print("N_before  P_ins/P_del  Theory   Error%   CavityBias")
    print("-" * 60)
    
    for trial in range(100):
        # Get initial state
        active_before = {i for i, r in enumerate(state.residues) if r.active}
        n_before = len(active_before)
        
        if n_before == 0:
            continue  # Skip empty system
        
        # Attempt insertion
        ins_result = mover.attemptInsertion(state)
        
        if not ins_result.accepted:
            continue
            
        # Find newly inserted residue
        active_after = {i for i, r in enumerate(state.residues) if r.active}
        new_residues = active_after - active_before
        
        if not new_residues:
            continue  # Skip if no new residue was actually added
            
        inserted_idx = list(new_residues)[0]
        
        # Immediately attempt deletion of the same particle
        del_result = mover.attemptDeletion(state, inserted_idx)
        
        if not del_result.accepted:
            # Clean up the insertion if deletion failed
            state.residues[inserted_idx].active = False
            continue
            
        p_ins = ins_result.acceptanceProbability
        p_del = del_result.acceptanceProbability
        cavity_bias = getattr(ins_result, 'cavityBiasFactor', 1.0)
        
        if p_del > 1e-10:  # Avoid division by zero
            ratio = p_ins / p_del
            
            # Theoretical ratio for ideal gas with Λ³:
            # P_ins/P_del = exp(βμ) · V/((N+1) · Λ³)
            # The cavity bias factors should cancel out in this ratio
            lambda_cubed = params.thermalLambdaNm**3
            theory = math.exp(beta * params.chemicalPotential) * V / ((n_before + 1) * lambda_cubed)
            
            error_pct = abs(ratio - theory) / max(theory, 1e-12) * 100
            ratios.append(ratio)
            
            print(f"{n_before:8d}  {ratio:10.6f}  {theory:8.6f}  {error_pct:6.1f}%   {cavity_bias:.3f}")
            
            # After fixing cavity bias sign, errors should be minimal for ideal gas
            if error_pct > 10.0:  # Much stricter after sign fix
                pytest.fail(f"Detailed balance violated with cavity+Λ: error {error_pct:.1f}% > 10%")
    
    print("-" * 60)
    
    # Overall statistics
    assert len(ratios) >= 20, f"Too few successful pairs: {len(ratios)} < 20"
    
    mean_ratio = np.mean(ratios)
    std_ratio = np.std(ratios)
    max_error = max(abs(r - mean_ratio) / max(mean_ratio, 1e-12) * 100 for r in ratios)
    
    print(f"Successful pairs: {len(ratios)}")
    print(f"Mean ratio: {mean_ratio:.6f}")
    print(f"Std ratio: {std_ratio:.6f}")
    print(f"Max deviation: {max_error:.1f}%")
    
    # The standard deviation should be reasonable (not indicating systematic bias)
    relative_std = std_ratio / max(mean_ratio, 1e-12) * 100
    assert relative_std < 30.0, f"Too much scatter in ratios: {relative_std:.1f}% > 30%"


def test_cavity_bias_lambda_consistency():
    """Test that cavity bias calculation is consistent between Λ and no-Λ paths."""
    
    # Setup two identical systems, one with Λ=1.0, one with Λ≠1.0
    state1 = pygcmc.MCState()
    state1.info.box = np.array([2.5, 2.5, 2.5])
    
    state2 = pygcmc.MCState()
    state2.info.box = np.array([2.5, 2.5, 2.5])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]  # Ideal gas
    ff.ljSigma = [0.1]
    state1.forcefield = ff
    state2.forcefield = ff
    
    # System 1: no thermal wavelength (default Λ=1.0)
    params1 = pygcmc.movement.MovementParams()
    params1.temperature = 298.15
    params1.chemicalPotential = -10.0
    params1.useCavityBias = True
    params1.thermalLambdaNm = 1.0  # Default value
    params1.seed = 12345
    params1.updateDerivedParameters()
    
    # System 2: with thermal wavelength
    params2 = pygcmc.movement.MovementParams()
    params2.temperature = 298.15
    params2.chemicalPotential = -10.0 + 3.0 * math.log(0.5)  # Adjusted for Λ³ difference
    params2.useCavityBias = True  
    params2.thermalLambdaNm = 0.5  # Different from default
    params2.seed = 12345
    params2.updateDerivedParameters()
    
    mover1 = pygcmc.movement.MovementModule()
    mover1.setParams(params1)
    
    mover2 = pygcmc.movement.MovementModule()
    mover2.setParams(params2)
    
    # Add a few particles to both systems
    np.random.seed(12345)
    for _ in range(50):
        if np.random.random() < 0.7:
            mover1.attemptInsertion(state1)
            mover2.attemptInsertion(state2)
    
    # Compare cavity bias factors - they should be similar for similar particle densities
    result1 = mover1.attemptInsertion(state1)  
    result2 = mover2.attemptInsertion(state2)
    
    bias1 = getattr(result1, 'cavityBiasFactor', 1.0)
    bias2 = getattr(result2, 'cavityBiasFactor', 1.0)
    
    # Cavity bias should be similar regardless of Λ (it's a geometric property)
    if bias1 > 0.1 and bias2 > 0.1:  # Both non-trivial
        relative_diff = abs(bias1 - bias2) / max(bias1, bias2) * 100
        assert relative_diff < 20.0, f"Cavity bias differs between Λ paths: {bias1:.3f} vs {bias2:.3f}"
    
    print(f"Cavity bias consistency test:")
    print(f"  Λ=1.0 path: bias={bias1:.3f}")
    print(f"  Λ≠1.0 path: bias={bias2:.3f}")
    print(f"  Difference: {abs(bias1-bias2)/max(bias1,bias2)*100:.1f}%")


if __name__ == "__main__":
    test_cavity_bias_with_lambda_detailed_balance()
    test_cavity_bias_lambda_consistency()
    print("\\n✓ All cavity bias + Λ³ tests passed")