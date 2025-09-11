#!/usr/bin/env python
"""
Performance benchmark tests for GCMC
Ensures no performance regression when adding features
"""

import pytest
import time
import os
import sys
import gc
import numpy as np

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

class PerformanceBenchmark:
    """Performance benchmark suite"""
    
    def __init__(self):
        self.baseline_rate = None  # Will be set by first run
        
    def create_test_system(self, seed=12345):
        """Create a standard test system"""
        state = pygcmc.MCState()
        state.info.box = (50.0, 50.0, 50.0)  # 50nm box
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljSigma = [3.15]
        ff.ljEps = [0.2]  # Moderate interaction
        state.forcefield = ff
        
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        template.atoms = [pygcmc.MCAtom()]
        template.atoms[0].type = 0
        template.atoms[0].charge = 0.0
        
        reservoir = pygcmc.movement.FragmentReservoir()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(seed)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(125000.0)  # 50^3 nm^3
        acceptance.setActivity(0, 1.0)  # Moderate activity
        engine.setAcceptanceCalculator(acceptance)
        
        return engine, state, reservoir
    
    def run_simulation(self, engine, steps=10000):
        """Run a standard simulation"""
        insertion_count = 0
        deletion_count = 0
        translation_count = 0
        rotation_count = 0
        
        for i in range(steps):
            move_type = i % 4
            
            if move_type == 0:
                engine.attemptInsertion(0)
                insertion_count += 1
            elif move_type == 1:
                # Try deletion
                for j in range(100):
                    instance = engine.reservoir_.getInstance(j) if hasattr(engine, 'reservoir_') else None
                    if instance and instance.isActive:
                        engine.attemptDeletion(0)
                        deletion_count += 1
                        break
            elif move_type == 2:
                # Try translation
                for j in range(100):
                    instance = engine.reservoir_.getInstance(j) if hasattr(engine, 'reservoir_') else None
                    if instance and instance.isActive:
                        engine.attemptTranslation(j)
                        translation_count += 1
                        break
            else:
                # Try rotation
                for j in range(100):
                    instance = engine.reservoir_.getInstance(j) if hasattr(engine, 'reservoir_') else None
                    if instance and instance.isActive:
                        engine.attemptRotation(j)
                        rotation_count += 1
                        break
        
        return {
            'insertions': insertion_count,
            'deletions': deletion_count,
            'translations': translation_count,
            'rotations': rotation_count,
            'total': insertion_count + deletion_count + translation_count + rotation_count
        }
    
    def measure_performance(self, config_name, env_vars, steps=10000):
        """Measure performance with given configuration"""
        
        # Set environment variables
        for key, val in env_vars.items():
            os.environ[key] = val
        
        # Force garbage collection before measurement
        gc.collect()
        
        # Create fresh system
        engine, state, reservoir = self.create_test_system()
        
        # Wire environment variables to engine settings
        if "GCMC_ENABLE_STATS" in os.environ:
            if os.environ["GCMC_ENABLE_STATS"] == "1":
                engine.enableStatistics(True)
        if "GCMC_STATS_INTERVAL" in os.environ:
            engine.setStatisticsInterval(int(os.environ["GCMC_STATS_INTERVAL"]))
        # Note: GCMC_STORE_PROB is handled via environment variable in C++
        # No need to wire it here as it's checked directly in shouldStoreProbability()
        
        # Warm up
        self.run_simulation(engine, steps=100)
        
        # Measure using CPU time and disable GC to reduce jitter
        gc_was_enabled = gc.isenabled()
        if gc_was_enabled:
            gc.disable()
        try:
            start = time.process_time()  # Use CPU time instead of wall clock time
            stats = self.run_simulation(engine, steps=steps)
            elapsed = time.process_time() - start
        finally:
            if gc_was_enabled:
                gc.enable()
        
        # Calculate rate
        rate = steps / elapsed if elapsed > 0 else float('inf')
        
        # Clean up environment variables
        for key in env_vars:
            if key in os.environ:
                del os.environ[key]
        
        return {
            'name': config_name,
            'steps': steps,
            'time': elapsed,
            'rate': rate,
            'stats': stats
        }
    
    def run_benchmark_suite(self, steps=10000):
        """Run complete benchmark suite"""
        
        configs = [
            ("baseline", {}),
            ("with_prob", {"GCMC_STORE_PROB": "1"}),
            ("with_stats", {"GCMC_ENABLE_STATS": "1"}),
            ("with_debug", {"GCMC_DEBUG": "1"}),
            ("full_features", {
                "GCMC_STORE_PROB": "1",
                "GCMC_ENABLE_STATS": "1",
                "GCMC_STATS_INTERVAL": "100"
            })
        ]
        
        results = []
        
        print(f"\n{'='*60}")
        print(f"Performance Benchmark - {steps} steps")
        print(f"{'='*60}")
        
        for name, env in configs:
            result = self.measure_performance(name, env, steps)
            results.append(result)
            
            print(f"\n{name:20s}: {result['rate']:8.0f} steps/s ({result['time']:6.2f}s)")
            print(f"  Moves: I={result['stats']['insertions']}, "
                  f"D={result['stats']['deletions']}, "
                  f"T={result['stats']['translations']}, "
                  f"R={result['stats']['rotations']}")
        
        # Set baseline if not set
        if self.baseline_rate is None:
            self.baseline_rate = results[0]['rate']
        
        # Check for regression
        print(f"\n{'='*60}")
        print("Performance Impact Analysis:")
        print(f"{'='*60}")
        
        baseline = results[0]['rate']
        for result in results[1:]:
            impact = (result['rate'] - baseline) / baseline * 100
            status = "✓" if impact > -5 else "⚠" if impact > -10 else "✗"
            print(f"{result['name']:20s}: {impact:+6.1f}% {status}")
        
        return results

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
class TestPerformance:
    """Performance regression tests"""
    
    def test_baseline_performance(self):
        """Test that baseline performance meets minimum requirements"""
        benchmark = PerformanceBenchmark()
        result = benchmark.measure_performance("baseline", {}, steps=5000)
        
        # Minimum acceptable rate (adjust based on hardware and build type)
        # Note: Debug builds and GitHub CI runners may be slower
        min_rate = 300  # steps/second
        
        assert result['rate'] >= min_rate, \
            f"Performance too low: {result['rate']:.0f} < {min_rate} steps/s"
        
        print(f"\n✓ Baseline performance: {result['rate']:.0f} steps/s")
    
    def test_probability_storage_impact(self):
        """Test that probability storage doesn't severely impact performance"""
        import numpy as np
        
        # Clean environment before test to avoid contamination
        orig_prob = os.environ.pop("GCMC_STORE_PROB", None)
        
        try:
            benchmark = PerformanceBenchmark()
            
            # Run multiple measurements alternating to reduce bias
            impacts = []
            for i in range(3):  # 3 rounds for median
                # Alternate order to avoid time-based bias
                if i % 2 == 0:
                    baseline = benchmark.measure_performance("baseline", {}, steps=5000)
                    with_prob = benchmark.measure_performance("with_prob", {"GCMC_STORE_PROB": "1"}, steps=5000)
                else:
                    with_prob = benchmark.measure_performance("with_prob", {"GCMC_STORE_PROB": "1"}, steps=5000)
                    baseline = benchmark.measure_performance("baseline", {}, steps=5000)
                
                impact = (with_prob['rate'] - baseline['rate']) / baseline['rate']
                impacts.append(impact)
            
            # Use median to reduce noise from system load variations
            median_impact = float(np.median(impacts))
            
            # Allow up to 40% performance degradation for probability storage
            # This is acceptable since it's only used for debugging/testing
            # Using median reduces false failures from load spikes
            assert median_impact > -0.40, \
                f"Probability storage impact too high: {median_impact*100:.1f}% (samples: {[f'{x*100:.1f}%' for x in impacts]})"
            
            print(f"\n✓ Probability storage impact: {median_impact*100:+.1f}% (median of {len(impacts)} runs)")
            
        finally:
            # Restore original environment
            if orig_prob is not None:
                os.environ["GCMC_STORE_PROB"] = orig_prob
    
    def test_stats_collection_impact(self):
        """Test that statistics collection has minimal impact"""
        import numpy as np
        
        # Clear any existing environment variables first
        stats_vars = ["GCMC_ENABLE_STATS", "GCMC_STATS_INTERVAL", "GCMC_STORE_PROB"]
        original_env = {}
        for var in stats_vars:
            if var in os.environ:
                original_env[var] = os.environ.pop(var)
        
        try:
            benchmark = PerformanceBenchmark()
            
            # Run multiple measurements alternating to reduce bias
            impacts = []
            for i in range(3):  # 3 rounds for median
                # Alternate order to avoid time-based bias
                if i % 2 == 0:
                    baseline = benchmark.measure_performance("baseline", {}, steps=5000)
                    with_stats = benchmark.measure_performance(
                        "with_stats", 
                        {"GCMC_ENABLE_STATS": "1", "GCMC_STATS_INTERVAL": "10000"},
                        steps=5000
                    )
                else:
                    with_stats = benchmark.measure_performance(
                        "with_stats", 
                        {"GCMC_ENABLE_STATS": "1", "GCMC_STATS_INTERVAL": "10000"},
                        steps=5000
                    )
                    baseline = benchmark.measure_performance("baseline", {}, steps=5000)
                
                impact = (with_stats['rate'] - baseline['rate']) / baseline['rate']
                impacts.append(impact)
            
            # Use median to reduce noise from system load variations
            median_impact = float(np.median(impacts))
            
            # Stats collection with large interval should have minimal impact
            # Allow up to 30% degradation (test environment has variability)
            # Using median reduces false failures from load spikes
            assert median_impact > -0.30, \
                f"Stats collection impact too high: {median_impact*100:.1f}% (samples: {[f'{x*100:.1f}%' for x in impacts]})"
            
            print(f"\n✓ Stats collection impact: {median_impact*100:+.1f}% (median of {len(impacts)} runs)")
            
        finally:
            # Restore original environment
            for var, val in original_env.items():
                os.environ[var] = val
    
    def test_memory_usage(self):
        """Test that memory usage is reasonable"""
        import tracemalloc
        
        # Start tracing
        tracemalloc.start()
        
        # Create and run system
        benchmark = PerformanceBenchmark()
        engine, state, reservoir = benchmark.create_test_system()
        benchmark.run_simulation(engine, steps=1000)
        
        # Get memory usage
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        # Convert to MB
        current_mb = current / 1024 / 1024
        peak_mb = peak / 1024 / 1024
        
        print(f"\n✓ Memory usage: current={current_mb:.1f}MB, peak={peak_mb:.1f}MB")
        
        # Check reasonable limits
        assert peak_mb < 100, f"Peak memory too high: {peak_mb:.1f}MB"

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running Performance Benchmarks...")
        
        # Run full benchmark suite
        benchmark = PerformanceBenchmark()
        results = benchmark.run_benchmark_suite(steps=10000)
        
        # Run regression tests
        test = TestPerformance()
        test.test_baseline_performance()
        test.test_probability_storage_impact()
        test.test_stats_collection_impact()
        test.test_memory_usage()
        
        print("\n✅ All performance tests passed!")
    else:
        print("PyGCMC not available, skipping benchmarks")