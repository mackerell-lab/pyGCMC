#!/bin/bash
# Script to run long performance tests in background

# Set up environment
export PYTHONPATH=$PYTHONPATH:./modules/bindings

# Create log directory
mkdir -p performance_logs

# Function to run test and save results
run_test() {
    local test_name=$1
    local script=$2
    local log_file="performance_logs/${test_name}_$(date +%Y%m%d_%H%M%S).log"
    
    echo "Starting $test_name in background..."
    echo "Log file: $log_file"
    
    # Run in background with nohup
    nohup python $script > $log_file 2>&1 &
    local pid=$!
    
    echo "Process ID: $pid"
    echo "To check progress: tail -f $log_file"
    echo "To check if running: ps -p $pid"
    echo ""
    
    # Save PID for later reference
    echo "$pid $test_name $log_file" >> performance_logs/running_tests.txt
}

# Main menu
echo "========================================"
echo "Background Performance Test Runner"
echo "========================================"
echo ""
echo "Available tests:"
echo "1. Quick test (8-216 waters, ~1 min)"
echo "2. Medium test (8-512 waters, ~5 min)"
echo "3. Full test (8-1000 waters, ~15 min)"
echo "4. Drude test (single system)"
echo "5. Check running tests"
echo "6. Kill all tests"
echo ""

read -p "Select option (1-6): " option

case $option in
    1)
        cat > performance_logs/quick_test.py << 'EOF'
import sys
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')
import time
import pygcmc

# Test only small systems
sizes = [(2,8), (3,27), (4,64), (5,125), (6,216)]

print("Quick Performance Test - Small Systems Only")
print("=" * 60)

# [Insert test code here - reuse from test_water_all_sizes.py]
# But only test the smaller systems
EOF
        run_test "quick_test" "performance_logs/quick_test.py"
        ;;
        
    5)
        echo "Currently running tests:"
        if [ -f performance_logs/running_tests.txt ]; then
            while IFS=' ' read -r pid name logfile; do
                if ps -p $pid > /dev/null 2>&1; then
                    echo "✓ PID $pid: $name (running)"
                    echo "  Log: $logfile"
                else
                    echo "✗ PID $pid: $name (finished)"
                    echo "  Log: $logfile"
                fi
            done < performance_logs/running_tests.txt
        else
            echo "No test records found."
        fi
        ;;
        
    6)
        echo "Killing all running tests..."
        if [ -f performance_logs/running_tests.txt ]; then
            while IFS=' ' read -r pid name logfile; do
                if ps -p $pid > /dev/null 2>&1; then
                    kill $pid
                    echo "Killed PID $pid ($name)"
                fi
            done < performance_logs/running_tests.txt
            rm performance_logs/running_tests.txt
        fi
        ;;
        
    *)
        echo "Invalid option"
        ;;
esac

echo ""
echo "Tip: You can always check results with:"
echo "  ls -la performance_logs/"
echo "  tail -f performance_logs/*.log"