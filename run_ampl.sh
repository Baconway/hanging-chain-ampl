#!/usr/bin/env bash
# ==============================================================================
# Helper Script: run_ampl.sh
# Runs the AMPL models for Exercise 4 or runs the Python companion suite
# ==============================================================================

set -e

echo "=================================================================="
echo "    EXERCISE 4: THE HANGING CHAIN PROBLEM EXECUTION SCRIPT        "
echo "=================================================================="

if command -v ampl &> /dev/null; then
    echo "AMPL binary detected: $(command -v ampl)"
    echo ""
    echo "--- [1/2] Running Question 2 (chain_basic.run) ---"
    ampl chain_basic.run
    echo ""
    echo "--- [2/2] Running Question 3 (chain_obstacle.run) ---"
    ampl chain_obstacle.run
    echo ""
    echo "AMPL executions finished successfully."
else
    echo "Note: 'ampl' binary not found in standard PATH."
    echo "To run with AMPL, install AMPL and ensure 'ampl' is in your PATH,"
    echo "or run inside VS Code with the official AMPL extension."
    echo ""
    echo "Running companion Python solver & visualization suite (visualize_chain.py)..."
    python3 visualize_chain.py
fi

