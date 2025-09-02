#!/bin/bash
# Test script for single QTL scan

echo "======================================"
echo "QTL Scan Test"
echo "======================================"
echo ""

# Request minimal resources for testing
echo "This test uses chromosome 19 only - should run in < 1 minute"
echo ""

# Activate environment
echo "Activating miQTL environment..."
source /proj/raulab/users/brian/claude-test/activate_miqtl.sh

# Run the test
echo "Running test scan..."
Rscript test_single_qtl_scan.R

echo ""
echo "Test complete!"
echo "Check for output file: test_scan_*_chr19.rds"