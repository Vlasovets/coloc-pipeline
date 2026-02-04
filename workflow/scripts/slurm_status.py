#!/usr/bin/env python3
"""
SLURM job status checker for Snakemake cluster execution
Returns job status: running, success, or failed
"""

import sys
import subprocess
import time

jobid = sys.argv[1]

# Try to get job status from sacct
try:
    # Query SLURM for job status
    result = subprocess.run(
        ["sacct", "-j", jobid, "--format=State", "--noheader", "--parsable2"],
        capture_output=True,
        text=True,
        timeout=10
    )
    
    status = result.stdout.strip().split("\n")[0]
    
    # Map SLURM states to Snakemake expectations
    running_states = ["PENDING", "CONFIGURING", "RUNNING", "COMPLETING"]
    success_states = ["COMPLETED"]
    failed_states = ["FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL"]
    
    if status in running_states:
        print("running")
    elif status in success_states:
        print("success")
    elif status in failed_states:
        print("failed")
    else:
        # Unknown state, assume running
        print("running")
        
except subprocess.TimeoutExpired:
    print("running")
except Exception as e:
    # If sacct fails, assume job is still running
    print("running")

sys.exit(0)
