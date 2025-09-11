#!/bin/bash

# Script to set up CASA in x86_64 mode on Apple Silicon Macs

echo "Setting up CASA environment for x86_64 architecture..."

# Create a conda environment with x86_64 architecture
CONDA_SUBDIR=osx-64 conda create -n casa_x86 python=3.10 -c conda-forge --yes

# Activate the environment and install CASA
conda activate casa_x86
CONDA_SUBDIR=osx-64 conda install casatools casatasks -c conda-forge --yes

echo "CASA x86_64 environment setup complete!"
echo "To use this environment, run: conda activate casa_x86"
