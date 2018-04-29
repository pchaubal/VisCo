#!/bin/bash
# Set job parameters
#BSUB -q serial
#BSUB -J Om
#BSUB -o Om.o%J
#BSUB -e Om.e%J

# Set number of CPUs
#BSUB -n 1


python analyze_Om.py