# pyWRF_main.py
import os

print("=== Setting environment variables ===")
# Set GFORTRAN environment variables
os.environ.update({
    'GFORTRAN_CONVERT_UNIT': 'big_endian',
    'GFORTRAN_FORMATTED_BUFFER_SIZE': '8192',
    'GFORTRAN_UNFORMATTED_BUFFER_SIZE': '8192'
})

# Critical: Limit OpenBLAS threads to avoid conflicts with MPI
os.environ.update({
    'OPENBLAS_NUM_THREADS': '1',          # Limit OpenBLAS to single thread
    'MKL_NUM_THREADS': '1',               # If using Intel MKL
    'NUMBA_NUM_THREADS': '1',             # Limit Numba threads
    'OMP_NUM_THREADS': '1'                # Limit OpenMP threads
})

print(f"GFORTRAN_CONVERT_UNIT: {os.environ.get('GFORTRAN_CONVERT_UNIT')}")
print(f"OPENBLAS_NUM_THREADS: {os.environ.get('OPENBLAS_NUM_THREADS')}")

print("=== Import and run ===")
import pywrf
from interface import transpond_class
from utils.logging_utils import create_logging_wrapper

transpond = transpond_class()

# Create transpond wrapper with logging functionality  
transpond_with_logging = create_logging_wrapper(transpond)

#print(pywrf.python_interface.python_da.__doc__)
# WRF handles MPI initialization internally
pywrf.plug.call_wrf_main(transpond_with_logging)
