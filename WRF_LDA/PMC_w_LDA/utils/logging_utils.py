"""
Logging utilities module for MPI environment Python logging
Author: Xianrui Zhu
Date: 2025/6/14
"""

import os
import sys
import logging
import numpy as np

def setup_mpi_logging():
    """
    Set up logging for MPI environment
    Must be called after MPI initialization
    Returns: (rank, size) tuple
    """
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    except ImportError:
        # If no mpi4py, assume single process
        rank = 0
        size = 1
    except:
        # MPI might not be initialized yet, use default values
        rank = 0
        size = 1
    
    # Create logs directory
    log_dir = "logs"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    # Set log filename
    log_filename = os.path.join(log_dir, f"python_rank_{rank:04d}.log")
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format=f'[Rank {rank:04d}] %(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename, mode='w'),
            logging.StreamHandler(sys.stdout)  # Also output to console
        ]
    )
    
    # Redirect print output to log file
    class PrintLogger:
        def __init__(self, logger):
            self.logger = logger
            self.terminal = sys.stdout
            
        def write(self, message):
            if message.strip():  # Avoid logging empty lines
                self.logger.info(message.strip())
            # If you want to show in terminal as well, uncomment below
            # self.terminal.write(message)
            
        def flush(self):
            pass
    
    # Redirect stdout to log
    logger = logging.getLogger()
    sys.stdout = PrintLogger(logger)
    
    print(f"=== MPI process {rank}/{size} logging setup complete ===")
    print(f"Log file: {log_filename}")
    
    return rank, size


def create_logging_wrapper(original_function):
    """
    Create a logging wrapper for functions
    The wrapper maintains exactly the same function signature and behavior as the original
    Key features:
    - No return value (consistent with original function)
    - Arrays passed by reference, modify original arrays directly
    - Fully compatible with f2py callback mechanism
    
    Args:
        original_function: Original function to be wrapped
        
    Returns:
        Wrapped function with identical properties to the original
    """
    def transpond_wrapper(
            time_str: str,
            u: np.array,
            v: np.array,
            t: np.array,
            q: np.array,
            xlat: np.array,
            xlon: np.array,
            x_size : int = None,
            y_size : int = None,
            z_size : int = None
    ):
        # Set up logging on first call
        if not hasattr(transpond_wrapper, '_logging_setup'):
            rank, size = setup_mpi_logging()
            transpond_wrapper._logging_setup = True
            print(f"=== MPI process {rank}/{size} starting data assimilation ===")
        
        # Call original function with exactly the same parameter passing
        # Note: No return value, arrays u,v,t passed by reference and modified directly
        original_function(time_str, u, v, t, q, xlat, xlon)
        
        # Consistent with original function: no return value
        return
    
    return transpond_wrapper
