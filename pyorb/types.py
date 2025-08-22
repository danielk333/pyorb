"""
This module contains convenient type information so that typing can be precise 
but not too verbose in the code itself.
"""

import numpy.typing as npt

NDArray = npt.NDArray

NDArray_N = npt.NDArray
"(n,) shaped ndarray"

NDArray_3 = npt.NDArray
"(3,) shaped ndarray"

NDArray_6 = npt.NDArray
"(6,) shaped ndarray"

NDArray_3xN = npt.NDArray
"(3,n) shaped ndarray"

NDArray_6xN = npt.NDArray
"(6,n) shaped ndarray"
