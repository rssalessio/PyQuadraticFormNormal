"""pybind11 QuadraticFormDistributions plugin"""
from __future__ import annotations
from typing import List, Tuple
import numpy as np
_Shape = Tuple[int, ...]

__all__ = [
    "davies_method"
]


def davies_method(x: List[np.float64], coeff: List[np.float64], nc: List[np.float64], df: List[np.uint64], sigma: float = 1, limit: int = 10000, accuracy: float = 1e-05) -> List[np.float64]:
    """
    Distribution function of quadratic forms in normal variables using Davies’s method.

    :param x: List of points to evaluate
    :param coeff: Coefficients of the chi^2 distributions.
    :param nc: Non-centrality parameters. Needs to be the same size as coeff.
    :param df: Degrees of freedom. Needs to be the same size as coeff.
    :param sigma: Std of the gaussian 
    :param limit: Maximum number of iterations
    :param accuracy: Desired accuracy
    :result results: Probability of the evaluted points    
    """
