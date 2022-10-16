"""pybind11 QuadraticFormDistributions plugin"""
from __future__ import annotations
import QuadraticFormDistributions
import typing
import numpy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "davies_method"
]


def davies_method(x: numpy.ndarray[numpy.float64], coeff: numpy.ndarray[numpy.float64], nc: numpy.ndarray[numpy.float64], df: numpy.ndarray[numpy.uint64], sigma: float = 1, limit: int = 10000, accuracy: float = 1e-05) -> typing.Tuple[numpy.ndarray[numpy.float64], numpy.ndarray[numpy.float64], int]:
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
    :result trace: Diagnostics
    :result fault: If 0, the algorithm has terminated successfully.        
    """
