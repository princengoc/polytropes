# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 00:37:17 2013

@author: tran

Collection of utilities functions. These are static methods that are shared by
several classes.
"""
import itertools as itt
from functools import reduce
import numpy as np


def el(i, j, n=4):
    """Returns the ij entry of a row vector of length n^2 - n"""
    if i < j:
        return i * (n - 1) + j - 1
    else:
        return i * (n - 1) + j


def toZeroDiagonalMatrix(point):
    r"""reshape point \in \R^{n^2 -n} to a zero diagonal matrix"""
    n = int(np.ceil(np.sqrt(len(point))))
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                mat[i, j] = point[el(i, j)]
    return mat


def isAcyclic(matr):
    """Returns true if matr is acyclic, False else"""
    n = matr.shape[0]
    power = 2
    previous = matr
    while power <= n:
        previous = reduce(lambda x, _: x.dot(x), [previous, matr])
        if sum(np.diagonal(previous) == 0) < n:
            return False
        power += 1
    return True
