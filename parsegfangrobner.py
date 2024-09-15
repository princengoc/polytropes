# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:24:30 2013

@author: tran

This program parses the output of gfan_topolyhedralcones
and count the number of combinatorial tropical types of full-dimensional polytropes.
If option count = False, it prints out a list of representatives,
one for each full-dimensional polytrope cone.

Input files: (see README of the input folder for descriptions)
input/rays.gfp4
input/cones.gfp4

"""

import numpy as np
import utilities as ut
import itertools as itt


def _rayline(line):
    if len(line) < 5:
        return None
    vchar = line.split('\t')[0].split()
    return map(lambda x: int(x), vchar)


def rayProcess(fn = 'input/rays.gfp4'):
    gf = open(fn,'r')
    rayList = []
    ray = _rayline(gf.readline())
    while ray is not None:
        rayList.append(ray)
        ray = _rayline(gf.readline())
    return np.vstack(rayList)


#now we process the cone representatives
def _coneline(line):
    if len(line) < 5:
        return None
    vchar = line.split('\t')[0]
    vchar = vchar.replace('{','')
    vchar = vchar.replace('}','')
    vchar = vchar.split()
    return map(lambda x: int(x), vchar)


def rayToPol(rays, raymat, count = True, tol = 1e-10):
    """Takes a list of rays and returns a pair (matrix, boolean).
    If the cone defined by the rays give a full-dimensional polytrope,
    returns a representative (point in the cone relative interior) and True
    Else, returns None and False
    If count, returns the boolean argument only.
    """
    m = len(rays)
    #find an interior point
    point = [0]*len(raymat[0])
    for i in rays:
        point += raymat[i]
    point = point*1./m
    #check if the point gives an acyclic matrix of shortest path.
    #if it is, then the cone is full-dimensional
    indicator = [1 if abs(i) < 1e-10 else 0 for i in point]
    mat = ut.toZeroDiagonalMatrix(indicator)
    if count:
        return ut.isAcyclic(mat)
    else:
        if(ut.isAcyclic(mat)):
            return (ut.toZeroDiagonalMatrix(point), True)
        else:
            return (None, False)


def coneProcess(raymat, fn = 'input/cones.gfp4', count = True):
    """Read in the list of cone representatives of GFP.
    Returns a list of pair (matrix, boolean), one for each cone.
    Boolean is
        True if the cone is a full-dimensional polytrope,
        False otherwise
    Matrix is
        some point in the cone interior if Boolean is True
        None otherwise
    If count is True
        Returns the total number of full-dimensional polytropes only"""
    gf = open(fn,'r')
    gf.readline() #ignore the empty cone
    poList = []
    rays = _coneline(gf.readline())
    if count:
        counter = 0
        while rays is not None:
            counter += rayToPol(rays, raymat)
            rays = _coneline(gf.readline())
        return counter
    else:
        while rays is not None:
            pol = rayToPol(rays, raymat)
            poList.append(pol)
            rays = _coneline(gf.readline())
        return poList


if __name__ == '__main__':
    raymat = rayProcess('input/rays.gfp4')
    ctr = coneProcess(raymat, fn = 'input/cones.gfp4')
    print(ctr)  #number of full-dimensional types of polytropes
