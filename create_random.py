#!/usr/bin/python3

import netCDF4 as nc4
import numpy as np
import os
import sys


def create_random_file(size1,*args):
    try:
        os.remove('/tmp/random.ncf')
    except FileNotFoundError:
        pass
    fd = nc4.Dataset('/tmp/random.ncf', 'w')
    rand_mat = np.random.rand(size1,*args)
    sizes = rand_mat.shape
    numdim = len(sizes)
    dimname = 'a'
    dims = []
    for i in range(numdim):
        dim = fd.createDimension(dimname)
        dimname = chr(ord(dimname) + 1)
        dims.append(dim)
    v = fd.createVariable('r', rand_mat.dtype, dimensions=dims)
    v[:] = rand_mat
    fd.title = "Random matrix"
    fd.close();


if __name__=='__main__':
    sizes = []
    for i in range(1,len(sys.argv)):
        sizes.append(int(sys.argv[i]))
    create_random_file(*sizes)
    pass
