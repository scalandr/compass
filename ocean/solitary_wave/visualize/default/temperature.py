#!/usr/bin/env python
'''
This script plots results from MPAS-Ocean convergence test.
'''
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import transforms
matplotlib.use('Agg')

grids = ['nonhydro', 'hydro']
nGrids = len(grids)
plt.figure(1, figsize=(12.0,6.0))

for j in range(nGrids):
    grid = grids[j]
    ncfileIC = Dataset('../../../' + grid + '/default/' + grid  + '/init.nc', 'r')
    ncfile = Dataset('../../../' + grid + '/default/' + grid + '/output.nc', 'r')
    temp = ncfile.variables['temperature'][27,0:1200,:]
    xCell = ncfileIC.variables['xCell'][0:1200]
    zMid = ncfile.variables['zMid'][27,0,:]
    L0 = 1436
    a0 = 220
    x = xCell/L0
    z = zMid/a0
    z1 = z[0:35]
    temp1 = temp[:,0:35]    
    plt.ylabel('z/a0')
    plt.subplot(2,1,j+1)
    plt.contour(x, z1, temp1.T)

plt.xlabel('x/L0')
plt.ylabel('z/a0')
ncfileIC.close()
ncfile.close()
plt.savefig('plotTemp.png')
