import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import xarray
from netCDF4 import Dataset
from numpy import array

#plt.figure(1, figsize=(12.0, 6.0))
#nx = 150
plottime = 12
#nCells = 15000
#nEdges = 45300
#nVertLevels = 68

ds = xarray.open_dataset('output.nc')
ds1 = xarray.open_dataset('init.nc')
#ds = ds.drop_vars(np.setdiff1d([j for j in ds.variables],
#                                              ['yEdge', 'streamFuncSubmeso']))
#ymean = ds.isel(Time=plottime).groupby('yEdge')
ymean = ds.isel(Time=plottime)
print(ymean.dims)
yEdge = ds1.yEdge
refZMid = ds1.refZMid
refZMidGrid, yEdgeGrid = np.meshgrid(refZMid, yEdge)
#mask = numpy.zeros((ymean.sizes['nEdges'], ymean.sizes['nVertLevels']))
mask = np.where(np.logical_and(refZMidGrid > -30, np.logical_and(yEdgeGrid < 125e3,
                                       yEdgeGrid > 75e3)),True,False)
#plt.imshow(mask)
print('mask', np.shape(mask))
print('mask type', type(mask))
print('refZMidGrid', np.shape(refZMidGrid))
print('yEdgeGrid', np.shape(yEdgeGrid))
print(np.max(yEdgeGrid[mask]))
print(np.min(yEdgeGrid[mask]))
print(np.max(refZMidGrid[mask]))
print(np.min(refZMidGrid[mask]))
#ymean = ymean.mean(dim='nVertLevels')
#print(ymean.dims)
y = ymean.streamFuncSubmeso.values
print('y', np.shape(y))
y = y[mask]
print('y mask', np.shape(y))
print('mean of y', np.mean(y), 'time (h)', plottime)
#plt.title('streamFunc plot')

#plt.savefig('streamFunc.png')
