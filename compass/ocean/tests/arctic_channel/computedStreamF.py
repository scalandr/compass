import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import xarray
from netCDF4 import Dataset
from numpy import array

#plt.figure(1, figsize=(12.0, 6.0))
nx = 600
ny = 462
plottime = 12
nCells = 277200
nVertLevels = 68

ds = xarray.open_dataset('output.nc')
ds1 = xarray.open_dataset('init.nc')

ds = ds.drop_vars(np.setdiff1d([j for j in ds.variables],
                                  ['velocityZonal', 'layerThickness', 'velZonTimesThick', 'yCell', 'nCells']))

x = np.zeros((nCells, nVertLevels))
y = np.zeros((nCells, nVertLevels))
z = np.zeros((nCells, nVertLevels))
count = 0
#for plottime in range (1,25):
count = count + 1
xmean = ds.isel(Time=plottime)#.groupby('yCell').mean(dim='nCells')
x = x + xmean.velocityZonal.values
print(np.shape(x))
y = y + xmean.layerThickness.values
print(np.shape(y))
z = z + xmean.velZonTimesThick.values
print(count)
x = x/24
y = y/24
z = z/24
print(np.shape(x*y))
print(np.shape(z))
streamFuncHighRes = z-(x*y)

#ymean = ds.isel(Time=plottime)
#print(ymean.dims)
yCell = ds1.yCell
refZMid = ds1.refZMid
refZMidGrid, yCellGrid = np.meshgrid(refZMid, yCell)
mask = np.where(np.logical_and(refZMidGrid > -30, np.logical_and(yCellGrid < 125e3,
                                       yCellGrid > 75e3)),True,False)
print('mask', np.shape(mask))
#print('mask type', type(mask))
print('refZMidGrid', np.shape(refZMidGrid))
print('yCellGrid', np.shape(yCellGrid))
print(np.max(yCellGrid[mask]))
print(np.min(yCellGrid[mask]))
print(np.max(refZMidGrid[mask]))
print(np.min(refZMidGrid[mask]))
#y = ymean.streamFuncSubmeso.values
y = streamFuncHighRes
print('y', np.shape(y))
y = y[mask]
print('y mask', np.shape(y))
print('mean of y', np.mean(y))
