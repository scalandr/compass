import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import xarray
from netCDF4 import Dataset
from numpy import array

nx = 600
ny = 462
plottime = 12
nCells = 277200
nVertLevels = 68

ds = xarray.open_dataset('output.nc')

ds = ds.drop_vars(numpy.setdiff1d([j for j in ds.variables],
                                  ['velocityZonal', 'layerThickness', 'velZonTimesThick', 'yCell', 'nCells']))

x = numpy.zeros((nCells, nVertLevels))
y = numpy.zeros((nCells, nVertLevels))
z = numpy.zeros((nCells, nVertLevels))
count = 0
count = count + 1
xmean = ds.isel(Time=plottime)#.groupby('yCell').mean(dim='nCells')
x = x + xmean.velocityZonal.values
print(numpy.shape(x))
y = y + xmean.layerThickness.values
print(numpy.shape(y))
print(count)
count = 0
for time in range (1,plottime+1):
    zmean = ds.isel(Time=time)
    z = z + zmean.velZonTimesThick.values
    count = count + 1
print(count)
z = z/plottime
print(numpy.shape(x*y))
print(numpy.shape(z))
streamFuncHighRes = z-(x*y)

figsize = [6.4, 4.8]
ds = xarray.open_dataset('output.nc')
ds1 = ds.isel(Time=plottime)
ds1 = ds1.sortby('yCell')

nCells = ds1.sizes['nCells']
nEdges = ds1.sizes['nEdges']
nVertLevels = ds1.sizes['nVertLevels']

xCell = numpy.zeros((nCells))
xCell = ds1.xCell
yCell = numpy.zeros((nCells))
yCell = ds1.yCell

#xCell_mid = numpy.median(xCell)
xCell_mid = 150000.0
cellMask_x = numpy.equal(xCell, xCell_mid)

zIndex = xarray.DataArray(data=numpy.arange(nVertLevels),
                                      dims='nVertLevels')
zInterface = numpy.zeros((nCells, nVertLevels + 1))
zInterface[:, 0] = ds1.ssh.values
for zIndex in range(nVertLevels):
    thickness = ds1.layerThickness.isel(nVertLevels=zIndex)
    thickness = thickness.fillna(0.)
    zInterface[:, zIndex + 1] = \
        zInterface[:, zIndex] - thickness.values

zMid = numpy.zeros((nCells, nVertLevels))
for zIndex in range(nVertLevels):
    zMid[:, zIndex] = (zInterface[:, zIndex] +
                       numpy.divide(zInterface[:, zIndex + 1] -
                                    zInterface[:, zIndex], 2.))

#yCells = cellMask_x.sortby('yCell')
yCells = yCell[cellMask_x]

zInterfaces_mesh, yCells_mesh = numpy.meshgrid(zInterface[0, :],
                                                           yCells)

#streamFuncSubmeso = numpy.zeros((nCells, nVertLevels))
streamFuncSubmeso = xarray.DataArray(data=streamFuncHighRes, dims=['nCells', 'nVertLevels'])
streamFuncSubmeso_xmesh = streamFuncSubmeso[cellMask_x, :]
print('streamFuncSubmeso_xmesh', numpy.shape(streamFuncSubmeso_xmesh))
print('zInterfaces_mesh', numpy.shape(zInterfaces_mesh[:,:-1]))
print(numpy.shape(numpy.divide(yCells_mesh, 1e3)[:,:-1]))

# Figures
plt.figure(figsize=figsize, dpi=100)
cmax = numpy.max(numpy.abs(streamFuncSubmeso_xmesh))
plt.pcolormesh(numpy.divide(yCells_mesh, 1e3)[:,:-1],
               zInterfaces_mesh[:,:-1],
               streamFuncSubmeso_xmesh.values,
               cmap='RdBu', vmin=-1. * cmax, vmax=cmax)
plt.xlabel('y (km)')
plt.ylabel('z (m)')
cbar = plt.colorbar()
#plt.clim([-0.5,0.5])
#cbar.ax.set_title('uNormal (m/s)')
#plt.savefig('uNormal_depth_section_t{}.png'.format(j),
#            bbox_inches='tight', dpi=200)
plt.savefig('streamFunc_section.png')
plt.close()
