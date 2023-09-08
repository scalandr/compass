import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import xarray
from netCDF4 import Dataset
from numpy import array

ds = xarray.open_dataset('output.nc')
figsize = [6.4, 4.8]
ds1 = ds.isel(Time=12)
ds1 = ds1.sortby('yEdge')
#ds1 = ds1.groupby('yEdge')

nCells = ds1.sizes['nCells']
nEdges = ds1.sizes['nEdges']
nVertLevels = ds1.sizes['nVertLevels']

xEdge = numpy.zeros((nEdges))
xEdge = ds1.xEdge
yCell = numpy.zeros((nCells))
yCell = ds1.yCell

xEdge_mid = numpy.median(xEdge)
edgeMask_x = numpy.equal(xEdge, xEdge_mid)

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

cellsOnEdge = ds1.cellsOnEdge
cellsOnEdge_x = cellsOnEdge[edgeMask_x, :]
yEdges = numpy.zeros((len(cellsOnEdge_x) + 1))
for i in range(len(cellsOnEdge_x)):
    if cellsOnEdge[i, 1] == 0:
        yEdges[i] = yCell[cellsOnEdge_x[i, 0] - 1]
        yEdges[i + 1] = yCell[cellsOnEdge_x[i, 0] - 1]
    elif cellsOnEdge[i, 1] == 0:
        yEdges[i] = yCell[cellsOnEdge_x[i, 1] - 1]
        yEdges[i + 1] = yCell[cellsOnEdge_x[i, 1] - 1]
    else:
        yEdges[i] = min(yCell[cellsOnEdge_x[i, 0] - 1],
                        yCell[cellsOnEdge_x[i, 1] - 1])
        yEdges[i + 1] = max(yCell[cellsOnEdge_x[i, 0] - 1],
                            yCell[cellsOnEdge_x[i, 1] - 1])

zInterfaces_mesh, yEdges_mesh = numpy.meshgrid(zInterface[0, :],
                                                           yEdges)

streamFuncSubmeso = numpy.zeros((nCells, nVertLevels))
streamFuncSubmeso = ds1.streamFuncSubmeso
streamFuncSubmeso_xmesh = streamFuncSubmeso[edgeMask_x, :]

# Figures
plt.figure(figsize=figsize, dpi=100)
cmax = numpy.max(numpy.abs(streamFuncSubmeso_xmesh.values))
plt.pcolormesh(numpy.divide(yEdges_mesh, 1e3),
               zInterfaces_mesh,
               streamFuncSubmeso_xmesh.values,
               cmap='RdBu', vmin=-1. * cmax, vmax=cmax)
plt.xlabel('y (km)')
plt.ylabel('z (m)')
cbar = plt.colorbar()
#cbar.ax.set_title('uNormal (m/s)')
#plt.savefig('uNormal_depth_section_t{}.png'.format(j),
#            bbox_inches='tight', dpi=200)
plt.savefig('streamFunc_section.png')
plt.close()
