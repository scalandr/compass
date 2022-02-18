import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

fig = plt.gcf()
fig.set_size_inches(8.0,10.0)
#iTime=[3,6]
#time=['3 hrs','6 hrs']

initfile = Dataset('../../../hydro/default/hydro/init.nc','r')
ncfileH = Dataset('../../../hydro/default/hydro/output.nc','r')
ncfileNH = Dataset('../../../nonhydro/default/nonhydro/output.nc','r')
normalVelocityH = ncfileH.variables['normalVelocity'] 
vertAleTransportTopH = ncfileH.variables['vertAleTransportTop']
zMidH = ncfileH.variables['zMid']  
normalVelocityNH = ncfileNH.variables['normalVelocity'] 
vertAleTransportTopNH = ncfileNH.variables['vertAleTransportTop']
zMidNH = ncfileNH.variables['zMid'] 
cellsOnEdge = initfile.variables['cellsOnEdge']
edgesOnCell = initfile.variables['edgesOnCell']

#horizontall velocity
zMidEdge = 0.5*(zMidH[12, 31, :] + zMidH[12, 32, :])
zMidEdge1 = zMidEdge/16
print(np.shape(zMidEdge1))
for i in range(0,6):
    iEdge = edgesOnCell[31,i] - 1
    for j in range(0,6):
        jEdge  = edgesOnCell[32,j] - 1
        if (iEdge == jEdge):
            midEdge = iEdge
normalVelocity1 = normalVelocityH[12,midEdge,:]/max(normalVelocityH[12,midEdge,:])
print(np.shape(normalVelocity1))
zMidEdge = 0.5*(zMidNH[12, 31, :] + zMidNH[12, 32, :])
zMidEdge2 = zMidEdge/16
print(np.shape(zMidEdge2))
for i in range(0,6):
    iEdge = edgesOnCell[31,i] - 1
    for j in range(0,6):
        jEdge  = edgesOnCell[32,j] - 1
        if (iEdge == jEdge):
            midEdge = iEdge
normalVelocity2 = normalVelocityNH[12,midEdge,:]/max(normalVelocityNH[12,midEdge,:])
print(np.shape(normalVelocity2))

#vertical velocity
zMid_origin1 = zMidH[12, 0, :]/16
print(np.shape(zMid_origin1))
vertAleTransportTop_origin1 = vertAleTransportTopH[12, 0, 0:40]/max(abs(vertAleTransportTopH[12, 0, 0:40]))
print(np.shape(vertAleTransportTop_origin1))
zMid_origin2 = zMidNH[12, 0, :]/16
print(np.shape(zMid_origin2))
vertAleTransportTop_origin2 = vertAleTransportTopNH[12, 0, 0:40]/max(abs(vertAleTransportTopNH[12, 0, 0:40]))
print(np.shape(vertAleTransportTop_origin2))

#plots
plt.figure(figsize=(8.4,4.2))
plt.subplot(1, 2, 1) 
plt.plot(normalVelocity1, zMidEdge1, 'r')
plt.plot(normalVelocity2, zMidEdge2, 'b')
plt.plot(normalVelocity2, zMidEdge2, '--', color='black')
plt.xlabel('u/u_max')
plt.ylabel('z/H')
plt.yticks([0, -0.2, -0.4, -0.6, -0.8, -1])
plt.title('Stratified Wave - hor profile')

plt.subplot(1, 2, 2)
plt.plot(vertAleTransportTop_origin1, zMid_origin1, 'r', label='H model')
plt.plot(vertAleTransportTop_origin2, zMid_origin2, 'b', label='NH model')
plt.plot(vertAleTransportTop_origin2, zMid_origin2, '--', color='black', label='eigenfunction analysis')
plt.xlim([-1.1, 1.1])
plt.xlabel('w/w_max')
plt.legend()
plt.title('Stratified Wave - vert profile')

ncfileH.close()
ncfileNH.close()
initfile.close()
plt.savefig('plotVertAndHor.png')

