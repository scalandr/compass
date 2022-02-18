import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

fig = plt.gcf()
fig.set_size_inches(8.0,10.0)
#iTime=[3,6]
#time=['3 hrs','6 hrs']

initfile = Dataset('../forward/init.nc','r')
ncfile = Dataset('../forward/output.nc','r')
normalVelocity = ncfile.variables['normalVelocity'] 
vertAleTransportTop = ncfile.variables['vertAleTransportTop']
zMid = ncfile.variables['zMid']  
cellsOnEdge = initfile.variables['cellsOnEdge']
edgesOnCell = initfile.variables['edgesOnCell']

#horizontall velocity
zMidEdge = 0.5*(zMid[12, 31, :] + zMid[12, 32, :])
zMidEdge1 = zMidEdge/16
print(np.shape(zMidEdge1))
for i in range(0,6):
    iEdge = edgesOnCell[31,i] - 1
    for j in range(0,6):
        jEdge  = edgesOnCell[32,j] - 1
        if (iEdge == jEdge):
            midEdge = iEdge
normalVelocity1 = normalVelocity[12,midEdge,:]/max(normalVelocity[12,midEdge,:])
print(np.shape(normalVelocity1))

#vertical velocity
zMid_origin1 = zMid[12, 0, :]/16
print(np.shape(zMid_origin1))
vertAleTransportTop_origin1 = vertAleTransportTop[12, 0, 0:40]/max(abs(vertAleTransportTop[12, 0, 0:40]))
print(np.shape(vertAleTransportTop_origin1))

#plots
plt.subplot(1, 2, 1) 
plt.plot(normalVelocity1, zMidEdge1)
plt.xlabel('u/u_max')
plt.ylabel('z/H')
plt.title('Stratified Internal Wave')

plt.subplot(1, 2, 2)
plt.plot(vertAleTransportTop_origin1, zMid_origin1)
plt.xlim([-1.1, 1.1])
plt.xlabel('w/w_max')
plt.title('Stratified Internal Wave')

ncfile.close()
initfile.close()
plt.savefig('plotVertAndHor.png')

