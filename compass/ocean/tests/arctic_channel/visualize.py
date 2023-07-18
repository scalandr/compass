import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from compass.ocean.tests.arctic_channel import horiz
from compass.step import Step

matplotlib.use('Agg')


class Visualize(Step):
    """
    A step for visualizing a cross-section through the density
    in the arctic channel test group
    """
    def __init__(self, test_case):
        """
        Create the step

        Parameters
        ----------
        test_case : compass.TestCase
            The test case this step belongs to
        """
        super().__init__(test_case=test_case, name='visualize')

        self.add_input_file(
            filename='output.nc',
            target='../forward/output.nc')
        self.add_input_file(
            filename='init.nc',
            target='../initial_state/initial_state.nc')
        self.add_output_file('plotTemp.png')

    def run(self):
        """
        Run this step of the test case
        """
        config = self.config

        section = config['visualize']
        # L0 = section.getfloat('L0')
        # a0 = section.getfloat('a0')
        time = section.getint('plotTime')

        ncfileIC = xr.open_dataset('init.nc')
        ds = xr.open_dataset('output.nc')

        # get grid variables
        nCells = ds.sizes['nCells']
        nVertLevels = ds.sizes['nVertLevels']
        section = config['horizontal_grid']
        nx = section.getint('nx')
        yCell = ncfileIC.variables['yCell'].values
        zMid = ds.variables['zMid'][time, 0, :]
        y = yCell[range(int(nx / 2), nCells, nx)] / 1e3
        z = zMid  # / a0

        # Prep variables for cell quantities
        xEdge = ncfileIC.xEdge
        yEdge = ncfileIC.yEdge
        xEdge_mid = np.median(xEdge)
        edgeMask_x = np.equal(xEdge, xEdge_mid)
        # cellsOnEdge = ncfileIC.cellsOnEdge
        # cellsOnEdge_x = cellsOnEdge[edgeMask_x, :]
        # cellIndex = np.subtract(cellsOnEdge_x[1:, 0], 1)
        yEdge_x = yEdge[edgeMask_x]
        layerThickness = ds.layerThickness.isel(Time=time)

        zInterface = np.zeros((nCells, nVertLevels + 1))
        zInterface[:, 0] = 0.  # TODO ds.ssh.values
        for zIndex in range(nVertLevels):
            thickness = layerThickness.isel(nVertLevels=zIndex)
            thickness = thickness.fillna(0.)
            zInterface[:, zIndex + 1] = \
                zInterface[:, zIndex] - thickness.values
        zInterfaces_mesh, yCells_mesh = np.meshgrid(zInterface[0, :],
                                                    yEdge_x)

        # Import cell quantities
        # temperature = ds.temperature.isel(Time=time)
        # temperature_x = temperature[cellIndex, :]
        salinity_i = ds.salinity.isel(Time=0)
        salinity = ds.salinity.isel(Time=time)
        # salinity_x = salinity[cellIndex, :]
        # w = ds.vertVelocityTop.isel(Time=time)
        # w_x = w[cellIndex, 1:]

        zlim = [-80, 0]
        # temp_levels = np.arange(-3, -1, 0.25)
        # density_levels = np.arange(23, 27, 0.25)
        sa_levels = np.arange(24, 32, 0.5)

        # Density
        # var = ds.variables['density'][time, range(int(nx/2), nCells, nx), :]
        # var = var - 1000
        # plt.contour(y, z, var.T, levels=density_levels, cmap='jet')
        # plt.xlabel('y (m)')
        # plt.ylabel('z (m)')
        # plt.colorbar(shrink=0.7)
        # plt.savefig(f'density_{time}.png')
        # plt.close()

        # Temperature
        # var = temperature[time, range(int(nx/2), nCells, nx), :]
        # plt.contour(y, z, var.T, levels=temp_levels, cmap='jet')
        # plt.xlabel('y (m)')
        # plt.ylabel('z (m)')
        # plt.colorbar(shrink=0.7)
        # plt.savefig(f'temperature_{time}.png')
        # plt.close()

        # Salinity
        plot_mean = False
        plt.figure()
        var = np.zeros((len(y), len(z)))
        var_i = np.zeros((len(y), len(z)))
        y_unique = np.unique(yCell)
        if plot_mean:
            for i in range(len(y_unique)):
                y_index = yCell == y_unique[i]
                var[i, :] = salinity[y_index, :].mean(axis=0)
                var_i[i, :] = salinity_i[y_index, :].mean(axis=0)
        else:
            var = salinity[range(int(nx / 2), nCells, nx), :]
        plt.contour(y, z, var_i.T, linestyles='--', levels=sa_levels,
                    cmap='jet')
        c1 = plt.contour(y, z, var.T, levels=sa_levels, cmap='jet')
        plt.ylim(zlim)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel('y (km)', fontsize=16)
        plt.ylabel('z (m)', fontsize=16)
        cbar = plt.colorbar(c1, shrink=0.7)
        cbar.ax.tick_params(labelsize=14)
        cbar.ax.set_ylabel('salinity (PSU)', fontsize=16)
        plt.tight_layout()
        plt.savefig(f'salinity_{time}.png')
        plt.close()

        # plt.figure()
        # plt.pcolormesh(np.divide(yCells_mesh, 1e3),
        #                zInterfaces_mesh,
        #                salinity_x.values, cmap='viridis')
        # plt.xticks(fontsize=14)
        # plt.yticks(fontsize=14)
        # plt.xlabel('y (km)', fontsize=16)
        # plt.ylabel('z (m)', fontsize=16)
        # cbar = plt.colorbar()
        # cbar.ax.set_title('salinity (PSU)', fontsize=16)
        # plt.savefig(f'sa_depth_section_t{time}.png',
        #             bbox_inches='tight', dpi=200)
        # plt.close()

        # # Vertical velocity
        # plt.figure()
        # plt.pcolormesh(np.divide(yCells_mesh, 1e3),
        #                zInterfaces_mesh,
        #                w_x.values, cmap='coolwarm', clim=[-1e-3,1e-3])
        # plt.xlabel('y (km)')
        # plt.ylabel('z (m)')
        # cbar = plt.colorbar()
        # cbar.ax.set_title('w (m/s)')
        # plt.savefig(f'w_depth_section_t{time}.png',
        #             bbox_inches='tight', dpi=200)
        # plt.close()

        for var in ['salinity', 'temperature', 'velocityZonal',
                    'velocityMeridional']:
            horiz.plot_horiz_field(ds, ncfileIC, f'{var}',
                                   out_file_name=f'{var}_surface_t{time}',
                                   t_index=time, z_index=0)
        ncfileIC.close()
        ds.close()
