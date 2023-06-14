from compass.ocean.tests.arctic_channel.forward import Forward
from compass.ocean.tests.arctic_channel.initial_state import InitialState
from compass.ocean.tests.arctic_channel.visualize import Visualize
from compass.testcase import TestCase


class Nonhydro(TestCase):
    """
    The nonhydro test case for the arctic channel group simply creates
    the mesh and initial condition, then performs a forward runs with
    the nonhydrostatic version of MPAS-O, and then plots the density
    profile.

    """

    def __init__(self, test_group):
        """
        Create the test case

        Parameters
        ----------
        test_group : compass.ocean.tests.arctic_channel.ArcticChannel
            The test group that this test case belongs to
        """
        name = 'nonhydro'
        super().__init__(test_group=test_group, name=name)

        self.add_step(
            InitialState(test_case=self))

        step = Forward(
            test_case=self, name='forward', ntasks=16,
            min_tasks=1, openmp_threads=1)
        step.add_namelist_file(
            'compass.ocean.tests.arctic_channel.nonhydro',
            'namelist.forward')
        step.add_streams_file(
            'compass.ocean.tests.arctic_channel.nonhydro',
            'streams.forward')
        self.add_step(step)

        self.add_step(
            Visualize(test_case=self))

    # no run() is needed because we're doing the default: running all steps
