from compass.ocean.tests.arctic_channel.hydro import Hydro
from compass.ocean.tests.arctic_channel.nonhydro import Nonhydro
from compass.testgroup import TestGroup


class ArcticChannel(TestGroup):
    """
    A test group for lock exchange type test cases
    """
    def __init__(self, mpas_core):
        """
        mpas_core : compass.MpasCore
            the MPAS core that this test group belongs to
        """
        super().__init__(mpas_core=mpas_core, name='arctic_channel')

        self.add_test_case(
            Hydro(test_group=self))
        self.add_test_case(
            Nonhydro(test_group=self))
