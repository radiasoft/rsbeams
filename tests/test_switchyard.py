import unittest
import os
import pathlib
import numpy as np
from scipy.constants import c
from rsbeams.rsdata import switchyard
from rsbeams.rsdata.SDDS import readSDDS

_ELEGANT_READ_FILE = 'test_resources/bunch_5001.sdds'
_OPAL_READ_FILE = 'test_resources/opal.h5'
_OPAL_MONITOR_FILE = 'test_resources/output_x_108.h5'


class TestReaders(unittest.TestCase):
    def setUp(self):
        pass

    def test_read_elegant(self):
        obj = switchyard.Switchyard()
        obj.read(_ELEGANT_READ_FILE, 'elegant')

        expected = readSDDS(_ELEGANT_READ_FILE)
        expected.read()

        coordinates = expected.columns.squeeze()
        coordinates['xp'] = coordinates['xp'] * coordinates['p']
        coordinates['yp'] = coordinates['yp'] * coordinates['p']
        coordinates['t'] = coordinates['t'] * c

        species_name = switchyard._DEFAULT_SPECIES_NAME
        self.assertTrue(np.all(np.isclose(coordinates['x'], obj.species[species_name].x, atol=1e-14)))
        self.assertTrue(np.all(np.isclose(coordinates['xp'], obj.species[species_name].ux, atol=1e-14)))
        self.assertTrue(np.all(np.isclose(coordinates['y'], obj.species[species_name].y, atol=1e-14)))
        self.assertTrue(np.all(np.isclose(coordinates['yp'], obj.species[species_name].uy, atol=1e-14)))
        self.assertTrue(np.all(np.isclose(coordinates['t'], obj.species[species_name].ct, atol=1e-14)))
        self.assertTrue(np.all(np.isclose(coordinates['p'], obj.species[species_name].pt, atol=1e-14)))

    def test_read_opal_phase_space(self):
        obj = switchyard.Switchyard()
        obj.read(_OPAL_READ_FILE, 'opal')

    def test_read_opal_monitor(self):
        obj = switchyard.Switchyard()
        obj.read(_OPAL_MONITOR_FILE, 'opal')

    def test_read_bad_format(self):
        # File exists but format name does not
        obj = switchyard.Switchyard()
        self.assertRaises(LookupError, obj.read, _ELEGANT_READ_FILE, 'not_supported_format')


class TestWriters(unittest.TestCase):
    def setUp(self):
        # Filenames to be created
        self._FILES = []

    def test_opal_to_elegant(self):
        fn = 'opal2elegant.sdds'
        obj = switchyard.Switchyard()
        obj.read(_OPAL_READ_FILE, 'opal')
        obj.write(fn, 'elegant')
        path = pathlib.Path(fn)

        self.assertTrue(path.is_file())
        self._FILES.append(fn)

    def test_elegant_to_genesis(self):
        fn = 'elegant2genesis.in'
        obj = switchyard.Switchyard()
        obj.read(_ELEGANT_READ_FILE, 'elegant')
        obj.write(fn, 'genesis')
        path = pathlib.Path(fn)

        self.assertTrue(path.is_file())
        self._FILES.append(fn)

    def test_opal_write(self):
        fn = 'opal_write.txt'
        obj = switchyard.Switchyard()
        obj.read(_ELEGANT_READ_FILE, 'elegant')
        obj.write(fn, 'opal')
        path = pathlib.Path(fn)

        self.assertTrue(path.is_file())
        self._FILES.append(fn)

    def tearDown(self):
        for filename in self._FILES:
            try:
                os.remove(filename)
            except OSError:
                pass
