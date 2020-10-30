import unittest
import os
import numpy as np
from scipy.constants import c
from rsbeams.rsdata import switchyard
from rsbeams.rsdata.SDDS import readSDDS

_ELEGANT_READ_FILE = 'bunch_5001.sdds'
_OPAL_READ_FILE = 'test_resources/opal.h5'

class TestReaders(unittest.TestCase):
    def setUp(self):
        pass

    def test_read_elegant(self):
        obj = switchyard.Switchyard()
        obj.read(_ELEGANT_READ_FILE, 'elegant')

        expected = readSDDS('bunch_5001.sdds')
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

    def test_read_opal(self):
        obj = switchyard.Switchyard()
        obj.read(_OPAL_READ_FILE, 'opal')

    def test_read_bad_format(self):
        # File exists but format name does not
        obj = switchyard.Switchyard()
        self.assertRaises(LookupError, obj.read, _ELEGANT_READ_FILE, 'not_supported_format')

class TestWriters(unittest.TestCase):
    def setUp(self):
        # Filenames to be created
        self._FILES = {
            'ELEGANT_TO_GENESIS': 'write1.in',
            'OPAL_TO_ELEGANT': 'write2.sdds'
        }

    def test_opal_to_elegant(self):
        obj = switchyard.Switchyard()
        obj.read(_OPAL_READ_FILE, 'opal')
        obj.write(self._FILES['OPAL_TO_ELEGANT'], 'elegant')

    def test_elegant_to_genesis(self):
        obj = switchyard.Switchyard()
        obj.read(_ELEGANT_READ_FILE, 'elegant')
        obj.write(self._FILES['ELEGANT_TO_GENESIS'], 'genesis')


    def tearDown(self):
        for filename in self._FILES.values():
            try:
                os.remove(filename)
            except OSError:
                pass
