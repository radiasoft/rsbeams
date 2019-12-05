import unittest
from rsbeams.rsdata.SDDS import writeSDDS, readSDDS
from subprocess import Popen, PIPE
import numpy as np

# Tests that use SDDS Tools Distribution


def read_file(filename, par=None, col=None):
    call1 = Popen('sddscheck {}'.format(filename), shell=True, stdout=PIPE, stderr=PIPE)
    check, err = call1.communicate()
    if err:
        return None, None, None
    else:
        check = check.decode('utf-8')
    if par:
        call2 = Popen('sdds2stream -par={} {}'.format(par, filename), shell=True, stdout=PIPE, stderr=PIPE)
        par, err = call2.communicate()
        if err:
            return None, None, None
        else:
            par = par.decode('utf-8')
    if col:
        call2 = Popen('sdds2stream -col={} {}'.format(col, filename), shell=True, stdout=PIPE, stderr=PIPE)
        col, err = call2.communicate()
        if err:
            return None, None, None
        else:
            col = col.decode('utf-8')

    return check, par, col


class TestWriteBinary(unittest.TestCase):

    def setUp(self):
        test1 = writeSDDS('file1.sdds')
        test1.create_parameter('par1', 2342.452, 'double')
        test1.create_parameter('par2', 42, 'long')
        test1.create_column('col1', np.array([1.23892e-3, 52452.2]).ravel(), 'double')
        test1.create_column('col2', np.array([135.252, 52452.2e4]).ravel(), 'double', colUnits='m', colSymbol='&n',
                            colDescription="A test description")
        test1.save_sdds('file1.sdds', 'binary')

        self.status, self.par1, self.col1 = read_file('file1.sdds', par='par1', col='col1')
        _, self.par2, self.col2 = read_file('file1.sdds', par='par2', col='col2')

    def test_status_one(self):
        self.assertEqual(self.status, 'ok\n')

    def test_parameter_one(self):
        self.assertAlmostEqual(float(self.par1),  2342.452, delta=1e8)

    def test_parameter_two(self):
        self.assertEqual(float(self.par2),  42)


class TestWriteAscii(unittest.TestCase):

    def setUp(self):
        test1 = writeSDDS('file1.sdds')
        test1.create_paremater('par1', 2342.452, 'double')
        test1.create_paremater('par2', 42, 'long')
        test1.create_column('col1', np.array([1.23892e-3, 52452.2]).ravel(), 'double')
        test1.create_column('col2', np.array([135.252, 52452.2e4]).ravel(), 'double', colUnits='m', colSymbol='&n',
                            colDescription="A test description")
        test1.save_sdds('file1.sdds', 'ascii')

        self.status, self.par1, self.col1 = read_file('file1.sdds', par='par1', col='col1')
        _, self.par2, self.col2 = read_file('file1.sdds', par='par2', col='col2')

    def test_status_one(self):
        self.assertEqual(self.status, 'ok\n')

    def test_parameter_one(self):
        self.assertAlmostEqual(float(self.par1),  2342.452, delta=1e-8)

    def test_parameter_two(self):
        self.assertEqual(float(self.par2),  42)


class TestReadBinary1(unittest.TestCase):
    filename = 'test_read_1_bunch.out'

    def setUp(self):
        test_read = readSDDS(self.filename)
        self.parameters, self.columns = test_read.read()

    def test_length(self):
        self.assertEqual(2, len(self.parameters))
        self.assertEqual(2, self.columns.shape[0])

    def test_parameters(self):
        targets = {'rowCount': 10000,
                   'pCentral': 0.14699293165778182,
                   'Particles': 10000}

        for key, target_value in targets.items():
            for page in self.parameters:
                self.assertAlmostEqual(page[key], target_value, delta=1e-8)


# Test for binary string write out in columns
# class TestStringWriteBinary(unittest.TestCase):
#
#     def setUp(self):
#         test1 = writeSDDS('file1.sdds')
#         test1.create_column('col1', [7, "Chicago".encode(), 8, "New York".encode()], 'string')
#         test1.save_sdds('file1.sdds', 'binary')
#
#         self.status, _, self.col1 = read_file('file1.sdds', col='col1')
#
#     def test_status_one(self):
#         self.assertEqual(self.status, 'ok\n')


if __name__ == '__main__':
    unittest.main()