import unittest
from rsbeams.rsdata.SDDS import readSDDS
from headers import header1, header2, header3


class TestHeaderRead(unittest.TestCase):

    def setUp(self):
        with open('header1.sdds', 'w') as f1:
            f1.write(header1)
        with open('header2.sdds', 'w') as f2:
            f2.write(header2)
        with open('header3.sdds', 'w') as f3:
            f3.write(header3)

    def test_header1(self):
        reader = object.__new__(readSDDS)
        reader.openf = open('header1.sdds', 'rb')
        reader.header = []

        reader._read_header()
        with open('header1.sdds', 'r') as f1:
            f1.readline()
            for file_line, sdds_reader_line in zip(f1.readlines()[2:], reader.header):
                self.assertEqual(file_line, sdds_reader_line)

    def test_header2(self):
        reader = object.__new__(readSDDS)
        reader.openf = open('header2.sdds', 'rb')
        reader.header = []

        with open('header2.sdds', 'r') as f2:
            f2.readline()
            for file_line, sdds_reader_line in zip(f2.readlines()[2:], reader.header):
                self.assertEqual(file_line, sdds_reader_line)

    def test_header3(self):
        # TODO: Less easy to test this one is right with many lines per namelist
        reader = object.__new__(readSDDS)
        reader.openf = open('header3.sdds', 'rb')
        reader.header = []

        reader._read_header()
        for it in reader.header:
            print('new namelist')
            print(it)

    def tearDown(self):
        import os
        os.remove('header1.sdds')
        os.remove('header2.sdds')
        os.remove('header3.sdds')

# class TestReadBinary1(unittest.TestCase):
#     filename = 'test_read_1_bunch.out'
#
#     def setUp(self):
#         test_read = readSDDS(self.filename)
#         self.parameters, self.columns = test_read.read()
#
#     def test_length(self):
#         self.assertEqual(2, len(self.parameters))
#         self.assertEqual(2, self.columns.shape[0])
#
#     def test_parameters(self):
#         targets = {'rowCount': 10000,
#                    'pCentral': 0.14699293165778182,
#                    'Particles': 10000}
#
#         for key, target_value in targets.items():
#             for page in self.parameters:
#                 self.assertAlmostEqual(page[key], target_value, delta=1e-8)