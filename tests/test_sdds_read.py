import unittest
import pandas as pd
from io import StringIO
import numpy as np
from subprocess import Popen, PIPE
from rsbeams.rsdata.SDDS import readSDDS, supported_namelists
from headers import header1, header2, header3, header4


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


# TODO: This test only verifies that _parse_header runs, it does not check the output
class TestHeaderParse(unittest.TestCase):
    # Skip read header method and extract directly into a list since we know our example structure
    def setUp(self):
        with open('header1.sdds', 'w') as f1:
            f1.write(header1)
        with open('header2.sdds', 'w') as f2:
            f2.write(header2)
        with open('header3.sdds', 'w') as f3:
            f3.write(header3)
        with open('header4.sdds', 'w') as f4:
            f4.write(header4)

    def test_header1(self):
        reader = object.__new__(readSDDS)
        reader.data = {key: [] for key in supported_namelists.keys()}
        reader.header = []
        reader.buffer = False
        with open('header1.sdds', 'r') as f1:
            f1.readline()
            for file_line in f1.readlines()[2:]:
                reader.header.append(file_line)
        reader._parse_header()
        for key, val in reader.data.items():
            for insets in val:
                print("Raw Header: {}".format(insets.namelist))
                print("Namelist type: {}".format(key))
                for k,v in insets.fields.items():
                    print("{}: {}".format(k, v))

    def test_header2(self):
        reader = object.__new__(readSDDS)
        reader.data = {key: [] for key in supported_namelists.keys()}
        reader.header = []
        reader.buffer = False
        with open('header2.sdds', 'r') as f2:
            f2.readline()
            for file_line in f2.readlines()[2:]:
                reader.header.append(file_line)
        reader._parse_header()
        for key, val in reader.data.items():
            for insets in val:
                print("Raw Header: {}".format(insets.namelist))
                print("Namelist type: {}".format(key))
                for k,v in insets.fields.items():
                    print("{}: {}".format(k, v))

    def test_header4(self):
        reader = object.__new__(readSDDS)
        reader.data = {key: [] for key in supported_namelists.keys()}
        reader.header = []
        reader.buffer = False
        with open('header4.sdds', 'r') as f4:
            f4.readline()
            for file_line in f4.readlines()[2:]:
                reader.header.append(file_line)
        reader._parse_header()
        for key, val in reader.data.items():
            for insets in val:
                print("Raw Header: {}".format(insets.namelist))
                print("Namelist type: {}".format(key))
                for k,v in insets.fields.items():
                    print("{}: {}".format(k, v))

    def tearDown(self):
        import os
        os.remove('header1.sdds')
        os.remove('header2.sdds')
        os.remove('header3.sdds')
        os.remove('header4.sdds')


# TODO: Generate verfication data when compose form is finalized
class TestColumnCompose(unittest.TestCase):
    def setUp(self):
        with open('header1.sdds', 'w') as f1:
            f1.write(header1)
        with open('header2.sdds', 'w') as f2:
            f2.write(header2)
        with open('header3.sdds', 'w') as f3:
            f3.write(header3)
        with open('header4.sdds', 'w') as f4:
            f4.write(header4)

    def test_header1(self):
        reader = object.__new__(readSDDS)
        reader.openf = open('header1.sdds', 'rb')
        reader.data = {key: [] for key in supported_namelists.keys()}
        reader.header = []
        reader._parameter_keys = []
        reader._column_keys = []
        reader.buffer = False
        reader._read_header()
        reader._parse_header()

        reader._compose_column_datatypes()
        print()
        for dt in reader._column_keys:
            print(dt)

    def test_header2(self):
        reader = object.__new__(readSDDS)
        reader.openf = open('header2.sdds', 'rb')
        reader.data = {key: [] for key in supported_namelists.keys()}
        reader.header = []
        reader._parameter_keys = []
        reader._column_keys = []
        reader._read_header()
        reader._parse_header()

        reader._compose_column_datatypes()
        print()
        for dt in reader._column_keys:
            print(dt)

    def tearDown(self):
        import os
        os.remove('header1.sdds')
        os.remove('header2.sdds')
        os.remove('header3.sdds')
        os.remove('header4.sdds')


# TODO: Generate verfication data when compose form is finalized
class TestParameterCompose(unittest.TestCase):
    def setUp(self):
        with open('header1.sdds', 'w') as f1:
            f1.write(header1)
        with open('header2.sdds', 'w') as f2:
            f2.write(header2)
        with open('header3.sdds', 'w') as f3:
            f3.write(header3)
        with open('header4.sdds', 'w') as f4:
            f4.write(header4)

    def test_header1(self):
        reader = object.__new__(readSDDS)
        reader.openf = open('header1.sdds', 'rb')
        reader.data = {key: [] for key in supported_namelists.keys()}
        reader.header = []
        reader.buffer = False
        reader._parameter_keys = []
        reader._column_keys = []
        reader._read_header()
        reader._parse_header()

        reader._compose_parameter_datatypes()
        print()
        for dt in reader._parameter_keys:
            print(dt)

    def test_header4(self):
        reader = object.__new__(readSDDS)
        reader.openf = open('header4.sdds', 'rb')
        reader.data = {key: [] for key in supported_namelists.keys()}
        reader.header = []
        reader._parameter_keys = []
        reader._column_keys = []
        reader._read_header()
        reader._parse_header()

        reader._compose_parameter_datatypes()
        print()
        for dt in reader._parameter_keys:
            print(dt)

    def tearDown(self):
        import os
        os.remove('header1.sdds')
        os.remove('header2.sdds')
        os.remove('header3.sdds')
        os.remove('header4.sdds')


class TestFirstPageParameterRead(unittest.TestCase):
    # Just a test of a single `_get_parameter_data` call
    # Does not test multipage or position calculation dependent read to avoid
    # Additional dependencies

    def test1(self):
        filename = 'file1.sdds'
        reader = readSDDS(filename, buffer=False)
        # Perform set from `read`
        position = 0
        arr, pos = reader._get_parameter_data(reader._parameter_keys, position)
        print(arr)

    def test2(self):
        filename = 'elegant_final.fin'
        use_buffer = False
        reader = readSDDS(filename, buffer=use_buffer)
        reader._compose_parameter_datatypes()
        print(reader.header_end_pointer)
        # Perform set from `read`
        if use_buffer:
            position = reader.header_end_pointer
        else:
            position = 0
        print(position)
        arr, pos = reader._get_parameter_data(reader._parameter_keys, position)
        print(arr)

    def test3(self):
        filename = 'elegant_final.fin'
        use_buffer = True
        reader = readSDDS(filename, buffer=use_buffer)
        reader._compose_parameter_datatypes()
        print(reader.header_end_pointer)
        # Perform set from `read`
        if use_buffer:
            position = reader.header_end_pointer
        else:
            position = 0
        print(position)
        arr, pos = reader._get_parameter_data(reader._parameter_keys, position)
        print(arr)


class TestParameterRead(unittest.TestCase):
    # Just a test of a single `_get_parameter_data` call
    # Does not test multipage or position calculation dependent read to avoid
    # Additional dependencies

    def test1(self):
        filename = 'elegant_final.fin'
        use_buffer = False
        reader = readSDDS(filename, buffer=use_buffer)
        reader.read2()
        print(reader.parameters)

    def test2(self):
        filename = 'elegant_final_ascii.fin'
        use_buffer = False
        reader = readSDDS(filename, buffer=use_buffer)
        reader.read2()

    def test3(self):
        filename = 'elegant_final.fin'
        use_buffer = False
        reader = readSDDS(filename, buffer=use_buffer)
        reader.read2(pages=[5, 10])
        print(reader.parameters)

    def test4(self):
        filename = 'elegant_final.fin'
        use_buffer = True
        reader = readSDDS(filename, buffer=use_buffer)
        reader.read2()
        print(reader.parameters)

    def test5(self):
        filename = 'elegant_final_ascii.fin'
        use_buffer = True
        reader = readSDDS(filename, buffer=use_buffer)
        reader.read2()
        print(reader.parameters)


class TestColumnRead(unittest.TestCase):

    def setUp(self):
        self.tests = {
            'test1': {'filename': 'linear_dipole_20_nopar_nostr.sig',
                      'columns': [],
                      'data': None},
            'test2': {'filename': 'linear_dipole_20_nopar.sig',
                      'columns': []},
            'test3': {'filename': 'exit.w0',
                      'columns': []}
        }
        use_buffer = True
        for testname, test in self.tests.items():
            filename = test['filename']
            reader = readSDDS(filename, buffer=use_buffer)
            for col in reader.data['&column']:
                test['columns'].append(col.fields['name'])
            run_string = 'sdds2stream  -delimiter=" " {filename} -col=' + ','.join(test['columns'])
            data_pipe = Popen(run_string.format(filename=filename), shell=True, stdout=PIPE, stderr=PIPE)
            cols, err = data_pipe.communicate()
            dat = pd.read_csv(StringIO(cols.decode().lstrip()), delim_whitespace=True, names=test['columns'])
            test['data'] = dat

    def test1(self):
        test = self.tests['test1']
        filename = test['filename']
        use_buffer = True
        reader = readSDDS(filename, buffer=use_buffer)
        reader.read2(pages=[0])
        for name in test['columns']:
            self.assertTrue(np.all(np.isclose(reader.columns[name].squeeze(),
                                              test['data'][name],
                                              rtol=1e-6)))

    def test2(self):
        test = self.tests['test2']
        filename = test['filename']
        use_buffer = True
        reader = readSDDS(filename, buffer=use_buffer)
        reader.read2(pages=[0])
        for name in test['columns']:
            print('DSFDF', reader.columns[name].shape)
            self.assertTrue(np.all(np.isclose(reader.columns[name].squeeze(),
                                              test['data'][name],
                                              rtol=1e-6)))

    def test3(self):
        test = self.tests['test3']
        filename = test['filename']
        page_length = test['data'].count()[0] // 20
        cols = [0, 1]
        use_buffer = True
        reader = readSDDS(filename, buffer=use_buffer)
        reader.read2(pages=cols)
        for name in test['columns']:
            for col in cols:
                self.assertTrue(np.all(np.isclose(reader.columns[name][:, col].squeeze(),
                                                  test['data'][name][col*page_length:(col+1)*page_length],
                                                  rtol=1e-6)))


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