import shlex
import numpy as np
# TODO: Find sdds use case with 'short' type
# SDDS defaults to 32 bit unsigned longs on all test systems while numpy uses 64 bits for the np.int_ in test cases
#  therefore int datatypes are hardcoded in size
data_types = {'double': np.float64, 'short': np.int16, 'long': np.int32, 'string': 'S{}', 'char': np.char}


def _return_type(string):
    target = 'type='
    type_field_position = string.find(target)
    if type_field_position < 0:
        return None
    for key in data_types.keys():
        if string[type_field_position:type_field_position+len(key)].find(key) > -1:
            return key
    return None


def _read_line(open_file):
    line = ''
    while line == '\n' or line == '' or line[0] == '!':
        line = str(open_file.readline(), 'latin-1')
    return line


def _shlex_split(line):
    # Cannot just use shlex.split
    # Need to ignore single-quotes due to their use as common characters in accelerator notation
    lex = shlex.shlex(line, posix=True)
    lex.quotes = '"'
    lex.whitespace_split = True
    lex.commenters = ''

    return list(lex)


def iter_always():
    i = 0
    while True:
        yield i
        i += 1