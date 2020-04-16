import numpy as np


class StructData:
    def __init__(self, data_type, max_string_length):
        self.data = None
        self._data = []
        self.max_string_length = max_string_length
        self.data_type = data_type
        self._oldarray = False

    @property
    def data_type(self):
        return self._data_type

    @data_type.setter
    def data_type(self, data_type):
        if type(data_type[0]) == list:
            self._data_type = []
            merge_hold = [b for a in data_type for b in a]
            for descr in merge_hold:
                if descr[0] == 'record_length':
                    continue
                elif type(descr[1]) == str:
                    self._data_type.append((descr[0], 'U{}'.format(self.max_string_length)))
                else:
                    self._data_type.append(descr)
        else:
            self._data_type = data_type

    def add(self, data):
        data_hold = None
        # Stack page data data together
        for datum in data:
            datum = self._merge(datum)

            if data_hold is None:
                data_hold = datum
            else:
                data_hold = np.concatenate([data_hold, datum])

        self._data.append(data_hold)

    def concat(self):
        self.data = np.array(self._data)

    def _merge(self, data):
        if len(data) == 1:
            return data[0]
        else:
            new_array = np.empty(1, dtype=self.data_type)
            for arr in data:
                names = [n for n in arr.dtype.names if n != 'record_length']
                if names:
                    new_array[names] = arr[names]
            return new_array

    def _check_type(self, data): #Next: Need to look at adjusting string size each time up to maximum? Probably just don't merge'
        if data.dtype != np.dtype(self.data_type):
            raise TypeError("Array Datatypes do not match")
