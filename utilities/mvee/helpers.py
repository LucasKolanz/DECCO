import json
import copy
import numpy as np



class NumpyEncoder(json.JSONEncoder):
    ''' Credit to StackOverflow post '''
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32,
                              np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)



def change_dict_to_integer_keys(d):
    '''
    d is a dictionary with strings as keys.

    We want a dictionary with integers as keys. This can happen because
    json stores dictionary keys as strings regardless of what type they
    actually are.

    Modifies d in place
    '''
    # json stores keys as strings -- we need ints
    oldkeys = copy.copy(list(d.keys()))
    for key in oldkeys:
        d[int(key)] = d[key]
        del d[key]

    return d
