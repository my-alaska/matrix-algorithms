import sys
import numpy as np

def total_size(obj, seen=None):
    """Recursively finds size of object and its children"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([total_size(v, seen) for v in obj.values()])
        size += sum([total_size(k, seen) for k in obj.keys()])
    elif isinstance(obj, np.ndarray):
        size += obj.nbytes
    elif hasattr(obj, '__dict__'):
        size += total_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([total_size(i, seen) for i in obj])
    return size
