import numpy as np

def add_fill_value(data, fill_value):

    data = np.ma.fix_invalid(data, fill_value=fill_value)
    np.ma.set_fill_value(data, fill_value)
  
    return data
