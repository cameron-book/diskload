import numpy as np
import pathlib 
from os.path import dirname
        
def read(filepath = None):
    if filepath == None:
        filepath = pathlib.Path(dirname(__file__)).parent.parent.parent.joinpath("REF_6371_loading_love_numbers_0_40000.txt")

    exponents = lambda s: float(s.decode('utf8').replace('D','E'))
    converters = {1:exponents,2:exponents,3:exponents}
    dtype = np.dtype([('n', int), ('h', float), ('l', float) , ('k', float)])
    result = np.loadtxt(filepath, dtype=dtype, comments='#', delimiter=None, converters=converters, skiprows=1)

    for i in range(len(result)):
        if result[i][0] != i:
            raise ValueError("The n-th Love numbers must appear on the n-th line in the file.")
        
    result = np.array( [result['h'], result['l'], result['k']] )
    
    return result
