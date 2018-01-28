import numpy as np

def extrapolate( love, nmax ):
    oo = love.shape[1] - 1
    h, l, k = love[0], love[1], love[2]
    
    h_oo = h[oo]
    l_oo = l[oo]
    k_oo = k[oo]

    more_h = np.ones( nmax + 1 - love.shape[1] ) * h_oo
    more_l = np.ones( nmax + 1 - love.shape[1] ) * l_oo * oo / np.arange( love.shape[1], nmax + 1 )
    more_k = np.ones( nmax + 1 - love.shape[1] ) * k_oo * oo / np.arange( love.shape[1], nmax + 1 )    

    rest = np.array( [more_h, more_l, more_k] )
    return np.concatenate((love,rest), axis=1)
