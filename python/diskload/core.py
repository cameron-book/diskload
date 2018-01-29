import diskload.earth as earth
from diskload.constants import *
from diskload.compensation import *
from math import sin, cos, tan, pi
import numba

@numba.jit(nopython = True)
def core( alpha, icomp, theta, w, core_value, core_derivative, nmax, love_numbers, earth_model=earth.default, computeU=True, computeV=True, computeG=True ):
    if nmax > len(love_numbers[0]):
        raise ValueError("Too few Love numbers available for the computation")

    h = love_numbers[0]
    l = love_numbers[1]
    k = love_numbers[2]

    if l[0] != 0:
        raise ValueError("Love number l_0 is not zero.")
    if k[0] != 0:
        raise ValueError("Love number k_0 is not zero.")

    x = cos(theta * DEGREES)
    sintheta = sin(theta * DEGREES)
    cosalpha = cos(alpha * DEGREES)

    p0 = 1.0
    p1 = x

    pp0 = 1.0
    pp1 = cosalpha
    pp2 = 0.0
    
    sigma = 0.0
    if icomp == Compensation.UNCOMPENSATED:
        sigma = (1-cosalpha)/2.0

    u = h[0] * sigma * p0
    v = 0.0
    g = sigma

    # keeping track of low-order bits
    Uc = 0.0
    Vc = 0.0
    Gc = 0.0
    
    coreFactor = 2.0
    if icomp == Compensation.COMPENSATED:
        coreFactor = 1.0 + cosalpha

    u = h[nmax] * core_value / coreFactor + (h[0] - h[nmax]) * sigma * p0

    l_oo = (nmax) * l[nmax]
    v = - sintheta * l_oo * core_derivative / coreFactor

    g = core_value / coreFactor
  
    if icomp == Compensation.COMPENSATED:
        # BADBAD: this differs from the paper -- we should be subtracting the first term?
        u -= h[nmax] * (1 - cosalpha) / coreFactor
        g -= (1 - cosalpha) / coreFactor
  
    for n in range(1, nmax+1):
        # Recursively compute the loading factor        
        pp2 = ((2*n+1)*cosalpha*pp1 - n*pp0)/(n+1)
        
        sigma = (pp0 - pp2)/(2*n+1)/coreFactor

        if computeU:
            z = ((h[n] - h[nmax]) * sigma * p1) - Uc         
            t = u + z
            Uc = (t - u) - z
            u = t

        # Compute the derivative of P_n(cos(x)) with respect to x
        dp1=(-n/sintheta)*(p0 - x*p1)

        if computeV:
            # this is /n and not /(n+1)
            z = ((l[n] - l_oo/n) * sigma * dp1) - Vc;            
            t = v + z
            Vc = (t-v) - z
            v = t

        if computeG:
            z = (k[n] * sigma * p1) - Gc;            
            t = g + z
            Vc = (t - g) - z
            g = t
    
        pp0 = pp1;
        pp1 = pp2;
    
        # Recursively compute the Legendre polynomial values
        p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
        p0 = p1;
        p1 = p2;

    radiusMeters = earth_model.radius * 1000
    rhoWater = 1000; # density of pure water(kg/m^3)
    rhoEarth = 3.0*earth_model.gravity/4.0/NEWTONS_G/pi/radiusMeters; # Average Earth density in (kg/m^3)
    metersToMillimeters = 1000;
    outputScale = (3*rhoWater/rhoEarth) * w * metersToMillimeters;

    u = u * outputScale
    v = v * outputScale
    g = g * outputScale

    return u, v, g
