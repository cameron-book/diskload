import diskload.earth as earth
from diskload.constants import *
from math import sin, cos, tan, pi
import numba

@numba.jit(nopython = True)
def point( theta, w, nmax, love_numbers, earth_model=earth.default ):
    if nmax > len(love_numbers[0]):
        raise ValueError("Too few Love numbers available for the computation")

    h = love_numbers[0]
    l = love_numbers[1]
    k = love_numbers[2]

    if l[0] != 0:
        raise ValueError("Love number l_0 is not zero.")
    if k[0] != 0:
        raise ValueError("Love number k_0 is not zero.")

    alpha = 0.0
    x = cos(theta * DEGREES)
    sintheta = sin(theta * DEGREES)  
    cosalpha = 1.0

    p0 = 1.0
    p1 = x
 
    pp0 = 1.0
    pp1 = cosalpha
    pp2 = 0.0

    coreValue = 1.0 / (2*sin(theta*DEGREES/2.0))
    # BADBAD: What whould this be?
    coreDerivative = 1
  
    sigma = 0.0
    coreFactor = 2.0

    u = (h[nmax] / (2.0*sin(theta*DEGREES/2.0)) + (h[0] - h[nmax]) * p0)/coreFactor
    l_oo = (nmax + 1) * l[nmax]
    v = - sintheta * l_oo * coreDerivative / coreFactor
    g = coreValue / coreFactor
   
    # keeping track of low-order bits
    Uc = 0.0
    Vc = 0.0
    Gc = 0.0
  
    for n in range(1,nmax+1):
        # Recursively compute the loading factor
        pp2 = ((2*n+1)*cosalpha*pp1 - n*pp0)/(n+1)

        sigma = (pp0 - pp2)/(2*n+1)/coreFactor
    
        # Kahan summation for U
        z = ((h[n] - h[nmax]) * p1 / coreFactor) - Uc
        t = u + z
        Uc = (t - u) - z
        u = t

        # Compute the derivative of P_n(cos(x)) with respect to x
        dp1=(-n/sintheta)*(p0 - x*p1)

        # Kahan summation for V
        z = ((l[n] - l_oo/(n+1)) * sigma * dp1) - Vc
        t = v + z
        Vc = (t - v) - z
        v = t
        
        # Kahan summation for the geoid term
        z = (k[n] * sigma * p1) - Gc
        t = g + z
        Vc = (t - g) - z
        g = t

        pp0 = pp1
        pp1 = pp2
    
        # Recursively compute the Legendre polynomial values
        p2 = ((2*n+1)*x*p1 - n*p0)/(n+1)
        p0 = p1
        p1 = p2

    radiusMeters = earth_model.radius * 1000
    rhoWater = 1000; # density of pure water(kg/m^3)
    rhoEarth = 3.0*earth_model.gravity/4.0/NEWTONS_G/pi/radiusMeters; # Average Earth density in (kg/m^3)
    metersToMillimeters = 1000;
    outputScale = (3*rhoWater/rhoEarth) * w * metersToMillimeters;

    u = u * outputScale
    v = v * outputScale
    g = g * outputScale

    return u, v, g

