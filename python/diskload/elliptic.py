import diskload.core
from diskload.constants import *
import mpmath
from mpmath import cos, sin, tan, asin, pi
from mpmath import re, sqrt
from mpmath import ellipk, ellippi, ellipf, ellipe
import diskload.earth as earth

mpmath.mp.dps = 22

# This is really G(cos(x),cos(y))
def elliptic_core_g(x,y):
  x = x/2
  y = y/2
  
  factor = (cos(x)/sin(y) + sin(y)/cos(x) - (cos(y)/tan(y)/cos(x) + sin(x)*tan(x)/sin(y)))/pi
  k = tan(x)/tan(y)
  m = k*k
  u = asin(tan(y)/tan(x))
  n = (sin(x)/sin(y))*(sin(x)/sin(y))

  complete = ellipk(m) - ellippi(n, m)
  incomplete = ellipf(u,m) - ellippi(n/k/k,1/m)/k

  return re(1.0 - factor*(incomplete + complete))

def one_side(x,y):
    m = ((-1 + x)*(1 + y))/((1 + x)*(-1 + y))
    k = sqrt(m)
    u = asin(1/k)
    EE = ellipe(m)
    EF = re(ellipf(u,m))
    n = (-1 + x)/(-1 + y)
    EPI = ellippi(n/m,1/m)/k
    return re(-(EE*(1 + x)*(-1 + y) + (x - y)*(EF + EPI*(x-y) + EF*y))/sqrt(((1 + x)*(1 - y))))

def elliptic_core_m(x,y):
    x = cos(x)
    y = cos(y)
    middle = abs(x-y)/2.0;
    return -middle + (one_side(x,y) + one_side(-x,-y)) / pi;

def elliptic_core_h(core_g,x,y):
    M = elliptic_core_m(x,y)
    factor = -sin(y)*tan(y)
    sine2 = sin(y)*sin(y)
    return (M - core_g*cos(y) + cos(y) - cos(x)*cos(y))/sine2
    
def elliptic( alpha, icomp, theta, w, nmax, love_numbers, earth_model=earth.default, computeU=True, computeV=True, computeG=True ):
    core_value = float(elliptic_core_g(alpha * DEGREES,theta * DEGREES))

    core_derivative = 0
    if computeV:
        core_derivative = float(elliptic_core_h(core_value, alpha * DEGREES,theta * DEGREES))
    
    return diskload.core( alpha, icomp, theta, w, core_value, core_derivative, nmax, love_numbers, earth_model, computeU, computeV, computeG )
