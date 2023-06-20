from pint import UnitRegistry
from sympy.physics.units.util import convert_to
from sympy import Eq, symbols, solve
from math import log, pi
from sympy.physics.units import newton, feet, second, pound, cm
from math import pi

import pandas as pd

from baseFunctions import CHE362

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

'''
Estimate diameter of top tray assuming:
    vapor is 75% flooding velocity
    active area is 80% 
'''

r_ = 2
rhoV, rhoL = q(1.154, 'kg/m**3'), q(792, 'kg/m**3')
f_LV = (r_ / (r_ + 1)) * (rhoV / rhoL)**.5

kV = symbols('kV')
kV = q(solve(Eq(
                10**(-.94506 - .70234 * log(f_LV) - .22618 * log(f_LV)**2),
                kV,
            ),
            kV
        )[0]*10, 'ft/s')

dyne = 0.00001 
sigma = 19.7 # dyne / cm
#uC = convert_to(kV * (sigma / 20)**.2 * ((rhoL - rhoV) / rhoV), feet / second)
uC = kV * (sigma / 20)**.2 * ((rhoL - rhoV) / rhoV)**.5
#uC = q(10.385, 'ft/s')

percFlood, activA = .75, .8
uO = q(uC * percFlood, 'ft/s')

d_, mW = q(39.796*10**3, 'mol/hr'), q(32, 'kg/mol')
v_ = d_ * (r_ + 1)
vDot = (v_ * mW / rhoV).to('ft**3/s')

area = (vDot / uO).to('ft**2')
actualArea = area / activA
diameter = (4 * actualArea / pi)**.5

print( 
      f"kV: {kV}",
      f"uC: {uC}", 
      f"uO: {uO}", 
      f"V: {v_}",
      f"vDot: {vDot}", 
      f"area: {area}", 
      f"actual area: {actualArea}", 
      f"diameter: {diameter}", 
    sep='\n')
