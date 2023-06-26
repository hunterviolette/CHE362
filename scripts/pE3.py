from pint import UnitRegistry
from sympy.physics.units import mol, hour
from math import log, modf

from baseFunctions import CHE362

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
uReg.define('kmol = 10**3 * mol')

q  = uReg.Quantity

class PE3(CHE362):

  @staticmethod
  def Two():
    f_ = 400*10**3 * mol/hour
    xF, xD, xB,  = .3, .99, .01
    
    print("==== Problem 1a ===")
    b_, d_ = PE3.Solve_DB(f_, xF, xD, xB)

    print("==== Part B/C ===")
    RminLiq = PE3.Solve_Rmin(xD, xF, .51) # vertical xF y-value
    RminVap = PE3.Solve_Rmin(xD, xF, .152, satLiquid=False) # horizontal xF x-value

    r_, minStages = 3, 11
    yInt = xD / (r_ + 1)
    
    stagesIdeal, eta = 23, .7
    stagesReal = PE3.Real_Stages(stagesIdeal, eta)

    print(
          f"==== Part D/E ====",
          f"Minimum stages: {minStages}",
          f"Ideal stages: {stagesIdeal}",
          f"Feed stage: {10}",
          f"Real stages: {stagesReal}",
          sep='\n'
        )

    mW = q(86.3, 'g/mol'),
    rhoV = q(3.07, 'kg/m**3').to('kg/m**3'),
    rhoL = q(615, 'kg/m**3').to('kg/m**3'),
    sigma = 13.3, # dyne / cm,
    percentFlood = .75,
    activeArea= .8,

    f_LV = (r_ / (r_ + 1)) * (rhoV / rhoL)**.5
    kV = q(10**(-.94506 - .70234 * log(f_LV, 10) - .22618 * log(f_LV, 10)**2), 'ft/s')

    uC = (kV * (sigma / 20)**.2 * ((rhoL - rhoV) / rhoV)**.5).to('ft/s')
    uO = uC * percentFlood

    v_ = d_ * (r_ + 1)
    vDot = (v_ * mW / rhoV).to('ft**3/s')

    area = (vDot / uO).to('ft**2')
    actualArea = area / activeArea
    diameter = (4 * actualArea / pi)**.5

    frac, whole = modf(diameter.magnitude)
    if frac >= 0.5:
      diameter = round(diameter, 0)
    elif (frac < 0.5) and (frac > .01):
      diameter = q(whole + .5, 'ft')
    else:
      diameter = round(diameter, 0)

PE3.Two()