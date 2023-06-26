from pint import UnitRegistry
from sympy.physics.units import mol, hour
from math import log, modf, pi

from baseFunctions import CHE362

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
uReg.define('lbmol = 1 * mol')

q  = uReg.Quantity

class E3(CHE362):

  @staticmethod
  def One():

    print(
      f"==== Problem 1 ====",
      f"Slope: {E3.Vaporization_Slope(.20)}",
      f"Verify Slope: {E3.Vaporization_SlopeCheck(.7, .3, .2, .3)}",
      f"minimum xWater liq: .062",
      f"maxmimum xWater: .8",
      f"20% vaporization: xWater (liq): .2",
      f"20% vaporization: xWater (vap): .3",
      f"10% liq, percent Vaporization: {((1 / (-1 * ((.42 - .3) / (.1 - .3)))) - 1)*100}",
      f"Flash temperatrue: 126 C",
      sep='\n'
    )

  @staticmethod
  def Two():
    
    f_ = 250 * mol/hour
    xF, xD, xB,  = .3, .99, .01
    
    b_, d_ = E3.Solve_DB(f_, xF, xD, xB)
    d_ = q(d_.magnitude, 'lbmol/h')

    print("==== Problem 2a ===",
          b_.to('lbmol/h'),
          d_.to('lbmol/h'),
          sep='\n'
        )
        
    print("==== Part B/C ===")
    RminLiq = E3.Solve_Rmin(xD, xF, .965) # vertical xF y-value
    RminVap = E3.Solve_Rmin(xD, xF, .71, satLiquid=False) # horizontal xF x-value

    r_, minStages = 2 * RminLiq, 6
    yInt = xD / (r_ + 1)
    print(f"y-intercept Sat liquid: {yInt}")
    
    stagesIdeal, eta = 11, .7

    print(
          f"==== Part D/E ====",
          f"Minimum stages: {minStages}",
          f"Ideal stages: {stagesIdeal}",
          f"Feed stage: {3}",
          sep='\n'
        )

    mW = q(106.4, 'lb/lbmol')
    rhoV = q(.197, 'lb/ft**3')
    rhoL = q(47.5, 'lb/ft**3')
    sigma = 16.7 # dyne / cm,
    percentFlood = .70
    activeArea= .8

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

    print( 
        f"==== Diameter Calculations ====",
        f"F_LV: {f_LV}",
        f"kV: {kV}",
        f"uC: {uC}", 
        f"uO: {uO}", 
        f"V: {v_}",
        f"Vdot: {vDot}", 
        f"area: {area}", 
        f"actual area: {actualArea}", 
        f"diameter: {diameter}", 
        sep='\n')
    


E3.One()
E3.Two()