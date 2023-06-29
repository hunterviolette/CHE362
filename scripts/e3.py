from pint import UnitRegistry
from sympy.physics.units import mol, hour

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
    xF, xD, xB,  = .9, .99, .01
    
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
          f"Real stages: {E3.Real_Stages(stagesIdeal, eta)}",
          sep='\n'
        )

    E3.Distillation_Diameter(
      mW = q(106.4, 'lb/lbmol'),
      rhoV = q(.197, 'lb/ft**3'),
      rhoL = q(47.5, 'lb/ft**3'),
      sigma = 16.7, # dyne / cm
      percentFlood = .70,
      activeArea= .8,
        r_ = r_,
        d_ = q(d_.to('mol/h').magnitude, 'mol/h'),
        q_ = q
      )


E3.One()
E3.Two()