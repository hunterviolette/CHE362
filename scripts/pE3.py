from pint import UnitRegistry
from sympy.physics.units import mol, hour

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

    PE3.Distillation_Diameter(
        mW = q(86.3, 'g/mol'),
        rhoV = q(3.07, 'kg/m**3').to('kg/m**3'),
        rhoL = q(615, 'kg/m**3').to('kg/m**3'),
        sigma = 13.3, # dyne / cm,
        percentFlood = .75,
        activeArea= .8,
        r_ = r_,
        d_ = d_
      )

PE3.Two()