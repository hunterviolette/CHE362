from pint import UnitRegistry
from sympy.physics.units.util import convert_to
from sympy import Eq, symbols
from sympy.physics.units import cm, atm, m, molar_gas_constant, \
    kelvin, second, mol, gram, kilogram, mmHg, mol, hour
from numpy import arange

import plotly_express as px
import pandas as pd

#px.defaults.template = "plotly_dark"
#px.defaults.color_continuous_scale = px.colors.sequential.Blackbody

from baseFunctions import CHE362



uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

class HW4(CHE362):
  
  @staticmethod
  def One():
    print('==== Problem 1 ====')
    Sh_y, Sh_x, p = 10, 36, 1*atm
    yB, xB, t = .01, .005, 298*kelvin
    vpA = q(126, 'mmHg').to('atm').magnitude * atm
    
    print('==== Liquid Phase ====')
    c = (1.012 * gram / cm**3) / (18 * gram / mol)
    kX = HW4.SherwoodSolver(Sh_x, # Sherwood Number
                1*cm, # characteristic length
                c, # molarity: mol/volume
                q(HW4.WilkeChangDiffusivity(
                          298, # temperature in K
                          15.999 + 2, # molecular weight of B in g/mol
                          0.8818, # viscosity of B in cp 
                          42.78, # molar volume of A at normal boiling in cm**3/mol 
                          2.6 # association factor ??
                        ), 'cm**2/s').magnitude * cm**2 / second, # diffusivity: area/time
                symbols('k'), # mtc: mol/(area*time)
                symbols('k'), # variable to solve for
                False # Low mass transfer
              )
    
    print('==== Gas Phase ====')
    c = convert_to(p / (molar_gas_constant * t), mol/m**3)
    print(f"Molarity: {c}")

    kY = HW4.SherwoodSolver(Sh_y, # Sherwood Number
                1*cm, # characteristic length
                c, # molarity: mol/volume
                q(HW4.FullerDiffusivity(
                          1, # pressure in atmospheres
                          298, # temperature in kelvin 
                          32.042, # molecular weight species A
                          14.007 * 2, # molecular weight species B 
                          15.9 + (2.31 * 4) + 6.11, # diffusion volume species A 
                          18.5, # diffusion volume species B
                          'Methanol', # name of species A
                          'N2' # name of species B
                        ), 'cm**2/s').magnitude * cm**2 / second, # diffusivity: area/time
                symbols('k'), # mtc: mol/(area*time)
                symbols('k'), # variable to solve for
                False # Low mass transfer
              )

    xI, yI = HW4.Solve_MoleFracInterface(
                              HW4.RaoultsLaw(
                                          xI = symbols('xI'),
                                          yI = symbols('yI'),
                                          vpA = vpA,
                                          p = p,
                                          gamma = 1.725,
                                          idealSolution = False
                                        ),
                              xB, # mole fraction at bulk for x phase
                              yB, # mole fraction at bulk for y phase
                              kX, # mtc of phase x
                              kY, # mtc of phase y
                            )
    
    nAX, nAY = HW4.FluxBetweenPhases(xB, xI, yB, yI, kX, kY) 

    oKX, oKY = HW4.SOLVE_OVERALL_MTC(
                          kX,
                          kY,
                          HW4.Slope_OMTC(vpA, p,
                                          lawUsed='raoults',
                                          idealSolution=False,
                                          gamma=1.725
                                        ),
                          symbols('oK'),
                          symbols('oK'),
                        )

    HW4.PhaseResistances(kX, oKX, kY, oKY)

  @staticmethod
  def Two():
    kmol = 10**3 * mol

    # cocurrent absorber 
    y0, y1 = .02, .001 # [2%, 0.1%] MeOH
    vp, p = 1 * atm, 2 * atm
    v = 100 * kmol / hour # Total flow rate of extract stream
    vS = v * (1 - y0)
    print(vS)

    df = pd.DataFrame()
    for x in arange(0, 1, .01):
      y = float((vp / p) * x)
      capX, capY = HW4.CapXY(x, f'x_{x}'), HW4.CapXY(y, f'y_{x}')
      df = pd.concat([df, pd.DataFrame({"x":[x], "y":[y], "X":[capX], "Y":[capY]}
                                    )], ignore_index=True)
      
    print(df)
    px.line(data_frame=df, x='X', y='Y').show()

  @staticmethod
  def Three():
    xB, yB = 1 * 10**-6, 0.2
    kX = 1.5 * (mol / (m**2 * second))
    kY = 0.8 * (mol / (m**2 * second))
    He, p = 43800*atm, 1*atm

    xI, yI = HW4.Solve_MoleFracInterface(
                              HW4.HenrysLaw(
                                          xI = symbols('xI'),
                                          yI = symbols('yI'),
                                          He = He,
                                          p = p
                                        ),
                              xB, # mole fraction at bulk for x phase
                              yB, # mole fraction at bulk for y phase
                              kX, # mtc of phase x
                              kY, # mtc of phase y
                            )
    
    nAX = HW4.FluxBetweenPhases(
                              xxI = xI,
                              xxB= xB,
                              k = kX,
                              phase = 'x'
                            ) 

    oKX, oKY = HW4.SOLVE_OVERALL_MTC(
                          kX,
                          kY,
                          He,
                          p,
                          symbols('oK'),
                          symbols('oK'),
                        )

    HW4.PhaseResistances(
                      kX,
                      oKX,
                      kY,
                      oKY
                    )


h = HW4()
h.Two()
