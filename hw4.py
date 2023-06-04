from pint import UnitRegistry
from sympy.physics.units.util import convert_to
from sympy.physics.units import cm, atm, m, molar_gas_constant, \
    kelvin, second, mol, gram, kilogram

from sympy import Eq, symbols
import sympy.abc as cnst
from baseFunctions import MassTransfer, Diffusion

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

class HW4(MassTransfer, Diffusion):
  
  @staticmethod
  def One():
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
h.One()