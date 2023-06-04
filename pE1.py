from pint import UnitRegistry
from sympy.physics.units.util import convert_to
from sympy.physics.units import cm, atm, m, molar_gas_constant, \
    kelvin, second, mol, gram, kilogram

from sympy import Eq, symbols
from baseFunctions import MassTransfer, Diffusion

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

class PE1(MassTransfer, Diffusion):
  
  @staticmethod
  def One():
    '''
    water is oxygenated by bubbling air through water.
    bubble of gas rising through liquid, mtc on liquid side is correlated to:
      Sh_L = 2 + (0.6 * Re**(1/2) * Sc**(1/3)
      
      Sh_G = 10
      
    characteristic length = 1*cm
    velocity_L = 0.2 * m/s
    t, p = 298*kelvin, 1*atm
    '''

    print('==== Problem 1 ====')
    print('==== GAS CALCULATIONS ====')
    # density and viscosity of N2 at 298K 
    d_AB = q(2.10 * 10**-5, 'm**2/s').to('cm**2/s') 


    Sh = 10
    print(f"Sherwood Number: {Sh}")

    length, p, t = 1*cm, 1*atm, 298*kelvin
    c = convert_to(p / (molar_gas_constant * t), mol/m**3)
    print(f"Molarity: {c}")

    kY = PE1.SherwoodSolver(Sh, # Sherwood Number
                    length, # characteristic length
                    c, # molarity: mol/volume
                    d_AB.magnitude * cm**2 / second, # diffusivity: area/time
                    symbols('k'), # mtc: mol/(area*time)
                    symbols('k'), # variable to solve for
                    False # Low mass transfer
                  )

    print('==== Liquid CALCULATIONS ====')
    # density and viscosity of water at 298K 
    viscosity = q(9.227 * 10**-4, 'kg/(m*s)') 
    density = q(1000, 'kg/m**3') 
    d_AB = q(2.10 * 10**-9, 'm**2/s').to('cm**2/s') 

    Re = PE1.ReynoldsNumber(L = q(1, 'cm').to('m'),
                            rho = density.to('kg/m**3'),
                            v = q(0.2, 'm/s').to('m/s'),
                            mu = viscosity.to('kg/(m*s)')
                          )

    Sc = PE1.SchmidtNumber(mu = viscosity.to('kg/(m*s)'),
                          rho = density.to('kg/m**3'),
                          d_AB = d_AB.to('m**2/s')
                        )

    Sh = 2 + (0.6 * Re**(1/2) * Sc**(1/3)) # CHECK COORELATION EQN
    print(f"Sherwood Number: {Sh}")

    c = (density.magnitude * kilogram / m**3) / (18 * gram / mol)
    print(f"Molarity: {c}")

    kX = PE1.SherwoodSolver(Sh, # Sherwood Number
                    1*cm, # characteristic length
                    c, # molarity: mol/volume
                    d_AB.magnitude * cm**2 / second, # diffusivity: area/time
                    symbols('k'), # mtc: mol/(area*time)
                    symbols('k'), # variable to solve for
                    False # Low mass transfer
                  )
    
  @staticmethod
  def Two():
    print('==== Problem 2 ====')
    xB, yB = 1 * 10**-6, 0.2
    kX = 3 * (mol / (m**2 * second))
    kY = 1 * (mol / (m**2 * second))
    He, p = 43800*atm, 1*atm

    xI, yI = PE1.Solve_MoleFracInterface(
                              PE1.HenrysLaw(
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
    
    nAX = PE1.FluxBetweenPhases(
                              xxI = xI,
                              xxB= xB,
                              k = kX,
                              phase = 'x'
                            ) 

    oKX, oKY = PE1.SOLVE_OVERALL_MTC(
                          kX,
                          kY,
                          He,
                          p,
                          symbols('oK'),
                          symbols('oK'),
                        )

    PE1.PhaseResistances(
                      kX,
                      oKX,
                      kY,
                      oKY
                    )


h = PE1()
h.One()
h.Two()