from pint import UnitRegistry
from sympy.physics.units.util import convert_to
from sympy.physics.units import cm, atm, m, molar_gas_constant, \
    kelvin, second, mol, gram, kilogram

import sympy.abc as cnst
from baseFunctions import MassTransfer, Diffusion

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

class HW3(MassTransfer, Diffusion):
  
  @staticmethod
  def One():
    '''
    #1 Absorb SO2 from air using water (bubble of SO2 gas rising through water), 
    
    the mass transfer coefficient on the liquid side has been correlated as:
        Sh = 2 + (0.6 * Re**(1/2) * Sc**(1/3)

    where the characteristic length is the bubble diameter.

    Estimate the liquid phase mass transfer coefficient 
      for a 1 cm diameter bubble of SO2 in air rising at a rate of 0.2 m/s 
      at a temperature of 298 K and pressure of 20 atm.
    '''
    print ('==== Problem 1 ====')
    # density and viscosity of water 
    viscosity = q(0.8848, 'centipoise') # LEAVE IN cp
    density = q(1.012, 'g/cm**3') # LEAVE IN g/mol

    Re = HW3.ReynoldsNumber(L = q(1, 'cm').to('m'),
                            rho = density.to('kg/m**3'),
                            v = q(.2, 'm/s'),
                            mu = viscosity.to('kg/(m*s)')
                          )
    
    d_AB = q(HW3.WilkeChangDiffusivity(298, # temperature in K
                          15.999 + 2, # molecular weight of B in g/mol
                          viscosity.magnitude, # viscosity of B in cp 
                          43.82, # molar volume of A at normal boiling in cm**3/mol 
                          2.6 # association factor ??
                        ), 'cm**2/s')

    Sc = HW3.SchmidtNumber(mu = viscosity.to('kg/(m*s)'),
                          rho = density.to('kg/m**3'),
                          d_AB = d_AB.to('m**2/s')
                        )

    Sh = 2 + (0.6 * Re**(1/2) * Sc**(1/3)) # CHECK COORELATION EQN
    print(f"Sherwood Number: {Sh}")
    
    c = (density.magnitude * gram / cm**3) / (18 * gram / mol)
    print(f"Molarity: {c}")

    kX = HW3.SherwoodSolver(Sh, # Sherwood Number
                    1*cm, # characteristic length
                    c, # molarity: mol/volume
                    d_AB.magnitude * cm**2 / second, # diffusivity: area/time
                    cnst.k, # mtc: mol/(area*time)
                    cnst.k, # variable to solve for
                    False # Low mass transfer
                  )


  @staticmethod
  def Two():
    '''
    #2 strip toluene (A) from liq water using N2 gas (B). 

    The gas phase correlations is: 
      Sh = 1.2 * Re**0.64 * Sc**0.33
        WHERE the characteristic length is the packing size.

    The liquid phase correlations is: 
      Sh = 0.1 * Re**0.3 * Sc**0.5
        WHERE the characteristic length is the packing size.

    Estimate the mass transfer coefficients (kx and ky) 

    T=298 K, P = 1 atm, Packing size = 5 cm, Gas velocity = 1 ft/s, Liquid velocity = 1 ft/s

    You can assume dilute solutions. For the liquid and gas, use the physical properties of
    pure water and nitrogen respectively. Use the methods we have covered in class to
    estimate the diffusivities. 
    '''

    print('==== Problem 2 ====')
    print('==== GAS CALCULATIONS ====')
    # density and viscosity of N2 at 298K 
    viscosity = q(.01724, 'centipoise') # LEAVE IN cp
    density = q(0.0012506, 'g/cm**3') # LEAVE IN g/mol

    Re = HW3.ReynoldsNumber(L = q(5, 'cm').to('m'),
                            rho = density.to('kg/m**3'),
                            v = q(1, 'foot/s').to('m/s'),
                            mu = viscosity.to('kg/(m*s)')
                          )
    
    # toulene is C7H8
    d_AB = q(HW3.FullerDiffusivity(
                          1, # pressure in atmospheres
                          298, # temperature in kelvin 
                          (12 * 7) + (1 * 8), # molecular weight species A
                          14.007 * 2, # molecular weight species B 
                          (15.9 * 7) + (2.31 * 8) - 18.3, # diffusion volume species A 
                          18.5, # diffusion volume species B
                          'Toluene', # name of species A
                          'N2' # name of species B
                        ), 'cm**2/s')
    

    Sc = HW3.SchmidtNumber(mu = viscosity.to('kg/(m*s)'),
                          rho = density.to('kg/m**3'),
                          d_AB = d_AB.to('m**2/s')
                        )

    Sh = 1.2 * Re**(.64) * Sc**(1/3) # CHECK COORELATION EQN
    print(f"Sherwood Number: {Sh}")

    length, p, t = 5*cm, 1*atm, 298*kelvin
    c = convert_to(p / (molar_gas_constant * t), mol/m**3)
    print(f"Molarity: {c}")

    kY = HW3.SherwoodSolver(Sh, # Sherwood Number
                    length, # characteristic length
                    c, # molarity: mol/volume
                    d_AB.magnitude * cm**2 / second, # diffusivity: area/time
                    cnst.k, # mtc: mol/(area*time)
                    cnst.k, # variable to solve for
                    False # Low mass transfer
                  )

    print('==== Liquid CALCULATIONS ====')
    # density and viscosity of water at 298K 
    viscosity = q(0.8848, 'centipoise') # LEAVE IN cp
    density = q(1.012, 'g/cm**3') # LEAVE IN g/mol

    Re = HW3.ReynoldsNumber(L = q(5, 'cm').to('m'),
                            rho = density.to('kg/m**3'),
                            v = q(1, 'foot/s').to('m/s'),
                            mu = viscosity.to('kg/(m*s)')
                          )
    
    d_AB = q(HW3.WilkeChangDiffusivity(298, # temperature in K
                          15.999 + 2, # molecular weight of B in g/mol
                          viscosity.magnitude, # viscosity of B in cp 
                          118.5, # molar volume of A at normal boiling in cm**3/mol 
                          2.6 # association factor ??
                        ), 'cm**2/s')

    Sc = HW3.SchmidtNumber(mu = viscosity.to('kg/(m*s)'),
                          rho = density.to('kg/m**3'),
                          d_AB = d_AB.to('m**2/s')
                        )

    Sh = 0.1 * Re**0.3 * Sc**0.5 # CHECK COORELATION EQN
    print(f"Sherwood Number: {Sh}")

    c = (density.magnitude * gram / cm**3) / (18 * gram / mol)
    print(f"Molarity: {c}")

    kX = HW3.SherwoodSolver(Sh, # Sherwood Number
                    5*cm, # characteristic length
                    c, # molarity: mol/volume
                    d_AB.magnitude * cm**2 / second, # diffusivity: area/time
                    cnst.k, # mtc: mol/(area*time)
                    cnst.k, # variable to solve for
                    False # Low mass transfer
                  )

h = HW3()
h.One()
h.Two()