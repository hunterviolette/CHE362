from pint import UnitRegistry
from sympy.physics.units.util import convert_to
from sympy.physics.units import cm, atm, m, molar_gas_constant, \
    kelvin, second, mol, gram, kilogram

from sympy import Eq, symbols
from baseFunctions import CHE362

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

class E1(CHE362):
    
    @staticmethod
    def One():
        print('==== Problem 1 ====')
        print('==== Gas CALCULATIONS ====')
        viscosity = q(1.86 * 10**-5, 'kg/(m*s)') 
        density = q(1.176, 'kg/m**3') 
        d_AB = q(1.98 * 10**-5, 'm**2/s').to('cm**2/s') 
        
        cL = q(0.5, 'm').to('m')
        gas_flow = q(0.1, 'm**3/s').to('m**3/s')
        area = q(0.5 / 2, 'm')**2 * 3.14159
  
        gas_velocity = gas_flow / area

        Re = E1.ReynoldsNumber(L = cL,
                                rho = density.to('kg/m**3'),
                                v = gas_velocity.to('m/s'),
                                mu = viscosity.to('kg/(m*s)')
                            )

        Sc = E1.SchmidtNumber(mu = viscosity.to('kg/(m*s)'),
                            rho = density.to('kg/m**3'),
                            d_AB = d_AB.to('m**2/s')
                            )

        Sh = 0.023 * Re**.8 * Sc**(1/3) # CHECK COORELATION EQN
        print(f"Sherwood Number: {Sh}")

        c = convert_to(1*atm / (molar_gas_constant * (298 * kelvin)), mol/m**3)
        print(f"Molarity: {c}")


        kX = E1.SherwoodSolver(Sh, # Sherwood Number
                        c*cm, # characteristic length
                        c, # molarity: mol/volume
                        d_AB.magnitude * cm**2 / second, # diffusivity: area/time
                        symbols('k'), # mtc: mol/(area*time)
                        symbols('k'), # variable to solve for
                        False # Low mass transfer
                    )
    
        print('==== Liquid Phase ====')
        d_AB = 2.0 * 10**-9
        c = (1.000 * gram / cm**3) / (18 * gram / mol)
        cL = q(0.001, 'm').to('m')
        kX = E1.SherwoodSolver(3.41, # Sherwood Number
                    cL.magnitude*cm, # characteristic length
                    c, # molarity: mol/volume
                    d_AB * cm**2 / second, # diffusivity: area/time
                    symbols('k'), # mtc: mol/(area*time)
                    symbols('k'), # variable to solve for
                    False # Low mass transfer
                )
        
    def Two():
        print('==== Problem 2 ====')
        xB, yB = .001, 0.1
        kX = .2 * (mol / (m**2 * second))
        kY = .1 * (mol / (m**2 * second))
        He, p = 1*atm, 2*atm

        xI, yI = E1.Solve_MoleFracInterface(
                                E1.HenrysLaw(
                                            xI = symbols('xI'), # xI variable is currently unknown
                                            yI = symbols('yI'), # yI variable is currently unknown
                                            He = He,
                                            p = p
                                            ),
                                xB, # mole fraction at bulk for x phase
                                yB, # mole fraction at bulk for y phase
                                kX, # mtc of phase x
                                kY, # mtc of phase y
                                )
        
        nAX, nAY = E1.FluxBetweenPhases(xB, xI, yB, yI, kX, kY) 

        oKX, oKY = E1.SOLVE_OVERALL_MTC(
                            kX,
                            kY,
                            E1.Slope_OMTC(He, p),
                            symbols('oMTC'), # oMTC variable is currently unknown
                            symbols('oMTC'), # solving equation for oMTC 
                            )

        E1.PhaseResistances(
                        kX,
                        oKX,
                        kY,
                        oKY
                        )

        
E1.One()
E1.Two()