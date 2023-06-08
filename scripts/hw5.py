from pint import UnitRegistry
from sympy.physics.units.util import convert_to
from sympy import Eq, symbols
from sympy.physics.units import cm, atm, m, molar_gas_constant, \
    kelvin, second, mol, gram, kilogram, mmHg, mol, hour

import pandas as pd

from baseFunctions import CHE362

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

class HW5(CHE362):

    @staticmethod
    def One():
        '''
        MeOH removed from N2 gas by water, countercurrent absorber
        '''

        y0, x0, y1 = .02, float(0), .001
        v_ = (100 * 10**3) * mol / hour
        Vs, gamma = v_ * (1 - y0), 1.725

        vp, p = q(126, 'mmHg').to('atm').magnitude * atm, 1 * atm

        x1 = HW5.Solve_EquilibRelationship(
                vp, p, symbols('x1'), y1, symbols('x1'), gamma
            )

        
        X1, X2 = HW5.CapXY(x0), HW5.CapXY(x1)
        Y1, Y2 = HW5.CapXY(y0), HW5.CapXY(y1)

        print('==== Part b ====')
        Ls = HW5.Solve_MaterialBal_Streams('countercurrent', 
                                    X1, X2, Y1, Y2,
                                    Ls = symbols('Ls'), 
                                    Vs = Vs, 
                                    solveFor = symbols('Ls')
                                )
        
        lines = []
        lines.append(pd.DataFrame({"X": [X1, X2], "Y":[Y1, Y2]}).astype(float))

        print('==== Part c ====')
        Xn_1_4Ls = HW5.Solve_MaterialBal_Streams('countercurrent', 
                                    X1 = HW5.CapXY(x0), 
                                    X2 = symbols('Xn'), 
                                    Y1 = HW5.CapXY(y0), 
                                    Y2 = HW5.CapXY(y1), 
                                    Ls = Ls*1.4, 
                                    Vs = Vs, 
                                    solveFor = symbols('Xn'),
                                )
        
        lines.append(pd.DataFrame({"X": [X1, Xn_1_4Ls], "Y":[Y1, Y2]}).astype(float))
        
        fig = HW5.EquilibriumDiagram(vp, p, # vapor pressure of A, total pressure
                xStep = .005, yStep = .01, # Formatting Figure tickmarks
                xRange = [0, .1], yRange = [0, .1],  # Formatting Figure axis range
                dataGen = [0, .5, .001], # Generate Equilibrium diagram data [startingValue, EndValue, StepValue]
                gamma = gamma
            )
        
        fig = HW5.PlotLines(fig, lines).show()
        
        print('==== Part d ====')
        stages = HW5.KremserEquation(
                                    xN = X2,
                                    yN = Y2,
                                    x0 = X1,
                                    y1 = Y2,
                                    m = HW5.Slope_OMTC(vp, p, 'raoults', 
                                            idealSolution=False, gamma=gamma
                                        ),
                                    Ls = 1.4*Ls, Vs= Vs, absorption=True
                                    
                                )

h = HW5().One()