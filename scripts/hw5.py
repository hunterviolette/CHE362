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

        X1, X2 = 0, .0752 # Found on graph
        Y1, Y2 = HW5.CapXY(y0), HW5.CapXY(y1)

        print('==== Part b ====')
        Ls = HW5.Solve_MaterialBal_Streams('countercurrent', 
                                    X1, X2, Y1, Y2,
                                    Ls = symbols('Ls'), 
                                    Vs = Vs, 
                                    solveFor = symbols('Ls')
                                )
        
        lines = []
        lines.append(pd.DataFrame({"X": [X1], "Y":[Y1]}).astype(float))
        lines.append(pd.DataFrame({"X": [X1], "Y":[Y2]}).astype(float))

        
        fig = HW5.EquilibriumDiagram(vp, p, # vapor pressure of A, total pressure
                xStep = .005, yStep = .01, # Formatting Figure tickmarks
                xRange = [0, .1], yRange = [0, .1],  # Formatting Figure axis range
                dataGen = [0, .5, .001], # Generate Equilibrium diagram data [startingValue, EndValue, StepValue]
                gamma = gamma
            )
        fig.add_hline(y=Y1, line_dash = "dot")

        # Line for part b
        #lines.append(pd.DataFrame({"X": [X1, .07526], "Y":[Y2, .02042]}).astype(float))

        print('==== Part c ====')
        Ls *= 1.4
        xN = HW5.Solve_MaterialBal_Streams('countercurrent', 
                                    X1, symbols('X2'), Y1, Y2,
                                    Ls = Ls, 
                                    Vs = Vs, 
                                    solveFor = symbols('X2')
                                )
        
        # Line for part c
        lines.append(pd.DataFrame({"X": [X1, xN], "Y":[Y2, Y1]}).astype(float))
        
        fig = HW5.PlotLines(fig, lines).show()

        print('==== Part d ====')
        HW5.KremserEquation(yN = Y2, # xN Only needed if stripping
                            x0 = X1, # y1 only needed if absorbing
                            m = HW5.Slope_OMTC(vp, p, 'raoults', gamma=gamma),
                            Ls = abs(Ls), Vs = abs(Vs), absorption=True, y1 = Y1,
                        )

    @staticmethod
    def Two():
        pass

h = HW5().Two()