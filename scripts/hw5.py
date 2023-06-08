from pint import UnitRegistry
from sympy.physics.units.util import convert_to
from sympy import Eq, symbols
from sympy.physics.units import cm, atm, m, molar_gas_constant, \
    kelvin, second, mol, gram, kilogram, mmHg, mol, hour, pounds, minute

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
        print('=== See Figure 1 ===')
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
        print('=== See Figure 2 ===')
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
                            x1 = X1, # y1 only needed if absorbing
                            m = HW5.Slope_OMTC(vp, p, 'raoults', gamma=gamma),
                            Ls = abs(Ls), Vs = abs(Vs), absorption=True, y1 = Y1,
                        )

    @staticmethod
    def Two():
        '''
        recover oil from safflower seeds. Crushed seeds contacted with solvent to recover oil
        The seeds contain 0.2 lbs-oil / lbs-seeds
        x0 = 0
        Recover 95% oil from seeds. 1000 lbs of seeds / min
        
    
        a) If a single cocurrent stage is used, what is the minimum solvent rate needed?  Show your 
        solution graphically.
        '''
        print('=== Problem 2 ====') 
        x1, y1 = .2, 0
        x2 = x1 * .05
        Ls = 1000 * (pounds / minute) * (1 - x1)
        print('=== Part a ====') 
        print('=== See Figure 1 ===')
        x1_, y1_, x2_ = HW5.CapXY(x1), HW5.CapXY(y1), HW5.CapXY(x2)
        print(f"[{x1_}, {y1_}], [{x2_}, Y2]")
        # Find y2 on graph .01
        x2_, y2_ = .005, .006

        Vs = HW5.Solve_MaterialBal_Streams('cocurrent',
                x1_, x2_, y1_, y2_, Ls, symbols('Vs'), symbols('Vs'))
        
        print('=== Part b ====') 
        print('=== See Figure 2 ===')
        Vs *= 2
        y2_ = HW5.Solve_MaterialBal_Streams('cocurrent',
                x1_, x2_, y1_, symbols('y2'), Ls, Vs, symbols('y2'))
        
        y2 = (y2_ / (1 + y2_))
        print(f"y2: {y2} lb-oil / lb-solvent")
        
        print('=== Part c ====') 
        x1_, y1_, x2_ = HW5.CapXY(x1), HW5.CapXY(y1), HW5.CapXY(x2)
        Ls = 1000 * (pounds / minute) * (1 - x1)

        x2_, y2_ = .01, .185
        Vs = HW5.Solve_MaterialBal_Streams('countercurrent',
                x1_, x2_, y1_, y2_, Ls, symbols('Vs'), symbols('Vs'))
                
        print('=== Part d ====') 
        print('=== See Figure 3 ===')
        y2_1_5 = y2_ / 1.5
        print(f"Find Yn for 1.5Vs: {y2_1_5}")

        print('=== Part e ====') 
        print('=== See Figure 4 ===')
        x1_, y1_, x2_ = HW5.CapXY(x1), HW5.CapXY(y1), HW5.CapXY(x2)

        print('=== Part f ====') 
        HW5.KremserEquation(x1 = x1_,
                            y1 = 0,
                            m = ((y2_ - y1_) / (x1_ - x2_)),
                            Ls = 800,
                            Vs = 1600,
                            absorption = False,
                            xN = x2_
                        )

h = HW5().One()