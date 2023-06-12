from pint import UnitRegistry
from sympy.physics.units.util import convert_to
from sympy import Eq, symbols, solve
from sympy.physics.units import atm, mol, mol, hour, meter, gram, feet, second, pound
from math import pi

import pandas as pd

from baseFunctions import CHE362

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

class HW6(CHE362):

    @staticmethod
    def One():
      '''
      Acetone removed from N2 stream with water
      99% of acetone removed from gas at yN
      x0 = 0, y0 = .04, yN = y0 * (1 - .01), V = 100 kmol / hr
      T = 30C, P = 1atm, vpA = 5.505 pisa, gamma = 6.1615
      Fp = 32, deltaP = .25 inches H20 / ft
      Ls_min = 1.4 * Ls
      '''
      x0, ynp1, = 0, .04
      y1 = ynp1 * .01
      
      kmol = 10**3 * mol
      v_ = 100 * (kmol / hour)
      Vs = v_ * (1 - ynp1)

      vp = q(5.505, 'psi').to('atm').magnitude * atm
      p, gamma = 1 * atm, 6.1615

      Ls_min = HW6.PackedColumn_MaterialBal(
              X0 = HW6.CapXY(x0),
              Ls = symbols('Ls'),
              Yn_1 = HW6.CapXY(ynp1),
              Vs = Vs, 
              Xn = HW6.Solve_EquilibRelationship(
                    vp, p, symbols('xn'), ynp1, symbols('xn'), gamma
                  ),
              Y1 = HW6.CapXY(y1),
              solveFor = symbols('Ls')
            )
            
      print('==== Problem 1 ====')
      print('==== Part A ====')
      Ls = 1.4 * Ls_min
      mW_Ls, rho_Ls = q(18, 'gram/mol').to('kg/mol'), q(997, 'kg/m**3')
      Ls_gpm = ((q(Ls * (hour / mol), 'mol/hour') * mW_Ls) / rho_Ls).to('gallon/hour')

      print(f'Ls_min: {Ls_min}, Ls:{Ls}, Ls_gpm:{Ls_gpm}')

      print('==== Part B ====')
      xn = HW6.Solve_EquilibRelationship(
                    vp, p, symbols('xn'), ynp1, symbols('xn'), gamma
                  )
      print(f'Mole fraction of acetone in water: {xn}')

      print('==== Part C ====')
      idealStages = HW6.PackedColumn_Kremser(
                    Ls, Vs, HW6.Slope_OMTC(vp, p, gamma), 
                    HW6.CapXY(ynp1), x0, HW6.CapXY(y1), HW6.CapXY(xn)
                  )
      print(f"Number of ideal stages: {idealStages}")

      print('==== Part D ====')
      term = (
                mW_Ls.to('g/mol').magnitude * 
                q(.8818, 'centipoise').magnitude
              ) / rho_Ls.to('lb/ft**3').magnitude
      
      _x_ = HW6.Slope_OMTC(vp, p, gamma) * term
      print(f"X_term for efficiency calculation:{_x_}")

      efficiency = HW6.PackedColumn_Efficiency(
                    _x_, symbols('e'), symbols('e')
                  )
      print(f"Column overall efficiency: {efficiency}")

      print('==== Part E ====')
      realStages = idealStages / efficiency
      print(f"Number of real stages: {realStages}")

      print('==== Part F ====')
      
      rho_A = ((q(1, 'atm') * q(58.08, 'g/mol')) / 
                (q(1,'R') * q(30, 'degC'))
              ).to('kg/m**3')
      
      mW_A = 58.08 * (gram / mol)
      mW_Ls = mW_Ls.to('g/mol').magnitude * (gram / mol)
      Mg = (ynp1 * mW_A) + ((1 - ynp1) * mW_Ls) 
      Ml = (xn * mW_A) + ((1 - xn) * mW_Ls) 
      
      Ln = solve(Eq(symbols('Ln') * (1 - xn), Ls), symbols('Ln'))[0]
      
      Gx, Gy = Ln * Ml, Vs * Mg
      print(f"Mg:{Mg}, Ml:{Ml}, Ln:{Ln}, Vs:{Vs}, Gx:{Gx}, Gy:{Gy}")
      xAxis = (Gx / Gy) * (rho_A / rho_Ls)**(1/2)
      print(f"With xAxis: {xAxis} and PressureDrop: .25, the Yaxis: 0.95")

      Cs = q(.95 / (32**.5 * .8818**(.05)), 'ft/s')
      print(f"Capacity Factor: {Cs}")

      rho_Ls, rho_A = rho_Ls.to('lb/ft**3').magnitude, rho_A.to('lb/ft**3').magnitude
      gasSuperficialV = (Cs / (rho_A / (rho_Ls - rho_A))**.5) * feet / second
      print(f"gasSuperficialVelocity: {gasSuperficialV}")
      

      Vflow = convert_to(Vs / (rho_A * (pound/feet**3)) * Mg, feet**3/second)

      print(Vflow / gasSuperficialV)

      diameter = solve(pi * (symbols('d') / 2)**2, Vflow / gasSuperficialV, symbols('d'))
      print(Vflow, diameter)

HW6.One()
