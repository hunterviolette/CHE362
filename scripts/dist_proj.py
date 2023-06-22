from pint import UnitRegistry

from sympy.physics.units import atm, mol, hour
from math import pi, log

import pandas as pd
import plotly_express as px
import plotly.graph_objects as go

from baseFunctions import CHE362

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity


class DP(CHE362):

  def __init__(self) -> None:
    pass
  
  def XY_Diagram(self, plot: bool = False):
    f_ = 500*10**3 * mol/hour
    xF, xD, xB,  = .2, .99, .01
    yF, yD, yB = .533, .99, .01

    self.b_, self.d_ = DP.Solve_DB(f_, xF, xD, xB)
    self.r_ = DP.Solve_Rmin(xD, yF, xF, 1.2)
    
    if plot:
      fig = DP.Generate_YX_Diagram(pd.read_csv('tables/proj_base.csv'))
      fig.add_vline(x=xF, line_dash="dot", annotation_text="Feed line", annotation_font_color="white")
      
      lines = []
      lines.append(pd.DataFrame({"X": [0, xD], "Y":[xD / (self.r_ + 1), yD]}).astype(float))
      lines.append(pd.DataFrame({"X": [xB, xF], "Y":[yB, .498]}).astype(float))

      fig = DP.PlotLines(fig, lines).show()

  def Area(self):
    DP.XY_Diagram(self)
    self.d_ = q(self.d_ * hour / mol, 'mol/h')

    t_CWI, t_CWO = q(30, 'degC'), q(45, 'degC')
    t_D, t_B, t_ST = q(80, 'degC'), q(136, 'degC'), q(185.477, 'degC')

    dT_Ln = (((t_D - t_CWI) - (t_D - t_CWO)) / log((t_D - t_CWI) / (t_D - t_CWO))).to('delta_degF')
    dT = (t_ST - t_B).to('delta_degF')
    if dT.magnitude > 50:
      dT = q(50, 'delta_degF') 

    q_C = (self.d_ * (self.r_ + 1) * q(30.75,'kJ/mol')).to('GJ/hr')
    q_B = q_C

    muOc = q(150, 'Btu/(h * ft**2 * degF)')
    muOb = q(200, 'Btu/(h * ft**2 * degF)')

    a_C = (q_C / (muOc * dT_Ln)).to('m**2')
    a_B = (q_B / (muOb * dT)).to('m**2')

    print(
          "==== Area Calculations ====",
          f"dT Ln, condenser: {dT_Ln}",
          f"dT, reboiler: {dT}",
          f"Heat duty, condenser: {q_C}",
          f"Heat duty, reboiler: {q_B}",
          f"Heat transfer coefficient, condenser: {muOc}",
          f"Heat transfer coefficient, reboilder: {muOb}",
          f"Heat transfer area, condenser: {a_C}",
          f"Heat transfer area, reboiler: {a_B}",
        sep='\n'
      )

  def Diameter(self):
    DP.Area(self)
    #### Input vars

    mW = q(78.11, 'g/mol') # Molecular weight of benzene
    rhoV= ( q(35.01, 'mol/m**3') * mW).to('kg/m**3') # rhoV Table 135 Yaws at Tb
    rhoL =( q(10450, 'mol/m**3') * mW).to('kg/m**3')  # Molecular weight of liquid benzene

    sigma = 21.11 # dyne / cm
    percFlood, activeArea = .75, .8
    ####

    f_LV = (self.r_ / (self.r_ + 1)) * (rhoV / rhoL)**.5
    kV = q(10**(-.94506 - .70234 * log(f_LV, 10) - .22618 * log(f_LV, 10)**2), 'ft/s')

    uC = (kV * (sigma / 20)**.2 * ((rhoL - rhoV) / rhoV)**.5).to('ft/s')
    uO = uC * percFlood

    v_ = self.d_ * (self.r_ + 1)
    vDot = (v_ * mW / rhoV).to('ft**3/s')

    area = (vDot / uO).to('ft**2')
    actualArea = area / activeArea
    self.d = (4 * actualArea / pi)**.5

    print( 
          f"==== Diameter Calculations ====",
          f"F_LV: {f_LV}",
          f"kV: {kV}",
          f"uC: {uC}", 
          f"uO: {uO}", 
          f"V: {v_}",
          f"Vdot: {vDot}", 
          f"area: {area}", 
          f"actual area: {actualArea}", 
          f"diameter: {self.d}", 
        sep='\n')

  def Economics(self):
    DP.Diameter(self)

    nTrays, traySpacing = 23, q(2, 'ft')
    fP, fM, vC = 1, 1, ((self.d / 2)**2).to('m**2').magnitude

DP().Diameter()