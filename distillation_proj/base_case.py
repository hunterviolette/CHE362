from sympy.physics.units import atm, mol, hour
from math import pi, log

import pint
import pandas as pd

from sys import path
path.append('..')
from scripts.baseFunctions import CHE362


uReg = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
uReg.define('dollars = 1')
q  = uReg.Quantity


class DP(CHE362):

  def __init__(self) -> None:
    pass
  
  def XY_Diagram(self, plot: bool = False):
    f_ = 500*10**3 * mol/hour
    xF, xD, xB,  = .2, .99, .01
    yF, yD, yB = .533, .99, .01

    self.b_, self.d_ = DP.Solve_DB(f_, xF, xD, xB)
    self.r_ = DP.Solve_Rmin(xD, yF, xF, False) * 1.2
    
    if plot:
      fig = DP.Generate_YX_Diagram(pd.read_csv('base_case_raw.csv'))
      fig.add_vline(x=xF, line_dash="dot", annotation_text="Feed line", annotation_font_color="white")
      
      lines = []
      lines.append(pd.DataFrame({"X": [0, xD], "Y":[xD / (self.r_ + 1), yD]}).astype(float))
      lines.append(pd.DataFrame({"X": [xB, xF], "Y":[yB, .498]}).astype(float))

      fig = DP.PlotLines(fig, lines).show()

  def Area(self, steamTemp: pint.Quantity = q(185.477, 'degC')):
    DP.XY_Diagram(self, False)
    self.d_ = q(self.d_.magnitude, 'mol/h')

    t_CWI, t_CWO = q(30, 'degC'), q(45, 'degC')
    t_D, t_B, t_ST = q(80, 'degC'), q(136, 'degC'), steamTemp
    hVap = q(30.75, 'kJ/mol') # heat of vaporization 

    dT_Ln = (((t_D - t_CWI) - (t_D - t_CWO)) / log((t_D - t_CWI) / (t_D - t_CWO))).to('delta_degC')
    dT = (t_ST - t_B).to('delta_degC')
    if dT.magnitude > 30:
      dT = q(30, 'delta_degC') 

    dT = dT.to('delta_degF')

    q_C = (self.d_ * (self.r_ + 1) * hVap).to('GJ/hr')
    q_B = (self.d_ * (self.r_ + 1) * hVap).to('GJ/hr')

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

    fP, fM, eff = 1, 1, .75
    nTrays, traySpacing = (14 - 1) / eff, q(2, 'ft')
    areaC, areaB = q(55.016, 'm**2'), q(110.164, 'm**2')
    qB, qC = q(4.399, 'GJ/h'), q(4.399, 'GJ/h')

    areaC, areaB = areaC.to('m**2').magnitude, areaB.to('m**2').magnitude
    h_ = q(10, 'ft') + nTrays * traySpacing
    vol = (pi * (self.d / 2)**2 * h_).to('m**3').magnitude
    area = (pi * (self.d / 2)**2).to('m**2').magnitude

    if nTrays > 20:
      fQ = 1
    else:
      fQ = 10**(.4771 + .08516 * log(nTrays, 10) - .3473 * log(nTrays, 10)**2)

    cShell = 10**(3.4974 + .4485 * log(vol, 10) + .1074 * log(vol, 10)**2) * \
                (2.25 + (1.82 * fP * fM))

    cTrays = 10**(3.3322 + .4838 * log(area, 10) + .3434 * log(area, 10)**2) * \
                fQ * nTrays
    
    cCondenser = 10**(4.3247 - .303 * log(areaC, 10) + .1634 * log(areaC, 10)**2) * \
                (1.63 + (1.66 * fM * fP))
    
    cBoiler = 10**(4.4646 - .5277 * log(areaB, 10) + .3955 * log(areaB, 10)**2) * \
                (1.63 + (1.66 * fM * fP))
    
    cSteam = qC * q(14.83, 'dollars/GJ') * q(350, 'day').to('h')
    cWater = qB * q(.354, 'dollars/GJ') * q(350, 'day').to('h')

    capCost = cShell + cTrays + cCondenser + cBoiler
    opCost = cSteam + cWater

    eAOC = opCost + (capCost / 5)
    steamAOC = cSteam / eAOC

    print( 
          f"==== Economic Calculations ====",
          f"cost: shell {cShell}",
          f"cost: trays {cTrays}",
          f"cost: condenser {cCondenser}", 
          f"cost: boiler {cBoiler}", 
          f"cost: steam {cSteam}",
          f"cost: water {cWater}", 
          f"cost: capital {capCost}", 
          f"annual: operating cost {opCost}", 
          f"5-yr annual: total cost: {eAOC}", 
          f"Steam cost / total cost: {steamAOC.magnitude}", 
        sep='\n')

DP().Economics()