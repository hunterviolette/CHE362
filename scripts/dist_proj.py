from sympy.physics.units import atm, mol, hour
from math import pi, log, modf

import pint
import pandas as pd
import plotly_express as px
import plotly.graph_objects as go

from baseFunctions import CHE362

uReg = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
uReg.define('dollars = 1')
q  = uReg.Quantity

class DP(CHE362):

  def __init__(self) -> None:
    pd.DataFrame().to_csv('tables/proj_export.csv')

    ### Heat Transfer Area Constants ###
    self.hVap = q(30.75, 'kJ/mol') # heat of vaporization 
    t_CWI = q(30, 'degC') # cooling water in: condenser
    t_CWO = q(45, 'degC') # cooling water out: condenser
    t_D = q(80, 'degC') # temperature of distillate stream boiling point
    t_B = q(136, 'degC') # temperature of bottom stream boiling point
    t_ST = q(185.477, 'degC') # temperature of steam

    self.hdUnit = 'GJ/hr' # Heat duty import units

    self.dT_Ln = (
              ((t_D - t_CWI) - (t_D - t_CWO)) 
              / log(
                    (t_D - t_CWI) / 
                    (t_D - t_CWO)
                  )
                  ).to('delta_degF')
    
    dT = (t_ST - t_B).to('delta_degC')
    if dT.magnitude > 30:
      dT = q(30, 'delta_degC') 

    self.dT = dT.to('delta_degF')

    ### Economic Constants ###
    self.fP, self.fM, self.eta = 1, 1, .75
    self.traySpacing = q(2, 'ft')

    self.df = pd.read_csv('tables/proj_import.csv')


  def HeatTransferArea(self, hD, condenser: bool = True):

    muOc = q(150, 'Btu/(h * ft**2 * degF)')
    muOb = q(200, 'Btu/(h * ft**2 * degF)')

    if condenser:
      return (q(hD, self.hdUnit) / (muOc * self.dT_Ln)).to('m**2')
    else:
      return (q(hD, self.hdUnit) / (muOb * self.dT)).to('m**2')
  
  def main(self):

    fdf = pd.DataFrame()
    for i, x in self.df.iterrows():
      nTrays = (x['idealTrays'] - 1) / self.eta

      areaC = DP.HeatTransferArea(self, hD=abs(x['Qc'])).to('m**2').magnitude
      areaB = DP.HeatTransferArea(self, hD=x['Qb'], condenser=False).to('m**2').magnitude
      area = (pi * (q(x["diameter"], 'ft') / 2)**2).to('m**2').magnitude

      h_ = (q(10, 'ft') + nTrays * self.traySpacing).to('ft')
      frac, whole = modf(h_.magnitude)
      if frac >= 0.5:
        h_ = round(h_, 0)
      elif frac < 0.5:
        h_ = q(whole + .5, 'ft')

      d = q(x["diameter"], 'ft')
      vol = (pi * (d / 2)**2 * h_).to('m**3').magnitude

      if nTrays > 20:
        fQ = 1
      else:
        fQ = 10**(.4771 + .08516 * log(nTrays, 10) - .3473 * log(nTrays, 10)**2)

      cShell = 10**(3.4974 + .4485 * log(vol, 10) + .1074 * log(vol, 10)**2) * \
                  (2.25 + (1.82 * self.fP * self.fM))

      cTrays = 10**(3.3322 + .4838 * log(area, 10) + .3434 * log(area, 10)**2) * \
                  fQ * nTrays
      
      cCondenser = 10**(4.3247 - .303 * log(areaC, 10) + .1634 * log(areaC, 10)**2) * \
                  (1.63 + (1.66 * self.fM * self.fP))
      
      cBoiler = 10**(4.4646 - .5277 * log(areaB, 10) + .3955 * log(areaB, 10)**2) * \
                  (1.63 + (1.66 * self.fM * self.fP))
      
      cSteam = q(abs(x["Qc"]), self.hdUnit) * q(14.83, 'dollars/GJ') * q(350, 'day').to('h')
      cWater = q(x["Qb"], self.hdUnit) * q(.354, 'dollars/GJ') * q(350, 'day').to('h')

      capCost = cShell + cTrays + cCondenser + cBoiler
      opCost = cSteam + cWater

      eAOC = opCost + (capCost / 5)
      steamAOC = cSteam / eAOC

      y = {}
      y["ideal stages"] = [x["idealTrays"]]
      y["real trays"] = [nTrays]
      y["reflux ratio"] = [x['reflux_ratio']]
      y["column height (ft)"] =  [h_.magnitude]
      y["condenser heat duty (GJ/h)"] = [x["Qc"]]
      y["boiler heat duty (GJ/h)"] = [x["Qb"]]
      y["condenser area (m**2)"] = [areaC]
      y["boiler area (m**2)"] = [areaB]
      y["shell cost"] = [cShell]
      y["tray cost"] = [cTrays]
      y["condenser cost"] = [cCondenser]
      y["reboiler cost"] = [cBoiler]
      y["annual steam cost"] = [cSteam.magnitude]
      y["annual cooling water cost"] = [cWater.magnitude]
      y["total capital cost"] = [capCost]
      y["annual utility cost"] = [opCost.magnitude]
      y["EAOC"] = [eAOC.magnitude]
      y["EAOC steam"] = [steamAOC.magnitude]

      fdf = pd.concat([pd.DataFrame(y), fdf], axis=0, ignore_index=True)
      
    fdf = fdf.round(2).sort_values()
    fdf.round(2).to_csv('tables/proj_export.csv')
    print(fdf)


DP().main()