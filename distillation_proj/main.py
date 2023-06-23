import pint
import pandas as pd
import plotly_express as px

from math import pi, log, modf
from plotly.offline import plot as save_plot
from plotly.subplots import make_subplots

from sys import path
path.append('..')
from scripts.baseFunctions import CHE362

uReg = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
uReg.define('dollars = 1')
q  = uReg.Quantity

class DP(CHE362):

  def __init__(self, steamTemp: float, steamCost: float, steamType: str) -> None:

    self.steamType = steamType
    ### Heat Transfer Area Constants ###
    self.hVap = q(30.75, 'kJ/mol') # heat of vaporization 
    t_CWI = q(30, 'degC') # cooling water in: condenser
    t_CWO = q(45, 'degC') # cooling water out: condenser
    t_D = q(80.1, 'degC') # temperature of distillate stream boiling point
    t_B = q(136, 'degC') # temperature of bottom stream boiling point
    t_ST = q(steamTemp, 'degC') # temperature of steam

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
    self.steamCost = steamCost
    self.fP, self.fM, self.eta = 1, 1, .75
    self.traySpacing = q(2, 'ft')

    df = pd.read_csv('main_raw.csv')
    self.df = df[df.converged == True] 

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
      
      frac, whole = modf(x['idealStages'] / self.eta)
      if frac >= 0.5:
        nTrays = whole + 1
      elif (frac < 0.5) and (frac > .0001) :
        nTrays = whole + 1
      else:
        nTrays = whole


      h_ = (q(10, 'ft') + (nTrays - 1) * self.traySpacing).to('ft')
      frac, whole = modf(h_.magnitude)
      if frac >= 0.5:
        h_ = round(h_, 0)
      elif (frac < 0.5) and (frac > .01) :
        h_ = q(whole + .5, 'ft')

      frac, whole = modf(x["diameter"])
      if frac >= 0.5:
        diameter = q(round(x["diameter"], 0), 'ft')
      elif (frac < 0.5) and (frac > .01) :
        diameter = q(whole + .5, 'ft')

      areaC = DP.HeatTransferArea(self, hD=abs(x['qC'])).to('m**2').magnitude
      areaB = DP.HeatTransferArea(self, hD=x['qB'], condenser=False).to('m**2').magnitude
      area = (pi * (diameter / 2)**2).to('m**2').magnitude

      vol = (pi * (diameter / 2)**2 * h_).to('m**3').magnitude

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
      
      cSteam = q(x["qB"], self.hdUnit) * q(self.steamCost, 'dollars/GJ') * q(350, 'day').to('h')
      cWater = q(abs(x["qC"]), self.hdUnit) * q(.354, 'dollars/GJ') * q(350, 'day').to('h')

      capCost = cShell + cTrays + cCondenser + cBoiler
      opCost = cSteam + cWater

      eAOC = opCost + (capCost / 5)
      steamAOC = cSteam / eAOC

      y = {}
      y["steam"] = [self.steamType]
      y["ideal stages"] = [x["idealStages"]]
      y["feed stage"] = [x["feedTray"]]
      y["real trays"] = [nTrays]
      y["diameter (ft)"] = [diameter.magnitude]
      y["reflux ratio"] = [x['refluxRatio']]
      y["column height (ft)"] =  [h_.magnitude]
      y["condenser heat duty (GJ/h)"] = [x["qC"]]
      y["boiler heat duty (GJ/h)"] = [x["qB"]]
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
    return fdf 
  
class DP_ANALYSIS():

  def __init__(self) -> None:
    pass
  
  def MultiPressure(self):
    fdf = pd.DataFrame()
    # (50 psig: 147.579 degC, 14.05), ( 150 psig 185.477, 14.83), (600 psig: 253.808 degC, 17.70)
    tasks = [[147.579, 14.05, '50psi'], [185.477, 14.83, '150psi'], [253.808, 17.7, '600psi']] 
    for task in tasks:
      fdf = pd.concat([DP(task[0], task[1], task[2]).main(), fdf], axis=0, ignore_index=True)

    fdf["annualized capX"] = fdf["total capital cost"] / 5
    self.df = fdf.round(3).sort_values(by='reflux ratio')
    self.df.to_csv('main_export.csv')
  
  def FigGenerator(self, xVal:str, yVal: str):
    fig = px.scatter(self.df, 
              self.df[xVal], 
              self.df[yVal],
              color=self.df["steam"],
              hover_data=["ideal stages", "feed stage", "real trays", "column height (ft)"]
            )

    fig.add_hline(y=self.df[yVal].min(), 
                  line_dash="dot", 
                  annotation_text=f"min {yVal}: {round(self.df[yVal].min(),0)}", 
                  annotation_font_color="white"
                )
    fig.add_hline(y=self.df[yVal].mean(), 
                  line_dash="dot", 
                  annotation_text=f"mean {yVal}: {round(self.df[yVal].mean(),0)}", 
                  annotation_font_color="white"
                )
    fig.add_hline(y=self.df[yVal].max(), 
                  line_dash="dot", 
                  annotation_text=f"max {yVal}: {round(self.df[yVal].max(),0)}", 
                  annotation_font_color="white"
                )
    return fig

  def Graphs(self):
    DP_ANALYSIS.MultiPressure(self)

    figures = [
      DP_ANALYSIS.FigGenerator(self, 
                              'reflux ratio',
                              'EAOC', 
                            ), 
      DP_ANALYSIS.FigGenerator(self, 
                              'reflux ratio', 
                              'annual utility cost',
                            ), 
      DP_ANALYSIS.FigGenerator(self, 
                              'reflux ratio', 
                              'annualized capX',
                            ), 
    ]

    fig = make_subplots(
        rows=len(figures), 
        cols=1, 
        shared_xaxes=True,
        vertical_spacing=.1, 
        subplot_titles=[
              'Reflux Ratio (X) Vs. EAOC (Y)',
              'Reflux Ratio (X) Vs. Annualized Utilites Cost (Y)',
              'Reflux Ratio (X) Vs. Annualized Capital Cost (Y)',
        ]) 
    
    for i, figure in enumerate(figures):
        for trace in range(len(figure["data"])):
            fig.append_trace(figure["data"][trace], row=i+1, col=1)

    fig.update_layout(hovermode="x unified", template='plotly_dark', height=2400)

    save_plot(fig, filename='main_figures.html')

DP_ANALYSIS().Graphs()
