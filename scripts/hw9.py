from pint import UnitRegistry
from sympy.physics.units.util import convert_to
from sympy import Eq, symbols, solve
from sympy.physics.units import atm, mol, hour
from math import pi

import pandas as pd
import plotly_express as px
import plotly.graph_objects as go

from baseFunctions import CHE362

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

class HW9(CHE362):

  @staticmethod
  def One():
    f_ = 500*10**3 * mol/hour
    xF, xD, xB,  = .3, .99, .01
    yF, yD, yB = .46, .99, .01

    b_, d_ = symbols('B'), symbols('D')
    soln = solve((
                    Eq( 
                        xD * d_ + xB * b_, 
                        xF * f_
                      ),
                    Eq( 
                        d_ + b_, 
                        f_ 
                      )
                  ),
                  (b_, d_)
                )
    print(f"Flow rates: {soln}")
    b_, d_ = soln[b_], soln[d_]

    Rmin = symbols('R')
    Rmin = solve( 
                  Eq(
                    Rmin / (Rmin + 1), 
                    (xD - yF) / (xD - xF)
                  ),
                  Rmin
                )[0]
    r_ = Rmin * 1.2
    print(f"Rmin: {Rmin}, R: {r_}")
    
    fig = HW9.Generate_YX_Diagram(pd.read_csv('tables/hw9a.csv'))
    fig.add_vline(x=xF, line_dash="dot", annotation_text="Feed line", annotation_font_color="white")
    fig.show()
    
    lines = []
    lines.append(pd.DataFrame({"X": [0, xD], "Y":[xD / (r_ + 1), yD]}).astype(float))
    lines.append(pd.DataFrame({"X": [xB, xF], "Y":[yB, .438]}).astype(float))

    fig = HW9.PlotLines(fig, lines).show()

    hL, hV = q(.1537, 'kJ/mol'), q(33.02, 'kJ/mol')
    d_ = q(d_ * hour / mol, 'mol/hour')
    qC = (d_ * (r_ + 1) * (hV - hL))

    print(f"Part D: {qC.to('kW')}")
    
    '''
    df["TOLy"] = (((yD - (xD / (r_ + 1))) / (xD - 0)) * df["xVal"] + (xD / (r_ + 1))).astype('float')
    for i, x in df.iterrows():
      if x["xVal"] < xB:
        df.loc[i,'BOLy'] = 0
      elif x["xVal"] == xB:
        df.loc[i,'BOLy'] = yB
      else:
        df.loc[i,'BOLy'] = (((yF - yB) / (xF - xB)) * x["xVal"]) - yB
    
    print('====')
    stage, stages = 0, [] 
    while True:
      
      print(stage, stages)      
      if stage == 0:
        x, y = xD, df.loc[df['xVal'] == xD, 'TOLy'].values[0]
        xn = df.loc[df['yVal'] == y, 'xVal'].values[0]
        stages.append([xn, y, x])
        stage += 1
      elif stage > 0:


      if (x < xB) or (stage > 5):
        break
    
    fig.add_trace(go.Scatter(x=df['xVal'], y=df['TOLy']))
    fig.add_trace(go.Scatter(x=df['xVal'], y=df['BOLy']))
    fig.add_vline(x=.3, line_dash="dot", annotation_text="Feed line", annotation_font_color="white")

    fig.show()
    '''

  @staticmethod
  def Two():
    xF, xD, xB  = .2976, .99, .01
    yF, yD, yB = .17, .99, .01

    Rmin = symbols('R')
    Rmin = solve( 
                  Eq(
                    Rmin / (Rmin + 1), 
                    (xD - xF) / (xD - yF)
                  ),
                  Rmin
                )[0]
    r_ = Rmin * 1.2
    print(f"Rmin: {Rmin}, R: {r_}")

    fig = HW9.Generate_YX_Diagram(pd.read_csv('tables/hw9a.csv'))
    fig.add_hline(y=xF, line_dash="dot", annotation_text="Feed line", annotation_font_color="white")
    fig.show()

    lines = []
    lines.append(pd.DataFrame({"X": [0, xD], "Y":[xD / (r_ + 1), yD]}).astype(float))
    lines.append(pd.DataFrame({"X": [xB, .192], "Y":[yB, .298]}).astype(float))
    fig = HW9.PlotLines(fig, lines).show()

  @staticmethod
  def Three():
    f_ = 100*10**3 * mol/hour
    xF, xD, xB = .2, .65, .01

    b_, d_ = symbols('B'), symbols('D')
    soln = solve((
                    Eq( 
                        xD * d_ + xB * b_, 
                        xF * f_
                      ),
                    Eq( 
                        d_ + b_, 
                        f_ 
                      )
                  ),
                  (b_, d_)
                )
    print(f"Flow rates: {soln}")
    b_, d_ = soln[b_], soln[d_]
  
    fig = HW9.Generate_YX_Diagram(pd.read_csv('tables/hw9b.csv'))
    fig.add_hline(y=xF, line_dash="dot", annotation_text="Feed line", annotation_font_color="white")
    fig.show()
    
    yF, yD, yB = .424, .65, .01

    Rmin = symbols('R')
    Rmin = solve( 
                  Eq(
                    xD / (Rmin + 1), .39),
                  Rmin
                )[0]
    
    r_ = Rmin * 1.5
    print(f"Rmin: {Rmin}, R: {r_}")

    lines = []
    lines.append(pd.DataFrame({"X": [0, xD], "Y":[xD / (r_ + 1), yD]}).astype(float))
    lines.append(pd.DataFrame({"X": [xB, xF], "Y":[yB, yF]}).astype(float))
    fig = HW9.PlotLines(fig, lines).show()

HW9.One()

