from pint import UnitRegistry

from baseFunctions import CHE362

import pandas as pd

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
uReg.define('kmol = 10**3 * mol')

q  = uReg.Quantity

class HW14(CHE362):

  @staticmethod
  def One():
    df = pd.read_csv('tables/hw13.csv')
    fig = HW14.Generate_YX_Diagram(df)

    xF, xN = .2, .01
    x_F, x_N = HW14.CapXY(xF), HW14.CapXY(xN)

    lines = []
    lines.append(pd.DataFrame({"X": [x_F], "Y":[0]}).astype(float))
    lines.append(pd.DataFrame({"X": [x_N], "Y":[0]}).astype(float))

    lines.append(pd.DataFrame({"X": [x_N, x_F], "Y":[0.00, .062 / 1.5]}).astype(float))

    fig = HW14.PlotLines(fig, lines).show()    

    Vs_min = ((.062 - 0) / (.25 - .01)) * q(100 * (1 - xF), 'kmol/h')
    print(f"Vs_min: {Vs_min}")

HW14.One()