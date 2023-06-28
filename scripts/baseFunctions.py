from sympy import Eq, solve, symbols
from sympy.physics.units.util import convert_to
from sympy.physics.units import mol, m, second, hour
import sympy.abc as cnst
import pint
from math import log, modf, pi

from numpy import arange
from numpy.random import random

import plotly.graph_objects as go
import plotly_express as px
import pandas as pd

uReg = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
uReg.define('lbmol = 1 * mol')

q  = uReg.Quantity

class Diffusion():

  @staticmethod
  def percentDifference(x:float, x2:float):
      value = ((x - x2)/x2)*100
      print(f"Percent difference is {round(value,14)} [x, x_ref] [{round(x,14)}, {round(x2,14)}]")
      return value

  @staticmethod
  def FullerDiffusivity(
                        p: float, # pressure in atmospheres
                        t: float, # temperature in kelvin 
                        molW1: float, # molecular weight species A 
                        molW2: float, # molecular weight species B 
                        dVol1: float, # diffusion volume species A 
                        dVol2: float, # diffusion volume species B
                        specA: str = 'A', # name of species A
                        specB: str = 'B', # name of species B
                      ): 

    m_AB = 2 / ((1 / molW1) + (1 / molW2))
    d_AB = (0.00143 * t**1.75) / (p * m_AB**(1/2) * (dVol1**(1/3) + dVol2**(1/3))**2 )
    
    print(f'mole_weight_ratios_{specA}_{specB}: {round(m_AB,2)}')
    print(f'FullerDiffusivity{specA}_{specB}: {q(round(d_AB, 6), "cm**2/s")}')

    return d_AB

  @staticmethod
  def SimplifiedFullerDiffusivity(
                        d_AB_1: float, 
                        t_1: float,
                        t_2: float,
                      ):
    d_AB = d_AB_1 * (t_2 / t_1)**(1.75)
    
    print(f'simple_fuller_diffusivity: {q(round(d_AB, 6), "cm**2/s")}')
    return d_AB
  
  @staticmethod
  def WilkeChangDiffusivity( # A:solute, B:solvent
        t: float, # temperature in K
        mW_B: float, # molecular weight of B in g/mol
        mu_B: float, # viscosity of B in cp
        mVol_A_normal: float, # molar volume of A at normal boiling in cm**3/mol 
        aF: float = 1.5 # association factor {2.6: water, 1.9: methanol, 1.5: ethanol, 1.0: unassociated}
      ):
    
    d_AB = (7.4 * 10**-8 * (aF*mW_B)**(1/2) * t) / (mu_B * mVol_A_normal**0.6)
    
    print(f'WilkeChangDiffusivity: {q(round(d_AB, 6), "cm**2/s")}')
    return d_AB

  @staticmethod
  def TynCalusDiffusivity( # A:solute, B:solvent
        t: float, # abs temperature in kelvin
        mVol_A_normal: float, # molar volume at normal boiling point in cm**3/mol
        mVol_B_normal: float, # molar volume at normal boiling point in cm**3/mol
        sTension_A_t: float, # surface tension of A at temp, t in dyn/cm
        sTension_B_t: float, # surface tension of B at temp, t in dyn/cm
        mu_B: float, # viscosity of B in cp
        solute: str # [water, organicAcid, nonpolar_into_monohydroxy_alcohols] 
      ):
    
    if mu_B > 30:
      raise Exception('mu_B > 30 cp')
    
    def Parachor(mVol_at_t: float, sTension_at_t: float):
      return (mVol_at_t * sTension_at_t**(1/4))
    
    if solute in ['water', 'organic_acid']:
      #2. If solute is water, double mVol_A and para_A
      #3. If solute is organic acid and solvent is other than water, methanol, or butanol, double mVol_A and para_A

      para_A = Parachor(mVol_A_normal, sTension_A_t) * 2
      para_B = Parachor(mVol_B_normal, sTension_B_t)

      d_AB = ( 
              8.93 * 10**-8 * 
              (mVol_A_normal * 2 / mVol_B_normal**2)**(1/6) * 
              (para_B / para_A)**.6 * t / mu_B
            )    

    elif solute == 'nonpolar_into_monohydroxy_alcohols':
      # For nonpolar solutes diffusing into monohydroxy alcohols, multiply mVol_B and para_B by 8*mu_B

      para_A = Parachor(mVol_A_normal, sTension_A_t)
      para_B = Parachor(mVol_B_normal, sTension_B_t) * 8 * mu_B
      
      d_AB = (8.93 * 10**-8 * 
              (mVol_A_normal / (mVol_B_normal**2 * 8 * mu_B))**(1/6) * 
              (para_B / para_A)**.6 * t / mu_B
            )
    
    else:
      para_A = Parachor(mVol_A_normal, sTension_A_t)
      para_B = Parachor(mVol_B_normal, sTension_B_t) 

      d_AB = ( 
              8.93 * 10**-8 * 
              (mVol_A_normal / mVol_B_normal**2)**(1/6) * 
              (Parachor(mVol_A_normal, sTension_A_t) / Parachor(mVol_B_normal, sTension_B_t)) *
              (para_B / para_A)**.6 * t / mu_B
            )
    
    print(f'para_A: {round(para_A, 4)}', f'para_B: {round(para_B, 4)}')
    print(f'TynCalusDiffusivity: {q(round(d_AB, 6), "cm**2/s")}')
    return d_AB

class MassTransfer():

  @staticmethod
  def ReynoldsNumber(L: pint.Quantity,
                    rho: pint.Quantity,
                    v: pint.Quantity,
                    mu: pint.Quantity
                  ):

    Re = (L * rho * v) / mu
    print(f'Reynolds Number: {Re}')
    return Re
  
  @staticmethod
  def SherwoodSolver(Sh: float, # Sherwood Number
                    L, # characteristic length
                    c, # molarity: mol/volume
                    d_AB, # diffusivity: area/time
                    k = cnst.k, # mtc: mol/(area*time)
                    solveFor = cnst.k, # variable to solve for
                    low_mass_transfer: bool = False, # low mass transfer rates
                  ):
    
    if low_mass_transfer:
      soln = solve(Eq((k  * L ) / d_AB, Sh), solveFor)
    else:
      soln = solve(Eq((k * L) / (c * d_AB), Sh), solveFor)

    print(f'Mass transfer coefficient: {convert_to(soln[0], mol / (m**2 * second))}')
    return soln[0]
  
  @staticmethod
  def SchmidtNumber(mu: pint.Quantity,
                    rho: pint.Quantity,
                    d_AB: pint.Quantity,
                  ):
    
    Sc = mu / (rho * d_AB)
    print(f'Schmidt Number: {Sc}')
    return Sc

  @staticmethod
  def LewisNumber(Pr: float,
                Sc: float
              ):
    
    return Pr / Sc
  
  @staticmethod
  def HenrysLaw(xI, # mole fraction of x at interface
                yI, # mole fraction of y at interface
                He, # Henrys law coefficient (pressure)
                p # Total pressure
              ):
    
    # ABSORPTION
    # Gasses dissolved in water and diulte, and non-ideal solutions
    print('====  Henrys law to solve mole fractions at interface ====')
    return Eq((He / p) * xI, yI) # cnst.[x,y] mole fractions at interface 

  @staticmethod
  def RaoultsLaw(xI: float, # mole fraction of x at interface
                yI: float, # mole fraction of y at interface
                vpA, # vapor pressure of A
                p, # Total Pressure
                gamma: float = 1, # Activity coefficient methods
                idealSolution: bool = False
              ):
    
    # DISTILLATION
    # Ideal solutions
    print('==== Raoults law to solve mole fractions at interface ====')
    if idealSolution:
      eq = Eq((vpA / p) * xI, yI)
    else:
      if gamma == 1:
        raise Exception('Solution is non-Ideal and gamma is default value')
      else:
        eq = Eq(((gamma * vpA) / p) * xI, yI)

    return eq

  @staticmethod
  def Solve_MoleFracInterface(eq2,
                              xB: float,
                              yB: float,
                              kX,
                              kY
                            ):
    
    xI, yI = symbols('xI'), symbols('yI')

    soln = solve((
                  Eq(kY * (yB - yI), kX * (xI - xB)),
                  eq2
                ),
                (xI, yI))
    
    xI, yI = soln[xI], soln[yI] 
    print(f'Mole fraction of [x,y] at the interface [{xI}, {yI}]')
    return [xI, yI]

  @staticmethod
  def FluxBetweenPhases(xB: float,
                        xI: float,
                        yB: float,
                        yI: float,
                        kX,
                        kY
                      ):
    res = []
    for phase in ['x', 'y']:

      if phase == 'x':
        nA = kX * (xI - xB)
        otherPhase = 'y'
      else:
        nA = kY * (yB - yI)
        otherPhase = 'x'

      print(f'Flux from {otherPhase} to {phase} phase: {nA}')
      res.append(nA)

    return res[0], res[1]

  @staticmethod
  def Slope_OMTC(
                He_or_vpA, # Henry coefficient (henrys) or vapor pressure of A (raoults)
                p, # total pressure
                gamma: float = 1.0, # Activity coefficient            
              ):

    m = (He_or_vpA / p) * gamma 
    print(f"Slope, m: {round(m, 4)} ")
    return m

  @staticmethod
  def SOLVE_OVERALL_MTC(kX, # mtc in x phase
                kY, # mtc in y phase
                m, # slope 
                oK, # Overall MTC for a phase
                solveFor, # input var being solved
              ):
    
    res = []
    for phase in ['x', 'y']:
      if phase == 'x':
        soln = solve(Eq((1 / kX) + (1 / (m * kY)), 1 / oK), solveFor)
      else:
        soln = solve(Eq((1 / kY) + (m / kX), 1 / oK), solveFor)
      
      print(f'Overall MTC for {phase} phase: {soln[0]} ')
      res.append(soln[0])

    return res[0], res[1]

  @staticmethod
  def PhaseResistances(
                    kX,
                    bigKX,
                    kY,
                    bigKY,
                    ):
    
    res = []
    for phase in ['x', 'y']:
      if phase == 'x':
        k, bigK = kX, bigKX
      else:
        k, bigK = kY, bigKY

      pR = (1 / k) / (1 / bigK)

      print(f'Phase resistance for {phase} phase: {round(pR, 8)}')
      res.append(pR)

    one = res[0] + res[1]
    if round(one, 4) != 1:
      raise Exception(f"Total resistance does not add to one: {one}")
    else:
      print(f"Total phase resistance sums to: {one}")

class SeparationProcesses():

  @staticmethod
  def CapXY(value: float):
    return value / (1 - value)
  
  @staticmethod
  def Solve_EquilibRelationship(vp, p, xn, yn, solveFor, gamma: float = 1):
    return solve(Eq(((gamma * vp) / p) * xn, yn), solveFor)[0]
    
  @staticmethod
  def Solve_MaterialBal_Streams(stream: str,
                X1, X2, Y1, Y2, Ls, Vs, solveFor, ndigits: int = 3):

    vars = []
    strNames = ["X0", "Xn", "Y0", "Yn", "Ls", "Vs"]
    for i, x in enumerate([X1, X2, Y1, Y2, Ls, Vs]):
      try: 
          vars.append(f"{strNames[i]}: {round(x, ndigits)}")
      except:
          vars.append(f"{strNames[i]}: {x}")

    print(vars)

    if stream == 'cocurrent':
      res = solve(
                  Eq((Y2 - Y1) / (X2 - X1), -Ls / Vs),
                  solveFor
                )[0]
    
    elif stream == 'countercurrent':
      res = solve(Eq((Y2 - Y1) / (X2 - X1), -Ls / Vs), solveFor)[0]
      #res = solve(Eq((X1 * Ls) + (Y2 * Vs), (X2 * Ls) + (Y1 * Vs)))[0]

    elif stream == 'packedcolumn':
      res = solve(
                Eq(X))

    else:
      raise Exception('Stream type not defined')
      
    print(f"{stream} stream solved for {solveFor}: {res}")
    return res
  
  @staticmethod
  def PackedColumn_MaterialBal(X0, Ls, Yn_1, Vs, Xn, Y1, solveFor):
    '''
    return solve(
              Eq((X0 * Ls) + (Yn_1 * Vs), (Xn * Ls) + (Y1 * Vs)),
              solveFor
            )[0]
    '''
  
    return solve(
                Eq((Yn_1 - Y1) / (Xn - X0), Ls / Vs),
                solveFor
              )[0]
  
  @staticmethod
  def PackedColumn_Kremser(Ls, Vs, m, ynp1, x0, y1, xn, absorption: bool = True):
          
    print(f"Kremer Inputs: Ls: {Ls}, Vs: {Vs}, x0: {x0}, xn: {xn}, ynp1: {ynp1}, y1: {y1}")

    a_ = (Ls / Vs) / m
    if absorption:
      return log(
              ((ynp1 - (m * x0)) / (y1 - (m * x0))) *
              (1 - (1 / a_)) + (1 / a_)
            ) / log(a_)
    else:
      s_ = 1 / a_
      return log(
                ((x0 - (ynp1 / m)) / (xn - ynp1)) *
                (1 - (1 / s_)) + (1 / s_)
              ) / log(s_)

  @staticmethod
  def PackedColumn_Efficiency(x, epsilon, solveFor):
    return solve(
            Eq(log(1.597 - (0.199 * log(x)) - (0.0896 * log(x)**2)), epsilon),
            solveFor
          )[0]

  @staticmethod
  def KremserEquation(x1: float, y1: float, # CAPXY Values
                      m: float, Ls, Vs, 
                      absorption: bool = True, # False == stripping
                      xN: float = 0, # xN needed for stripping
                      yN: float = 0, # y1 only needed if absorbing
                      # Absorption (V -> L), Stripping (L -> V)
                      ndig: int = 4
                    ):
    
    
    a_ = (Ls / Vs) / m
    if absorption:
      print(f"Kremer Inputs: x1: {x1}, y1: {y1}, yN: {yN}, m: {m}, A: {a_}")
      process = 'absorption'
      numStages = log(
                      (y1 - (m * x1)) / (yN - (m * x1)) *
                      (1 - (1 / a_)) + (1 / a_)
                    ) / log(a_)
      
    else:
      a_ = 1 / a_
      print(f"Kremer Inputs: x1: {x1}, xN: {xN}, yN: {yN}, m: {m}, S: {a_}")
      process = 'stripping'
      

      numStages = log(
                      (x1 - (yN / m)) / (xN - (yN / m)) *
                      (1 - (1 / a_)) + (1 / a_)
                    ) / log(a_)
      
    print(f"Kremser Equation - [number_stages, process]: [{numStages}, {process}]")
    return numStages

class Util(SeparationProcesses):

  @staticmethod
  def EquilibriumDiagram(vp, p, 
                        xStep: float, yStep: float, # axis Formatting Figure tickmarks
                        xRange: list, yRange: list,  # axis Formatting Figure axis range
                        dataGen: list, # [startingValue, EndValue, StepValue]
                        gamma: float = 1
                      ):
    
    df = pd.DataFrame()
    for x in arange(dataGen[0], dataGen[1], dataGen[2]): 
      y = float(((vp * gamma) / p) * x)
      df = pd.concat([df, pd.DataFrame({"x":[x], "y":[y], 
                                        "X":[Util.CapXY(x)], "Y":[Util.CapXY(y)
                                      ]}
                                    )], ignore_index=True)
      
    fig = px.line(data_frame=df, x='X', y='Y')
    fig.update_layout(
              xaxis=dict(
                  tickmode='linear',
                  dtick=xStep,
                  range=xRange
                ),
              yaxis=dict(
                  tickmode='linear',
                  dtick=yStep,
                  range=yRange,
                ),
              template='plotly_dark'
            )
    
    return fig
  
  @staticmethod
  def PlotLines(fig, 
                lines: list = [], 
              ):
    if lines:
      for i, x in enumerate(lines):
        if len(x) == 1:
          mode = 'markers'
        else:
          mode = 'lines'
        
        fig.add_trace(go.Scatter(x=x['X'], y=x['Y'], name=f"{mode}{i}", mode=mode,
                      marker=dict(color=tuple(random(size=3) * 256) , size=12, 
                                  line=dict(width=1))))

    return fig

  @staticmethod
  def FortyFiveLine():
    df = pd.DataFrame()
    for x in arange(0, 1.02, .02): 
      y = float(x)
      df = pd.concat([df, pd.DataFrame({"x":[x], "y":[y]}
                                    )], ignore_index=True)

    return go.Scatter(
            x = df['x'], y = df['y'], 
            name=f"Forty Five Degree Line",
            marker = dict(color = tuple(random(size = 3) * 256), 
            size = 12, line = dict(width = 1))
          )

  @staticmethod
  def Generate_YX_Diagram(df: pd.DataFrame):
    fig = px.line(df, df['xVal'], df['yVal'])
    fig.add_trace(Util.FortyFiveLine())

    return fig.update_layout(
            xaxis=dict(
                tickmode='linear',
                dtick=.01,
                range=[0,1]
              ),
            yaxis=dict(
                tickmode='linear',
                dtick=.01,
                range=[0,1],
              ),
            template='plotly_dark'
          )
  
class Distillation():

  @staticmethod
  def Solve_DB(f_, xF, xD, xB):
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
    b_ = q(soln[b_]  * hour / mol, 'mol/h')
    d_ = q(soln[d_]  * hour / mol, 'mol/h')
    
    print(
          f"Bottoms flow rate: {b_}",
          f"Distillate flow rate: {d_}",
          sep='\n'
        )
    return b_, d_
  
  @staticmethod
  def Solve_Rmin(xD, zF, xx, satLiquid: bool = True):
    Rmin = symbols('R')
    if satLiquid:
      Rmin = solve( 
                    Eq(
                      Rmin / (Rmin + 1), 
                      (xD - xx) / (xD - zF)
                    ),
                    Rmin
                  )[0]
      print(f"Rmin satLiq: {Rmin}")

    else:
      Rmin = solve( 
                    Eq(
                      Rmin / (Rmin + 1), 
                      (xD - zF) / (xD - xx)
                    ),
                    Rmin
                  )[0] 
      print(f"Rmin satVap: {Rmin}")
    return Rmin
  
  @staticmethod
  def Distillation_Diameter(
                    mW: pint.Quantity,
                    rhoV: pint.Quantity,
                    rhoL: pint.Quantity,
                    sigma: float, # surface tension dyne/cm
                    percentFlood: float,
                    activeArea: float,
                    r_: float, 
                    d_: pint.Quantity
                  ):

    f_LV = (r_ / (r_ + 1)) * (rhoV / rhoL)**.5
    kV = q(10**(-.94506 - .70234 * log(f_LV, 10) - .22618 * log(f_LV, 10)**2), 'ft/s')

    uC = (kV * (sigma / 20)**.2 * ((rhoL - rhoV) / rhoV)**.5).to('ft/s')
    uO = uC * percentFlood

    v_ = d_ * (r_ + 1)
    vDot = (v_ * mW / rhoV).to('ft**3/s')

    area = (vDot / uO).to('ft**2')
    actualArea = area / activeArea
    diameter = (4 * actualArea / pi)**.5

    frac, whole = modf(diameter.magnitude)
    if frac >= 0.5:
      diameter = round(diameter, 0)
    elif (frac < 0.5) and (frac > .01):
      diameter = q(whole + .5, 'ft')
    else:
      diameter = round(diameter, 0)

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
        f"diameter: {diameter}", 
        sep='\n')
    
    return diameter

  @staticmethod
  def Real_Stages(idealStages: float, eta: float):
    frac, whole = modf(idealStages / eta)
    if frac >= 0.5:
      stagesReal = whole + 1
    elif (frac < 0.5) and (frac > .01) :
      stagesReal = whole + 1
    else:
      stagesReal = whole  
    return stagesReal
  
  @staticmethod
  def Vaporization_Slope(vapPercent):
    return (1 / vapPercent) -1
  
  @staticmethod
  def Vaporization_SlopeCheck(y2: float, y1: float, x2: float, x1: float):
    return (y2 - y1) / (x2 - x1)
  
class CHE362(Diffusion, MassTransfer, Util, Distillation):
  pass

if __name__ == "__main__":
  s = CHE362()
  '''
  s.Solve_MaterialBal_Streams(
              stream='countercurrent',
              X1 = 0,
              X2 = .074, 
              Y1 = s.CapXY(.001), 
              Y2 = s.CapXY(.02), 
              Ls = symbols('Ls'),
              Vs = 98 * 10**3 * mol / hour,
              solveFor = symbols('Ls')
            )
  
  '''

  s.KremserEquation(0, .02, 0, .001, .286, 9.785, 27.222)
