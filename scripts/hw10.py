from pint import UnitRegistry
from math import pi, log

from baseFunctions import CHE362

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

class HW10(CHE362):

  @staticmethod
  def C():
    
    t_CWI, t_CWO = q(30, 'degC'), q(45, 'degC')
    t_D, t_B, t_ST = q(110.6, 'degC'), q(136, 'degC'), q(185.477, 'degC')

    dT_Ln = (((t_D - t_CWI) - (t_D - t_CWO)) / log((t_D - t_CWI) / (t_D - t_CWO))).to('delta_degF')
    dT = (t_ST - t_B).to('delta_degF')
    if dT.magnitude > 50:
      dT = q(50, 'delta_degF') 

    q_C = q(6720.22, 'kW').to('GJ/hr') # from HW9
    q_B = q_C 

    muOc = q(150, 'Btu/(h * ft**2 * degF)')
    muOb = q(200, 'Btu/(h * ft**2 * degF)')

    a_C = (q_C / (muOc * dT_Ln)).to('m**2')
    a_B = (q_B / (muOb * dT)).to('m**2')

    print(
          "==== Part C ====",
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
  
  @staticmethod
  def D():
    #### Input vars
    r_, d_= 3.975, q(147.959*10**3, 'mol/hr') # from HW9.One

    mW = q(92.14, 'g/mol') # Molecular weight of toulene
    rhoV= (q(33.18, 'mol/m**3') * mW).to('kg/m**3') # rhoV Table 135 Yaws at Tb
    rhoL =q(.777, 'g/cm**3').to('kg/m**3')  # Molecular weight of liquid toulene

    sigma = 18.18 # dyne / cm
    percFlood, activeArea = .75, .8
    ####

    f_LV = (r_ / (r_ + 1)) * (rhoV / rhoL)**.5
    kV = q(10**(-.94506 - .70234 * log(f_LV, 10) - .22618 * log(f_LV, 10)**2), 'ft/s')

    uC = (kV * (sigma / 20)**.2 * ((rhoL - rhoV) / rhoV)**.5).to('ft/s')
    uO = uC * percFlood

    v_ = d_ * (r_ + 1)
    vDot = (v_ * mW / rhoV).to('ft**3/s')

    area = (vDot / uO).to('ft**2')
    actualArea = area / activeArea
    diameter = (4 * actualArea / pi)**.5

    print( 
        f"==== Part D ====",
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

  @staticmethod
  def One():
    HW10.C()
    HW10.D()

HW10.One()