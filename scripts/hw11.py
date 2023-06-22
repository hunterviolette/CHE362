from pint import UnitRegistry
from math import pi, log

uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
q  = uReg.Quantity

class HW11():

  @staticmethod
  def Two():
    t_CWI, t_CWO = q(30, 'degC'), q(45, 'degC')
    t_D, t_B, t_ST = q(110.6, 'degC'), q(136, 'degC'), q(185.477, 'degC')

    dT_Ln = (((t_D - t_CWI) - (t_D - t_CWO)) / log((t_D - t_CWI) / (t_D - t_CWO))).to('delta_degF')
    dT = (t_ST - t_B).to('delta_degF')
    if dT.magnitude > 50:
      dT = q(50, 'delta_degF') 

    q_C = q(13.38695324118, 'GJ/hr')
    q_B = q(13.614113476152, 'GJ/hr')  

    muOc = q(150, 'Btu/(h * ft**2 * degF)')
    muOb = q(200, 'Btu/(h * ft**2 * degF)')

    a_C = (q_C / (muOc * dT_Ln)).to('m**2')
    a_B = (q_B / (muOb * dT)).to('m**2')

    print(
          "==== Problem 2 ====",
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
    
HW11.Two()