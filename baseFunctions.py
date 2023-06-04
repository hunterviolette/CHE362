from sympy import Eq, solve, symbols
import sympy.abc as cnst
import pint

uReg = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
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
  '''
  mass flux = mass transfer coefficent * concentration driving force
  mf = mtc * c

  low mass transfer rates (dilute mixture of A):
    n_A = k_x * delta_xA
      WHERE k_x units are length/time

  high mass transfer rates (diffusion of A through non-diffusing B)
    n_A = k_x * log((1 - x_A2) / (1 - x_A1))
      where k_x units are mols/cm**2

  other forms of mtc:
    Gases: 
      n_A = k_p * delta_pA = k_y * delta_yA = k_c * delta_cA = k_c * c * delta_yA
        pA: partial pressure of A

    Liquids:
      n_A = k_x * delta_xA = k_c * delta_cA

  converting b/w forms of mtc:
    Gases:
      low mass transfer:
        k_y = k_c * (p / (R * t)) = k_p * P
          k_y: gas mtc
          p: pressure
      
      high mass transfer:
        k_yy = k_G * (delta_pB / log(pB2 / pB1))
        k_yy = (k_y * (delta_pB / log(pB2 / pB1))) / p
        k_yy = k_c * ((delta_pB / log(pB2 / pB1) / (R * t))
          
    Liquids:
      low mass transfer:
        k_x = k_c * c
          # k_x: liquid mtc
          # c: molar concentration (mols/volume)

      high mass transfer:
        k_xx = k_x * (delta_xB / log(xB2 / xB1)
        k_xx = k_L * (delta_xB / log(xB2 / xB1) * c

  Re = (L * rho * v) / mu
    # L: length, meter
    # rho: density, kg/m**3
    # v: velocity, m/s
    # mu: viscosity, kg/(m*s)
  
  Sh_prime = (k_prime * L) / (c * d_AB)
    # d_AB: diffusivity m**2/s

  Sh = (k * L) / (c * d_AB) = (k_c * L) / d_AB

  Sc = mu / (rho * d_AB)
  
  Le = Pr / Sc
  
  '''

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

    print(f'Mass transfer coefficient: {soln[0]}')
    return soln
  
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
  def HenrysLaw(xI: float, # mole fraction of x at interface
                yI: float, # mole fraction of y at interface
                He, # Henrys law coefficient (pressure)
                p # Total pressure
              ):
    
    # ABSORPTION
    # Gasses dissolved in water and diulte, and non-ideal solutions
    return Eq((He / p) * xI, yI) # cnst.[x,y] mole fractions at interface 

  def RaoultsLaw(xI: float, # mole fraction of x at interface
                yI: float, # mole fraction of y at interface
                vpA, # vapor pressure of A
                p, # Total Pressure
                gamma = 1, # Activity coefficient methods
                idealSolution: bool = True
              ):
    
    # DISTILLATION
    # Ideal solutions
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
    
    if type(eq2) != 'sympy.core.relational.Equality':
      raise Exception('Equation 2 is not the right class') 
    
    xI, yI = symbols('xI'), symbols('yI')
    soln = solve((
                  Eq(kY * (yB - yI), kX * (xI - xB)),
                  eq2),
                (xI, yI))
    print(soln)
    return soln

  @staticmethod
  def FluxBetweenPhases(xxI: float, # mole frac of phase at interface
                        xxB: float, # mole frac of phase at bulk
                        k, # mtc of phase
                        phase: str = 'x' 
                      ):
    
    if phase not in ['x', 'y']:
      raise Exception("Invalid phase option: x or y")

    if phase == 'x':
      nA = k * (xxI - xxB)
      otherPhase = 'y'
    else:
      nA = k * (xxB - xxI)
      otherPhase = 'x'

    print(f'Flux from {otherPhase} to {phase}: {nA}')
    return nA

  @staticmethod
  def PhaseResistance(k,
                      bigK,
                      phase: str = 'x'
                    ):
    
    pR = (1 / k) / (1 / bigK)

    print(f'Phase resistance for {phase} phase: {round(pR, 3)}')
    return pR
  
  @staticmethod
  def SOLVE_OMTC(kX, # mtc in x phase
                kY, # mtc in y phase
                He, # Henrys law coefficient
                p, # Total pressure
                oK, # Overall MTC for a phase
                solveFor, # input var being solved
                phase: str = 'x'):
    
    if phase not in ['x', 'y']:
      raise Exception("Invalid phase option: x or y")
    
    m = He / p

    if phase == 'x':
      soln = solve(Eq((1 / kX) + (1 / (m * kY)), 1 / oK), solveFor)
    else:
      soln = solve(Eq((1 / kY) + (m / kX), 1 / oK), solveFor)
    
    print(f'Overall MTC for {phase} phase: {soln} ')
    return soln

  @staticmethod
  def VerifyPhaseResistances(xPR, yPR):
    
    one = xPR + yPR
    if one != 1:
      raise Exception("Resistaance between phases doesn't balance")
    else:
      print(f"Phase resistance sum is 1")

if __name__ == "__main__":
  b = MassTransfer()
  
  reg = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)
  reg.default_format = "~P"
  q  = reg.Quantity 


  re = b.ReynoldsNumber(
                    L = q(5, 'cm').to('m'),
                    rho = q(1000, 'kg/m**3'),
                    v = q(10, 'm/s'),
                    mu = q(.95, 'kg/(m*s)')
                  )
  print(re)

  '''
  d_AB = b.FullerDiffusivity(1, # pressure in atmospheres
                  299, # temperature in kelvin 
                  112.56, # molecular weight species A
                  28.97, # molecular weight species B 
                  109.65, # diffusion volume species A 
                  19.7, # diffusion volume species B
                  'Chlorobenzene', # name of species A
                  'Air' # name of species B
                )
  
  b.percentDifference(d_AB, .074)
  '''
  
  '''
  d_AB_1, t_1 = 0.783, 298
  d_AB_2_table, t_2 = 2.149, 533
  
  d_AB = b.SimplifiedFullerDiffusivity(d_AB_1, t_1, t_2)
  b.percentDifference(d_AB, d_AB_2_table)
  '''

  '''
  mVol_at_t = 59.47
  sTension_at_t = 19.22
  
  print(mVol_at_t * sTension_at_t**(1/4))
  '''