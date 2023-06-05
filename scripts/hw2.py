
from baseFunctions import Diffusion

from pint import UnitRegistry
uReg = UnitRegistry(autoconvert_offset_to_baseunit = True)
uReg.default_format = "~P"
Q  = uReg.Quantity

class HomeworkTwo(Diffusion):
  
  @staticmethod
  def One():
    '''
    #1 a) Use Fuller’s method to estimate the diffusivity of benzene in nitrogen at 311.3 K and 1
    atm. Compare your estimate to the value reported in Table 3.2 in the text. Report the %
    difference (% difference = (estimated value – table value)/table value * 100%) between your
    estimated value and the reported value.
    '''
    print('=============================')
    print('Question #1a')
    d_AB = HomeworkTwo.FullerDiffusivity(1, # pressure in atmospheres
                      311.3, # temperature in kelvin 
                      78.114, # molecular weight species A
                      28.013, # molecular weight species B 
                      15.9 * 6 + 2.31 * 6 - 18.3, # diffusion volume species A 
                      18.5, # diffusion volume species B
                      'Benzene', # name of species A
                      'Nitrogen' # name of species B
                    )

    HomeworkTwo.percentDifference(d_AB, 0.102)

    '''
    b) Use Fuller’s method to estimate the diffusivity of cyclohexane in hydrogen at 288.6 and 1
    atm. Compare your estimate to the value reported in Table 3.2 in the text. Report the %
    difference (% difference = (estimated value – table value)/table value * 100%) between your
    estimated value and the reported value. % difference.
    '''
    print('=============================')
    print('Question #1b')
    
    d_AB = HomeworkTwo.FullerDiffusivity(1, # pressure in atmospheres
                      288.6, # temperature in kelvin 
                      84.15, # molecular weight species A
                      2.016, # molecular weight species B 
                      15.9 * 6 + 2.31 * 12, # diffusion volume species A 
                      6.2, # diffusion volume species B
                      'Cyclohexane', # name of species A
                      'Hydrogen' # name of species B
                    )

    HomeworkTwo.percentDifference(d_AB, 0.319)

    '''
    c) Use the table value of the diffusivity of ammonia in hydrogen at 298K and 1 atm to estimate
    the diffusivity at 533K. Report the % difference (% difference = (estimated value – table
    value)/table value * 100%) between your estimated value and the reported value. % difference.
    '''
    print('=============================')
    print('Question #1c')
    
    d_AB_1, t_1 = .783, 298
    d_AB_2_table, t_2 = 2.149, 533
    
    d_AB = HomeworkTwo.SimplifiedFullerDiffusivity(d_AB_1, t_1, t_2)
    
    HomeworkTwo.percentDifference(d_AB, d_AB_2_table)

  @staticmethod
  def Two():
    '''
    #2 a) Estimate the liquid diffusivity of acetic acid (solute) in a dilute solution of acetone at 288K
    using the Wilke-Chang method and the Tyn and Calus method. Compare your estimate to the
    value reported in Table 3.4 in the text. Report the % difference (% difference = (estimated value
    – table value)/table value * 100%) between your estimated value and the reported value.
    '''
    
    # A = acetic acid (t_b 118 C), B = acetone (t_b 56.29 C)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    t, mu_b, mVol_A_normal = 288, 0.2379, 87.32
    
    print('=============================')
    print('Question #2a WilkeChangDiffusivity')
    
    d_AB = HomeworkTwo.WilkeChangDiffusivity(t, # temperature in K
                          58.080, # molecular weight of B in g/mol
                          mu_b, # viscosity of B in cp 
                          mVol_A_normal, # molar volume of A at normal boiling in cm**3/mol 
                          1.0 # association factor ??
                        )
    
    d_AB_Table = 2.92 * 10**-5
    HomeworkTwo.percentDifference(d_AB, d_AB_Table)

    print('=============================')
    print('Question #2a TynCalusDiffusivity')
    
    d_AB = HomeworkTwo.TynCalusDiffusivity(t, # abs temperature in kelvin
                              mVol_A_normal, # molar volume of A at normal boiling point in cm**3/mol
                              59.47, # molar volume of B at normal boiling point in cm**3/mol
                              17.98, # surface tension of A at temp, t in dyn/cm
                              19.22, # surface tension of B at temp, t in dyn/cm
                              mu_b, # viscosity of B in cp
                              'organic_acid', # [water, organic_acid, nonpolar_into_monohydroxy_alcohols] 
                            )
                            
    HomeworkTwo.percentDifference(d_AB, d_AB_Table)

    '''
    b) Estimate the liquid diffusivity of water (solute) in a dilute solution of ethanol (solvent) at
    298K using the Wilke-Chang method and the Tyn and Calus method. Compare your estimate to
    the value reported in Table 3.4 in the text. Report the % difference (% difference = (estimated
    value – table value)/table value * 100%) between your estimated value and the reported value. 
    '''
    print('=============================')
    print('Question #2b WilkeChangDiffusivity')
    
    # A = water (t_b 100 C), B = ethanol (t_b 78.29 C)
    t, mu_b, mVol_A_normal = 298, 0.4488, 18.93 
    
    
    d_AB = HomeworkTwo.WilkeChangDiffusivity(t, # temperature in K
                          46.08, # molecular weight of B in g/mol
                          mu_b, # viscosity of B in cp 
                          mVol_A_normal, # molar volume of A at normal boiling in cm**3/mol 
                          1.5 # association factor
                        )
    d_AB_Table = 1.26 * 10**-5
    HomeworkTwo.percentDifference(d_AB, d_AB_Table)
    
    print('=============================')
    print('Question #2b TynCalusDiffusivity')
    d_AB = HomeworkTwo.TynCalusDiffusivity(t, # abs temperature in kelvin
                              mVol_A_normal, # molar volume of A at normal boiling point in cm**3/mol
                              62.62, # molar volume of B at normal boiling point in cm**3/mol
                              57.09, # surface tension of A at temp, t in dyn/cm
                              17.57, # surface tension of B at temp, t in dyn/cm
                              mu_b, # viscosity of B in cp
                              'water' # [water, organic_acid, nonpolar_into_monohydroxy_alcohols] 
                            )
    
    HomeworkTwo.percentDifference(d_AB, d_AB_Table)

d = HomeworkTwo()
d.One()
d.Two()