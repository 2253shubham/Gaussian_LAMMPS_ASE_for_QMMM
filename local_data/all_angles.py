from __future__ import unicode_literals

######################################################################################
# k value in kcal/mol, eqlb_ang in degrees

def C_C_H():
    return dict(
        k_value = 37.50,
        eqlb_ang = 110.7,
    )


def H_C_C():
    return dict(
        k_value = 37.50,
        eqlb_ang = 110.7,
    )


def H_C_H():
    return dict(
        k_value = 58.35,
        eqlb_ang = 112.7,
    )


def H_Si_O(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.80,
    )


def O_Si_H(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 110,
    )


def H_Si_H(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.80,
    )


def H_Si_HL(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 110.0,
    )


def HL_Si_H(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 110.0,
    )


def Si_O_Si(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 148.0,
    )


def Si_HL_Si(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 148.0,
    )


def Si_HL_O(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 148.0,
    )


def O_HL_Si(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 148.0,
    )


def Si_O_HL(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 148.0,
    )


def HL_O_Si(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 148.0,
    )


def O_Si_O(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def O_Si_HL(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def HL_Si_O(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def HL_Si_HL(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def Os_Si_Os(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0, # 33.00,
        eqlb_ang = 107.80,
    )


def Os_Si_Oa(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0, # 33.00,
        eqlb_ang = 107.80,
    )


def Oa_Si_Os(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0, # 33.00,
        eqlb_ang = 107.80,
    )


def Si_Oa_C2(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def C2_Oa_Si(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def Al_Oa_C2(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def C2_Oa_Al(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def Si_HL_Os(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 148.0,
    )


def Os_HL_Si(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 148.0,
    )


def Os_Si_HL(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def HL_Si_Os(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def Oa_Al_Oa(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0,
        eqlb_ang = 107.80,
    )


def Oa_Si_HL(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def HL_Si_Oa(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 107.5,
    )


def Oa_H_C1(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0, # 33.00,
        eqlb_ang = 107.80,
    )


def C1_H_Oa(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0, # 33.00,
        eqlb_ang = 107.80,

    )


def Oa_H_C2(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0, # 33.00,
        eqlb_ang = 107.80,
    )


def C2_H_Oa(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0, # 33.00,
        eqlb_ang = 107.80,
    )


def Si_Os_Si(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0,
        eqlb_ang = 107.80,
    )


def Si_Oa_Al(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0,
        eqlb_ang = 107.80,
    )


def Al_Oa_Si(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 0,
        eqlb_ang = 107.80,
    )


def Si_Os_HL(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 148.0,
    )


def HL_Os_Si(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 148.0,
    )


def Si_Oa_H(): # framework atoms fixed
    return dict(
        k_value = 35,
        eqlb_ang = 115.5,
    )


def H_Oa_Si(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_ang = 104.6,
    )


def Al_Oa_H(): # framework atoms fixed (OPLSAA) #trying
    return dict(
        k_value = 35,
        eqlb_ang = 102.2,
    )


def H_Oa_Al(): # framework atoms fixed (OPLSAA)
    return dict(
        k_value = 35,
        eqlb_ang = 102.20,
    )


def C1_C2_C2(): # OPLSAA
    return dict(
        k_value = 58.35,
        eqlb_ang = 112.7,
    )


def C2_C2_C1(): # OPLSAA
    return dict(
        k_value = 58.35,
        eqlb_ang = 112.7,
    )


def Oa_C2_H(): # OPLSAA/AMBER
    return dict(
        k_value = 35,
        eqlb_ang = 109.5,
    )


def H_C2_Oa(): # OPLSAA/AMBER
    return dict(
        k_value = 35,
        eqlb_ang = 109.5,
    )


def Oa_C2_C2(): # OPLSAA/AMBER, TraPPE-UA (of ether)
    return dict(
        k_value = 50, # 99.96,
        eqlb_ang = 109.5, # 112.0,
    )


def C2_C2_Oa(): # OPLSAA/AMBER, TraPPE-UA (of ether)
    return dict(
        k_value = 50, # 99.96,
        eqlb_ang = 109.5, # 112.0,
    )


def Oa_C2_C1(): # OPLSAA/AMBER, TraPPE-UA (of ether)
    return dict(
        k_value = 50, # 99.96,
        eqlb_ang = 109.5, # 112.0,
    )


def C1_C2_Oa(): # OPLSAA/AMBER, TraPPE-UA (of ether)
    return dict(
        k_value = 50, # 99.96,
        eqlb_ang = 109.5, # 112.0,
    )


def C1_C2_H(): # OPLSAA
    return dict(
        k_value = 37.50,
        eqlb_ang = 110.7,
    )


def H_C2_C1(): # OPLSAA
    return dict(
        k_value = 37.50,
        eqlb_ang = 110.7,
    )


def C2_H_C2(): # TraPPE-EH (rigid), change according to cracking
    return dict(
        k_value = 0,
        eqlb_ang = 110.7,
    )


def C2_C2_C2(): # OPLSAA
    return dict(
        k_value = 58.35,
        eqlb_ang = 112.7,
    )


def C2_C1_H(): # OPLSAA
    return dict(
        k_value = 37.50,
        eqlb_ang = 110.7,
    )


def H_C1_C2(): # OPLSAA
    return dict(
        k_value = 37.50,
        eqlb_ang = 110.7,
    )


def C2_C2_H(): # OPLSAA
    return dict(
        k_value = 37.50,
        eqlb_ang = 110.7,
    )


def H_C2_C2(): # OPLSAA
    return dict(
        k_value = 37.50,
        eqlb_ang = 110.7,
    )


def C2_H_HL(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 110.7,
    )


def HL_H_C2(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 110.7,
    )


def C2_HL_H(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 110.7,
    )


def H_HL_C2(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 110.7,
    )


def C2_C2_HL(): # TraPPE-EH 
    return dict(
        k_value = 0,
        eqlb_ang = 110.7,
    )


def HL_C2_C2(): # TraPPE-EH 
    return dict(
        k_value = 0,
        eqlb_ang = 110.7,
    )


def C2_HL_C2(): # TraPPE-EH 
    return dict(
        k_value = 0,
        eqlb_ang = 110.7,
    )


def H_C1_H(): # OPLSAA
    return dict(
        k_value = 33.00,
        eqlb_ang = 107.8,
    )


def H_C2_H(): # OPLSAA
    return dict(
        k_value = 33.00,
        eqlb_ang = 107.8,
    )


def H_C2_HL(): # TraPPE-EH 
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def HL_C2_H(): # TraPPE-EH 
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Si_H_C2(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def C2_H_Si(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Os_H_C2(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def C2_H_Os(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Si_H_C1(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def C1_H_Si(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Os_H_C1(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def C1_H_Os(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Os_Si_H(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def H_Si_Os(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Si_Os_H(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def H_Os_Si(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Os_H_Si(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Si_H_Os(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Oa_Al_H(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def H_Al_Oa(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Al_H_C1(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def C1_H_Al(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Al_H_C2(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def C2_H_Al(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Oa_H_Si(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def Si_H_Oa(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def C1_H_C2(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )


def C2_H_C1(): # not feasible
    return dict(
        k_value = 0,
        eqlb_ang = 107.8,
    )