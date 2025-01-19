from __future__ import unicode_literals

######################################################################################
# k value in kcal/mol, eqlb_dist in angstroms

def C_C(): # OPLS AA
    return dict(
        k_value = 268,
        eqlb_dist = 1.53,
    )


def C1_C2(): # OPLS AA
    return dict(
        k_value = 268,
        eqlb_dist = 1.53,
    )


def C2_C1(): # OPLS AA
    return dict(
        k_value = 268,
        eqlb_dist = 1.53,
    )


def C2_C2(): # OPLS AA
    return dict(
        k_value = 268,
        eqlb_dist = 1.53,
    )


def C_H(): # OPLS AA
    return dict(
        k_value = 340,
        eqlb_dist = 1.09,
    )


def H_C(): # OPLS AA
    return dict(
        k_value = 340,
        eqlb_dist = 1.09,
    )


def C1_H(): # OPLS AA
    return dict(
        k_value = 340,
        eqlb_dist = 1.09,
    )


def C2_H(): # OPLS AA
    return dict(
        k_value = 340,
        eqlb_dist = 1.09,
    )


def H_C1(): # OPLS AA
    return dict(
        k_value = 340,
        eqlb_dist = 1.09,
    )


def H_C2(): # OPLS AA
    return dict(
        k_value = 340,
        eqlb_dist = 1.09,
    )


def H_Si(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_dist = -1.47,
    ) 


def Si_H(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_dist = -1.47,
    ) 


def Si_O(): # framework atoms fixed
    return dict(
        k_value = 340,
        eqlb_dist = 1.59,
    )    


def O_Si(): # framework atoms fixed
    return dict(
        k_value = 340,
        eqlb_dist = 1.59,
    )


def Si_HL(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_dist = 1.47,
    )    


def HL_Si(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_dist = 1.47,
    ) 


def O_HL(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_dist = 0.12,
    )


def Os_H(): # framework atom fixed
    return dict(
        k_value = 0,
        eqlb_dist = -0.960,
    )


def H_Os(): # framework atom fixed
    return dict(
        k_value = 0,
        eqlb_dist = -0.960,
    )


def HL_O(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_dist = 0.12,
    )   


def Os_Si(): # framework atoms fixed
    return dict(
        k_value = 0, # 340,
        eqlb_dist = 1.59,
    ) 


def Si_Os(): # framework atoms fixed
    return dict(
        k_value = 0, # 340,
        eqlb_dist = 1.59,
    ) 


def HL_Os(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_dist = 0.12,
    )   


def Os_HL(): # framework atoms fixed
    return dict(
        k_value = 0,
        eqlb_dist = 0.12,
    )   


def Si_Oa(): # framework atoms fixed
    return dict(
        k_value = 0, # 340,
        eqlb_dist = 1.59,
    )   


def Oa_Si(): # framework atoms fixed
    return dict(
        k_value = 0, # 340,
        eqlb_dist = 1.59,
    )  


def Al_Oa(): # framework atoms fixed
    return dict(
        k_value = 0, # 340,
        eqlb_dist = 1.59,
    )   


def Oa_Al(): # framework atoms fixed
    return dict(
        k_value = 0, # 340,
        eqlb_dist = 1.59,
    ) 


def H_Oa(): # AMBER/OPLSAA
    return dict(
        k_value = 450,
        eqlb_dist = 1.1,
    )   


def Oa_H(): # AMBER/OPLSAA
    return dict(
        k_value = 450,
        eqlb_dist = 1.1,
    ) 


def C2_Oa(): # AMBER/OPLSAA
    return dict(
        k_value = 320,
        eqlb_dist = 1.41,
    )


def Oa_C2(): # AMBER/OPLSAA
    return dict(
        k_value = 320,
        eqlb_dist = 1.41,
    )


def C2_HL(): # Link atoms 
    return dict(
        k_value = 0,
        eqlb_dist = 0.55,
    )


def HL_C2(): # Link atoms
    return dict(
        k_value = 0,
        eqlb_dist = 0.55,
    )


def H_HL(): # not feasible
    return dict(
        k_value = 0,
        eqlb_dist = -0.55,
    )


def HL_H(): # not feasible
    return dict(
        k_value = 0,
        eqlb_dist = -0.55,
    )


def H_Al(): # not feasible
    return dict(
        k_value = 0,
        eqlb_dist = -0.55,
    )


def Al_H(): # not feasible
    return dict(
        k_value = 0,
        eqlb_dist = -0.55,
    )


def H_H(): # not feasible
    return dict(
        k_value = 0,
        eqlb_dist = -0.55,
    )