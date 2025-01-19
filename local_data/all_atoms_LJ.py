from __future__ import unicode_literals

######################################################################################
# epsln in kcal/mol, sigma in angstroms

def C(): # OPLSAA
    return dict(
        epsln = 0.066,
        sigma = 3.5,
        q = -0.12
    )


def H(): # OPLSAA
    return dict(
        epsln = 0.030,
        sigma = 2.5,
        q = 0.06,
    )


def HL(): # dummy atom
    return dict(
        epsln = 0,
        sigma = 0,
        q = 0,
    )


#def O(): # 2011 Zimmermann
#    return dict(
#        epsln = 0.075,
#        sigma = 1.77,
#        q = 0,
#    )   


#def Si(): # 2011 Zimmermann
#    return dict(
#        epsln = 0.20,
#        sigma = 2.2,
#        q = 0,
#    )   


def Al(): # TraPPE-zeo  
    return dict(
        epsln = 0.044,
        sigma = 2.3,
        q = 0.99,
    )


def Oa(): # TraPPE-zeo [Si-[Oa]-Al]  
    return dict(
        epsln = 0.115,
        sigma = 3.6,
        q = -0.58
    )       


def Os(): # TraPPE-zeo [Si-[Os]-Si]  
    return dict(
        epsln = 0.105,
        sigma = 3.3,
        q = -0.75
    )   


def Si(): # TraPPE-zeo 
    return dict(
        epsln = 0.044,
        sigma = 2.3,
        q = 1.5,
    )     


#def H(): # TraPPE AA (changing to OPLSAA for sigma)
#    return dict(
#        epsln = 0.030,
#        sigma = 2.5,
#        q = 0,
#    )


#def C1(): # TraPPE AA [C]H3
#    return dict(
#        epsln = 0.008,
#        sigma = 3.3,
#        q = 0,
#    )


#def C2(): # TraPPE AA [C]H2
#    return dict(
#        epsln = 0.01,
#        sigma = 3.65,
#        q = 0,
#    )


def C1(): # OPLSAA [C]H3
    return dict(
        epsln = 0.066,
        sigma = 3.5,
        q = -0.18,
    )


def C2(): # OPLSAA [C]H2
    return dict(
        epsln = 0.066,
        sigma = 3.5,
        q = -0.12,
    )