def frac_flow(Sw, corey_exp=2, muw=1, muo=1):
    '''
    Calculates fractional flow of water
    '''
    krw = (Sw)**corey_exp
    kro = (1-Sw)**corey_exp
    F = 1/(1+(kro/muo)*(muw/krw))
    
    return F