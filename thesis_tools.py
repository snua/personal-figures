def long(conc,t):
    '''
    Takes 1D numpy array and backcalculates the longitudinal diffusion
    based on Taylor 1953
    '''
    upper_limit=0.9
    lower_limit=0.1
    taylor_constant = 3.625

    x_upper = np.min(np.where(conc<upper_limit)) / len(conc)
    x_lower = np.min(np.where(conc<lower_limit)) / len(conc)

    diff = (1/t)*((x_lower - x_upper)/taylor_constant)**2

    return diff

def frac_flow(Sw, corey_exp=2, muw=1, muo=1):
    '''
    Calculates fractional flow of water
    '''
    krw = (Sw)**corey_exp
    kro = (1-Sw)**corey_exp
    F = 1/(1+(kro/muo)*(muw/krw))

    return F