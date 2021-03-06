import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import optimize
from scipy.special import erf, erfc
from scipy.integrate import simps

#plotting section starts---------------------------------------------------------------------
def plt_config(use_tex=False):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from cycler import cycler

    #mpl.style.use('classic') # change everything to the old matplotlib look
    mpl.rcParams['image.cmap'] = 'jet'

    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}

    plt.rc('font', **font)
    plt.rc('legend',fontsize='large') # using a named size
    plt.rc('lines', linewidth=2.5)
    plt.rc('axes', prop_cycle=(cycler('color', ['b', 'r', 'g', 'darkorange', 'y']) +
                           cycler('linestyle', ['-', '--', ':', '-.','-'])))
    plt.rc('axes', labelsize='large')

    #the option to use tex formating,
    #example:
    # plt.title(r'$\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!', fontsize=13)
    # plt.legend([r'\textit{voltage}'], loc=2, fontsize=13)
    if use_tex:
        plt.rc('text', usetex=True)
        plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'],
                                'monospace': ['Computer Modern Typewriter']})
    else:
        plt.rc('font',family='serif')


def plt_default_plot(style=1, nrows=1, ncols=1, figsize=(8,4)):
    '''
    Creates general charts. The parameters specified here are the ones thought to be useful as default parameters.
    All parameters can be overridden once called.

    Example usage:
        fig, ax = tt.plt_default_plot(nrows=2, figsize=(8,8))
        ax0, ax1 = ax.flatten()
        ax0.plot(np.linspace(0,100,100),np.linspace(0,100,100))
        ax1.scatter(np.linspace(0,100,100),np.linspace(0,100,100))
    '''

    import string
    enumerate_alphabet = list(string.ascii_lowercase)

    if style == 1:
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
        # ax0, ax1 = axes.flatten() # use this to assign each axis to separate variable
        if nrows == 1 and ncols == 1:
            axes.set_xlabel('xlabel')
            axes.set_ylabel('ylabel')

        elif (nrows == 1 and ncols > 1) or \
             (nrows > 1 and ncols == 1):
            for i, ax in enumerate(axes):
                ax.set_title('('+enumerate_alphabet[i]+')')
                ax.set_xlabel('xlabel')
                ax.set_ylabel('ylabel')

        else:
            for i, row in enumerate(axes):
                for j, ax in enumerate(row):
                #pass
                    ax.set_title('('+enumerate_alphabet[i*ncols+j]+')')
                    ax.set_xlabel('xlabel')
                    ax.set_ylabel('ylabel')

        plt.subplots_adjust(right=0.95,left=0.05,top=0.9,bottom=0.1,hspace=0.0,wspace=0.01)

    elif style == 2:
        # This is just an example on how to overlap the axes.
        # If overlapped axes required, use style 1 instead and use plt.subplot2grid as shown below
        fig, axes = plt.subplots(nrows=nrows*3, ncols=ncols*3, figsize=figsize)

        ax1 = plt.subplot2grid((nrows*3, nrows*3), (0, 0))
        ax2 = plt.subplot2grid((nrows*3, nrows*3), (0, 1), colspan=2)
        ax3 = plt.subplot2grid((nrows*3, nrows*3), (1, 0), colspan=2, rowspan=2)
        ax4 = plt.subplot2grid((nrows*3, nrows*3), (1, 2), rowspan=2)

        ax1.set_title('(a)'), ax1.set_xlabel('xlabel'), ax1.set_ylabel('ylabel')
        ax2.set_title('(b)'), ax2.set_xlabel('xlabel'), ax2.set_ylabel('ylabel')
        ax3.set_title('(c)'), ax3.set_xlabel('xlabel'), ax3.set_ylabel('ylabel')
        ax4.set_title('(d)'), ax4.set_xlabel('xlabel'), ax4.set_ylabel('ylabel')


    elif style == 11:
        # This a sample template to plot multiple imshow figures
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
        for i, ax in enumerate(axes):
            ax.axis('off')
            ax.set_title('('+enumerate_alphabet[i]+')')
            ax.imshow(np.random.randint(0,100,1000).reshape(10,100))

        plt.subplots_adjust(right=0.95,left=0.05,top=0.9,bottom=0.1,hspace=0.0,wspace=0.01)

    return fig, axes

def plt_savefig(default_name=True, file_name='default_name.png'):
    if default_name:
        file_name = "".join([os.path.splitext(os.path.basename(__file__))[0], ".png"])
        #splittext required to get rid of extension name

    plt.tight_layout()
    plt.savefig(file_name)

def plt_color(type='Tab',id=0):
    '''
    Allows different plotting color styles to be specified from a pre-defined dictionary.
    Example:
        plt.scatter(np.linspace(0,100,100),np.linspace(0,100,100), c=tt.plt_color(id=1))
    '''
    if type == 'Tab':
        cdict={0: 'tab:blue', 1: 'tab:orange', 2: 'tab:green', 3: 'tab:red', 4: 'tab:purple',
           5: 'tab:brown', 6: 'tab:pink', 7: 'tab:gray', 8: 'tab:olive', 9: 'tab:cyan'}
    else:
        pass

    return cdict[id]
#plotting section ends---------------------------------------------------------------------

#def diff(conc,td,diff_type='long',upper_limit=0.9,lower_limit=0.1,taylor_constant=3.625):
#    '''
#    Takes conc (1D numpy array)) and backcalculates the longitudinal diffusion
#    based on Taylor 1953
#    '''
#    if diff_type == 'long': # calculate longitudinal diffusion
#        x_upper = np.min(np.where(conc<upper_limit)) / len(conc)
#        x_lower = np.min(np.where(conc<lower_limit)) / len(conc)
#        diff_val = (1/td)*((x_lower - x_upper)/taylor_constant)**2
#
#    elif diff_type == 'trans': # TO DO: calculate longitudinal diffusion
#        pass
#
#    return diff_val

def diff(distance,conc,var=0.5,diff_type='long'):
    '''
    Backcalculates the longitudinal or transverse diffusion (Taylor 1953, Hiby 1962)
    given
        1. distance = xd for longitudinal, yd for transverse
        2. conc = concentration to be fitted (as a function of distance)
        3. var = td for longitudinal and xd for transverse
    
    Example usage:
        
        # C below is the tracer profile to be fitted with error function
        xd = np.linspace(0.5/50,1,50)
        DL, C_calc = diff(xd, C, var=0.5, diff_type='long')
        plt.plot(xd,C)
        plt.plot(xd,C_calc)
    '''   
    if diff_type == 'long': # TO DO: calculate longitudinal diffusion
        def long_diff_func(td_xd, DL):
                td = td_xd['td'].iloc[0] 
                xd = td_xd['xd'].values
                return 0.5 * (erfc((xd-td)/(2*np.sqrt(DL*xd))))
        
        td_xd = pd.DataFrame()
        td_xd['xd'] = distance
        td_xd['td'] = var
        best_fit_DL, DL_covariance = optimize.curve_fit(long_diff_func, td_xd, conc, p0=0.0001) 
        best_fit_conc = long_diff_func(td_xd, best_fit_DL)
        
        return best_fit_DL, best_fit_conc

    elif diff_type == 'trans': # calculate transverse diffusion
        def trans_diff_func(xd_yd, DT):
                xd = xd_yd['xd'].iloc[0]
                yd = xd_yd['yd'].values
                return 1 - 0.5 * (erf((yd+0.5)/(2*np.sqrt(DT*xd)))+erf((yd-0.5)/(2*np.sqrt(DT*xd))))          
        
        xd_yd = pd.DataFrame()
        xd_yd['yd'] = distance
        xd_yd['xd'] = var
        best_fit_DT, DT_covariance = optimize.curve_fit(trans_diff_func, xd_yd, conc, p0=0.0001)
        best_fit_conc = trans_diff_func(xd_yd, best_fit_DT)
        
        return best_fit_DT, best_fit_conc
    
#def frac_flow(Sw,n=2,mw=1,mo=1, Swc=0.0, Sor=0.0):
#    '''
#    Calculates fractional flow of water
#    Corey's exponent n is assumed the same for water and oil phases
#    '''
#    F = (1/(1+((((1-Sw-Sor)/(1-Swc-Sor))**n*mw)/(((Sw-Swc)/(1-Swc-Sor))**n*mo))))
#
#    return F
#
#def frac_flow_grad(Sw,n=2,mw=1,mo=1, Swc=0.0, Sor=0.0):
#    '''
#    Calculates gradient of fractional flow (dF/dS)
#    Corey's exponent is assumed the same for water and oil phases.
#    The expression used here is fron symbolic differentiation on Matlab.
#    It probably can be simplify further using SciPy 'simplify' function.
#    '''
#    dF = (-((mw*n*((Sor + Sw - 1)/(Sor + Swc - 1))**(n - 1))/
#        (mo*(-(Sw - Swc)/(Sor + Swc - 1))**n*(Sor + Swc - 1)) +
#        (mw*n*((Sor + Sw - 1)/(Sor + Swc - 1))**n)/(mo*(-(Sw - Swc)/(Sor + Swc - 1))**(n + 1)*
#        (Sor + Swc - 1)))/((mw*((Sor + Sw - 1)/(Sor + Swc - 1))**n)/
#        (mo*(-(Sw - Swc)/(Sor + Swc - 1))**n) + 1)**2);
#
#    return dF
#
#def solve_bl(initial=0.75,n=2,mw=1,mo=1, Swc=0.0, Sor=0.0):
#    '''
#    Solves Buckley-Leverett's shock front. fsolve seems not to converge 
#    to correct solution for 'hard problems'
#    '''
#    from scipy.optimize import fsolve
#
#    def buckley_leverett(Sw):
#        return frac_flow_grad(Sw)*(Sw-Swc) - frac_flow(Sw)
#
#    Swf = fsolve(buckley_leverett, initial)
#
#    return Swf

def frac_flow(Sw,nw=2,no=2,mw=1,mo=1, Swc=0.0, Sor=0.0):
    '''
    Calculates fractional flow of water
    nw = Corey's exponent (water), no = Corey's exponent (oil)
    '''
    F = (1/(1+((((1-Sw-Sor)/(1-Swc-Sor))**no*mw)/(((Sw-Swc)/(1-Swc-Sor))**nw*mo))))

    return F

def frac_flow_grad(Sw,nw=2,no=2,mw=1,mo=1, Swc=0.0, Sor=0.0):
    '''
    Calculates gradient of fractional flow (dF/dS).
    nw = Corey's exponent (water), no = Corey's exponent (oil)
    The expression found using the following commands:
    
    import sympy as sp

    F, Sw, Swc, Sor, nw, no, mw, mo = sp.symbols('F Sw Swc Sor nw no mw mo')
    F = (1/(1+((((1-Sw-Sor)/(1-Swc-Sor))**no*mw)/(((Sw-Swc)/(1-Swc-Sor))**nw*mo))))
    dF = sp.simplify(sp.diff(F,Sw))
    
    '''

    dF = mo*mw*((-Sw + Swc)/(Sor + Swc - 1))**nw*((Sor + Sw - 1)/(Sor + Swc - 1))**no*(-no*(Sw - Swc) + nw*(Sor + Sw - 1))/((Sw - Swc)*(mo*((-Sw + Swc)/(Sor + Swc - 1))**nw + mw*((Sor + Sw - 1)/(Sor + Swc - 1))**no)**2*(Sor + Sw - 1))
    return dF

def solve_bl(nw=2,no=2,mw=1,mo=1,Swc=0.0,Sor=0.0):
    '''
    Solves BL using Welge's method.
    Returns shock front saturation and velocity (Swf and vwf),
    as well as the S and dF/dS up to shock front.
    
    '''
    
    res = 500 #resolution
    
    S = np.linspace(Swc+1/res,1-Sor,res)
    dFdS = frac_flow_grad(S,nw=nw,no=no,mw=mw,mo=mo,Swc=Swc,Sor=Sor)
    
    welge_area = np.full((res), np.nan)
    start_index = np.nanargmax(dFdS)
    
    # Iteratively going through from Sat with highest dF/dS to max Sat
    # and compare the Welge areas
    for i in range(start_index+1,res-1):
        Swf = S[i]
        vwf = frac_flow_grad(Swf,nw=nw,no=no,mw=mw,mo=mo,Swc=Swc,Sor=Sor)
        Swf_2_index = np.min(np.where(dFdS>vwf))
        Swf_2 = S[Swf_2_index]
        
        # Compute Welge areas using simple geometries 
        y1 = dFdS[np.where(dFdS>vwf)]
        area_1a = simps(y1, dx=(1-Swc-Sor)/res)
        area_1b = (Swf-Swf_2)*vwf
        
        y2= dFdS[:Swf_2_index]
        area_2a = (Swf_2-Swc)*vwf
        area_2b = simps(y2, dx=(1-Swc-Sor)/res)
        
        # Shock front sat is when welge_area=0
        welge_area[i] = np.abs((area_1a-area_1b)-(area_2a-area_2b))
        # Activate this plots for troubleshooting
#        plt.plot(dFdS,S)  
#        plt.plot([vwf,vwf], [0,Swf], '--k')
#        plt.scatter(vwf,Swf,c='b')
#        plt.scatter(vwf,Swf_2,c='r')
#        plt.title('Ignore this plot!!')   

    # Select Sat that meets the requirement  welge_area~=0
    selected_index = np.nanargmin(welge_area)
    Swf = S[selected_index]
    vwf = dFdS[selected_index]
    S_solution = S[selected_index:]
    v_solution = dFdS[selected_index:]
    
    return Swf, vwf, S_solution, v_solution 

def tl_frac_flow(C , M=1, omega=2/3):
    '''
    Calculates the Todd-Longstaff fractional flow curves
    '''
    Me = M**(1-omega)
    F = 1 / (1+(1/Me)*(1-C)/C)

    return F

def tl_frac_flow_grad(C, M=1, omega=2/3):
    '''
    to find diff using symbolic function:
    M_eff = Symbol('M_eff')
    C = Symbol('C')
    F = Symbol('F')
    F = 1 / (1+(1/M_eff(M))*(1-C)/C)
    dF = diff(F,C)
    simplify(dF)
    '''
    Me = M**(1-omega)
    dF = Me/(C*Me - C + 1)**2

    return dF

def tl_prod(C, M=1, omega=2/3, td=1):
    '''
    Computes the oil recovery in Todd-Longstaff model at time td.
    To find the curves, just call this function repeteadly with different values of td.
    '''
    #first we compute the dF curve
    dF = tl_frac_flow_grad(C, M, omega)
    #at time td, we have
    conc = dF*td
    # find the portion of the curve that goes beyond production well
    v1 = np.argwhere(conc>1)
    v1_loc = len(v1)
    # find area from injection well up to x=(1/Me)*td
    area_1 = conc.min()*1
    # find area from x=(1/Me)*td up to production well
    area_2 = np.trapz(C[v1_loc:], -conc[v1_loc:])
    # find area from production well up to finger tips
    area_3 = np.trapz(C[:v1_loc], -conc[:v1_loc])
    sanity_check = area_1 + area_2 + area_3 #must equal to 1
    Np = area_1 + area_2
    return Np
#------------------------------------------------------------------------------------------------
# Mistress-specific section

def delete_files(all_files=False,ext_name=''):
    if all_files:
        for f in glob.glob('*.rst1'):
            os.remove(f)
        for f in glob.glob('*.rst2'):
            os.remove(f)
        for f in glob.glob('*.vtk'):
            os.remove(f)
        for f in glob.glob('*.hist'):
            os.remove(f)
        for f in glob.glob('*.mout'):
            os.remove(f)
        for f in glob.glob('*.debug'):
            os.remove(f)

    else:

        for f in glob.glob('*.'+ext_name):
            os.remove(f)

def generate_ms_dat(from_template=False,template_name=' ',mistress_version='2D', modify_key=[], out_name='default'):
    '''
    Generates .dat file for Mistress 2D or 3D.

    Example:
        tt.generate_ms_dat(modify_key=['TITLE xxxx']) #all default values, except for TITLE
        tt.generate_ms_dat(from_template=True,template_name='test.dat',out_name='test2')
    #TO DO: add capability of adding new keyword without deleting existing one
            for e.g. multiple WELL or TOUT
    '''
    if from_template:
        file = open(template_name, 'r')
        template = file.read()
        file.close()

    else:
        if mistress_version == '2D':
            template = '''

TITLE Example data set for MISTRESS pseudos
NGRID   10   10
GSIZE   1.0  1.0
SOLVER ICCGS
*SOLVER BANDP
VISCW   1.0
VISCO   1.0
VISCS   1.0
DENSITYS 1.0
DENSITYO 1.0
THETA 1.0 ! IMPLICITNESS PARAM = 1.0 FOR PC AND DIFFUSION TERMS
SINIT   0.00
CINIT   0.00
SWCRIT  0.00
SORSDL  0.00
KROSWC  1.0
KRWSOR  1.0
*        NW   NO
RELPERM  2.0  2.0
DIFF 0.00 0.00
*ALPHA 0.0018 0.000036
*READTRAN
*MODERANX 1 555555
*READPERM
*PCMAX    10.0
*            BLX   TRX   BLY   TRY    SWI
*MODSINIT      1     10    1     1     0.8
*     GX    GY
GRAV  0.0   0.0
TOUT 1.5
TOUT 1.6
*TOUT 1.0
*TOUT 1.1
*TOUT 3.0
*TOUT 4.0
*TOUT 5.0
QINJ  1.00
*      VAL     TIME
*FWINJ  0.0     2.0
CINJ  0.0     0.0
FRQOUT 1
FRQDBG 10000
FRQRST 5000
*    BLX BLY  TRX TRY   TYPE          BHP    PI       CINJ
WELL  1   1    1   50   INJN                          1.0
WELL  1   51   1   NY   INJN                          0.0
*    BLX BLY  TRX TRY   TYPE          BHP    PI
WELL  NX  1    NX  NY   PROD          0.0    100000.0
COUR  0.4
* CHANGEVT 0.05
FCTS
FCTC
* Ask for pseudos for flow in x-direction.
* ----------------------------------------
* Direction No of.       Output frequency
*           grid blocks  (ie every 50 timesteps)
* ----------------------------------------------
*XPSEUDO     1              50
*XPSEUDO 2 50
*XPSEUDO 4 50
* Ask for output so we can generate effective relative permeabilities
* at a later date.
* -------------------------------------------------------------------
* Direction  Output frequency  Fluid (WATER or SOLV)
* -------------------------------------------------------------------
*EFFREPX          50            WATER
*READREG
OUTLEVEL 1
FULLSIZE
*DTMOVIE 1.0E-03 WATER
END

        '''
        elif mistress_version == '3D':
            template = '''
TITLE UNTITLE
*
SOLVER ICCG2S
*
NGRID  1000 100 1
GSIZE  10.0 1.0 1.0
*
GRAV 0.0 0.0 0.0
*
VISCW   1.0
VISCO   16.0
VISCS   1.0
*
VISCPOWR 4.0
*
DENSITYS 1.00
DENSITYO 1.00
*
SINIT 0.0
*CINIT CON 0.0
*
SWCRIT 0.0
SORSDL 0.0
KROSWC 1.0
KRWSOR 1.0
*
MODERANX 1 999 1
*
*       NW  NO
RELPERM  2   2
*
*
DIFF 0.0 0.0005 0.000
QINJ 1.0
FWC_INJ  50.0  0.0   1.0
*
*    BLX BLY BLZ TRX TRY  TRZ   TYPE    CINJ
WELL  1   1   1   1   NY    NZ   INJN
*    BLX BLY BLZ TRX TRY  TRZ   TYPE    BHP    PI
WELL  NX  1  1   NX  NY   NZ   PROD    1.0    10000.0
*
TOUT 1.1
*
COUR 0.4
THETA 1.0
*
FRQOUT 20
FRQRST 2000
*
FCTS
FCTC
OUTLEVEL 2
FULLSIZE
*
END
    '''
    textfile = open('temp_file.csv', 'w')
    textfile.write(template)
    textfile.close()
    df = pd.read_csv('temp_file.csv', header=None)
    #os.remove('temp_file.csv')
    for key in modify_key:
        a = df[0].str.contains(key[:6]) #match the first 4 letters in the keyword
        b = np.flatnonzero(a)
        df.iloc[int(b)] = key

    # for key in modify_key:
    #     for i in len(df):
    #     a = df[0].str.contains(key[:4]) #match the first 4 letters in the keyword
    #     b = np.flatnonzero(a)
    #     df.iloc[int(b)] = key

    out_name = out_name + '.dat'
    # Use this function to search for any files which match your filename
    files_present = glob.glob(out_name)

# if no matching files, write to csv, if there are matching files, print statement
    if not files_present:
        df.to_csv(out_name, header=None, index=None, mode='a')
        print('Sucessfully created new input file: '+out_name)
    else:
        print('WARNING: This file already exists!')

def read_hist(file_name):
    '''
    Loads .hist file into separate columns
    '''
    hist = pd.read_table(file_name, sep=r"\s*")
    return hist

#--------------------------------------------------------------------------
def generate_vtk(threeD_grid, variable_name, base_vtk_name, time_step):
    '''
    Takes 3D grid and turns it into legacy VTK format
    '''

    zi, yi, xi = threeD_grid.shape
    vtk_Name = base_vtk_name+"."+str(time_step)+".vtk"

    myFile = open(vtk_Name, "w")
    myFile.write("# vtk DataFile Version 2.0 \n")
    myFile.write(("Time Step:  %s \n") %(time_step))
    myFile.write("ASCII \n")
    myFile.write("DATASET STRUCTURED_POINTS \n")
    myFile.write(("DIMENSIONS  %s %s %s \n") %(xi+1,yi+1,zi+1))
    myFile.write("ORIGIN 0 0 0 \n")
    myFile.write("SPACING  1  1  1 \n")
    myFile.write(("CELL_DATA %s \n") %(xi*yi*zi))
    myFile.write(("SCALARS %s float \n") %(variable_name))
    myFile.write("LOOKUP_TABLE default\n")

    for z in range(zi):
        for y in range(yi):
            for x in range(xi):
                myFile.write('%0.3f \n' %(threeD_grid[z,y,x]) )

    myFile.close()

def mix_length(C,upper_limit=0.9,lower_limit=0.1):
    '''
    Returns the mixing length in 1D numpy array
    '''
    C1 = [x for x in C if x >= upper_limit]
    C2 = [x for x in C if x >= lower_limit]
    L_mix = float(len(C2)-len(C1))/len(C)

    return L_mix

def vtk_dim(file_name):
    '''
    Returns the grid dimension of a given vtk file
    '''
    df = pd.read_csv(file_name, header=None);
    dim_loc = np.flatnonzero(df[0].str.contains('DIMENSIONS'))
    dim = df.iloc[int(dim_loc)].str.split(expand=True)
    x_grid, y_grid, z_grid = \
        int(dim.iloc[0,1])-1, int(dim.iloc[0,2])-1, int(dim.iloc[0,3])-1

    return x_grid, y_grid, z_grid

def extract_2d_vtk(file_name, variable='Concentration'):
    '''
    Extracts the map of a specified variable (e.g. concentration)
    for a given vtk file
    '''

    x_grid, y_grid, z_grid = vtk_dim(file_name)

    df = pd.read_csv(file_name, header=None);
    start = np.flatnonzero(df[0].str.contains(variable))+2
    end = start+(x_grid*y_grid);
    map_2d = np.asfarray(df[int(start):int(end)])
    map_2d = np.reshape(map_2d,(y_grid,x_grid)) #saturation/concentration map

    return map_2d

def contour(var_map=0, type='vtk', file_name='', variable='Concentration', contour_val=0.5):
    '''
    Takes vtk file or concentration/saturation map and
    returns contour line and its peak-to-peak distance
    '''
    # we first load the concentration/saturation map


    if type == "vtk":
        x_grid, y_grid, z_grid = vtk_dim(file_name)
        map_2d = extract_2d_vtk(file_name=file_name, variable=variable)
    elif type == "map":
        x_grid, y_grid, z_grid = var_map.shape
        map_2d = var_map

    # then find the contour line
    con = np.zeros((y_grid,1))
    for j in range (0,y_grid):
        con_temp = [x for x in map_2d[j,:] if x >= contour_val] # location of contour line
        con[j,:] = (np.asarray(con_temp)).shape # contour_line

    L_mix = ((np.amax(con)-np.amin(con))/x_grid) # peak-to-peak distance of the contour

    return con, L_mix

def calc_finger(vtk_file):
    '''
    Calculates number of times where we have change in direction (in y- direction),
    based on Fortran code by Ann. To get the number of finger, divide nfinmx with 2
    '''
    
    C = extract_2d_vtk(vtk_file)
    C = np.transpose(C)
    nx, ny = C.shape
    
    ifxmax = 0
    nfinmx = 0
    small = 0.00005
    one = 0.9999
    
    for i in range(nx):
        nfingr = 0
        
        dcdy = C[i,1]-C[i,0]
        if (abs(dcdy)<small): #if dc/dy too small, see whether we are in C=1 or C=0 area..
            if(C[i,0]<small): #C=0 area
                gdsgn1 = 1.0 #1.0 means not the area with concentration
            elif (C[i,0]>one): #C=1 area
                gdsgn1 = -1.0
            else:
                gdsgn1 = np.sign(dcdy) #0<C<1 area
        else:
            gdsgn1 = np.sign(dcdy)
        
        for j in range(1,ny-1):
            dcdy = C[i,j+1]-C[i,j]    
            if (abs(dcdy)<small):
                if(C[i,0]<small):
                    gdsgn2 = 1.0
                elif (C[i,0]>one):
                    gdsgn2 = -1
                else:
                    gdsgn2 = gdsgn1
            else:
                gdsgn2 = np.sign(dcdy)
                
            gdchng = gdsgn1*gdsgn2 #indicate positive to negative dc/dy, or vice verca
            
            if (gdchng<small): #negative indicates change in direction
                nfingr = nfingr + 1.0
                
            gdsgn1 = gdsgn2
            
        if (nfingr>nfinmx):
            nfinmx = nfingr
            ifxmax = i
    
    return C, nfinmx/2

def extract_conc(file_name,upper_limit=0.9,lower_limit=0.1):
    '''
    Extracts the average concentration computed in Mistress (i.e. .conc file).
    The first column of the table is the time step. For each row (hence for 
    each time step), the average concentration across Y-axis is stored from 
    second column to the last column.
    '''   
    df = pd.read_csv(file_name, header=None, skiprows=1, delim_whitespace=True)
    C = df.values
    mix_len = np.full((len(df),2), np.nan)
    mix_len[:,0] = C[:,0]
    for j in range(len(df)):
        mix_len[j,1] = mix_length(C[j,1:]) # use the existing function to calculate mixing length
    
    return mix_len

#------------------------------------------------------------------------------------------------
# formatting section

def jptr_hide_code():
    '''
    Returns a string of HTML code to hide the cells with code in Jupyter.
    # http://blog.nextgenetics.net/?e=102
    # http://chris-said.io/2016/02/13/...
      ...how-to-make-polished-jupyter-presentations-with-optional-code-visibility/

    Usage example:
        from IPython.display import HTML
        HTML(tt.jupyter_hide_code())
    '''
    html_hide = '''
        <script>
          function code_toggle() {
            if (code_shown){
              $('div.input').hide('500');
              $('#toggleButton').val('Show Code')
            } else {
              $('div.input').show('500');
              $('#toggleButton').val('Hide Code')
            }
            code_shown = !code_shown
          }

          $( document ).ready(function(){
            code_shown=false;
            $('div.input').hide()
          });
        </script>
        <form action="javascript:code_toggle()"><input type="submit" id="toggleButton" value="Show Code"></form>
        '''
    return html_hide

def jptr_pdf_template():
    '''
    Creates template file that get rids of code cell when converting Jupyter to pdf file.

    Usage example:
        tt.jptr_pdf_template()
        !jupyter nbconvert test3.ipynb --to pdf --template hidecode.tplx
    '''
    template = '''
    ((*- extends 'article.tplx' -*))

    ((* block input_group *))
        ((*- if cell.metadata.get('nbconvert', {}).get('show_code', False) -*))
            ((( super() )))
        ((*- endif -*))
    ((* endblock input_group *))
    '''

    textfile = open('hidecode.tplx', 'w')
    textfile.write(template)
    textfile.close()

def test():
    return 'Hello'