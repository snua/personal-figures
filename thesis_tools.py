import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def diff(conc,td,diff_type='long',upper_limit=0.9,lower_limit=0.1,taylor_constant=3.625):
    '''
    Takes conc (1D numpy array)) and backcalculates the longitudinal diffusion
    based on Taylor 1953
    '''
    if diff_type == 'long': # calculate longitudinal diffusion
        x_upper = np.min(np.where(conc<upper_limit)) / len(conc)
        x_lower = np.min(np.where(conc<lower_limit)) / len(conc)
        diff_val = (1/td)*((x_lower - x_upper)/taylor_constant)**2

    elif diff_type == 'trans': # TO DO: calculate longitudinal diffusion
        pass

    return diff_val

def frac_flow(Sw,n=2,mw=1,mo=1, Swc=0.0, Sor=0.0):
    '''
    Calculates fractional flow of water
    Corey's exponent n is assumed the same for water and oil phases
    '''
    F = (1/(1+((((1-Sw-Sor)/(1-Swc-Sor))**n*mw)/(((Sw-Swc)/(1-Swc-Sor))**n*mo))))

    return F

def frac_flow_grad(Sw,n=2,mw=1,mo=1, Swc=0.0, Sor=0.0):
    '''
    Calculates gradient of fractional flow (dF/dS)
    Corey's exponent is assumed the same for water and oil phases
    '''
    dF = (-((mw*n*((Sor + Sw - 1)/(Sor + Swc - 1))**(n - 1))/
        (mo*(-(Sw - Swc)/(Sor + Swc - 1))**n*(Sor + Swc - 1)) +
        (mw*n*((Sor + Sw - 1)/(Sor + Swc - 1))**n)/(mo*(-(Sw - Swc)/(Sor + Swc - 1))**(n + 1)*
        (Sor + Swc - 1)))/((mw*((Sor + Sw - 1)/(Sor + Swc - 1))**n)/
        (mo*(-(Sw - Swc)/(Sor + Swc - 1))**n) + 1)**2);

    return dF

def solve_bl(initial=0.75,n=2,mw=1,mo=1, Swc=0.0, Sor=0.0):
    '''
    Solves Buckley-Leverett's shock front
    '''
    from scipy.optimize import fsolve

    def buckley_leverett(Sw):
        return frac_flow_grad(Sw)*(Sw-Swc) - frac_flow(Sw)

    Swf = fsolve(buckley_leverett, initial)

    return Swf

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
def generate_ms_dat(mistress_version='2D', modify_key=[], file_name='default'):
    '''
    Generates .dat file for Mistress 2D or 3D.
    #TO DO: add capability of adding new keyword without deleting existing one
            for e.g. multiple WELL or TOUT
    '''
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
        a = df[0].str.contains(key[:5]) #match the first 4 letters in the keyword
        b = np.flatnonzero(a)
        df.iloc[int(b)] = key

    # for key in modify_key:
    #     for i in len(df):
    #     a = df[0].str.contains(key[:4]) #match the first 4 letters in the keyword
    #     b = np.flatnonzero(a)
    #     df.iloc[int(b)] = key

    file_name = file_name + '.dat'
    # Use this function to search for any files which match your filename
    files_present = glob.glob(file_name)

# if no matching files, write to csv, if there are matching files, print statement
    if not files_present:
        df.to_csv(file_name, header=None, index=None, mode='a')
        print('Sucessfully created new .dat file')
    else:
        print('WARNING: This file already exists!')


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
