import numpy as np
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
