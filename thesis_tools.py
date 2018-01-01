def diff(conc,td,diff_type=long,upper_limit=0.9,lower_limit=0.1,taylor_constant = 3.625):
    '''
    Takes 1D numpy array and backcalculates the longitudinal diffusion
    based on Taylor 1953
    '''
    if diff_type == long: # calculate longitudinal diffusion
        x_upper = np.min(np.where(conc<upper_limit)) / len(conc)
        x_lower = np.min(np.where(conc<lower_limit)) / len(conc)
        diff_val = (1/td)*((x_lower - x_upper)/taylor_constant)**2

    elif diff_type == trans: # TO DO: calculate longitudinal diffusion
        pass

    return diff_val

def frac_flow(Sw, corey_exp=2, muw=1, muo=1):
    '''
    Calculates fractional flow of water
    '''
    krw = (Sw)**corey_exp
    kro = (1-Sw)**corey_exp
    F = 1/(1+(kro/muo)*(muw/krw))

    return F

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