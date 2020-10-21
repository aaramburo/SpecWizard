# A merge between initialize_ionization_tables and load_ionbal
import numpy as np
import h5py as h5
import sys
from phs
import physical_data as const
parameters = read_params('dummy.par')

ions_to_do = which_ions_to_do = [parameters.doH1,parameters.doHe2,parameters.doC2,
parameters.doC3,parameters.doC4,parameters.doN2,parameters.doN3,parameters.doN4,
parameters.doN5,parameters.doO1,parameters.doO3,parameters.doO4,parameters.doO5,
parameters.doO6,parameters.doO7,parameters.doMg2,parameters.doNe8,parameters.doAl2,
parameters.doAl3,parameters.doSi2,parameters.doSi3,parameters.doSi4,parameters.doS5,
parameters.doFe2,parameters.doFe3,parameters.do21cm]

number_ions = np.sum(ions_to_do) # was nion in the original F90 subroutine

number_species =  number_ions# was nspecies in the original F90

if (number_ions==0):
    print("ERROR: No ions were selected!")
    sys.exit()
if (parameters.doH1 and parameters.do_long_spectrum and parameters.nlyman > const.nlyman_all):
    print('ERROR: nlyman > n_lines_max')
    sys.exit()
if (conts.Lambda_H1[0] !=  const.lyalpha):    # a bit superflus
    print( 'ERROR: Lambda_H1[0] and lyalpha  are not equal ')

ionbal_names = np.array(["h1","he2","c2","c3","c4","n2","n3","n4","n5","o1","o3","o4","o5","o6",  #This elements of this list must be in the same order as
                "o7","mg2","ne8","al2","al3","si2","si3","si4","s5","fe2","fe3","21cm"]) # ions_to_do

ionbal_to use = ionbal_names[ions_to_do]

# something like

ionss = []
first_iteration = True
for ion_file in ionbal_to_use:
    with h5.File(parameters.ibdir+ion_file+".hdf5","r") as f:
        if first_iteration:
            first_iteration = False
            z_ranges_table = f["/redshift"][...]
            logt_table = f["/logt"][...]
            logd_talbe = f["/logd"][...]
        print(ion_file)
        ionss.append([ion_file,f["/ionbal"][...]])