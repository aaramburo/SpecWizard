import numpy as np
import h5py as h5
import sys
import physical_data as const
from read_parameter_file import read_params
import reading_data_for_los as RD
import matplotlib.pyplot as plt
from computeibb import computeibb as computeib
from prodata import projectdata
from projectdataf import project_data
import time
from intergrid import Intergrid   
import read_ionizationtables as RIT

def los_long_axis(Position,HubbleParam,ExpansionFactor):
    '''
    Returns the long axis that the LOS belongs to
    0-x
    1-y
    2-z
    '''
    return np.argmax([(Position[:,0].max()-Position[:,0].min()) /HubbleParam * ExpansionFactor,(Position[:,1].max()-Position[:,1].min()) /HubbleParam * ExpansionFactor,(Position[:,2].max()-Position[:,2].min()) /HubbleParam * ExpansionFactor])



def align_positions(los_long_axis,Position):
    '''
    Changes the frame of reference of the coordinates, depending on the long axis.
    '''
    Position_2 = np.zeros(np.shape(Position),dtype=np.float128)
    if los_long_axis == 0:
        Position_2[:,0] = Position[:,1]
        Position_2[:,1] = Position[:,2]
        Position_2[:,2] = Position[:,0]
    elif los_long_axis == 1:
        Position_2[:,0] = Position[:,2]
        Position_2[:,1] = Position[:,0]
        Position_2[:,2] = Position[:,1]

    else:
        Position_2[:,0] = Position[:,0]
        Position_2[:,1] = Position[:,1]
        Position_2[:,2] = Position[:,2]
    
    return Position_2

times = []
#Load dummy parameters from a "typical" specwizard run 
parameters = read_params('dummy.par')

ionpar = RIT.read_iontables(parameters)

for lm in range(40):

    #Load a small EAGLE box for testing
    fname = "/home/andres/small_snap/snap_028_z000p000"

    #Load a example of a postprocess LOS file
    fname2 = "/home/andres/part_los_z0.010"

    #Load class that contianes the particle information i,e position, velocities etc.
    loss = RD.read_particle_data_from_los(fname2,lm)        

    #Load the class the contains the header information
    simdata = RD.read_header(fname2)

    #Conversion factors 
    convert = RD.units_and_factors_for_los(fname2)

    # This outputs which is what is the long axis from the LOS  (0=x,1=y,z=2)

    loss.los_long_axis = los_long_axis(loss.Position, simdata.HubbleParam, simdata.ExpansionFactor)


    acurrent = simdata.ExpansionFactor
    zcurrent = simdata.Redshift

    simdata.CurrentHubbleCt = 100. * simdata.HubbleParam *  np.sqrt(1. + simdata.Omega0*(1./acurrent-1.) + simdata.OmegaLambda* (acurrent**2-1.)) /acurrent
    simdata.boxkms = simdata.BoxSize / simdata.HubbleParam * acurrent * simdata.CurrentHubbleCt

    simdata.vpixsizekms = 1
    simdata.nveloc = int(simdata.boxkms / simdata.vpixsizekms) + 1


    Coordinates_conv     = convert.Pos_aexp_exp * np.log10(simdata.ExpansionFactor) + convert.Pos_h_exp * np.log10(simdata.HubbleParam) + np.log10(convert.Pos_cgs_unit) - np.log10(convert.cm_per_mpc) # in physical Mpc
    Coordinates_conv     = 10.0**Coordinates_conv
    loss.Position = loss.Position * Coordinates_conv
    BoxPhys = simdata.BoxSize * Coordinates_conv
    loss.ParticleSmoothingLength = loss.ParticleSmoothingLength * Coordinates_conv


    Velocity_conv = convert.Vel_aexp_exp * np.log10(simdata.ExpansionFactor) + convert.Vel_h_exp * np.log10(simdata.HubbleParam) + np.log10(convert.Vel_cgs_unit) - np.log10(1e5) #$ physical km/s
    Velocity_conv = 10.0**Velocity_conv
    loss.Velocity = loss.Velocity * Velocity_conv


    Density_conv = convert.Dens_aexp_exp * np.log10(simdata.ExpansionFactor) + convert.Dens_h_exp * np.log10(simdata.HubbleParam) + np.log10(convert.Dens_cgs_unit)- np.log10(convert.proton_mass)
    Density_conv = 10**Density_conv

    loss.ParticleDensity = loss.ParticleDensity * Density_conv
    #Convert from total density to *Hydrogen* number density
    loss.ParticleDensity = loss.ParticleDensity * loss.MassFractions[:,0]
    loss.Metallicity = loss.Metallicity[...] / const.Zmass_solar  #metallicity in solar units


    Mass_conv = convert.Mass_aexp_exp * np.log10(simdata.ExpansionFactor) + convert.Mass_h_exp * np.log10(simdata.HubbleParam) + np.log10(convert.Mass_cgs_unit) - np.log10(convert.solar_mass)
    Mass_conv = 10.0**Mass_conv
    loss.Mass   = loss.Mass * Mass_conv


    Temp_Conv           = convert.Temp_aexp_exp * np.log10(simdata.ExpansionFactor) + convert.Temp_h_exp * np.log10(simdata.HubbleParam) + np.log10(convert.Temp_cgs_unit) #! cgs
    Temp_Conv           = 10.0**Temp_Conv
    loss.ParticleTemperature = loss.ParticleTemperature * Temp_Conv 


    loss.Position = align_positions(loss.los_long_axis,loss.Position)
    loss.x_physical = loss.x_position * Coordinates_conv  #LOS x-starting point
    loss.y_physical = loss.y_position * Coordinates_conv  #LOS y-starting point

    

    loss,tf = project_data(loss,simdata,parameters,ionpar)
    times.append(tf)

plt.hist(times,10)
plt.savefig("histtime.png",dpi=300)
print(np.average(times))