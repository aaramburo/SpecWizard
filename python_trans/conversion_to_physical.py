import numpy as np

def conversion_to_physical(loss,simdata,convert,const):
    """
    For simulations that need to transform simulation data into physical quantities 
    INPUT: line of sight (los._) class containing the simulation data
           simulation data (simdata._) class containing the header information
           convertion data (convert._) class
           Physical constants (const._) class 
           
    OUTPUT: line of sight (los._) with the quantities transformed 
    """
    
    # Positions physical -> cm
    Coordinates_conv              = convert.Pos_aexp_exp * np.log10(simdata.ExpansionFactor) + convert.Pos_h_exp *                                                   np.log10(simdata.HubbleParam) + np.log10(convert.Pos_cgs_unit) -                                                                 np.log10(convert.cm_per_mpc)
    Coordinates_conv              = 10.0**Coordinates_conv
    loss.Position                 = loss.Position * Coordinates_conv
    loss.BoxPhys                  = simdata.BoxSize * Coordinates_conv
    loss.ParticleSmoothingLength  = loss.ParticleSmoothingLength * Coordinates_conv
    
    # Velocity physical  -> km/s
    Velocity_conv                 = convert.Vel_aexp_exp * np.log10(simdata.ExpansionFactor) + convert.Vel_h_exp *                                                     np.log10(simdata.HubbleParam) + np.log10(convert.Vel_cgs_unit) - np.log10(1e5) 
    Velocity_conv                 = 10.0**Velocity_conv
    loss.Velocity                 = loss.Velocity * Velocity_conv
    
    # Density physica -> g/cm³
    Density_conv                  = convert.Dens_aexp_exp * np.log10(simdata.ExpansionFactor) + convert.Dens_h_exp *                                                 np.log10(simdata.HubbleParam) + np.log10(convert.Dens_cgs_unit)-                                                                 np.log10(convert.proton_mass)
    Density_conv                  = 10**Density_conv
    loss.ParticleDensity          = loss.ParticleDensity * Density_conv
    #Convert from total density to *Hydrogen* number density
    loss.ParticleDensity          = loss.ParticleDensity * loss.MassFractions[:,0]
    loss.Zmetal                   = loss.Metallicity[...] / const.Zmass_solar  #metallicity in solar units
         
    #Mass
    Mass_conv                     = convert.Mass_aexp_exp * np.log10(simdata.ExpansionFactor) + convert.Mass_h_exp *                                                 np.log10(simdata.HubbleParam) + np.log10(convert.Mass_cgs_unit) -                                                                 np.log10(convert.solar_mass)
    Mass_conv                     = 10.0**Mass_conv
    loss.Mass                     = loss.Mass * Mass_conv
    
    #Temperature
    Temp_Conv                     = convert.Temp_aexp_exp * np.log10(simdata.ExpansionFactor) + convert.Temp_h_exp *                                                 np.log10(simdata.HubbleParam) + np.log10(convert.Temp_cgs_unit) #! cgs
    Temp_Conv                     = 10.0**Temp_Conv
    loss.ParticleTemperature      = loss.ParticleTemperature * Temp_Conv 
    
    
    loss.x_physical = loss.x_position * Coordinates_conv  #LOS x-starting point
    loss.y_physical = loss.y_position * Coordinates_conv  #LOS y-starting point
    
    
    return loss