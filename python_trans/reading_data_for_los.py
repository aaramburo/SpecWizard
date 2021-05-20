import numpy as np
import h5py as h5
import sys
import physical_data as const
import numpy as np
import read_eagle
import random


class read_particle_data():

    def __init__(self,parameters,header,LOS_num=None):
        fname = parameters.datadir + parameters.snap_base
        f       = h5.File(fname+'.0.hdf5', 'r')

        def read_datasets(itype, att, nfiles=16):
            """ Read a selected dataset, itype is the PartType and att is the attribute name. """
            # Output array.
            data = []
            # Loop over each file and extract the data.
            for i in range(nfiles):
                f = h5.File(fname+'.%i.hdf5'%i, 'r')
                tmp = f['PartType%i/%s'%(itype, att)][...]
                data.append(tmp)

                f.close()
            # Combine to a single array.
            if len(tmp.shape) > 1:
                data = np.vstack(data)
            else:
                data = np.concatenate(data)
            return data

        def get_coordinates_for_los(header,LOS_num=None):
            """
            if use_random_los = True, it generates a random (x,y) starting point for the LOS
            otherwise it loads from  a coordinate file a fixed pair of points. 
            Outputs the x and y widths to define the region for read_eagle
            """

            if parameters.use_random_los:
                x_frac = random.random() * header.BoxSize
                y_frac = random.random() * header.BoxSize
            else:        
                if LOS_num==None: print("ERROR you must provide the LOS Number for this mode")

                cord_frac = h5.File(parameters.los_coordinates_file,'r')
                x_frac = cord_frac['Projection/x_fraction_array'][LOS_num] * header.BoxSize
                y_frac = cord_frac['Projection/y_fraction_array'][LOS_num] * header.BoxSize

            N1D    =  int(float(header.NumPart_Total[0])**1/3)
            extent = 1./N1D

            RegionExtentX = [x_frac - 8*extent, x_frac + 8*extent]
            RegionExtentY = [y_frac - 8*extent, y_frac + 8*extent]

            return RegionExtentX,RegionExtentY




        self.NumFilesPerSnapshot     = f['Header'].attrs.get("NumFilesPerSnapshot")

        if parameters.COLIBRE:
            Key = {'Cords':'Coordinates','Vel':'Velocities','Dens':'Densities','Mass':'Masses','Temp':'Temperatures','SFR':'StarFormationRates','SML':'SmoothingLengths','EMF':'ElementMassFractions'}

        else:
            Key = {'Cords':'Coordinates','Vel':'Velocity','Dens':'Density','Mass':'Mass','Temp':'Temperature','SFR':'StarFormationRate','SML':'SmoothingLength','EMF':'ElementAbundance','Met':'Metallicity'}

        # if parameters.use_smoothed_abundance:
        #     Key['EMF'] = 'SmoothedElemenAbundance'
        #     Key['Met'] = 'SmoothedMetallicity'

        if parameters.read_eagle: #Read_region

	   

	    fname = parameters.datadir+parameters.snap_base
	    x_range, y_range = get_coordinates_for_los(header,LOS_num)
	    snap = read_eagle.EagleSnapshot(fname+".0.hdf5")
	    snap.select_region(x_range[0], x_range[1],y_range[0], y_range[1],0, header.BoxSize)
	
            self.Position                       = snap.read_dataset(0, Key['Cords'])
            self.Velocity                       = snap.read_dataset(0, Key['Vel'])
            self.ParticleDensity                = snap.read_dataset(0, Key['Dens'])
            self.StarFormationRate              = snap.read_dataset(0, Key['SFR'])
            self.ParticleSmoothingLength        = snap.read_dataset(0, Key['SML'])
    #        self.PartID = snap.read_dataset(0, Key['Cords'])
            self.ParticleNeutralHFraction       = 0
            self.ParticleMolecularHFraction     = 0
            self.ParticleTemperature            = snap.read_dataset(0, Key['Temp'])
            self.Metallicity                    = snap.read_dataset(0, Key['Met'])
            self.Mass                           = snap.read_dataset(0, Key['Mass'])
            self.MassFractions                  = np.zeros((len(self.Mass),9))
            self.MassFractions[:,0]             = snap.read_dataset(0, Key['EMF']+'/Hydrogen')
            self.MassFractions[:,1]             = snap.read_dataset(0, Key['EMF']+'/Helium')
            self.MassFractions[:,2]             = snap.read_dataset(0, Key['EMF']+'/Carbon')
            self.MassFractions[:,3]             = snap.read_dataset(0, Key['EMF']+'/Silicon')
            self.MassFractions[:,4]             = snap.read_dataset(0, Key['EMF']+'/Iron')
            self.MassFractions[:,5]             = snap.read_dataset(0, Key['EMF']+'/Magnesium')
            self.MassFractions[:,6]             = snap.read_dataset(0, Key['EMF']+'/Nitrogen')
            self.MassFractions[:,7]             = snap.read_dataset(0, Key['EMF']+'/Neon')
            self.MassFractions[:,8]             = snap.read_dataset(0, Key['EMF']+'/Oxygen')



        elif  self.NumFilesPerSnapshot > 1:

            self.Position                       = read_datasets(0, Key['Cords'],self.NumFilesPerSnapshot)
            self.Velocity                       = read_datasets(0, Key['Vel'],self.NumFilesPerSnapshot)
            self.ParticleDensity                = read_datasets(0, Key['Dens'],self.NumFilesPerSnapshot)
            self.StarFormationRate              = read_datasets(0, Key['SFR'],self.NumFilesPerSnapshot)
            self.ParticleSmoothingLength        = read_datasets(0, Key['SML'],self.NumFilesPerSnapshot)
    #        self.PartID = snap.read_dataset(0, Key['Cords'])
            self.ParticleNeutralHFraction       = 0
            self.ParticleMolecularHFraction     = 0
            self.ParticleTemperature            = read_datasets(0, Key['Temp'],self.NumFilesPerSnapshot)
            self.Metallicity                    = read_datasets(0, Key['Met'],self.NumFilesPerSnapshot)
            self.Mass                           = read_datasets(0, Key['Mass'],self.NumFilesPerSnapshot)
            self.MassFractions                  = np.zeros((len(self.Mass),9))
            self.MassFractions[:,0]             = read_datasets(0, Key['EMF']+'/Hydrogen')
            self.MassFractions[:,1]             = read_datasets(0, Key['EMF']+'/Helium')
            self.MassFractions[:,2]             = read_datasets(0, Key['EMF']+'/Carbon')
            self.MassFractions[:,3]             = read_datasets(0, Key['EMF']+'/Silicon')
            self.MassFractions[:,4]             = read_datasets(0, Key['EMF']+'/Iron')
            self.MassFractions[:,5]             = read_datasets(0, Key['EMF']+'/Magnesium')
            self.MassFractions[:,6]             = read_datasets(0, Key['EMF']+'/Nitrogen')
            self.MassFractions[:,7]             = read_datasets(0, Key['EMF']+'/Neon')
            self.MassFractions[:,8]             = read_datasets(0, Key['EMF']+'/Oxygen')




class units_and_factors():
    """ Read the units and conversion factors for COLIBRE, Aurora and Eagle """
    def __init__(self, parameters):
        fname = parameters.datadir + parameters.snap_base
        if parameters.use_los_file:
            f       = h5.File(fname+'.hdf5', 'r')

        else:
            f       = h5.File(fname+'.0.hdf5', 'r')

        if parameters.COLIBRE:
            self.cm_per_mpc          = f['Units'].attrs.get("Unit length in cgs (U_L)")
            self.proton_mass         = f['PhysicalConstants/CGS'].attrs.get("proton_mass")
            self.solar_mass          = f['PhysicalConstants/CGS'].attrs.get("solar_mass")
            Key = {'Cords':'Coordinates','Vel':'Velocities','Dens':'Densities','Mass':'Masses','Temp':'Temperatures','SFR':'StarFormationRates','SML':'SmoothingLengths'}
            Atrr = {'hexp':'h-scale exponent','aexp':'a-scale exponent','CGSfac':'Conversion factor to CGS (not including cosmological corrections)'}
        else:

            self.cm_per_mpc          = f['Constants'].attrs.get("CM_PER_MPC")
            self.proton_mass         = f['Constants'].attrs.get("PROTONMASS")
            self.solar_mass          = f['Constants'].attrs.get("SOLAR_MASS")
            Key = {'Cords':'Coordinates','Vel':'Velocity','Dens':'Density','Mass':'Mass','Temp':'Temperature','SFR':'StarFormationRate','SML':'SmoothingLength'}
            Atrr = {'hexp':'h-scale-exponent','aexp':'a-scale-exponent','CGSfac':'CGSConversionFactor'}



        self.Pos_h_exp           = f['PartType0/'+Key['Cords']].attrs.get(Atrr['hexp'])
        self.Pos_aexp_exp        = f['PartType0/'+Key['Cords']].attrs.get(Atrr['aexp'])
        self.Pos_cgs_unit        = f['PartType0/'+Key['Cords']].attrs.get(Atrr['CGSfac'])
        self.Vel_h_exp           = f['PartType0/'+Key['Vel']].attrs.get(Atrr['hexp'])
        self.Vel_aexp_exp        = f['PartType0/'+Key['Vel']].attrs.get(Atrr['aexp'])
        self.Vel_cgs_unit        = f['PartType0/'+Key['Vel']].attrs.get(Atrr['CGSfac'])
        self.Dens_h_exp          = f['PartType0/'+Key['Dens']].attrs.get(Atrr['hexp'])
        self.Dens_aexp_exp       = f['PartType0/'+Key['Dens']].attrs.get(Atrr['aexp'])
        self.Dens_cgs_unit       = f['PartType0/'+Key['Dens']].attrs.get(Atrr['CGSfac'])
        self.Temp_h_exp          = f['PartType0/'+Key['Temp']].attrs.get(Atrr['hexp'])
        self.Temp_aexp_exp       = f['PartType0/'+Key['Temp']].attrs.get(Atrr['aexp'])
        self.Temp_cgs_unit       = f['PartType0/'+Key['Temp']].attrs.get(Atrr['CGSfac'])
        self.Mass_h_exp          = f['PartType0/'+Key['Mass']].attrs.get(Atrr['hexp'])
        self.Mass_aexp_exp       = f['PartType0/'+Key['Mass']].attrs.get(Atrr['aexp'])
        self.Mass_cgs_unit       = f['PartType0/'+Key['Mass']].attrs.get(Atrr['CGSfac'])
        f.close()

class read_header():
    """
    Read the header of the snapshot files, right now works for COLIBRE, Aurora and Eagle
    """
    def __init__(self, parameters):

        fname = parameters.datadir + parameters.snap_base

        if parameters.use_los_file:
            f       = h5.File(fname+'.hdf5', 'r')

        else:
            f       = h5.File(fname+'.0.hdf5', 'r')

        self.NumPart_ThisFile        = f['Header'].attrs.get('NumPart_ThisFile')
        self.NumPart_Total           = f['Header'].attrs.get('NumPart_Total')
        self.NumPart_Total_HighWord  = f['Header'].attrs.get("NumPart_Total_HighWord")
        self.NumFilesPerSnapshot     = f['Header'].attrs.get("NumFilesPerSnapshot")
        self.MassTable               = f['Header'].attrs.get("MassTable")

        if not parameters.COLIBRE:

            self.Flag_Sfr            = f['Header'].attrs.get("Flag_Sfr")
            self.Flag_Cooling        = f['Header'].attrs.get("Flag_Cooling")
            self.Flag_StellarAge     = f['Header'].attrs.get("Flag_StellarAge")
            self.Flag_Metals         = f['Header'].attrs.get("Flag_Metals")
            self.Flag_Feedback       = f['Header'].attrs.get("Flag_Feedback")
            self.ExpansionFactor     = f['Header'].attrs.get("ExpansionFactor")
            self.Redshift            = f['Header'].attrs.get("Redshift")
            self.BoxSize             = f['Header'].attrs.get("BoxSize")
            self.Omega0              = f['Header'].attrs.get("Omega0")
            self.OmegaLambda         = f['Header'].attrs.get("OmegaLambda")
            self.OmegaBaryon         = f['Header'].attrs.get("OmegaBaryon")
            self.HubbleParam         = f['Header'].attrs.get("HubbleParam")

        else:

            self.ExpansionFactor     = f['Header'].attrs.get("Scale-factor")
            self.Redshift            = f['Header'].attrs.get("Redshift")
            self.BoxSize             = f['Header'].attrs.get("BoxSize")
            self.Omega0              = f['Cosmology'].attrs.get("Omega0")
            self.OmegaBaryon         = f['Cosmology'].attrs.get("Omega_b")
            self.OmegaLambda         = f['Cosmology'].attrs.get("Omega_lambda")
            self.HubbleParam         = f['Cosmology'].attrs.get("h")

        self.Time                    = -1
        self.rhoc                    = 3 * (const.H0*self.HubbleParam)**2 / (8. * np.pi * const.G)
        self.rhocb                   = self.rhoc / self.ExpansionFactor**3 * self.OmegaBaryon

        f.close()



class read_particle_data_from_los():
    def __init__(self,parameters,LOSN=0):
        Key = {'Cords':'Positions',
              'Vel':'Velocity',
              'Dens':'Density',
              'Mass':'Mass',
              'Temp':'Temperature',
              'SFR':'StarFormationRate',
              'SML':'SmoothingLength',
              'EMF':'ElementAbundance'
              ,'Met':'Metallicity'}

        fname = parameters.datadir + parameters.snap_base

        f       = h5.File(fname+'.hdf5', 'r')
        losn = "LOS"+str(LOSN)+"/"
        print(losn+'PartType0/'+Key['Cords'])
        self.Position                       = f[losn+Key['Cords']]
        self.Velocity                       = f[losn+Key['Vel']]
        self.ParticleDensity                = f[losn+Key['Dens']]
        self.StarFormationRate              = f[losn+Key['SFR']]
        self.ParticleSmoothingLength        = f[losn+Key['SML']]
#        self.PartID = snap.read_dataset(0, Key['Cords'])
        self.ParticleNeutralHFraction       = 0
        self.ParticleMolecularHFraction     = 0
        self.ParticleTemperature            = f[losn+Key['Temp']]
        self.Metallicity                    = f[losn+Key['Met']]
        self.Mass                           = f[losn+Key['Mass']]
        self.MassFractions                  = np.zeros((len(self.Mass),9))
        self.MassFractions[:,0]             = f[losn+Key['EMF']+'/Hydrogen']
        self.MassFractions[:,1]             = f[losn+Key['EMF']+'/Helium']
        self.MassFractions[:,2]             = f[losn+Key['EMF']+'/Carbon']
        self.MassFractions[:,3]             = f[losn+Key['EMF']+'/Silicon']
        self.MassFractions[:,4]             = f[losn+Key['EMF']+'/Iron']
        self.MassFractions[:,5]             = f[losn+Key['EMF']+'/Magnesium']
        self.MassFractions[:,6]             = f[losn+Key['EMF']+'/Nitrogen']
        self.MassFractions[:,7]             = f[losn+Key['EMF']+'/Neon']
        self.MassFractions[:,8]             = f[losn+Key['EMF']+'/Oxygen']
        self.x_position                     = f[losn].attrs.get('x-position')
        self.y_position                     = f[losn].attrs.get('y-position')



class units_and_factors_for_los():
    """ Read the units and conversion factors for COLIBRE, Aurora and Eagle """
    def __init__(self,parameters,LOSN=0):

        fname = parameters.datadir + parameters.snap_base

        f       = h5.File(fname+'.hdf5', 'r')

        self.cm_per_mpc          = f['Constants'].attrs.get("CM_PER_MPC")
        self.proton_mass         = f['Constants'].attrs.get("PROTONMASS")
        self.solar_mass          = f['Constants'].attrs.get("SOLAR_MASS")
        Key = {'Cords':'Positions','Vel':'Velocity','Dens':'Density','Mass':'Mass','Temp':'Temperature','SFR':'StarFormationRate','SML':'SmoothingLength'}
        Atrr = {'hexp':'h-scale-exponent','aexp':'aexp-scale-exponent','CGSfac':'CGSConversionFactor'}
        losn = "LOS"+str(LOSN)+"/"

        self.Pos_h_exp           = f[losn+Key['Cords']].attrs.get(Atrr['hexp'])
        self.Pos_aexp_exp        = f[losn+Key['Cords']].attrs.get(Atrr['aexp'])
        self.Pos_cgs_unit        = f[losn+Key['Cords']].attrs.get(Atrr['CGSfac'])
        self.Vel_h_exp           = f[losn+Key['Vel']].attrs.get(Atrr['hexp'])
        self.Vel_aexp_exp        = f[losn+Key['Vel']].attrs.get(Atrr['aexp'])
        self.Vel_cgs_unit        = f[losn+Key['Vel']].attrs.get(Atrr['CGSfac'])
        self.Dens_h_exp          = f[losn+Key['Dens']].attrs.get(Atrr['hexp'])
        self.Dens_aexp_exp       = f[losn+Key['Dens']].attrs.get(Atrr['aexp'])
        self.Dens_cgs_unit       = f[losn+Key['Dens']].attrs.get(Atrr['CGSfac'])
        self.Temp_h_exp          = f[losn+Key['Temp']].attrs.get(Atrr['hexp'])
        self.Temp_aexp_exp       = f[losn+Key['Temp']].attrs.get(Atrr['aexp'])
        self.Temp_cgs_unit       = f[losn+Key['Temp']].attrs.get(Atrr['CGSfac'])
        self.Mass_h_exp          = f[losn+Key['Mass']].attrs.get(Atrr['hexp'])
        self.Mass_aexp_exp       = f[losn+Key['Mass']].attrs.get(Atrr['aexp'])
        self.Mass_cgs_unit       = f[losn+Key['Mass']].attrs.get(Atrr['CGSfac'])
        f.close()


# simzmin                    = 1e12
# simzmax                    = -1.
#
#
#
# ''''
# Read header
#
# ''''
# if (parameters.use_snapshot_file):
#
#     nlosfiles               = 1
#
#     if (parameters.EAGLE or parameters.AURORA): #must add proper flags for this
#         simfile             = parameters.snap_base + '.0.hdf5'
#
#     elif (parameters.COLIBRE):
#         simfile             = parameters.snap_base + '.hdf5'
#
#     else:
#         simfile             = '/snapshot_' # check this line SW_Sub.F90 968
#
#     head = read_header(parameters.datadir + simfile) #probably is better to manage it as a class
#
#     simz                    = head.Redshift
#     simzmin                 = head.Redshift
#     simzmax                 = head.Redshift
#
# else:
#     print("verificar en que casos sucede esto")
#
#
# if (parameters.do_long_spectrum):
#     print("something to do")
#
#
#
# ''''
# read or create coordiante values by doing fractions of the boxsize
# ''''
#
#
# if parameters.use_random_los:
#     x_fracction_array       = np.random.rand(parameters.nspec)
#     y_fracction_array       = np.random.rand(parameters.nspec)
#     z_fracction_array       = np.random.rand(parameters.nspec) #not sure this is right
#     phi_array               = np.zeros(parameters.nspec)
#     theta_array             = np.zeros(parameters.nspec)
#
# else:
#     print("jk")
#
#
# numspec                     = parameters.numspec
#
# ''''
# Get the fractions into co-moving and physical coordinates
# ''''
#
# #cMpc
# x_comoving                  = x_fraction_array * head.BoxSize
# y_comoving                  = y_fraction_array * head.BoxSize
# z_comoving                  = z_fraction_array * head.BoxSize
# #Mpc
# x_physical                  = x_comoving * head.ExpansionFactor / head.HubbleParam
# y_physical                  = y_comoving * head.ExpansionFactor / head.HubbleParam
# z_physical                  = z_comoving * head.ExpansionFactor / head.HubbleParam
#
#
# for i in numspec:
#
#
# if(parameters.readregion): # must add this flag to parameterfile
#     N1D                     = np.int(round(head.NumPart_Total[0]**(1./3.)))
#     extent                  = 1. / np.float(N1D)
#
#     RegionExtentX           = np.array([x_fraction_array[i] - (8. * extent),
#                               x_fraction_array[i] + (8. * extent)])
#     RegionExtentY           = np.array([y_fraction_array[i] - (8. * extent),
#                               y_fraction_array[i] + (8. * extent)])
#     RegionExtentZ           = np.array([0 ,1.])
#
#     # in comoving
#     RegionExtentX          = RegionExtentX * head.BoxSize
#     RegionExtentY          = RegionExtentX * head.BoxSize
#     RegionExtentZ          = RegionExtentZ * head.BoxSize
#
#     if(parameters.COLIBRE):
#         longfile = parameters.datadir+'/'+parameters.snap_base+'.hdf5'
#     else:
#         longfile = parameters.datadir+'/'+parameters.snap_base+'.0.hdf5'
#
#     snap = read_eagle.EagleSnapshot(longfile)
#     snap.select_region(RegionExtentX[0], RegionExtentX[1], RegionExtentY[0], RegionExtentY[1], RegionExtentZ[0], RegionExtentZ[1])
#
#     if (firstrun):
#         print("reading position")
#     Position =  snap.read_dataset(0, "Coordinates")
#
#     if (firstrun):
#         print("reading velocity")
#     Velocity = snap.read_dataset(0, "Velocity")
#
#     if (firstrun):
#         print("reading density")
#     ParticleDensity = snap.read_dataset(0, "Density")
#
#     if (firstrun):
#         print("reading SFR")
#     StarFormationRate = snap.read_dataset(0,"StarFormationRate")
#
#     if (firstrun):
#         print("Reading paricle smoothing length")
#     ParticleSmoothingLength = snap.read_dataset(0,"SmoothingLength")
#
#     ParticleNeutralHFraction = 0
#     ParticleMolecularHFraction = 0
#
#     if(firstrun):
#         print("reading temperatures")
#     ParticleTemperature = snap.read_dataset(0,"Temperature")
