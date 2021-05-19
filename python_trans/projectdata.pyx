cimport numpy as npy
import numpy as np
import physical_data as const
import sys
from read_parameter_file import read_params
import computeib

def projectdata(npy.ndarray Position, npy.ndarray Velocity, npy.ndarray ParticleSmoothingLength, npy.ndarray ParticleDensity,
npy.ndarray ParticleTemperature,npy.ndarray ParticleNeutralHFraction, npy.ndarray ParticleMolecularHFraction,
npy.ndarray Metallicity ,npy.ndarray ib_redshift, npy.ndarray ionizbal, npy.ndarray MassFractions , npy.ndarray Mass,
npy.ndarray ElementAtomicMass , npy.ndarray ion_elnr, float x_physical, float y_physical,
float zcurrent, float BoxSize, float ExpansionFactor, npy.ndarray ib_logd, npy.ndarray ib_logt, float HubbleParam,
int nveloc, int NGas, long nion, int nz, int nznt):

    cdef int i,ioff,j,ii,ion,iz1,iz2,h1_index,si2_index,iiz,iz
    cdef double xx,yy,zz,hh,h2,b,b2,hinv2,hinv3,vr,zmingrid,zmaxgrid,dzinv,dzgrid
    cdef double box,box_2,dzmax,zgrid,zf,zedge,ztrans,dr2,qi,qf,deltaz,dvbin,dz1
    cdef double dz2,Density,log_dens,log_temp, DensCon, impactparameter,q,z,kernel_factor
    cdef double dx,dy,Q1,Q2,temp_val
    cdef double RotationMatrix[3][3]
    cdef double[:] totnr_ion , rho_tot, veloc_tot, temp_tot, met_tot
    cdef double[:,:] n_ion,veloc_ion, temp_ion, rho_ion
    parameters = read_params('dummy.par')

    if (nz > 1):
        if (not parameters.use_maxdens_above_zmax):
            if (zcurrent < ib_redshift[0] or zcurrent > ib_redshift[-1]):
                print("ERROR: Z out of bounds of the ion table!")
                sys.exit()
        if(parameters.use_maxdens_above_zmax > ib_redshift[-1] and zcurrent > ib_redshift[-1]):
            iz1 = nz -1
            iz2 = nz
            dz1 = 0
            dz2 = 1.0
        else:
            iz2=1
            while(z > ib_redshift[iz2]):
                iz2 = iz2 =1
            iz1 = iz2-1
            dz1 = (ib_redshift[iz2] - zcurrent) / (ib_redshift[iz2] - ib_redshift[iz1])
            dz2 = 1. -dz1
    else:
        print("WARNING: Ion tables only have  1 redshift!")
        iz1 = 1
        iz2 = iz1
        dz1 = 0.5
        dz2 = 1. - dz1

    zmingrid = 0
    zmaxgrid = (BoxSize*ExpansionFactor) / HubbleParam
    box      = zmaxgrid
    dzgrid   = (zmaxgrid - zmingrid) / nveloc
    box_2    = 0.5 * box
    densscale = (ExpansionFactor/ExpansionFactor)**3
    for i in range(NGas):
      ##    cdef double[:] p       /* (how to define this in cpython) */

        xx      = Position[1,i]
        yy      = Position[2,i]
        zz      = Position[3,i]

        hh      = ParticleSmoothingLength[i]
        h2      = hh*hh
        hinv2   = 1. / h2
        hinv3   = hinv2 / hh

        dx      = abs(xx - x_physical)
        dy      = abs(yy - y_physical)

        if(dx > box_2):
            dx  = box - dx

        if(dy > box_2):
            dy = box - dy

        b2    = dx*dx + dy*dy
        b     = (b2)**0.5
        impactparameter = b

        if( impactparameter <= hh):

            ncontr = ncontr + 1

            #Density

            Density = ParticleDensity[i]

            if(Density < 0):
                print("ERROR:  Negative Density")

            Density = Density * densscale
            log_dens = np.log10(Density)

            log_temp = np.log10(ParticleTemperature[i])

            ionfrac = computeib(ib_redshift[iz1],ib_redshift[iz2],iz1, iz2, dz1, dz2, log_temp, log_dens, ib_logd, ionizbal, ib_logt, nz, nznt)

            if(parameters.urchin or parameters.aurora):
                  if(parameters.subtract_Hmol):
                      if(parameters.doH1):
                          ionfrac[h1_index] = ParticleNeutralHFraction[i] * (1. - ParticleMolecularHFraction[i])
                      if(parameters.doSi2):
                          ionfrac[si2_index] = ParticleNeutralHFraction[i] * (1. - ParticleMolecularHFraction[i])
                  else:
                      if(parameters.doH1):
                          ionfrac[h1_index] = ParticleNeutralHFraction[i]
                      if(parameters.doSi2):
                          ionfrac[si2_index] = ParticleNeutralHFraction[i]

            if(parameters.ionfracone):
                totnr_ion[:] = MassFractions[ion_elnr[:],i] * Mass[i] / ElementAtomicMass[ion_elnr[:]]
            else:
                totnr_ion[:] = ionfrac[:] * MassFractions[ion_elnr[:],i] * Mass[i] / ElementAtomicMass[ion_elnr[:]]


            if(parameters.NoPecVel):
                vr = 0.0
            else:
                vr = Velocity[2,i] # peculiar velocity in km/s !!!

    # central pixel to contribute to
            iz = np.int((zz - zmingrid) * dzinv + 1.0)
            #
            # if(parameters.gimic):
            #     if(boundary[i] == 1):
            #         spectrum_boundary[iz]=spectrum_boundary[iz]+1.0

            dzmax = np.sqrt(np.abs(h2 - b2))
            ioff =  np.int(dzmax * dzinv) + 2.0

  # if (ioff*2 > nveloc) stop 'ioff*2 > nveloc'

   # segment pixel loop
            for iiz in range(iz-ioff, iz+ioff+1):
                j = iiz
                j = np.mod(j-1+10*nveloc,nveloc)+1
                zgrid = zmingrid + (j-0.5)*dzgrid
                deltaz = np.abs(zgrid - zz)
                if (deltaz > box_2):
                    deltaz = box - deltaz
                dr2 = b2 + deltaz**2
                zf = deltaz + dzgrid*0.5
                zi = deltaz - dzgrid*0.5

                if(parameters.integrate_kernel):
                    if(parameters.use_gaussian_kernel):
                        kernel_factor   = 4#integrate_G3_line() # function to be defined
                    else:
                        print("lol")
                else:
                    if(parameters.use_gaussian_kernel):
                        print("vamos a baialr algo que esta perron")#kernel_factor  = gaussian_kernel() # Fuctnion TBD
                    else: #use cubic sppline
                        q = (dr2 *hinv2)**(0.5)
                        if (q < 0.5):
                            kernel_factor = (1.6+6.*q**2 *(q-1.))
                        elif(q < 1):
                            kernel_factor = 2.*(1.-q)**3
                        else:
                            kernel_facor = 0
                        kernel_factor = (kernel_factor * 8. * hinv3) / np.pi

                    # if(parameters.gimic):
                    #     if(boundary[i] == 1):
                    #         kernel_factor = 0

                    if(kernel_factor > 0):
                        n_ion[:,j]        = n_ion[:,j] + np.multiply(kernel_factor, totnr_ion[:])
                        veloc_ion[:,j]    = veloc_ion[:,j] + np.multiply(kernel_factor* vr, totnr_ion[:])
                        temp_ion[:,j]     = temp_ion[:,j]  + np.muliply(kernel_factor* ParticleTemperature[i],totnr_ion[:])
                        # Particle Density n_H = cgs, rescaled for long spectra -> Particle cgs rho
                        temp_val =  kernel_factor * Density * const.atomi_munit / MassFractions[0,i]
                        rho_ion[:,j]      = rho_ion[:,j]   + np.multiply(temp_val,totnr_ion[:] )
                        # .... weighted by mass
                        rho_tot[j]     = rho_tot[j]     + kernel_factor * Mass[i]
                        veloc_tot[j]   = veloc_tot[j]   + kernel_factor * Mass[i] * vr
                        temp_tot[j]    = temp_tot[j]    + kernel_factor * Mass[i] * ParticleTemperature[i]
                        met_tot[j]     = met_tot[j]     + kernel_factor * Mass[i] * Metallicity[i]

    DensCon = const.Msun / const.Mpc_to_cm**3
    DensCon = DensCon * densscale

    for ii in range(nion):
        for i in range(nveloc):
            if(n_ion[ii,i] > 0):
                veloc_ion[ii,i] = veloc_ion[ii,i] / n_ion[ii,i]
                temp_ion[ii,i]  = temp_ion[ii,i] /  n_ion[ii,i]
                rho_ion[ii,i]   = rho_ion[ii,i] / n_ion[ii,i]
            n_ion[ii,i] = n_ion[ii,i] * DensCon
    for i in range(nveloc):
        if(rho_tot[i]>0):
            veloc_tot[i] = veloc_tot[i] / rho_tot[i]
            temp_tot[i] = temp_tot[i] / rho_tot[i]
            met_tot[i]  = met_tot[i]  / rho_tot[i]
            rho_tot[i]  = rho_tot[i] * DensCon
