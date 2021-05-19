import numpy as np
import physical_data as const
import sys
from read_parameter_file import read_params
from computeibb import computeibb as computeib

def projectdata(Position,Velocity,ParticleSmoothingLength, ParticleDensity, ParticleTemperature, ParticleNeutralHFraction,ParticleMolecularHFraction,
Metallicity ,ib_redshift, ionizbal,  MassFractions , Mass, ElementAtomicMass, ion_elnr, x_physical, y_physical, zcurrent, BoxSize, ExpansionFactor, ib_logd, ib_logt,  HubbleParam,
nveloc, NGas, nion, nz, nznt):

    parameters = read_params('dummy.par')
    n_ion     =  np.zeros((nion,nveloc))
    veloc_ion =  np.zeros((nion,nveloc))
    temp_ion  =  np.zeros((nion,nveloc))
    rho_ion   =  np.zeros((nion,nveloc))
    veloc_tot =  np.zeros(nveloc)
    rho_tot   =  np.zeros(nveloc)
    met_tot   =  np.zeros(nveloc)
    met_tot   =  np.zeros(nveloc)
    temp_tot  =  np.zeros(nveloc)

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
            while(zcurrent > ib_redshift[iz2]):
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
    dzinv    = 1. / dzgrid
    box_2    = 0.5 * box
    densscale = (ExpansionFactor/ExpansionFactor)**3
    ncontr = 0
    for i in range(NGas):
      ##    cdef double[:] p       /* (how to define this in cpython) */

        xx      = Position[i,0]
        yy      = Position[i,1]
        zz      = Position[i,2]

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
            urchin = False
            aurora = False
            if(urchin or aurora):
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

            if(False):#parameters.ionfracone):
                totnr_ion[:] =  MassFractions[i,ion_elnr] * Mass[i] / ElementAtomicMass
            else:
        
                totnr_ion = ionfrac[:] * MassFractions[i,ion_elnr] * Mass[i] / ElementAtomicMass


            if(parameters.NoPecVel):
                vr = 0.0
            else:
                vr = Velocity[i,2] # peculiar velocity in km/s !!!

    # central pixel to contribute to
            iz = np.int((zz - zmingrid) * dzinv + 1.0)
            #
            # if(parameters.gimic):
            #     if(boundary[i] == 1):
            #         spectrum_boundary[iz]=spectrum_boundary[iz]+1.0

            dzmax = np.sqrt(np.abs(h2 - b2))
            ioff =  np.int(dzmax * dzinv) + 2

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

                if(False):#parameters.integrate_kernel):
                    if(parameters.use_gaussian_kernel):
                        kernel_factor   = 4#integrate_G3_line() # function to be defined
                    else:
                        print("lol")
                else:
                    if(False):#parameters.use_gaussian_kernel):
                        print("vamos a baialr algo que esta perron")#kernel_factor  = gaussian_kernel() # Fuctnion TBD
                    else: #use cubic sppline
                        q = (dr2 *hinv2)**(0.5)
                        if (q < 0.5):
                            kernel_factor = (1.6+6.*q**2 *(q-1.))
                        elif(q < 1):
                            kernel_factor = 2.*(1.-q)**3
                        else:
                            kernel_factor = 0
                        kernel_factor = (kernel_factor * 8. * hinv3) / np.pi

                    # if(parameters.gimic):
                    #     if(boundary[i] == 1):
                    #         kernel_factor = 0

                    if(kernel_factor > 0):
                        n_ion[:,j]        = n_ion[:,j] + np.multiply(kernel_factor, totnr_ion[:])
                        veloc_ion[:,j]    = veloc_ion[:,j] + np.multiply(kernel_factor* vr, totnr_ion[:])
                        temp_ion[:,j]    = temp_ion[:,j]  + np.multiply(kernel_factor* ParticleTemperature[i],totnr_ion[:])
                        # Particle Density n_H = cgs, rescaled for long spectra -> Particle cgs rho
                        temp_val =  kernel_factor * Density * const.atomi_munit / MassFractions[i,0]
                        rho_ion[:,j]     = rho_ion[:,j]   + np.multiply(temp_val,totnr_ion[:] )
                        # .... weighted by mass
                        rho_tot[j]     = rho_tot[j]     + kernel_factor * Mass[i]
                        veloc_tot[j]   = veloc_tot[j]   + kernel_factor * Mass[i] * vr
                        temp_tot[j]    = temp_tot[j]    + kernel_factor * Mass[i] * ParticleTemperature[i]
                        met_tot[j]     = met_tot[j]     + kernel_factor * Mass[i] * Metallicity[i]

    DensCon = const.M_sun / const.Mpc**3
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
    return veloc_tot,temp_tot,met_tot,rho_tot
