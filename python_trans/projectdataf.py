import physical_data as const
import numpy as np
import time
from computeibb import computeibb as computeib


def redshift_interp_param(nz,ib_redshift,zcurrent,parameters):
    
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
            dz2 = 1. - dz1

    return iz1,iz2,dz1,dz2 




def project_data(los,simdata,parameters,ionpar):
    t1 = time.time()

    ionbal_to_use  = ionpar.ionbal_to_use
    
    ElementAtomicMass = [const.atom_dic[item][0] for item in ionbal_to_use]
    MassFraction_indx = [const.atom_dic[item][-1] for item in ionbal_to_use]
    
    ParticleNeutralHFraction = los.ParticleNeutralHFraction
    ParticleMolecularHFraction = los.ParticleMolecularHFraction
    
    
    NGas = len(los.Mass)
    ib_logt = ionpar.logt_table
    ib_logd = ionpar.logd_table
    ib_redshift = ionpar.z_ranges_table
    ionizbal = np.log10(ionpar.ionss)
    ion_elnr = MassFraction_indx 
    nion = len(ionbal_to_use)
    nz = len(ionpar.z_ranges_table)
    nznt = nz*len(ib_logt)
    nveloc = simdata.nveloc
    
    n_ion, veloc_ion, temp_ion, rho_ion = [np.zeros((nion,nveloc)) for i in range(4)]
 
    rho_tot, temp_tot, met_tot, veloc_tot = [np.zeros(nveloc) for i in range(4)] 
        
    iz1,iz2,dz1,dz2 = redshift_interp_param(nz,ib_redshift,simdata.Redshift,parameters)
    
    #Usefull quantities   
    zmingrid = 0
    zmaxgrid = (simdata.BoxSize*simdata.ExpansionFactor) / simdata.HubbleParam
    box      = zmaxgrid
    dzgrid   = (zmaxgrid - zmingrid) / nveloc
    dzinv    = 1. / dzgrid 
    box_2    = 0.5 * box
    densscale = (simdata.ExpansionFactor/simdata.ExpansionFactor)**3
    ncontr= 0
    
    h2 = los.ParticleSmoothingLength**2
    hinv2 = np.power(h2,-1)
    hh      = los.ParticleSmoothingLength
    hinv3 = hinv2 / los.ParticleSmoothingLength

    dx = abs(los.Position[:,0] - los.x_physical)
    dy = abs(los.Position[:,1] - los.y_physical)

    zz = los.Position[:,2]
    dx[dx > box_2] = box - dx[dx > box_2]
    dy[dy > box_2] = box - dy[dy > box_2]

    b2 = dx*dx + dy*dy
    b = impactparameter = np.power(b2,0.5)

    particle_mask = np.where(impactparameter <= hh)

    #if parameters.NoPecVel:
    vr = 0
    #else:
    #    vr = Velocity[i,los_long_axis]

    Density = los.ParticleDensity[particle_mask] * densscale

    log_dens = np.log10(Density)
    log_temp = np.log10(los.ParticleTemperature[particle_mask])

    ionfrac = [computeib(ib_redshift[iz1],ib_redshift[iz2],iz1, iz2, dz1, dz2, log_temp[i], log_dens[i], ib_logd, ionizbal,                    ib_logt, nz, nznt,parameters.use_fitted_ibfactor,parameters.ibfactor) for i in range(len(log_temp))]


    iz = np.array((zz - zmingrid) * dzinv + 1).astype(int)

    #contribute to projection segment

    dzmax = np.sqrt(np.abs(h2 - b2))

    ioff  = np.array(dzmax * dzinv).astype(int) + 1


    iizg  = iz+ioff
    iizi  = iz-ioff



    totnr_ion = ionfrac*los.MassFractions[:,ion_elnr[:]][particle_mask,:][0] * los.Mass[particle_mask][:, None] / ElementAtomicMass[:]

    chunks = [list(range(iizi[i],iizg[i])) for i in range(len(iizi))]

    for i,chunk in enumerate(chunks):

        k = np.mod(np.array(chunk)+10*nveloc,nveloc)

        zzgrid = zmingrid + (k - 0.5)*dzgrid
        deltaz = abs(zzgrid - zz[i])

        deltaz[deltaz > box_2] = box - deltaz[deltaz > box_2] 
        dr2 = b2[i] + deltaz**2
        zf = deltaz + dzgrid*0.5
        zi = deltaz - dzgrid*0.5
        q = np.sqrt(dr2 * hinv2[i])

        kernel_factor = np.zeros(np.shape(q))
        kernel_factor[q<0.5] = (1.+6*q[q<0.5]**2*(q[q<0.5]-1))
        kernel_factor[(0.5 < q) & (q < 1)] = 2.*(1-q[(0.5 < q) & (q < 1)])**3

        kernel_factor = (kernel_factor * 8. * hinv3[i]) / np.pi
        krn_cond = kernel_factor>0
        kernel_chunks = k[krn_cond]
        
        n_ion[:,kernel_chunks]      = n_ion[:,kernel_chunks]     + kernel_factor[krn_cond] * totnr_ion[i,:,None]
        veloc_ion[:,kernel_chunks]  = veloc_ion[:,kernel_chunks] + kernel_factor[krn_cond] * totnr_ion[i,:,None] * vr
        temp_ion[:,kernel_chunks]   = temp_ion[:,kernel_chunks]  + kernel_factor[krn_cond] * totnr_ion[i,:,None] *                       los.ParticleTemperature[i]

        rho_ion[:,kernel_chunks]    = rho_ion[:,kernel_chunks]   + kernel_factor[krn_cond] * Density[i] * const.atomi_munit /       los.MassFractions[i,0]
        
        rho_tot[kernel_chunks]      = rho_tot[kernel_chunks]     + kernel_factor[krn_cond] * los.Mass[i]
        veloc_tot[kernel_chunks]    = veloc_tot[kernel_chunks]   + kernel_factor[krn_cond] * los.Mass[i] * vr
        temp_tot[kernel_chunks]     = temp_tot[kernel_chunks]    + kernel_factor[krn_cond] * los.Mass[i] *                               los.ParticleTemperature[i]
        met_tot[kernel_chunks]      = met_tot[kernel_chunks]     + kernel_factor[krn_cond] * los.Mass[i] * los.Metallicity[i]


    t2 = time.time()
    tf = t2-t1
    print("It took %f seconds" % (tf))
    
    los.n_ion      = nion
    los.veloc_ion  = veloc_ion
    los.temp_ion   = temp_ion
    los.rho_ion    = rho_ion
    los.rho_tot    = rho_tot
    los.veloc_tot  = veloc_tot
    los.temp_tot   = temp_tot
    los.met_tot    = met_tot
    
    return los,tf
    
    