def Interp3D(indata, z,t,dens, value_array, verbose = 0, order=1 ):                                                                                                                                                                                                 
    lo  = np.array([dens[0],t[0],z[0]])            
    hi  = np.array([dens[-1],t[-1],z[-1]])                                                                                                                                                 
    theMaps = [ None, None,z[:]];                
    
    interfunc = Intergrid( indata, lo=lo, hi=hi, maps= theMaps, prefilter=True, verbose = verbose,order=order );  
    
    return interfunc.at( value_array );


def projectdata_scypi():
    n_ion     =  np.zeros((nion,nveloc))
    veloc_ion =  np.zeros((nion,nveloc))
    temp_ion  =  np.zeros((nion,nveloc))
    rho_ion   =  np.zeros((nion,nveloc))
    rho_tot   =  np.zeros(nveloc)
    temp_tot   =  np.zeros(nveloc)
    met_tot   =  np.zeros(nveloc)
    veloc_tot   =  np.zeros(nveloc)

    ionss = []
    first_iteration = True

    for ion_file in ionbal_to_use:
        with h5.File(parameters.ibdir+ion_file+".hdf5","r") as f:
            if first_iteration:
                z_ranges_table = f["/redshift"][...]
                logt_table = f["/logt"][...]
                logd_table = f["/logd"][...]
                first_iteration = False

            ionss.append(f["/ionbal"][...])
            print(ion_file+' Ionization table loaded')


    zarray = zcurrent*np.ones(len(log_temp))                           
    iinput = np.array([log_dens,log_temp,zarray]).T

    h2 = ParticleSmoothingLength**2
    hinv2 = np.power(h2,-1)
    hh      = ParticleSmoothingLength
    hinv3 = hinv2 / ParticleSmoothingLength

    dx = abs(Position_2[:,0] - x_physical)
    dy = abs(Position_2[:,1] - y_physical)

    dx[dx > box_2] = box - dx[dx > box_2]
    dy[dy > box_2] = box - dy[dy > box_2]

    b2 = dx*dx + dy*dy
    b = impactparameter = np.power(b2,0.5)

    particle_mask = np.where(impactparameter <= hh)

    Density = ParticleDensity[particle_mask] * densscale

    log_dens = np.log10(Density)
    log_temp = np.log10(ParticleTemperature[particle_mask])


    des = [Interp3D(ionss[ll],z_ranges_table,logt_table,logd_table,iinput,verbose=0, order=1) for ll in range(len(ionss))]
    ionfrac = np.transpose(des)
    totnr_ion = ionfrac*loss.MassFractions[:,ion_elnr[:]][particle_mask,:][0] * Mass[particle_mask][:, None] / ElementAtomicMass[:]

    j = [list(range(iizi[i],iizg[i])) for i in range(len(iizi))]

    for i,jj in enumerate(j):
        chunks = np.mod(np.array(jj)+10*nveloc,nveloc)
        zzgrid = zmingrid + (chunks - 0.5)*dzgrid
        deltaz = abs(zzgrid - zz[i])

        deltaz[deltaz > box_2] = box - deltaz[deltaz > box_2] 

        dr2 = b2[i] + deltaz**2
        zf = deltaz + dzgrid*0.5
        zi = deltaz - dzgrid*0.5
        q = np.sqrt(dr2 * hinv2[i])

        kernel_factor = np.zeros(np.shape(q))
        q_cond0p5 = q<0.5
        q_cond1   = q<1
        kernel_factor[q_cond0p5] = (1.+6*q[q_cond0p5]**2*(q[q_cond0p5]-1))
        kernel_factor[(q_cond0p5) & (q_cond1)] = 2.*(1-q[(q_cond0p5) & (q_cond1)])**3

        kernel_factor = (kernel_factor * 8. * hinv3[i]) / np.pi
        krnl_cond = kernel_factor>0

        n_ion[:,chunks[krnl_cond]]      = n_ion[:,chunks[krnl_cond]]     + kernel_factor[krnl_cond] * totnr_ion[i,:,None]
        veloc_ion[:,chunks[krnl_cond]]  = veloc_ion[:,chunks[krnl_cond]] + kernel_factor[krnl_cond] * totnr_ion[i,:,None] * vr
        temp_ion[:,chunks[krnl_cond]]   = temp_ion[:,chunks[krnl_cond]]  + kernel_factor[krnl_cond] * totnr_ion[i,:,None] * loss.ParticleTemperature[i]

        rho_ion[:,chunks[krnl_cond]]    = rho_ion[:,chunks[krnl_cond]]   + kernel_factor[krnl_cond] * Density[i] * const.atomi_munit/ loss.MassFractions[i,0]

        rho_tot[chunks[krnl_cond]]      = rho_tot[chunks[krnl_cond]]     + kernel_factor[krnl_cond] * Mass[i]
        veloc_tot[chunks[krnl_cond]]    = veloc_tot[chunks[krnl_cond]]   + kernel_factor[krnl_cond] * Mass[i] * vr
        temp_tot[chunks[krnl_cond]]     = temp_tot[chunks[krnl_cond]]    + kernel_factor[krnl_cond] * Mass[i] * loss.ParticleTemperature[i]

        met_tot[chunks[krnl_cond]]      = met_tot[chunks[krnl_cond]]     + kernel_factor[krnl_cond] * Mass[i] * loss.Metallicity[i]
