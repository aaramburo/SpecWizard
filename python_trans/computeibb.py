from read_parameter_file import read_params
import numpy as np

def computeibb(z1, z2, iz1, iz2, dz1, dz2, log_temp, log_dens, ib_logd, ionizbal, ib_logt,nz, nznt):


    parameters = read_params('dummy.par')

    if (parameters.use_fitted_ibfactor):
        ibfactor_use = 1./(get_fitted_ibfactor[z1]*dz1)
    else:
        ibfactor_use = parameters.ibfactor

    logd = log_dens - np.log10(ibfactor_use)
    
#    print("logd:"+str(logd))
    

    #Ensure that the temperature is in the range of the ionization tables

    if (log_temp < ib_logt[0]):  #ib_logt is the one in the tables
        logt = ib_logt[0]

    elif (log_temp > ib_logt[-1]):
        logt = ib_logt[-1]

    else:
        logt = log_temp
#    print("logt:"+str(logt))

        
   # Same for density

#    if (logd < ib_logd[0]):
#        logd = ib_logd[0]
#    elif (logt > ib_logd[-1]):
#        logd = ib_logd[-1]
#    else:
#        logd = log_dens
        
#    print("logd:"+str(logd))


    it2 = 1
    while(logt > ib_logt[it2]):
        it2 = it2 + 1
    it1 = it2 -1
    dt1 = (ib_logt[it2] - logt) / (ib_logt[it2] - ib_logt[it1])
    dt2 = 1. - dt1

#    print("dt1:"+str(dt1))
#    print("dt2:"+str(dt2))

    id2 = 1
    while(logd > ib_logd[id2]):
        id2 = id2 + 1
    id1 = id2 - 1
    dd1 = (ib_logd[id2]-logd) / (ib_logd[id2] - ib_logd[id1])
    dd2 = 1. - dd1
    
#    print("dd1:"+str(dd1))
#    print("dd2:"+str(dd2))
    

    # Weights
    w111 = dz1 * dt1 * dd1
    w211 = dz2 * dt1 * dd1
    w121 = dz1 * dt2 * dd1
    w221 = dz2 * dt2 * dd1
    w112 = dz1 * dt1 * dd2
    w212 = dz2 * dt1 * dd2
    w122 = dz1 * dt2 * dd2
    w222 = dz2 * dt2 * dd2

#    print("w111:"+str(w111))
#    print("w211:"+str(w211))
#    print("w121:"+str(w121))
#    print("w221:"+str(w221))
#    print("w112:"+str(w112))
#    print("w212:"+str(w212))
#    print("w122:"+str(w122))
#    print("w222:"+str(w222))
    
    #indx
    fxd = 4018
    i111 = iz1 + (it1-1)*nz + (id1-1)*nznt +fxd
    i211 = iz2 + (it1-1)*nz + (id1-1)*nznt +fxd 
    i121 = iz1 + (it2-1)*nz + (id1-1)*nznt +fxd 
    i221 = iz2 + (it2-1)*nz + (id1-1)*nznt +fxd  
    i112 = iz1 + (it1-1)*nz + (id2-1)*nznt +fxd  
    i212 = iz2 + (it1-1)*nz + (id2-1)*nznt +fxd 
    i122 = iz1 + (it2-1)*nz + (id2-1)*nznt +fxd 
    i222 = iz2 + (it2-1)*nz + (id2-1)*nznt +fxd   
    
    
#    print("i111:"+str(i111))
#    print("i211:"+str(i211))
#    print("i121:"+str(i121))
#    print("i221:"+str(i221))
#    print("i112:"+str(i112))
#    print("i212:"+str(i212))
#    print("i122:"+str(i122))
#    print("i222:"+str(i222))
    

    fraction =  w111*ionizbal[:,i111] + w211*ionizbal[:,i211] + w121*ionizbal[:,i121] + w221*ionizbal[:,i221] + w112*ionizbal[:,i112] + w212*ionizbal[:,i212] + w122*ionizbal[:,i122] + w222*ionizbal[:,i222]

    return np.power(10,fraction)
