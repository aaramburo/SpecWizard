cimport numpy as npy
from read_parameter_file import read_params
from libc.math cimport log10
import numpy as np

def computeib(int z1,int z2,int iz1, int iz2, double dz1, double dz2, double logtemp, double logdens, npy.ndarray ib_logd,npy.ndarray  ionizbal, npy.ndarray ib_logt, int nz,int nznt):

    cdef int i, it1, it2, id1, id2
    cdef int i111, i211, i121, i112, i122, i212, i221, i222
    cdef double w111, w211, w121, w112, w122, w212, w221, w222
    cdef double dd1, dd2, dt1, dt2, logd, logt, ibfactor_use
    cdef double[:] get_fitted_ibfactor, fraction

    parameters = read_params('dummy.par')

    if (parameters.use_fitted_ibfactor):
      ibfactor_use = 1./(get_fitted_ibfactor[z1]*dz1)
    else:
      ibfactor_use = parameters.ibfactor

    logd = logdens - log10(ibfactor_use)

    #Ensure that the temperature is in the range of the ionization tables

    if (logtemp < ib_logt[0]):  #ib_logt is the one in the tables
      logt = ib_logt[0]

    elif (logtemp > ib_logt[-1]):
      logt = ib_logt[-1]

    else:
      logt = logtemp

   # Same for density

    if (logd < ib_logd[0]):
      logd = ib_logd[0]
    elif (logt > ib_logd[-1]):
      logd = ib_logd[-1]
    else:
      logd = logdens

    it2 = 1
    while(logt > ib_logt[it2]):
      it2 = it2 + 1
    it1 = it2 -1
    dt1 = (ib_logt[it2] - logt) / (ib_logt[it2] - ib_logt[it1])
    dt2 = 1. - dt1

    id2 = 1
    while(logd > ib_logd[id2]):
      id2 = id2 + 1
    id1 = id2 - 1
    dd1 = (ib_logd[id2]-logd) / (ib_logd(id2) - ib_logd(id1))
    dd2 = 1. - dd1

    # Weights
    w111 = dz1 * dt1 * dd1
    w211 = dz2 * dt1 * dd1
    w121 = dz1 * dt2 * dd1
    w221 = dz2 * dt2 * dd1
    w112 = dz1 * dt1 * dd2
    w212 = dz2 * dt1 * dd2
    w122 = dz1 * dt2 * dd2
    w222 = dz2 * dt2 * dd2

    #indx
    i111 = iz1 + (it1-1)*nz + (id1-1)*nznt
    i211 = iz2 + (it1-1)*nz + (id1-1)*nznt
    i121 = iz1 + (it2-1)*nz + (id1-1)*nznt
    i221 = iz2 + (it2-1)*nz + (id1-1)*nznt
    i112 = iz1 + (it1-1)*nz + (id2-1)*nznt
    i212 = iz2 + (it1-1)*nz + (id2-1)*nznt
    i122 = iz1 + (it2-1)*nz + (id2-1)*nznt
    i222 = iz2 + (it2-1)*nz + (id2-1)*nznt

    fraction[:] =  w111*ionizbal[:,i111] + w211*ionizbal[:,i211] + w121*ionizbal[:,i121] + w221*ionizbal[:,i221] + w112*ionizbal[:,i112] + w212*ionizbal[:,i212] + w122*ionizbal[:,i122] + w222*ionizbal[:,i222]
    return np.power(10,fraction)
