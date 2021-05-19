import numpy as np



#Physical Cosntants

Planck           = 6.6260693e-27    # [ergs*s]
LightSpeed       = 2.99792458e10    # [cm / s]
Boltz            = 13.3806505e-16   # [erg/K]
ThomsonCross     = 6.65245873e-25   # [cm^2]
G                = 6.6742e-8        # [cm^3/g/s^2]
yr               = 3.1558e7         # [s]
M_sun            = 1.989e33         # [g]
Mpc              = 3.085678e24      # [cm]
H0               = 1e7/Mpc          # [s^-1 (100 km/s/Mpc)]
lyalpha          = 1215.6701        #Rest wavelength source: Morton 2003, ApJS, 149, 205



#Atomic data

#21 cm
A21             = 2.85e-15 # inverse of lifetime of 21 cm level [s^-1]
nu21            = 1.4204e9 # 21cm frequency [s^-1]
atomi_munit     = 1.6605388e-24 # Atomic mass unit [g]

#Hydrogen atom
massH          = 1.00794 * atomi_munit # Mass of hidrogen [g]
nlyman_all     = 31 # probably superflus in python implementation
Lambda_H1      = np.array([1215.6701, 1025.7223, 972.5368, 949.7431,937.8035,
                 930.7483, 926.2257, 923.1504, 920.9631, 919.3514, 918.1294,
                 917.1806, 916.429, 915.824, 915.329, 914.919, 914.576, 914.286,
                 914.039, 913.826,913.641, 913.480, 913.339, 913.215,  913.104,
                 913.006, 912.918, 912.839, 912.768, 912.703, 912.645 ])

f_H1           = np.array([0.416400, 0.079120, 0.029000, 0.013940,0.007799,
                0.004814, 0.003183, 0.002216, 0.001605, 0.00120,  0.000921,
                7.226e-4, 0.000577, 0.000469, 0.000386, 0.000321, 0.000270,
                0.000230, 0.000197, 0.000170, 0.000148, 0.000129, 0.000114,
                0.000101, 0.000089, 0.000080, 0.000071, 0.000064, 0.000058,
                0.000053, 0.000048 ])
#Helium
massHe        = 4.002602 * atomi_munit
Lambda_He2    = np.array([303.7822])
f_He2         = np.array([0.416])

#21cm
Lambda_21cm   = np.array([2.1e9])
f_21cm        = np.array(0.416)

#Carbon     [C]
massC         = 12.0107 * atomi_munit
Lambda_C2     = np.array([1334.5323,1036.3367])
f_C2          = np.array([0.127800,0.118000])
Lambda_C3     = np.array([977.0201])
f_C3          = np.array([0.7570])
Lambda_C4     = np.array([1548.2041,1550.7812])
f_C4          = np.array([0.189900,0.094750])
Lambda_C5     = np.array([40.2678,34.9728])
f_C5          = np.array([0.648,0.141])
Lambda_C6     = np.array([33.7342,33.7396,28.4652,28.4663])
f_C6          = np.array([0.277,0.139,0.0527,0.0263])

#Nitrogen     [N]
massN         = 14.0067 * atomi_munit
Lambda_N2     = np.array([1083.9937])
f_N2          = np.array([0.111 ])
Lambda_N3     = np.array([989.799])
f_N3          = np.array([0.12287])
Lambda_N4     = np.array([ 765.148])
f_N4          = np.array([0.632])
Lambda_N5     = np.array([1238.821, 1242.804])
f_N5          = np.array([ 0.156000, 0.0770])
Lambda_N6     = np.array([28.7875])
f_N6          = np.array([0.675])
Lambda_N7     = np.array([24.7792, 24.7846])
f_N7          = np.array([0.277, 0.139])

#Oxygen       [O]
massO         = 15.9994 * atomi_munit
Lambda_O1     = np.array([1302.1685, 988.7734, 971.7382])
f_O1          = np.array([0.048000, 0.0465, 0.0116])
Lambda_O3     = np.array([702.332, 832.927])
f_O3          = np.array([ 0.126, 0.0998])
Lambda_O4     = np.array([787.711])
f_O4          = np.array([0.110])
Lambda_O5     = np.array([629.730])
f_O5          = np.array([ 0.499])
Lambda_O6     = np.array([ 1031.9261, 1037.6167])
f_O6          = np.array([0.13250, 0.06580])
Lambda_O7     = np.array([21.6019, 18.6284])
f_O7          = np.array([0.696, 0.146])
Lambda_O8     = np.array([18.9671, 18.9725, 16.0055, 16.0067])
f_O8          = np.array([0.277, 0.139, 0.0527, 0.0263])

#Neon     [Ne]
massNe        = 20.1797 * atomi_munit
Lambda_Ne8    = np.array([770.409, 780.324])
f_Ne8         = np.array([0.103, 0.0505])
Lambda_Ne9    = np.array([13.4471])
f_Ne9         = np.array([0.724])

#Magnesium  [Mg]
massMg        = 24.3050 * atomi_munit
Lambda_Mg2    = np.array([2796.3543, 2803.5315])
f_Mg2         = np.array([0.6155, 0.3058])

#Aluminium [Al]
massAl        = 26.981538 * atomi_munit
Lambda_Al2    = np.array([1670.79])
f_Al2         = np.array([1.73812])
Lambda_Al3    = np.array([1854.72, 1862.79])
f_Al3         = np.array([0.559399, 0.277866])

#Silicon [Si]
massSi        = 28.0855 * atomi_munit
Lambda_Si2    = np.array([1260.420])
f_Si2         = np.array([1.17621])
Lambda_Si3    = np.array([1206.500])
f_Si3         = np.array([1.63])
Lambda_Si4    = np.array([1393.76018, 1402.77291])
f_Si4         = np.array([0.513, 0.254])

#Sulfur [S]
massS         = 32.065 * atomi_munit
Lambda_S5     = np.array([786.48])
f_S5          = np.array([1.46])

#Iron   [Fe]
massFe        = 55.8452 * atomi_munit
Lambda_Fe2    = np.array([1144.9379, 1608.45085, 1063.1764, 1096.8769,1260.533,
                1121.9748, 1081.8748, 1143.2260,1125.4477])
f_Fe2         = np.array([0.083, 0.0577, 0.0547, 0.032700,0.024000, 0.0290,
                0.012600, 0.0192,0.0156])
Lambda_Fe3    = np.array([1122.52])
f_Fe3         = np.array([0.0544257])
Lambda_Fe17   = np.array([15.0140, 15.2610])
f_Fe17        = np.array([2.72, 0.614])
Lambda_Fe19   = np.array([13.5180, 13.5146])
f_Fe19        = np.array([0.717, 0.0199])

Lambda_Fe21   = np.array([12.2840])
f_Fe21        = np.array([1.24])


#Quick fix, maybe theres a smarther as always, the last element is the index in the hdf5
atom_dic = {"h1":  [massH,Lambda_H1,f_H1,0],
            "he2": [massHe,Lambda_He2,f_He2,1],
            "c2":  [massC,Lambda_C2,f_C2,2],
            "c3":  [massC,Lambda_C3,f_C3,2],
            "c4":  [massC,Lambda_C4,f_C4,2],
            "n2":  [massN,Lambda_N2,f_N2,6],
            "n3":  [massN,Lambda_N3,f_N3,6],
            "n4":  [massN,Lambda_N4,f_N4,6],
            "n5":  [massN,Lambda_N5,f_N5,6],
            "o1":  [massO,Lambda_O1,f_O1,8],
            "o3":  [massO,Lambda_O3,f_O3,8],
            "o4":  [massO,Lambda_O4,f_O4,8],
            "o5":  [massO,Lambda_O5,f_O5,8],
            "o6":  [massO,Lambda_O6,f_O6,8],
            "o7":  [massO,Lambda_O7,f_O7,8],
            "mg2": [massMg,Lambda_Mg2,f_Mg2,5],
            "ne8": [massNe,Lambda_Ne8,f_Ne8,7],
            "al2": [massAl,Lambda_Al2,f_Al2],
            "al3": [massAl,Lambda_Al3,f_Al3],
            "si2": [massSi,Lambda_Si2,f_Si2,3],
            "si3": [massSi,Lambda_Si3,f_Si3,3],
            "si4": [massSi,Lambda_Si4,f_Si4,3],
            "s5":  [massS,Lambda_S5,f_S5],
            "fe2": [massFe,Lambda_Fe2,f_Fe2,4],
            "fe3": [massFe,Lambda_Fe3,f_Fe3,4],
            "21cm":[massH,Lambda_21cm,f_21cm,0]}




# Solar Data assumed values for the Sun,
  #Solar abundances source: Hazy I, CLOUDY 94 (Anders & Grevesse 1989; Grevesse & Noels 1993)


  #Total metal mass fraction
Zmass_solar   = 0.0126637     # M_Metal/M_tot
Ymass_solar   = 0.2466        # M_Helium/M_total
Xmass_solar   = 1.-Ymass_solar-Zmass_solar # M_Hydrogen/M_total
YNumber_solar = (Ymass_solar/Xmass_solar)*(massH/massHe) # N_Helium/N_Hydrogen

#Metallicities in number relative to hydrogen.
ZC_solar      = 3.55e-4
ZN_solar      = 9.33e-5
ZO_solar      = 7.41e-4
ZNe_solar     = 1.23e-4
ZMg_solar     = 3.80e-5
ZAl_solar     = 2.95e-6
ZSi_solar     = 3.55e-5
ZS_solar      = 1.62e-5
ZFe_solar     = 3.24e-5
