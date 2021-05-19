class read_params(object):
    """
    Read the parameter file.
    """



    def __init__(self, par_file):
        def parsebool(str):
            str = str.strip()
            if str in ['true', 'True', 't', 'T', '1']:
                out = True
            elif str in ['false', 'False', 'f', 'F', '0']:
                out = False
            else:
                raise ValueError('Could not parse string {} as a boolean'.format(str))
            return out

        with open(par_file, 'r') as input_file:
            for line in input_file:
                try:
                    row = line.split()
                    if len(row)==0 or list(row)[0][0] == "#":
                        continue
                    label = row[0]
                    data = row[1:]  # rest of row is data list
                    attr = label_attr_map[label][0]
                    datatypes = label_attr_map[label][1:]

                    values = [parsebool(data[i]) if datatypes[i] == bool else\
                                    (datatypes[i](data[i])) for i in range(len(data))]
                    self.__dict__[attr] = values if len(values) > 1 else values[0]
                except:
                    print('Error at parameter file in : '+line)

label_attr_map = {
        "ibdir:": ["ibdir",str],
        "datadir:": [ "datadir",str],
        "file_list:": ["file_list",str],
        "outputdir:": ["outputdir",str],
        "gimic:": ["gimic",bool],
        "aurora:": ["aurora",bool],
        "COLIBRE:": ["COLIBRE",bool],
        "use_snapshot_file:": ["use_snapshot_file",bool],
        "snap_base:": ["snap_base",str],
        "use_random_los:": ["use_random_los",bool],
        "use_los_file:": ["use_los_file",bool],
        "los_coordinates_file:": ["los_coordinates_file",str],
        "do_long_spectrum:": ["do_long_spectrum",bool],
        "NoPecVel:": ["NoPecVel",bool],
        "overwrite:": ["overwrite",bool],
        "nspec:": ["nspec",int],
        "output_zspaceopticaldepthweighted_values:": ["output_zspaceopticaldepthweighted_values",bool],
        "output_realspacemassweighted_values:": ["output_realspacemassweighted_values",bool],
        "output_realspacenionweighted_values:": ["output_realspacenionweighted_values",bool],
        "output_frequency:": ["output_frequency",bool],
        "impose_eos:": ["impose_eos",bool],
        "imposed_eos_T0:": ["imposed_eos_T0",float],
        "imposed_eos_gamma:": ["imposed_eos_gamma",float],
        "imposed_eos_maxod:": ["imposed_eos_maxod",float],
        "doH1:": ["doH1",bool],
        "doHe2:": ["doHe2",bool],
        "doC2:": ["doC2",bool],
        "doC3:": ["doC3",bool],
        "doC4:": ["doC4",bool],
        "doN2:": ["doN2",bool],
        "doN3:": ["doN3",bool],
        "doN4:": ["doN4",bool],
        "doN5:": ["doN5",bool],
        "doO1:": ["doO1",bool],
        "doO3:": ["doO3",bool],
        "doO4:": ["doO4",bool],
        "doO5:": ["doO5",bool],
        "doO6:": ["doO6",bool],
        "doO7:": ["doO7",bool],
        "doMg2:": ["doMg2",bool],
        "doNe8:": ["doNe8",bool],
        "doAl2:": ["doAl2",bool],
        "doAl3:": ["doAl3",bool],
        "doSi2:": ["doSi2",bool],
        "doSi3:": ["doSi3",bool],
        "doSi4:": ["doSi4",bool],
        "doS5:": ["doS5",bool],
        "doFe2:": ["doFe2",bool],
        "doFe3:": ["doFe3",bool],
        "do21cm:": ["do21cm",bool],
        "doall:": ["doall",bool],
        "ibfactor:": ["ibfactor",float],
        "use_fitted_ibfactor:": ["use_fitted_ibfactor",bool],
        "ibfactor_he_reionization:": ["ibfactor_he_reionization",bool],
        "use_maxdens_above_zmax:": ["use_maxdens_above_zmax",bool],
        "modify_metallicity:": ["modify_metallicity",bool],
        "maxz_rel:": ["maxz_rel",float],
        "scale_simulation_abundances:": ["scale_simulation_abundances",bool],
        "z_rel:": ["z_rel",float],
        "impose_z_rho_relation:": ["impose_z_rho_relation",bool],
        "z_index:": ["z_index",float],
        "z_mean:": ["z_mean",float],
        "log_normal_scatter:": ["log_normal_scatter",bool],
        "z_sig_bin:": ["z_sig_bin",int],
        "z_sig_dex:": ["z_sig_dex",float],
        "ZC_rel:": ["ZC_rel",float],
        "ZN_rel:": ["ZN_rel",float],
        "ZO_rel:": ["ZO_rel",float],
        "ZNe_rel:": ["ZNe_rel",float],
        "ZMg_rel:": ["ZMg_rel",float],
        "ZAl_rel:": ["ZAl_rel",float],
        "ZSi_rel:": ["ZSi_rel",float],
        "ZS_rel:": ["ZS_rel",float],
        "ZFe_rel:": ["ZFe_rel",float],
        "read_part_ids_from_file:": ["read_part_ids_from_file",bool],
        "flagged_particle_metallicity:": ["flagged_particle_metallicity",float],
        "generate_noise:": ["generate_noise",bool],
        "use_noise_file:": ["use_noise_file",float],
        "noisefile:": ["noisefile",str],
        "sigtonoise:": ["sigtonoise",float],
        "minnoise:": ["minnoise",float],
        "minbother_blue:" : ["minbother_blue",float],
        "minbother_red:" : ["minbother_red",float],
        "vpixsizekms:" : ["vpixsizekms",float],
        "do_convolve_spectrum:" : ["do_convolve_spectrum",bool],
        "fwhm:": ["fwhm",float],
        "zqso:": ["zqso",float],
        "minlambda:": ["minlambda",float],
        "maxlambda:": ["minlambda",float],
        "zabsmin:": ["zabsmin",float],
        "zabsmax:": ["zabsmax",float],
        "nlyman:" : ["nlyman",int],
        "fzresol:": ["fzresol",float],
        "pixsize:": ["pixsize",float]
        }

parameters = read_params('dummy.par')


if parameters.use_snapshot_file:
    print("This is saying the use snapshot file is True")
