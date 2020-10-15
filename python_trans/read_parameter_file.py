class read_params(object):
    """
    Read the parameter file.
    """
    def __init__(self, par_file):
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

                    values = [(datatypes[i](data[i])) for i in range(len(data))]
                    self.__dict__[attr] = values if len(values) > 1 else values[0]
                except:
                    print(line)
                    
label_attr_map = {
        "ibdir:": ["ibdir",str],
        "datadir:": [ "datadir",str],
        "file_list:": ["file_list",str],
        "outputdir:": ["outputdir",str],
        "gimic:": ["gimic",int],
        "aurora:": ["aurora",int],
        "use_snapshot_file:": ["use_snapshot_file",int],
        "snap_base:": ["snap_base",str],
        "use_random_los:": ["use_random_los",int],
        "los_coordinates_file:": ["los_coordinates_file",str],
        "do_long_spectrum:": ["do_long_spectrum",int],
        "NoPecVel:": ["NoPecVel",int],
        "overwrite:": ["overwrite",int],
        "nspec:": ["nspec",int],
        "output_zspaceopticaldepthweighted_values:": ["output_zspaceopticaldepthweighted_values",int],
        "output_realspacemassweighted_values:": ["output_realspacemassweighted_values",int],
        "output_realspacenionweighted_values:": ["output_realspacenionweighted_values",int],
        "output_frequency:": ["output_frequency",int],
        "impose_eos:": ["impose_eos",int],
        "imposed_eos_T0:": ["imposed_eos_T0",float],
        "imposed_eos_gamma:": ["imposed_eos_gamma",float],
        "imposed_eos_maxod:": ["imposed_eos_maxod",float],
        "doH1:": ["doH1",int],
        "doHe2:": ["doHe2",int],
        "doC2:": ["doC2",int],
        "doC3:": ["doC3",int],
        "doC4:": ["doC4",int],
        "doN2:": ["doN2",int],
        "doN3:": ["doN3",int],
        "doN4:": ["doN4",int],
        "doN5:": ["doN5",int],
        "doO1:": ["doO1",int],
        "doO3:": ["doO3",int],
        "doO4:": ["doO4",int],
        "doO5:": ["doO5",int],
        "doO6:": ["doO6",int],
        "doO7:": ["doO7",int],
        "doMg2:": ["doMg2",int],
        "doNe8:": ["doNe8",int],
        "doAl2:": ["doAl2",int],
        "doAl3:": ["doAl3",int],
        "doSi2:": ["doSi2",int],
        "doSi3:": ["doSi3",int],
        "doSi4:": ["doSi4",int],
        "doS5:": ["doS5",int],
        "doFe2:": ["doFe2",int],
        "doFe3:": ["doFe3",int],
        "do21cm:": ["do21cm",int],
        "doall:": ["doall",int],
        "ibfactor:": ["ibfactor",float],
        "ibfactor_he_reionization:": ["ibfactor_he_reionization",int],
        "use_fitted_ibfactor:": ["use_fitted_ibfactor",int],
        "ibfactor_he_reionization:": ["ibfactor_he_reionization",int],
        "use_maxdens_above_zmax:": ["use_maxdens_above_zmax",int],
        "modify_metallicity:": ["modify_metallicity",int],
        "maxz_rel:": ["maxz_rel",float],
        "scale_simulation_abundances:": ["scale_simulation_abundances",int],
        "z_rel:": ["z_rel",float],
        "impose_z_rho_relation:": ["impose_z_rho_relation",int],
        "z_index:": ["z_index",float],
        "z_mean:": ["z_mean",float],
        "log_normal_scatter:": ["log_normal_scatter",int],
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
        "read_part_ids_from_file:": ["read_part_ids_from_file",int],
        "flagged_particle_metallicity:": ["flagged_particle_metallicity",float],
        "generate_noise:": ["generate_noise",int],
        "use_noise_file:": ["use_noise_file",float],
        "noisefile:": ["noisefile",str],
        "sigtonoise:": ["sigtonoise",float],
        "minnoise:": ["minnoise",float],
        "minbother_blue:" : ["minbother_blue",float],
        "minbother_red:" : ["minbother_red",float],
        "vpixsizekms:" : ["vpixsizekms",float],
        "do_convolve_spectrum:" : ["do_convolve_spectrum",int],
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
