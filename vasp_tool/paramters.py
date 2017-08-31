# The ugly python2-3.4 way
def merge_dict(*presets, **keywords):
    key_presets = []
    for p in presets:
        key_presets += list(p.items())
    return dict(key_presets + list(keywords.items()))
    
# Parameters that almost not changed
def_par = {"amix": 0.1,   # mixing parameters
           "bmix": 0.01,
           "prec": "Accurate",
           "encut": 800,
           "nelm": 200,  # max SC steps
           "nelmin": 5,  # min SC steps
}

new_start = {"istart": 0,
             "icharg": 2}

restart = {"istart": 1,
           "icharg": 2}

g_smear = {"ismear": 0,
           "sigma": 0.01}

fine_conv = {"ediff": 1e-8,
             "ediffg": 1e-7, }

rough_conv = {"ediff": 1e-5,
              "ediffg": 1e-4}

write_wave = {"lwave": True,
              "lcharg": True,
              "lvtot": False,}

nowrite_wave = {"lwave": False,
              "lcharg": False,
              "lvtot": False,}

relax_all = {"ibrion": 2,
             "isif": 3,
             "nsw": 500,
}

norelax = {"ibrion": -1,
           "isif": 1,           # still calculate extern pressure
           "nsw": 1,
}


# par_relax = {**def_par, **new_start,
             # **g_smear, **fine_conv,
             # **nowrite_wave,
             # **relax_all,
# }

par_relax = merge_dict(def_par, new_start,
                       g_smear, fine_conv,
                       nowrite_wave,
                       relax_all,
                       **{"npar": 4})


par_ground = merge_dict(def_par, new_start,
                        g_smear, fine_conv,
                        write_wave,
                        norelax,
                        **{"npar": 4})

# par_ground = {**def_par, **new_start,
              # **g_smear, **fine_conv,
              # **write_wave,
              # **norelax}

par_hybrid = merge_dict(def_par, restart,
              g_smear, fine_conv,
              write_wave,
              norelax,
              **{"algo": "Damped",
                 "precfock": "Normal",
                 "npar": 4})

# par_hybrid = {**def_par, **restart,
#               **g_smear, **fine_conv,
#               **write_wave,
#               **norelax,
#               "algo": "Damped",
#               "precfock": "Normal"}

par_rpa = merge_dict(def_par, restart,
                     g_smear, fine_conv,
                     write_wave,
                     norelax,
                     **{"loptics": True,
                        "nbands": 56,
                        "npar": 4})

par_diag = merge_dict(def_par, restart,
           g_smear, fine_conv,
           write_wave,
           norelax,
            **{"loptics": True,
               "algo": "Exact",
               "lpead": True,
               "nbands": 56,
               "nedos": 10**4,
               "npar": 4})

par_gw0 = merge_dict(def_par, restart,
                     g_smear, fine_conv,
                     write_wave,
                     norelax,
                     **{"algo": "GW0",
                        "nelm": 1,           # single shot, no scGW
                        "nelmin": 1,
                        "nomega": 64,
                        "nbands": 56,
                        "nedos": 10**4,
                        "encutgw": 300,
                        "lpead": True,
                        "lrpa": True,
                        "kpar": 2})

# only calculate WAVEDER
par_gw0_none = merge_dict(def_par, restart,
                          g_smear, fine_conv,
                          nowrite_wave,
                          norelax,
                          **{"algo": "None",
                             "nelm": 1,           # single shot, no scGW
                             "nelmin": 1,
                             "nomega": 64,
                             "nbands": 56,
                             "nedos": 10**4,
                             "encutgw": 300,
                             "lpead": True,
                             "loptics": True,
                             "lrpa": True,
                             "kpar": 2})

par_bse = merge_dict(def_par, restart,
                     g_smear, fine_conv,
                     write_wave,      # maybe need to recalc
                     norelax,
                     **{"algo": "BSE",
                        "nelm": 1,            # single shot
                        "nomega": 64,
                        "omegamax": 64,
                        "nbands": 56,
                        "nbandsgw": 56,      # needs to be tweaked
                        "nbandso": 4,        # needs to change to systems!
                        "nbandsv": 8,
                        "lpead": True,
                        "encutgw": 300,
                        "ladder": True,
                        "lhartree": True,
                        "lusew": True,
                        "lrpa": True,  # No kpar used for BSE calculation!
                     })

par_bandstruct_DFT = merge_dict(def_par, restart,
                                g_smear, fine_conv,
                                write_wave,
                                norelax,
                                **{"icharg": 11,  # keep previous charge density calculated
                                   "lorbit": 11,  # DOCAR + PROCAR
                                })


# par_rpa = {**def_par, **restart,
#            **g_smear, **fine_conv,
#            **write_wave,
#            **norelax,
#            "loptics": True,
#            "nbands": 56,
# }

# par_diag = {**def_par, **restart,
#            **g_smear, **fine_conv,
#            **write_wave,
#            **norelax,
#             "loptics": True,
#             "algo": "Exact",
#             "lpead": True,
#             "nbands": 56,
#             "nedos": 10**4,
# }

# par_gw0 = {**def_par, **restart,
#            **g_smear, **fine_conv,
#            **write_wave,
#            **norelax,
#            "algo": "GW0",
#            "nelm": 1,           # single shot, no scGW
#            "nelmin": 1,
#            "nomega": 64,
#            "omegamax": 64,
#            "nbands": 56,
#            "nedos": 10**4,
#            "encutgw": 300,
#            "lusew": True,
# }

# # only calculate WAVEDER
# par_gw0_none = {**def_par, **restart,
#                 **g_smear, **fine_conv,
#                 **nowrite_wave,
#                 **norelax,
#                 "algo": "None",
#                 "nelm": 1,           # single shot, no scGW
#                 "nelmin": 1,
#                 "nomega": 64,
#                 "omegamax": 64,
#                 "nbands": 56,
#                 "nedos": 10**4,
#                 "encutgw": 300,
#                 "lusew": True,
#                 "loptics": True,
# }

# par_bse = {**def_par, **restart,
#            **g_smear, **fine_conv,
#            **write_wave,      # maybe need to recalc
#            **norelax,
#            "algo": "BSE",
#            "nelm": 1,            # single shot
#            "nomega": 64,
#            "omegamax": 64,
#            "nbands": 56,
#            "nbandsgw": 56,      # needs to be tweaked
#            "nbandso": 4,        # needs to change to systems!
#            "nbandsv": 8,
#            "encutgw": 300,
#            "ladder": True,
#            "lhartree": True,
#            "lusew": True,
#            "antires": 1,
# }

default_parameters = {"relax": par_relax,
                      "ground": par_ground,
                      "rpa": par_rpa,
                      "diag": par_diag,
                      "gw0": par_gw0,
                      "hybrid": par_hybrid,
                      "gw0_none": par_gw0_none,
                      "bse": par_bse,
                      "bs_DFT": par_bandstruct_DFT,
                      "bs_hybrid": par_hybrid}

