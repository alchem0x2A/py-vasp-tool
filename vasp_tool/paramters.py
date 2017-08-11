# Parameters that almost not changed
def_par = {"amix": 0.1,   # mixing parameters
           "bmix": 0.01,
           "prec": "Accurate",
           "encut": 800,
           "npar": 4,
           "kpar": 2,
           "nelm": 200,  # max SC steps
           "nelmin": 4,  # min SC steps
}

new_start = {"istart": 0,
             "icharg": 2}

restart = {"istart": 1,
           "icharg": 2}

g_smear = {"ismear": 0,
           "sigma": 0.01}

fine_conv = {"ediff": 1e-8,
             "ediffg": 1e-7}

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


par_relax = {**def_par, **new_start,
             **g_smear, **fine_conv,
             **nowrite_wave,
             **relax_all,
}

par_ground = {**def_par, **new_start,
              **g_smear, **fine_conv,
              **write_wave,
              **norelax}

par_hybrid = {**def_par, **restart,
              **g_smear, **fine_conv,
              **write_wave,
              **norelax
              "algo": "Damped",
              "precfock": "Normal"}

par_rpa = {**def_par, **restart,
           **g_smear, **fine_conv,
           **write_wave,
           **norelax,
           "loptics": True,
           "nbands": 56,
}

par_diag = {**def_par, **restart,
           **g_smear, **fine_conv,
           **write_wave,
           **norelax,
            "loptics": True,
            "algo": "Exact",
            "lpead": True,
            "nbands": 56,
            "nedos": 10**4,
}

par_gw0 = {**def_par, **restart,
           **g_smear, **fine_conv,
           **write_wave,
           **norelax,
           "algo": "GW0",
           "nelm": 1,           # single shot, no scGW
           "nelmin": 1,
           "noemga": 64,
           "omegamax": 64,
           "nbands": 56,
           "nedos": 10**4,
           "encutgw": 300,
           "lusew": True,
}

# only calculate WAVEDER
par_gw0_none = {**def_par, **restart,
                **g_smear, **fine_conv,
                **nowrite_wave,
                **norelax,
                "algo": "None",
                "nelm": 1,           # single shot, no scGW
                "nelmin": 1,
                "noemga": 64,
                "omegamax": 64,
                "nbands": 56,
                "nedos": 10**4,
                "encutgw": 300,
                "lusew": True,
                "loptics": True,
}

par_bse = {**def_par, **restart
           **g_smear, **fine_conv,
           **write_wave,      # maybe need to recalc
           **norelax,
           "algo": "BSE",
           "nelm": 1,            # single shot
           "noemga": 64,
           "omegamax": 64,
           "nbands": 56,
           "nbandsgw": 56,      # needs to be tweaked
           "nbandso": 4,        # needs to change to systems!
           "nbandsv": 8,
           "encutgw": 300,
           "ladder": True,
           "lhartree": True,
           "lusew": True,
           "antires": 1,
}

default_parameters = {"relax": par_relax,
                      "ground": par_ground,
                      "rpa": par_rpa,
                      "diag": par_diag,
                      "gw0": par_gw0,
                      "hybrid": par_hybrid,
                      "gw0_none": par_gw0_none,
                      "bse": par_bse,}
