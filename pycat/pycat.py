import _pycat

def run_ihc(sound, cf, fs, cohc=1, cihc=1):
    # TODO: test input parameters
    vihc = _pycat.run_ihc(sound, cf, fs, cohc, cihc)

    return vihc



def run_synapse(vihc, cf, nrep, fs, anf_type='hsr', implnt='actual'):
    # TODO: test input pars

    anf_map = {'hsr': 3,
               'msr': 2,
               'lsr': 1}

    implnt_map = {'actual': 1,
                  'approx': 0}

    synout, psth = _pycat.run_synapse(vihc, cf, nrep, fs,
                                      anf_map[anf_type], implnt_map[implnt]);


    return synout, psth


