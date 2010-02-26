import numpy as np
import os

import thorns as th
import dsam

def par_dir(par_file):
    """
    Add directory path to par file name that is refered to the
    modules location.
    """
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, 'pars', par_file)



class AuditoryPeriphery(object):

    _train_type = [('typ', 'S3'),
                   ('cf', float),
                   ('id', int),
                   ('spikes', np.ndarray)]

    def __init__(self):
        pass


    def set_freq(self, freq):
        pass


    def _generate_anf(self, sg_type=None, anf_num=0):
        """
        Returns ANF ear module.
        """
        assert sg_type

        if sg_type == 'carney':
            anf = dsam.EarModule("An_SG_Carney")
            anf.read_pars(par_dir("anf_carney.par"))
        elif sg_type == 'binomial':
            anf = dsam.EarModule("An_SG_Binomial")
            anf.read_pars(par_dir("anf_binomial.par"))
        else:
            assert False

        anf.set_par("NUM_FIBRES", anf_num)

        return anf




    def _run_anf(self, anf_type, anf, fs, anf_num):
        """ Run spike generator several times and format the output. """

        # Explicit set of dt is required by some modules
        if anf.module_type.lower() == 'an_sg_binomial':
            anf.set_par("PULSE_DURATION", 1.1/fs)

        freq_map = self.get_freq_map()
        anf_trains = []
        for anf_id in range(anf_num):
            anf.run()
            anf_signal = anf.get_signal()
            anf_spikes = th.signal_to_spikes(fs, anf_signal)

            for freq, train in zip(freq_map, anf_spikes):
                anf_trains.append( (anf_type, freq, anf_id, train) )

        return anf_trains



    def run(self):
        """ Run the model """
        pass
