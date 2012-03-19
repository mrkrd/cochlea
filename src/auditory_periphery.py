from __future__ import division

import numpy as np
import os
import scipy.signal as dsp

import thorns as th

def par_dir(par_file):
    """
    Add directory path to par file name that is refered to the
    modules location.
    """
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, 'pars', par_file)



class AuditoryPeriphery(object):

    _anf_dtype = [('spikes', np.ndarray),
                  ('duration', float),
                  ('cf', float),
                  ('type', '|S3'),
                  ('index', int)]

    def __init__(self):
        pass


    def set_freq(self, freq):
        pass




    def _run_anf(self, anf_type, sg_module, fs, anf_num, accumulate):
        """ Run spike generator several times and format the output. """

        sg_module.set_par('pulse_duration', 1.1/fs)

        assert not accumulate, "Deprecated, `accumulate' is now disabled"

        if accumulate:
            run_num = 1
            sg_module.set_par('num_fibres', anf_num)
        else:
            run_num = anf_num
            sg_module.set_par('num_fibres', 1)

        freq_map = self.get_freq_map()
        trains = []
        for anf_idx in range(run_num):
            sg_module.run()

            anf_signal = sg_module.get_signal()
            anf_spikes = th.signal_to_trains(anf_signal, fs)

            for cf, train in zip(freq_map, anf_spikes):
                trains.append( (train['spikes'],
                                train['duration'],
                                cf,
                                anf_type,
                                anf_idx) )

        return trains



    def run(self):
        """ Run the model """
        pass




def run_human_me_filter_for_zilany2009(signal, fs):
    assert fs > 40e3

    signal_fft = np.fft.fft(signal)
    freqs = np.fft.fftfreq(len(signal), d=1/fs)


    me_filter_response = np.loadtxt(data_dir("me_filter_response.txt"))
    me_filter_freqs = np.loadtxt(data_dir("me_filter_freqs.txt"))
    fmin = me_filter_freqs.min()
    fmax = me_filter_freqs.max()


    # Convert dB to amplitudae ratio
    response_ratio = 10 ** (me_filter_response / 20)


    # Interpolate the filter to fit signal's FFT
    band = ((np.abs(freqs)>=fmin) & (np.abs(freqs)<=fmax))
    band_len = np.sum( band )
    ratio_interp = dsp.resample(response_ratio, band_len/2)
    ratio_interp = np.concatenate( (ratio_interp, ratio_interp[::-1]) )


    # Apply the filter
    signal_fft[band] *= ratio_interp

    signal_fft[ np.logical_not(band) ] = 0


    filtered = np.fft.ifft( signal_fft )
    filtered = np.array( filtered.real )

    return filtered


def data_dir(par_file):
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, 'data', par_file)
