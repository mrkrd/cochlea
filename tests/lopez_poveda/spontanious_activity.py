# Author: Marek Rudnicki
# Time-stamp: <2009-09-16 20:37:02 marek>
#
# Description: Scanning parameters space of IHC module to fit
# spontanious activity of HSR fibers

from collections import namedtuple


import numpy as np
import matplotlib.pyplot as plt

import dsam
from cochlea.auditory_periphery import par_dir


Pars = namedtuple('Pars', 'beta_ca gamma_ca gmax_ca conc_thresh_ca')

class IHCRP_IHC(object):
    def __init__(self):
        self.ihcrp = dsam.EarModule("IHCRP_LopezPoveda")
        self.ihc = dsam.EarModule("IHC_Meddis2000")
        self.ihc.read_pars(par_dir("ihc_hsr_Meddis2002.par"))

        dsam.connect(self.ihcrp, self.ihc)

        self.fs = 100000.0
        self.signal = np.zeros( np.floor(0.001 * self.fs) )


    def run(self):
        self.ihcrp.run(self.fs, self.signal)
        self.ihc.run()

        return self.ihc.get_signal().mean()


    def set_pars(self, pars):
        # self.ihc.set_par("DIAG_MODE", 'screen')

        self.ihc.set_par("BETA_CA", pars.beta_ca)
        self.ihc.set_par("GAMMA_CA", pars.gamma_ca)
        self.ihc.set_par("GMAX_CA", pars.gmax_ca)
        self.ihc.set_par("CONC_THRESH_CA", pars.conc_thresh_ca)



def main():

    ihc = IHCRP_IHC()

    beta_ca_list = [400] #np.arange(1, 2000, 10)
    gamma_ca_list = [130] #np.arange(1, 200, 1)
    gmax_ca_list = np.arange(8e-11, 8e-9, 8e-11)
    conc_thresh_ca_list = np.arange(4e-14, 4e-12, 4e-14)

    pars = ( Pars(beta_ca, gamma_ca, gmax_ca, conc_thresh_ca)
             for beta_ca in beta_ca_list
             for gamma_ca in gamma_ca_list
             for gmax_ca in gmax_ca_list
             for conc_thresh_ca in conc_thresh_ca_list )


    psp_list = []
    for par_set in pars:
        ihc.set_pars(par_set)
        psp = ihc.run()
        psp_list.append(psp)

    psp_list = np.asarray(psp_list)
    psp_list = psp_list.reshape( (len(gmax_ca_list),
                                  len(conc_thresh_ca_list)) )

    plt.imshow(psp_list, aspect='auto')
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    main()
