from sumner2002 import Sumner2002
from holmberg2008 import Holmberg2008
from sumner2002_vesicles import Sumner2002_Vesicles
from lopez_poveda2006 import LopezPoveda2006
from dsam import set_dB_SPL
from pycat import Carney2009


# TODO: (re)move main
def main():
    import cochlea

    import cProfile
    import matplotlib.pyplot as plt
    import stuff

    print "start:"
    ear = cochlea.LopezPoveda2006(freq=(tw.real_freq_map[0],
                                   tw.real_freq_map[99],
                                   100),
                             animal='human')
    # earH = cochlea.Holmberg2008()

    fs = 48000.0
    fstim = tw.bm_pars.real_freq_map[38]
    print "Fstim =", fstim
    t = np.arange(0, 0.1, 1.0/fs)
    s = np.sin(2 * np.pi * t * fstim)
    s = stuff.set_dB_SPL(30, s)
    s = s * np.hanning(len(s))

    print "Sumner2002 running..."
    hsr, msr, lsr = ear.run(fs, s, output_format='signals')
    print "done"
    # print "Holmberg2008 running..."
    # hsrH, msrH, lsrH = earH.run(fs, s, output_format='signals')
    # print "done"

    plt.imshow(hsr.T, aspect='auto', interpolation=None)
    plt.colorbar()
    plt.show()

    # plt.imshow(hsrH.T, aspect='auto', interpolation=None)
    # plt.colorbar()
    # plt.show()

    print "done."


if __name__ == "__main__":
    import cProfile

    # cProfile.run('main()')
    main()
