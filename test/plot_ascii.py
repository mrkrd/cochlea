#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plot_tab(tab):
    if tab.ndim == 1:
        plt.plot(tab);
        plt.show()

    else:
        plt.imshow(tab.T, aspect="auto") #, cmap=cm.gray)
        plt.show()


def main():
    tab = np.loadtxt(sys.argv[1])
    plot_tab(tab)



if __name__ == "__main__":
    main()
